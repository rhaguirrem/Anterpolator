#!/usr/bin/env python3
"""Standalone grid-to-BMF utility for Anterpolator.

This tool is intentionally standalone so it can later be wired into the
Anterpolator Operations panel without changing the core interpolation flow.

Subcommands:
- export: CSV regular grid -> Vulcan BMF, TBMS config-text BMF, or experimental TBMS2.0 container
- inspect: quick binary inspection for TBMS2.0-style BMF files
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
import os
import re
import struct
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd


TBMS_SIGNATURE = b"TBMS2.0\x00"
TBMS_HEADER_SIZE = 128
TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET = 2056
TBMS_EXPERIMENTAL_DIR_MAGIC = b"ATBXDIR\x00"
TBMS_EXPERIMENTAL_SECTION_ENTRY = struct.Struct("<16sQQII8x")
TBMS_EXPERIMENTAL_SECTION_FLAGS = {
    "json_utf8": 1,
    "int32_le": 2,
    "float64_le": 3,
}
MAX_DENSE_EXPORT_BYTES = int(os.environ.get("ANTERPOLATOR_BMF_MAX_DENSE_BYTES", str(8 * 1024 ** 3)))
SUPPORTED_EXPORT_BACKENDS = {"tbms-config-text", "tbms-experimental", "vulcan"}
DENSE_EXPORT_BACKENDS = {"tbms-config-text", "vulcan"}
MAX_SELECTED_CSV_OBJECT_BYTES = int(os.environ.get("ANTERPOLATOR_BMF_MAX_SELECTED_CSV_BYTES", str(8 * 1024 ** 3)))
MAX_SOURCE_ROW_PREP_BYTES = int(os.environ.get("ANTERPOLATOR_BMF_MAX_SOURCE_ROW_PREP_BYTES", str(1 * 1024 ** 3)))
_LOGGED_LEAPFROG_METADATA_SIGNATURES: set[tuple[str, str, tuple[str, ...]]] = set()
GRID_INDEX_COLUMNS = ("grid_i", "grid_j", "grid_k")
GEOMETRY_EXTENT_COLUMNS = ("__lower_x", "__upper_x", "__lower_y", "__upper_y", "__lower_z", "__upper_z")
GEOMETRY_SIZE_COLUMN_ROLES = ("dx", "dy", "dz")


def _normalize_column_token(value: object) -> str:
    return str(value or "").strip()


def _find_header_column(header_names: Sequence[str], candidates: Sequence[str]) -> str | None:
    exact = {str(name): str(name) for name in header_names}
    for candidate in candidates:
        token = _normalize_column_token(candidate)
        if token in exact:
            return exact[token]
    lowered: dict[str, str | None] = {}
    for name in header_names:
        key = str(name).strip().lower()
        if key in lowered:
            lowered[key] = None
        else:
            lowered[key] = str(name)
    for candidate in candidates:
        match = lowered.get(_normalize_column_token(candidate).lower())
        if match:
            return match
    return None


def _normalize_geometry_column_spec(
    spec: Sequence[str] | Mapping[str, str] | str | None,
    roles: Sequence[str],
    label: str,
) -> dict[str, str]:
    if spec is None:
        return {}
    if isinstance(spec, str):
        values = [part.strip() for part in spec.split(",") if part.strip()]
        if not values:
            return {}
        spec = values
    if isinstance(spec, Mapping):
        resolved = {role: _normalize_column_token(spec.get(role)) for role in roles}
    else:
        values = [_normalize_column_token(value) for value in spec]
        if not values:
            return {}
        if len(values) != len(roles):
            raise ValueError(f"{label} must provide exactly {len(roles)} column names.")
        resolved = dict(zip(roles, values))
    missing = [role for role in roles if not resolved.get(role)]
    if missing and len(missing) != len(roles):
        raise ValueError(f"{label} must provide all required columns; missing: {missing}")
    if len(missing) == len(roles):
        return {}
    return resolved


def _resolve_geometry_column_map(
    header_names: Sequence[str],
    spec: Sequence[str] | Mapping[str, str] | str | None,
    roles: Sequence[str],
    label: str,
) -> dict[str, str]:
    requested = _normalize_geometry_column_spec(spec, roles, label)
    if not requested:
        return {}
    resolved: dict[str, str] = {}
    missing = []
    for role, column_name in requested.items():
        match = _find_header_column(header_names, [column_name])
        if match is None:
            missing.append(column_name)
        else:
            resolved[role] = match
    if missing:
        raise ValueError(f"{label} column(s) not found in CSV: {missing}")
    return resolved


def _detect_geometry_size_columns(header_names: Sequence[str]) -> dict[str, str]:
    candidates = {
        "dx": ("dx", "dX", "DX", "size_x", "block_size_x", "width_x"),
        "dy": ("dy", "dY", "DY", "size_y", "block_size_y", "width_y"),
        "dz": ("dz", "dZ", "DZ", "size_z", "block_size_z", "width_z", "height_z"),
    }
    resolved: dict[str, str] = {}
    for role in GEOMETRY_SIZE_COLUMN_ROLES:
        match = _find_header_column(header_names, candidates[role])
        if match is None:
            return {}
        resolved[role] = match
    return resolved


def _detect_geometry_extent_columns(header_names: Sequence[str]) -> dict[str, str]:
    canonical: dict[str, str] = {}
    for role in GEOMETRY_EXTENT_COLUMNS:
        match = _find_header_column(header_names, [role])
        if match is None:
            canonical = {}
            break
        canonical[role] = match
    if canonical:
        return canonical
    candidates = {
        "__lower_x": ("lower_x", "x_lower", "min_x", "x_min", "from_x"),
        "__upper_x": ("upper_x", "x_upper", "max_x", "x_max", "to_x"),
        "__lower_y": ("lower_y", "y_lower", "min_y", "y_min", "from_y"),
        "__upper_y": ("upper_y", "y_upper", "max_y", "y_max", "to_y"),
        "__lower_z": ("lower_z", "z_lower", "min_z", "z_min", "from_z"),
        "__upper_z": ("upper_z", "z_upper", "max_z", "z_max", "to_z"),
    }
    resolved: dict[str, str] = {}
    for role, role_candidates in candidates.items():
        match = _find_header_column(header_names, role_candidates)
        if match is None:
            return {}
        resolved[role] = match
    return resolved


def _emit_progress(progress_callback, value: int, maximum: int = 100, message: str = "") -> None:
    if progress_callback is None:
        return
    try:
        progress_callback(max(0, min(int(value), int(maximum))), int(maximum), str(message or ""))
    except Exception:
        pass


def _scale_progress(progress_start: int, progress_end: int, fraction: float) -> int:
    fraction = max(0.0, min(float(fraction), 1.0))
    return int(round(progress_start + (progress_end - progress_start) * fraction))


def _format_byte_size(size_bytes: int | float) -> str:
    value = float(size_bytes)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(value) < 1024.0 or unit == "TiB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0
    return f"{value:.1f} TiB"


def detect_csv_delimiter(path: str | Path) -> str:
    path = Path(path)
    if not path.is_file():
        return ","
    try:
        lines = []
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.strip() and not line.lstrip().startswith(("#", "//")):
                    lines.append(line)
                if len(lines) >= 3:
                    break
    except OSError:
        return ","
    sample = "".join(lines)
    counts = {delimiter: sample.count(delimiter) for delimiter in [",", ";", "\t", "|"]}
    delimiter = max(counts, key=counts.get)
    return delimiter if counts[delimiter] > 0 else ","


def parse_header_line(path: str | Path, delimiter: str, line_number: int) -> list[str]:
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"File not found: {path}")
    line_number = max(int(line_number or 1), 1)
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for current, line in enumerate(handle, start=1):
            if current == line_number:
                tokens = [token.strip() for token in line.rstrip("\r\n").split(delimiter)]
                tokens = [token for token in tokens if token]
                if not tokens:
                    raise ValueError(f"Parsed header line {line_number} in '{path.name}' produced no tokens.")
                return tokens
    raise ValueError(f"Header line {line_number} exceeds total lines in file '{path.name}'.")


def resolve_effective_csv_header_line(path: str | Path, configured_line: int = 1) -> int:
    path = Path(path)
    line_number = max(int(configured_line or 1), 1)
    if not path.is_file():
        raise ValueError(f"File not found: {path}")
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for current_line_number, line_text in enumerate(handle, start=1):
            if current_line_number < line_number:
                continue
            stripped = line_text.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if current_line_number != line_number:
                print(
                    f"Header line {line_number} in '{path.name}' is metadata/comment; "
                    f"using first data header line {current_line_number}."
                )
            return current_line_number
    raise ValueError(f"Could not find a non-comment CSV header line in '{path.name}'.")


def parse_effective_header_line(path: str | Path, delimiter: str, line_number: int) -> list[str]:
    metadata = parse_leapfrog_block_metadata(path)
    log_leapfrog_metadata_summary(path, metadata, context="header scan")
    return parse_header_line(path, delimiter, resolve_effective_csv_header_line(path, line_number))


def _parse_metadata_numeric_values(text: object, numeric_type=float, stop_at_equals: bool = False) -> list[object]:
    value_text = str(text or "")
    if stop_at_equals:
        value_text = value_text.split("=", 1)[0]
    values = re.findall(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?", value_text)
    return [numeric_type(value) for value in values]


def parse_leapfrog_block_metadata(path: str | Path, max_lines: int = 100) -> dict[str, object]:
    path = Path(path)
    if not path.is_file():
        return {}

    metadata: dict[str, object] = {}
    raw_lines = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line_number, line_text in enumerate(handle, start=1):
            if line_number > max_lines:
                break
            stripped = line_text.strip()
            if not stripped:
                continue
            if not stripped.startswith("#"):
                break

            content = stripped.lstrip("#").strip()
            raw_lines.append(content)
            if ":" not in content:
                if "title" not in metadata and content:
                    metadata["title"] = content
                continue

            key, raw_value = content.split(":", 1)
            normalized_key = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")
            value = raw_value.strip()

            if normalized_key == "encoding":
                metadata["encoding"] = value
            elif normalized_key == "rotation_type":
                metadata["rotation_type"] = value
            elif normalized_key in {"azimuth", "dip", "pitch"}:
                numbers = _parse_metadata_numeric_values(value)
                if numbers:
                    metadata[f"{normalized_key}_degrees"] = float(numbers[0])
            elif normalized_key == "parent_block_size":
                numbers = _parse_metadata_numeric_values(value)
                if len(numbers) >= 3:
                    metadata["parent_block_size"] = [float(number) for number in numbers[:3]]
            elif normalized_key == "size_in_parent_blocks":
                numbers = _parse_metadata_numeric_values(value, numeric_type=int, stop_at_equals=True)
                if len(numbers) >= 3:
                    metadata["size_in_parent_blocks"] = [int(number) for number in numbers[:3]]
                total_values = _parse_metadata_numeric_values(value.split("=", 1)[1], numeric_type=int) if "=" in value else []
                if total_values:
                    metadata["parent_block_count"] = int(total_values[0])
            elif normalized_key in {"minimum_parent_centroid", "maximum_parent_centroid", "minimum_corner", "maximum_corner"}:
                numbers = _parse_metadata_numeric_values(value)
                if len(numbers) >= 3:
                    metadata[normalized_key] = [float(number) for number in numbers[:3]]
            elif normalized_key in {"sub_blocks", "subblocks"}:
                parts = value.split()
                if parts:
                    metadata["subblock_scheme"] = parts[0]
                numbers = _parse_metadata_numeric_values(value, numeric_type=int)
                if len(numbers) >= 3:
                    metadata["subblock_factors"] = [int(number) for number in numbers[:3]]

    if raw_lines:
        metadata["raw_lines"] = raw_lines
    return metadata


def _format_metadata_vector(metadata: Mapping[str, object], key: str) -> str | None:
    values = metadata.get(key)
    if not values:
        return None
    return "(" + ", ".join(f"{float(value):g}" for value in values) + ")"


def log_leapfrog_metadata_summary(path: str | Path, metadata: Mapping[str, object], context: str = "") -> None:
    if not metadata:
        return

    recognized_fields = [
        key for key in (
            "rotation_type",
            "azimuth_degrees",
            "dip_degrees",
            "pitch_degrees",
            "parent_block_size",
            "size_in_parent_blocks",
            "minimum_parent_centroid",
            "maximum_parent_centroid",
            "minimum_corner",
            "maximum_corner",
            "subblock_scheme",
            "subblock_factors",
        )
        if key in metadata
    ]
    if not recognized_fields:
        return

    path = Path(path)
    signature = (str(path.resolve()), str(context or "default"), tuple(recognized_fields))
    if signature in _LOGGED_LEAPFROG_METADATA_SIGNATURES:
        return
    _LOGGED_LEAPFROG_METADATA_SIGNATURES.add(signature)

    summary_parts = []
    parent_size = _format_metadata_vector(metadata, "parent_block_size")
    if parent_size:
        summary_parts.append(f"parent block size={parent_size}")
    grid_size = metadata.get("size_in_parent_blocks")
    if grid_size:
        summary_parts.append(f"parent grid={tuple(int(value) for value in grid_size)}")
    min_corner = _format_metadata_vector(metadata, "minimum_corner")
    if min_corner:
        summary_parts.append(f"minimum corner={min_corner}")
    min_centroid = _format_metadata_vector(metadata, "minimum_parent_centroid")
    if min_centroid:
        summary_parts.append(f"minimum parent centroid={min_centroid}")
    if any(key in metadata for key in ("azimuth_degrees", "dip_degrees", "pitch_degrees")):
        summary_parts.append(
            "rotation=(azimuth {azimuth:g}, dip {dip:g}, pitch {pitch:g})".format(
                azimuth=float(metadata.get("azimuth_degrees", 0.0) or 0.0),
                dip=float(metadata.get("dip_degrees", 0.0) or 0.0),
                pitch=float(metadata.get("pitch_degrees", 0.0) or 0.0),
            )
        )
    subblock_factors = metadata.get("subblock_factors")
    if subblock_factors:
        scheme = metadata.get("subblock_scheme") or "sub-blocks"
        summary_parts.append(f"{scheme}={tuple(int(value) for value in subblock_factors)}")

    label = f" ({context})" if context else ""
    print(f"Recognized Leapfrog block metadata{label} in '{path.name}': " + "; ".join(summary_parts))


def infer_bmf_export_geometry_from_metadata(
    metadata: Mapping[str, object],
    regularize_to_base_block: bool = False,
    dense_regular_grid: bool = False,
) -> dict[str, object]:
    if not metadata:
        return {}

    parent_block_size = metadata.get("parent_block_size")
    subblock_factors = metadata.get("subblock_factors")
    cell_size = None
    source = ""
    if (regularize_to_base_block or dense_regular_grid) and parent_block_size:
        cell_size = [float(value) for value in parent_block_size]
        source = "parent_block_size"
    elif parent_block_size and subblock_factors:
        factors = [int(value) for value in subblock_factors]
        if len(factors) == 3 and all(value > 0 for value in factors):
            cell_size = [float(size) / float(factor) for size, factor in zip(parent_block_size, factors)]
            source = "parent_block_size/subblock_factors"
    elif parent_block_size:
        cell_size = [float(value) for value in parent_block_size]
        source = "parent_block_size"

    if not cell_size or len(cell_size) != 3 or any(value <= 0 for value in cell_size):
        return {}

    origin = None
    origin_source = ""
    if metadata.get("minimum_corner"):
        origin = [float(value) for value in metadata["minimum_corner"]]
        origin_source = "minimum_corner"
    elif metadata.get("minimum_parent_centroid") and parent_block_size:
        origin = [
            float(centroid) - 0.5 * float(size)
            for centroid, size in zip(metadata["minimum_parent_centroid"], parent_block_size)
        ]
        origin_source = "minimum_parent_centroid-parent_block_size/2"

    result: dict[str, object] = {
        "cell_size": cell_size,
        "cell_size_source": source,
    }
    if dense_regular_grid and parent_block_size and subblock_factors and not regularize_to_base_block:
        result["dense_metadata_note"] = (
            "Detected sub-block metadata, but dense tbms-config-text/Vulcan output uses the parent block size by default. "
            "Using the minimum sub-block size would create the full dense sub-cell lattice for the entire model extent."
        )
    if origin is not None and len(origin) == 3:
        result["origin"] = origin
        result["origin_source"] = origin_source
    if any(key in metadata for key in ("azimuth_degrees", "dip_degrees", "pitch_degrees", "rotation_type")):
        result["rotation_metadata_detected"] = True
        result["rotation_note"] = "Rotation metadata is preserved in the export summary, but this BMF writer emits an axis-aligned grid."
    return result


def _count_csv_data_rows(
    path: Path,
    header_line: int = 1,
    progress_callback=None,
    progress_start: int = 5,
    progress_end: int = 8,
) -> int | None:
    try:
        file_size = max(int(path.stat().st_size), 1)
        total_lines = 0
        bytes_read = 0
        last_byte = b""
        report_step = max(file_size // 100, 1)
        next_report = report_step
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(8 * 1024 * 1024)
                if not chunk:
                    break
                bytes_read += len(chunk)
                total_lines += chunk.count(b"\n")
                last_byte = chunk[-1:]
                if bytes_read >= next_report or bytes_read >= file_size:
                    _emit_progress(
                        progress_callback,
                        _scale_progress(progress_start, progress_end, bytes_read / file_size),
                        100,
                        f"Scanning CSV rows {bytes_read / file_size:.0%}...",
                    )
                    next_report = bytes_read + report_step
        if file_size > 0 and last_byte not in (b"", b"\n", b"\r"):
            total_lines += 1
        return max(total_lines - max(int(header_line or 1), 1), 0)
    except Exception:
        return None


def _read_csv_with_chunk_progress(
    path: Path,
    total_rows: int | None,
    progress_callback=None,
    progress_start: int = 8,
    progress_end: int = 30,
    chunksize: int = 250_000,
    **read_kwargs,
) -> pd.DataFrame:
    chunks: List[pd.DataFrame] = []
    rows_read = 0
    reader = pd.read_csv(path, chunksize=chunksize, **read_kwargs)
    for chunk in reader:
        chunks.append(chunk)
        rows_read += len(chunk)
        if total_rows and total_rows > 0:
            fraction = min(rows_read / total_rows, 1.0)
            message = f"Reading CSV rows {rows_read:,}/{total_rows:,} ({fraction:.0%})..."
        else:
            fraction = 0.5
            message = f"Reading CSV rows {rows_read:,}..."
        _emit_progress(progress_callback, _scale_progress(progress_start, progress_end, fraction), 100, message)
    if chunks:
        _emit_progress(progress_callback, progress_end, 100, f"Combining {len(chunks):,} CSV chunks into a table...")
        return pd.concat(chunks, ignore_index=True)
    return pd.DataFrame()


def _build_csv_read_names(header_names: Sequence[str]) -> list[str]:
    seen_names: dict[str, int] = {}
    final_names = []
    for index, name in enumerate(header_names):
        base_name = str(name or "").strip() or f"Column_{index + 1}"
        count = seen_names.get(base_name, 0)
        seen_names[base_name] = count + 1
        final_names.append(base_name if count == 0 else f"{base_name}_{count + 1}")
    return final_names


def _iter_csv_chunks_for_export(
    path: Path,
    delimiter: str | None,
    header_line: int,
    usecols: Sequence[str] | None = None,
    chunksize: int = 250_000,
    progress_callback=None,
    progress_start: int = 5,
    progress_end: int = 30,
    progress_label: str = "Reading CSV rows",
    total_rows: int | None = None,
) -> Iterable[pd.DataFrame]:
    effective_header_line = resolve_effective_csv_header_line(path, header_line)
    csv_delimiter = delimiter or detect_csv_delimiter(path)
    final_names = _build_csv_read_names(parse_header_line(path, csv_delimiter, effective_header_line))
    read_kwargs = {
        "sep": csv_delimiter,
        "header": None,
        "names": final_names,
        "skiprows": effective_header_line,
        "usecols": list(usecols) if usecols else None,
        "dtype": str,
        "comment": "#",
        "chunksize": chunksize,
        "low_memory": True,
        "memory_map": True,
    }
    rows_read = 0
    try:
        reader = pd.read_csv(path, **read_kwargs)
        for chunk in reader:
            rows_read += len(chunk)
            if total_rows and total_rows > 0:
                fraction = min(rows_read / total_rows, 1.0)
                message = f"{progress_label} {rows_read:,}/{total_rows:,} ({fraction:.0%})..."
            else:
                fraction = 0.5
                message = f"{progress_label} {rows_read:,}..."
            _emit_progress(progress_callback, _scale_progress(progress_start, progress_end, fraction), 100, message)
            yield chunk
    except (MemoryError, pd.errors.ParserError) as exc:
        if "out of memory" not in str(exc).lower() and not isinstance(exc, MemoryError):
            raise
        read_kwargs.pop("low_memory", None)
        read_kwargs.pop("memory_map", None)
        read_kwargs["engine"] = "python"
        rows_read = 0
        reader = pd.read_csv(path, **read_kwargs)
        for chunk in reader:
            rows_read += len(chunk)
            if total_rows and total_rows > 0:
                fraction = min(rows_read / total_rows, 1.0)
                message = f"{progress_label} {rows_read:,}/{total_rows:,} ({fraction:.0%})..."
            else:
                fraction = 0.5
                message = f"{progress_label} {rows_read:,}..."
            _emit_progress(progress_callback, _scale_progress(progress_start, progress_end, fraction), 100, message)
            yield chunk


def _auto_read_csv(
    path: Path,
    delimiter: str | None = None,
    header_line: int = 1,
    usecols: Sequence[str] | None = None,
    progress_callback=None,
    progress_start: int = 5,
    progress_end: int = 30,
) -> pd.DataFrame:
    """Read CSV with optional explicit delimiter and header-line selection."""
    effective_header_line = resolve_effective_csv_header_line(path, header_line)
    csv_delimiter = delimiter or detect_csv_delimiter(path)
    header_names = parse_header_line(path, csv_delimiter, effective_header_line)
    final_names = _build_csv_read_names(header_names)
    read_kwargs = {
        "header": None,
        "names": final_names,
        "skiprows": effective_header_line,
        "usecols": list(usecols) if usecols else None,
        "dtype": str,
        "comment": "#",
    }
    total_rows = None
    if progress_callback is not None:
        total_rows = _count_csv_data_rows(path, effective_header_line, progress_callback, progress_start, min(progress_start + 3, progress_end))
        read_progress_start = min(progress_start + 3, progress_end)
    else:
        read_progress_start = progress_start
    try:
        if progress_callback is not None:
            return _read_csv_with_chunk_progress(
                path,
                total_rows,
                progress_callback=progress_callback,
                progress_start=read_progress_start,
                progress_end=progress_end,
                sep=csv_delimiter,
                low_memory=True,
                memory_map=True,
                **read_kwargs,
            )
        return pd.read_csv(path, sep=csv_delimiter, low_memory=True, memory_map=True, **read_kwargs)
    except (MemoryError, pd.errors.ParserError) as exc:
        if "out of memory" not in str(exc).lower() and not isinstance(exc, MemoryError):
            raise
        if progress_callback is not None:
            return _read_csv_with_chunk_progress(
                path,
                total_rows,
                progress_callback=progress_callback,
                progress_start=read_progress_start,
                progress_end=progress_end,
                sep=csv_delimiter,
                engine="python",
                **read_kwargs,
            )
        return pd.read_csv(path, sep=csv_delimiter, engine="python", **read_kwargs)


def _normalize_export_field_type(field_type: object) -> str | None:
    text = str(field_type or "").strip().lower()
    if not text:
        return None
    aliases = {
        "bool": "boolean",
        "boolean": "boolean",
        "logical": "boolean",
        "int": "int",
        "integer": "int",
        "long": "int",
        "double": "double",
        "float": "double",
        "float64": "double",
        "number": "double",
        "string": "string",
        "str": "string",
        "text": "string",
        "name": "string",
        "categorical": "string",
        "category": "string",
        "namedshort": "string",
    }
    normalized = aliases.get(text)
    if normalized is None:
        raise ValueError(
            f"Unsupported export field type {field_type!r}. Supported types: boolean, int, double, string."
        )
    return normalized


def _normalize_export_type_overrides(column_types: Mapping[str, object] | None) -> Dict[str, str]:
    normalized: Dict[str, str] = {}
    for raw_name, raw_type in (column_types or {}).items():
        column_name = str(raw_name or "").strip()
        if not column_name:
            continue
        normalized[column_name] = _normalize_export_field_type(raw_type) or ""
    return {key: value for key, value in normalized.items() if value}


def _coerce_bool(value, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", ""}:
        return False
    return default


def _normalize_value_exceptions(value_exceptions: Mapping[str, Mapping[str, object]] | None) -> Dict[str, Dict[str, Dict[str, object]]]:
    normalized: Dict[str, Dict[str, Dict[str, object]]] = {}
    for raw_column, raw_rules in (value_exceptions or {}).items():
        column_name = str(raw_column or "").strip()
        if not column_name or not isinstance(raw_rules, Mapping):
            continue
        rules: Dict[str, Dict[str, object]] = {}
        for raw_value, raw_rule in raw_rules.items():
            value_text = str(raw_value)
            if isinstance(raw_rule, Mapping):
                replacement = raw_rule.get("replacement", "")
                include_in_regularization = _coerce_bool(
                    raw_rule.get("include_in_regularization", raw_rule.get("include_regularization", False))
                )
            else:
                replacement = raw_rule
                include_in_regularization = False
            rules[value_text] = {
                "replacement": "" if replacement is None else str(replacement),
                "include_in_regularization": bool(include_in_regularization),
            }
        if rules:
            normalized[column_name] = rules
    return normalized


def _filter_value_exceptions_for_regularization(
    value_exceptions: Mapping[str, Mapping[str, object]] | None,
    include_in_regularization: bool | None,
) -> Dict[str, Dict[str, Dict[str, object]]]:
    normalized = _normalize_value_exceptions(value_exceptions)
    if include_in_regularization is None:
        return normalized

    filtered: Dict[str, Dict[str, Dict[str, object]]] = {}
    for column_name, column_rules in normalized.items():
        rules = {
            bad_value: rule
            for bad_value, rule in column_rules.items()
            if bool(rule.get("include_in_regularization", False)) == bool(include_in_regularization)
        }
        if rules:
            filtered[column_name] = rules
    return filtered


def _apply_value_exceptions(
    df: pd.DataFrame,
    value_exceptions: Mapping[str, Mapping[str, object]] | None,
    include_in_regularization: bool | None = None,
) -> pd.DataFrame:
    normalized = _filter_value_exceptions_for_regularization(value_exceptions, include_in_regularization)
    if not normalized:
        return df

    out = df.copy()
    for column_name, replacements in normalized.items():
        if column_name not in out.columns or not replacements:
            continue
        series = out[column_name]
        text_series = series.astype(str)
        stripped_series = text_series.str.strip()
        updated = series.copy()
        for bad_value, rule in replacements.items():
            replacement = str(rule.get("replacement", ""))
            mask = series.notna() & ((text_series == bad_value) | (stripped_series == bad_value))
            if mask.any():
                updated.loc[mask] = replacement
        out[column_name] = updated
    return out


def _infer_axis_cell_size(values: np.ndarray) -> float:
    finite_values = values[np.isfinite(values)]
    uniq = np.unique(finite_values)
    if uniq.size < 2:
        return 1.0

    diffs = np.diff(np.sort(uniq))
    diffs = diffs[diffs > 1e-12]
    if diffs.size == 0:
        return 1.0

    rounded = np.round(diffs, decimals=6)
    rounded = rounded[rounded > 0]
    if rounded.size == 0:
        return float(np.min(diffs))

    values, counts = np.unique(rounded, return_counts=True)
    best_count = int(np.max(counts))
    if best_count <= 1:
        return float(np.median(diffs))
    best_values = values[counts == best_count]
    return float(np.min(best_values))


def _infer_cell_size(df: pd.DataFrame, xyz_cols: Sequence[str]) -> np.ndarray:
    sizes: List[float] = []
    for col in xyz_cols:
        vals = pd.to_numeric(df[col], errors="coerce").dropna().to_numpy(dtype=float)
        sizes.append(_infer_axis_cell_size(vals))
    return np.asarray(sizes, dtype=float)


def _next_power_of_two_multiple(value: float, quantum: float) -> float:
    value = float(value)
    quantum = float(quantum)
    if value <= 0 or quantum <= 0 or not math.isfinite(value) or not math.isfinite(quantum):
        return float("nan")
    ratio = max(value / quantum, 1.0)
    level = int(math.ceil(math.log(ratio, 2))) if ratio > 1.0 else 0
    return quantum * (2.0 ** level)


def _hierarchy_level_for_multiple(value: float, quantum: float, tolerance: float) -> int | None:
    value = float(value)
    quantum = float(quantum)
    if value <= 0 or quantum <= 0 or not math.isfinite(value) or not math.isfinite(quantum):
        return None
    ratio = value / quantum
    if ratio <= 0:
        return None
    level = int(round(math.log(ratio, 2)))
    expected = quantum * (2.0 ** level)
    if level < 0 or abs(value - expected) > max(float(tolerance), abs(expected) * 1e-6, 1e-9):
        return None
    return level


def _hierarchy_spacing_masks(values: np.ndarray, quantum: float, tolerance: float) -> tuple[np.ndarray, np.ndarray]:
    spacings = np.asarray(values, dtype=float)
    valid_mask = np.asarray([
        _hierarchy_level_for_multiple(float(value), quantum, tolerance) is not None
        for value in spacings
    ], dtype=bool)
    return spacings[valid_mask], spacings[~valid_mask]


def _spacing_histogram(values: np.ndarray, precision: int = 6) -> list[dict[str, object]]:
    if len(values) == 0:
        return []
    rounded = np.round(np.asarray(values, dtype=float), decimals=precision)
    rounded = rounded[np.isfinite(rounded) & (rounded > 0)]
    if len(rounded) == 0:
        return []
    unique_values, counts = np.unique(rounded, return_counts=True)
    order = np.argsort(unique_values)
    return [
        {"spacing": float(unique_values[index]), "count": int(counts[index])}
        for index in order
    ]


def _infer_base_size_from_hierarchy_spacing_frequencies(
    values: np.ndarray,
    quantum: float,
    tolerance: float,
) -> dict[str, object]:
    levels: List[int] = []
    for value in np.asarray(values, dtype=float):
        level = _hierarchy_level_for_multiple(float(value), quantum, tolerance)
        if level is not None:
            levels.append(int(level))
    if not levels:
        return {"status": "insufficient", "confidence": "none", "frequency_table": []}

    level_counts = {level: int(levels.count(level)) for level in sorted(set(levels))}
    max_level = max(level_counts)
    counts = [int(level_counts.get(level, 0)) for level in range(max_level + 1)]
    frequency_table = [
        {
            "level": int(level),
            "spacing": float(quantum * (2.0 ** level)),
            "count": int(counts[level]),
        }
        for level in range(max_level + 1)
    ]

    best_candidate: dict[str, object] | None = None
    for level in range(max_level):
        current_count = int(counts[level])
        next_count = int(counts[level + 1])
        if current_count < 3:
            continue
        if next_count <= 0:
            continue
        if next_count > current_count * 0.5 or current_count - next_count < 2:
            continue
        previous_counts = [count for count in counts[:level + 1] if count > 0]
        rising_steps = sum(
            1
            for left, right in zip(previous_counts, previous_counts[1:])
            if right >= left
        )
        rising_score = rising_steps / max(len(previous_counts) - 1, 1)
        drop_ratio = next_count / current_count
        score = (1.0 - drop_ratio) + rising_score + math.log1p(current_count) * 0.05
        if best_candidate is None or score > float(best_candidate["score"]):
            best_candidate = {
                "level": int(level),
                "size": float(quantum * (2.0 ** level)),
                "count": current_count,
                "next_level": int(level + 1),
                "next_size": float(quantum * (2.0 ** (level + 1))),
                "next_count": next_count,
                "drop_ratio": float(drop_ratio),
                "rising_score": float(rising_score),
                "score": float(score),
            }

    if best_candidate is None:
        return {
            "status": "ambiguous",
            "confidence": "none",
            "frequency_table": frequency_table,
        }

    confidence = "medium"
    if (
        float(best_candidate["drop_ratio"]) <= 0.33
        and float(best_candidate["rising_score"]) >= 0.75
        and int(best_candidate["count"]) >= 5
    ):
        confidence = "high"
    return {
        "status": "inferred",
        "confidence": confidence,
        "frequency_table": frequency_table,
        "candidate_level": int(best_candidate["level"]),
        "candidate_size": float(best_candidate["size"]),
        "drop_from_count": int(best_candidate["count"]),
        "drop_to_level": int(best_candidate["next_level"]),
        "drop_to_size": float(best_candidate["next_size"]),
        "drop_to_count": int(best_candidate["next_count"]),
        "drop_ratio": float(best_candidate["drop_ratio"]),
        "rising_score": float(best_candidate["rising_score"]),
    }


def _infer_axis_base_block_from_centroids(values: np.ndarray, index_tolerance: float) -> dict[str, object]:
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    unique_values = np.unique(finite_values)
    if len(unique_values) == 0:
        return {"status": "empty", "confidence": "none"}

    if len(unique_values) < 2:
        return {
            "status": "insufficient_axis_variation",
            "confidence": "none",
            "unique_centers": int(len(unique_values)),
            "min_center": float(unique_values[0]),
            "max_center": float(unique_values[0]),
        }

    diffs = np.diff(np.sort(unique_values))
    diffs = diffs[np.isfinite(diffs) & (diffs > max(float(index_tolerance), 1e-9))]
    if len(diffs) == 0:
        return {
            "status": "insufficient_axis_variation",
            "confidence": "none",
            "unique_centers": int(len(unique_values)),
            "min_center": float(unique_values.min()),
            "max_center": float(unique_values.max()),
        }

    rounded_diffs = np.round(diffs, decimals=6)
    quantum = float(np.min(rounded_diffs[rounded_diffs > 0]))
    hierarchy_diffs, discarded_diffs = _hierarchy_spacing_masks(rounded_diffs, quantum, index_tolerance)
    frequency_break = _infer_base_size_from_hierarchy_spacing_frequencies(
        hierarchy_diffs,
        quantum,
        index_tolerance,
    )
    observed_span = float(unique_values.max() - unique_values.min() + quantum)
    enclosing_power_of_two_size = float(_next_power_of_two_multiple(observed_span, quantum))
    base_size = float(frequency_break.get("candidate_size", quantum)) if frequency_break.get("status") == "inferred" else quantum
    origin = float(unique_values.min() - 0.5 * quantum)
    base_level = _hierarchy_level_for_multiple(base_size, quantum, index_tolerance)

    confidence = "low"
    if frequency_break.get("status") == "inferred":
        confidence = str(frequency_break.get("confidence", "medium"))
    elif len(unique_values) >= 4 and len(discarded_diffs) == 0:
        confidence = "medium"

    return {
        "status": "ok",
        "confidence": confidence,
        "unique_centers": int(len(unique_values)),
        "quantum": float(quantum),
        "observed_span": float(observed_span),
        "base_size": float(base_size),
        "base_level": None if base_level is None else int(base_level),
        "parent_base_size_ambiguous": frequency_break.get("status") != "inferred",
        "base_size_inference": frequency_break,
        "enclosing_power_of_two_size": float(enclosing_power_of_two_size),
        "origin": float(origin),
        "min_center": float(unique_values.min()),
        "max_center": float(unique_values.max()),
        "spacing_histogram": _spacing_histogram(diffs),
        "hierarchy_spacing_histogram": _spacing_histogram(hierarchy_diffs),
        "discarded_spacing_histogram": _spacing_histogram(discarded_diffs),
    }


def infer_base_block_size_from_centroid_statistics(
    df: pd.DataFrame,
    xyz_cols: Sequence[str],
    index_tolerance: float = 1e-3,
) -> dict[str, object]:
    clean = _to_numeric_xyz(df, xyz_cols)
    if clean.empty:
        return {"status": "empty", "confidence": "none", "axes": {}}

    axis_names = ("x", "y", "z")
    axes: dict[str, object] = {}
    cell_size = np.full(3, np.nan, dtype=float)
    origin = np.full(3, np.nan, dtype=float)

    for axis, (axis_name, coord_col) in enumerate(zip(axis_names, xyz_cols)):
        axis_report = _infer_axis_base_block_from_centroids(
            pd.to_numeric(clean[coord_col], errors="coerce").to_numpy(dtype=float),
            index_tolerance,
        )
        axes[axis_name] = axis_report
        if axis_report.get("status") == "ok":
            cell_size[axis] = float(axis_report["base_size"])
            origin[axis] = float(axis_report["origin"])

    good_axes = np.isfinite(cell_size) & (cell_size > 0)
    if not np.any(good_axes):
        return {"status": "insufficient_axis_variation", "confidence": "none", "axes": axes}

    peer_size = float(np.median(cell_size[good_axes]))
    for axis, axis_name in enumerate(axis_names):
        if good_axes[axis]:
            continue
        values = pd.to_numeric(clean[xyz_cols[axis]], errors="coerce").dropna().to_numpy(dtype=float)
        if len(values) == 0:
            continue
        cell_size[axis] = peer_size
        origin[axis] = float(np.min(values) - 0.5 * peer_size)
        axis_report = dict(axes.get(axis_name, {}))
        axis_report.update({
            "status": "filled_from_peer_axes",
            "confidence": "low",
            "base_size": float(peer_size),
            "origin": float(origin[axis]),
        })
        axes[axis_name] = axis_report

    complete = bool(np.all(np.isfinite(cell_size) & (cell_size > 0) & np.isfinite(origin)))
    axis_confidences = [str(report.get("confidence", "none")) for report in axes.values() if isinstance(report, Mapping)]
    if not complete:
        confidence = "none"
        status = "partial"
    elif any(value == "low" for value in axis_confidences):
        confidence = "low"
        status = "ok"
    elif any(value == "medium" for value in axis_confidences):
        confidence = "medium"
        status = "ok"
    else:
        confidence = "high"
        status = "ok"

    return {
        "status": status,
        "confidence": confidence,
        "cell_size": [float(value) for value in cell_size] if complete else None,
        "origin": [float(value) for value in origin] if complete else None,
        "source": "centroid_spacing_statistics",
        "axes": axes,
    }


def _format_numeric_vector(values: Sequence[float] | np.ndarray, precision: int = 6) -> str:
    formatted = []
    for value in np.asarray(values, dtype=float):
        if not math.isfinite(float(value)):
            formatted.append(str(float(value)))
            continue
        text = f"{float(value):.{precision}f}".rstrip("0").rstrip(".")
        formatted.append(text or "0")
    return "[" + ", ".join(formatted) + "]"


def _summarize_statistical_base_size_inference(statistical_base_geometry: Mapping[str, object]) -> dict[str, object]:
    axes = statistical_base_geometry.get("axes") if isinstance(statistical_base_geometry, Mapping) else None
    if not isinstance(axes, Mapping):
        return {"source": "none", "message": "No centroid-spacing block-size statistics were available."}

    inferred_axes = []
    ambiguous_axes = []
    axis_details = {}
    for axis_name, raw_report in axes.items():
        if not isinstance(raw_report, Mapping):
            continue
        inference = raw_report.get("base_size_inference")
        inference_status = inference.get("status") if isinstance(inference, Mapping) else None
        detail = {
            "status": raw_report.get("status"),
            "quantum": raw_report.get("quantum"),
            "base_size": raw_report.get("base_size"),
            "base_level": raw_report.get("base_level"),
            "parent_base_size_ambiguous": bool(raw_report.get("parent_base_size_ambiguous", False)),
            "inference_status": inference_status,
        }
        if isinstance(inference, Mapping):
            detail.update({
                "candidate_size": inference.get("candidate_size"),
                "drop_to_size": inference.get("drop_to_size"),
                "drop_from_count": inference.get("drop_from_count"),
                "drop_to_count": inference.get("drop_to_count"),
                "drop_ratio": inference.get("drop_ratio"),
            })
        axis_details[str(axis_name)] = detail
        if inference_status == "inferred":
            inferred_axes.append(str(axis_name).upper())
        if raw_report.get("parent_base_size_ambiguous", False):
            ambiguous_axes.append(str(axis_name).upper())

    if inferred_axes:
        source = "centroid_spacing_frequency_break"
        message = (
            "Base block size inferred from hierarchy-valid centroid spacing frequencies; "
            f"frequency break detected on axis/axes {', '.join(inferred_axes)}."
        )
    elif ambiguous_axes:
        source = "centroid_spacing_quantum_fallback"
        message = (
            "Base block size remains ambiguous from centroid spacing statistics; "
            "using the smallest hierarchy quantum as a conservative fallback."
        )
    else:
        source = "centroid_spacing_statistics"
        message = "Base block size inferred from centroid spacing statistics."
    return {
        "source": source,
        "message": message,
        "axis_details": axis_details,
        "inferred_axes": inferred_axes,
        "ambiguous_axes": ambiguous_axes,
    }


def _build_block_size_determination_report(
    prepared: Mapping[str, object],
    metadata_geometry: Mapping[str, object],
    statistical_base_geometry: Mapping[str, object],
    user_cell_size_provided: bool,
    user_origin_provided: bool,
    use_grid_index_columns: bool,
    regularize_to_base_block: bool,
) -> dict[str, object]:
    cell_size = [float(x) for x in np.asarray(prepared["cell"], dtype=float)]
    origin = [float(x) for x in np.asarray(prepared["origin"], dtype=float)]
    index_source = str(prepared.get("grid_index_source") or "xyz")
    width_source = prepared.get("irregular_width_source")
    statistical_source = prepared.get("statistical_base_block_source")

    if width_source in {"lower_upper_columns", "mapped_lower_upper_columns"}:
        source = "csv_lower_upper_extent_columns"
        method = "explicit_row_extents"
        message = "Block row sizes read directly from CSV lower/upper extent columns."
        if width_source == "mapped_lower_upper_columns":
            source = "mapped_csv_lower_upper_extent_columns"
            message = "Block row sizes read directly from user-mapped CSV lower/upper extent columns."
    elif width_source == "dimension_columns":
        source = "csv_dimension_columns"
        method = "explicit_dimension_columns"
        message = "Block row sizes read directly from CSV dimension columns and converted to lower/upper extents."
    elif user_cell_size_provided:
        source = "user_cell_size"
        method = "manual_parent_size"
        message = "Parent/base block size supplied by the user; centroid hierarchy inference used only for missing row extents."
    elif metadata_geometry.get("cell_size"):
        source = str(metadata_geometry.get("cell_size_source") or "csv_metadata")
        method = "metadata_parent_size"
        message = f"Parent/base block size read from CSV metadata ({source})."
    elif use_grid_index_columns:
        source = "grid_index_columns"
        method = "grid_index_geometry"
        message = "Block geometry derived from CSV grid_i/grid_j/grid_k columns."
    elif statistical_source:
        statistical_report = _summarize_statistical_base_size_inference(statistical_base_geometry)
        source = str(statistical_report.get("source"))
        method = "statistical_centroid_spacing"
        message = str(statistical_report.get("message"))
    elif regularize_to_base_block:
        source = "regularized_xyz_grid"
        method = "regularized_grid"
        message = "Rows regularized to a base-cell grid before BMF export."
    elif not prepared.get("is_irregular"):
        source = "xyz_regular_grid"
        method = "regular_grid_from_xyz"
        message = "Regular grid block size inferred from XYZ coordinate spacing."
    else:
        source = index_source
        method = "centroid_hierarchy_inference"
        message = "Row extents inferred from centroid hierarchy using the selected parent/base block size."

    report = {
        "method": method,
        "source": source,
        "message": message,
        "cell_size": cell_size,
        "origin": origin,
        "index_source": index_source,
        "irregular_width_source": width_source,
        "user_cell_size_provided": bool(user_cell_size_provided),
        "user_origin_provided": bool(user_origin_provided),
        "metadata_cell_size_source": metadata_geometry.get("cell_size_source"),
        "metadata_origin_source": metadata_geometry.get("origin_source"),
    }
    geometry_mapping = prepared.get("geometry_column_mapping")
    if geometry_mapping:
        report["geometry_column_mapping"] = geometry_mapping
    if statistical_source:
        report["statistical_base_block_source"] = statistical_source
        report["statistical_inference"] = _summarize_statistical_base_size_inference(statistical_base_geometry)
    return report


def _build_grid_alignment_error_message(
    clean: pd.DataFrame,
    xyz_cols: Sequence[str],
    cell: np.ndarray,
    raw_idx: np.ndarray,
    rounded: np.ndarray,
    max_err: float,
    index_tolerance: float,
    explicit_cell_size: bool,
) -> str:
    axis_errors = np.max(np.abs(raw_idx - rounded), axis=0)
    inferred_cell = _infer_cell_size(clean, xyz_cols)
    message = (
        "Coordinates do not align to a regular grid within tolerance. "
        f"max_err={max_err:.6g}, tolerance={index_tolerance:.6g}, "
        f"axis_errors={dict(zip(xyz_cols, _format_numeric_vector(axis_errors).strip('[]').split(', ')))}."
    )
    if explicit_cell_size:
        ratios = np.divide(
            cell,
            inferred_cell,
            out=np.full_like(cell, np.nan, dtype=float),
            where=inferred_cell > 0,
        )
        message += (
            f" Explicit cell_size={_format_numeric_vector(cell)}; smallest CSV coordinate/centroid increment is "
            f"{_format_numeric_vector(inferred_cell)}; ratio={_format_numeric_vector(ratios, precision=3)}."
        )
        coarse_axes = [
            str(axis_name)
            for axis_name, ratio, axis_error in zip(xyz_cols, ratios, axis_errors)
            if math.isfinite(float(ratio)) and ratio > 1.001 and axis_error > index_tolerance
        ]
        if coarse_axes:
            message += (
                " The selected cell size is coarser than the actual coordinate spacing on "
                f"axis/axes {coarse_axes}, so the CSV appears to contain sub-block rows rather than one row per "
                "requested BMF cell. For a dense tbms-config-text/Vulcan BMF, aggregate or regularize the CSV to "
                "one row per exported cell first. For row-indexed output that preserves the CSV rows without dense "
                "regular-grid allocation, use the tbms-experimental backend."
            )
    else:
        message += (
            f" Inferred cell_size={_format_numeric_vector(cell)}. Set explicit Cell Size X/Y/Z if these coordinates "
            "represent a known regular grid."
        )
    return message

def _build_centroid_hierarchy_error_message(
    error: Exception,
    df: pd.DataFrame,
    xyz_cols: Sequence[str],
    parent_size: Sequence[float] | np.ndarray,
) -> str:
    message = str(error)
    try:
        clean = _to_numeric_xyz(df, xyz_cols)
        inferred_cell = _infer_cell_size(clean, xyz_cols)
        parent_arr = np.asarray(parent_size, dtype=float)
        ratios = np.divide(
            parent_arr,
            inferred_cell,
            out=np.full_like(parent_arr, np.nan, dtype=float),
            where=inferred_cell > 0,
        )
        message += (
            f" Explicit parent cell_size={_format_numeric_vector(parent_arr)}; smallest CSV coordinate/centroid increment is "
            f"{_format_numeric_vector(inferred_cell)}; ratio={_format_numeric_vector(ratios, precision=3)}."
        )
        coarse_axes = [
            str(axis_name)
            for axis_name, ratio in zip(xyz_cols, ratios)
            if math.isfinite(float(ratio)) and ratio > 1.001
        ]
        if coarse_axes:
            message += (
                " The selected parent block size is coarser than the actual coordinate spacing on "
                f"axis/axes {coarse_axes}, so the CSV appears to contain sub-block rows rather than one row per "
                "parent block. The row extents must form a parent_size, parent_size/2, parent_size/4, ... hierarchy; "
                "otherwise aggregate or regularize the CSV to one row per exported cell first."
            )
    except Exception:
        pass
    return message


def _to_numeric_xyz(df: pd.DataFrame, xyz_cols: Sequence[str]) -> pd.DataFrame:
    out = df.copy()
    for col in xyz_cols:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=list(xyz_cols))
    return out


def _prepare_grid(
    df: pd.DataFrame,
    xyz_cols: Sequence[str],
    cell_size: Sequence[float] | None,
    origin: Sequence[float] | None,
    index_tolerance: float,
) -> Dict[str, object]:
    clean = _to_numeric_xyz(df, xyz_cols)
    if clean.empty:
        raise ValueError("No valid XYZ rows found after numeric conversion.")

    explicit_cell_size = cell_size is not None
    if cell_size is None:
        cell = _infer_cell_size(clean, xyz_cols)
    else:
        cell = np.asarray(cell_size, dtype=float)
    if np.any(cell <= 0):
        raise ValueError(f"Invalid cell size values: {cell}")

    xyz = clean[list(xyz_cols)].to_numpy(dtype=float)

    if origin is None:
        min_xyz = xyz.min(axis=0)
        origin_arr = min_xyz - 0.5 * cell
    else:
        origin_arr = np.asarray(origin, dtype=float)

    center_reference = origin_arr + 0.5 * cell
    raw_idx = (xyz - center_reference) / cell
    rounded = np.rint(raw_idx)
    max_err = float(np.max(np.abs(raw_idx - rounded)))
    if max_err > index_tolerance:
        raise ValueError(
            _build_grid_alignment_error_message(
                clean,
                xyz_cols,
                cell,
                raw_idx,
                rounded,
                max_err,
                index_tolerance,
                explicit_cell_size,
            )
        )

    idx = np.floor(raw_idx + 1e-9).astype(int)
    if np.any(idx < 0):
        shift = idx.min(axis=0)
        idx = idx - shift
        origin_arr = origin_arr + shift * cell
        center_reference = origin_arr + 0.5 * cell

    dims = idx.max(axis=0) + 1
    extents = dims * cell

    grid_tuples = [tuple(map(int, r)) for r in idx]
    dup_count = int(len(grid_tuples) - len(set(grid_tuples)))

    return {
        "df": clean,
        "idx": idx,
        "grid_tuples": grid_tuples,
        "origin": origin_arr,
        "cell": cell,
        "dims": dims.astype(int),
        "extents": extents,
        "duplicates": dup_count,
        "max_index_error": max_err,
    }


def _infer_cell_size_from_grid_indices(clean: pd.DataFrame, xyz_cols: Sequence[str], idx: np.ndarray) -> np.ndarray:
    cell = np.full(3, np.nan, dtype=float)
    for axis, coord_col in enumerate(xyz_cols):
        pairs = pd.DataFrame({
            "grid_index": idx[:, axis],
            "coord": clean[coord_col].to_numpy(dtype=float),
        })
        grouped = pairs.groupby("grid_index", sort=True)["coord"].median()
        if len(grouped) >= 2:
            index_delta = np.diff(grouped.index.to_numpy(dtype=float))
            coord_delta = np.diff(grouped.to_numpy(dtype=float))
            valid = index_delta != 0
            ratios = coord_delta[valid] / index_delta[valid]
            ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
            if len(ratios):
                cell[axis] = float(np.median(ratios))
    missing = ~np.isfinite(cell) | (cell <= 0)
    if np.any(missing):
        inferred = _infer_cell_size(clean, xyz_cols)
        cell[missing] = inferred[missing]
    if np.any(~np.isfinite(cell)) or np.any(cell <= 0):
        raise ValueError(
            "Could not infer valid BMF cell size from grid_i/grid_j/grid_k and XYZ columns. "
            "Set explicit Cell Size X/Y/Z values."
        )
    return cell


def _prepare_grid_from_index_columns(
    df: pd.DataFrame,
    xyz_cols: Sequence[str],
    cell_size: Sequence[float] | None,
    origin: Sequence[float] | None,
) -> Dict[str, object]:
    required_cols = [*xyz_cols, *GRID_INDEX_COLUMNS]
    clean = df.copy()
    for col in required_cols:
        clean[col] = pd.to_numeric(clean[col], errors="coerce")
    clean = clean.dropna(subset=required_cols)
    if clean.empty:
        raise ValueError("No valid XYZ/grid_i/grid_j/grid_k rows found after numeric conversion.")

    raw_idx = clean[list(GRID_INDEX_COLUMNS)].to_numpy(dtype=float)
    rounded = np.rint(raw_idx)
    max_index_error = float(np.max(np.abs(raw_idx - rounded)))
    if max_index_error > 1e-6:
        raise ValueError(
            "CSV grid_i/grid_j/grid_k columns must contain integer grid indices. "
            f"Maximum index error from an integer is {max_index_error:g}."
        )
    idx = rounded.astype(np.int64)

    if cell_size is None:
        cell = _infer_cell_size_from_grid_indices(clean, xyz_cols, idx)
    else:
        cell = np.asarray(cell_size, dtype=float)
    if cell.shape != (3,) or np.any(cell <= 0):
        raise ValueError(f"Invalid cell size values: {cell_size}")

    if origin is None:
        xyz = clean[list(xyz_cols)].to_numpy(dtype=float)
        origin_candidates = xyz - (idx.astype(float) + 0.5) * cell
        origin_arr = np.median(origin_candidates, axis=0)
    else:
        origin_arr = np.asarray(origin, dtype=float)
    if origin_arr.shape != (3,):
        raise ValueError(f"Invalid origin values: {origin}")

    if np.any(idx < 0):
        shift = idx.min(axis=0)
        idx = idx - shift
        origin_arr = origin_arr + shift * cell

    dims = idx.max(axis=0) + 1
    extents = dims * cell
    grid_tuples = [tuple(map(int, row)) for row in idx]
    dup_count = int(len(grid_tuples) - len(set(grid_tuples)))

    return {
        "df": clean,
        "idx": idx,
        "grid_tuples": grid_tuples,
        "origin": origin_arr,
        "cell": cell,
        "dims": dims.astype(int),
        "extents": extents,
        "duplicates": dup_count,
        "max_index_error": max_index_error,
        "grid_index_source": "grid_i/grid_j/grid_k",
    }


def _cluster_axis_centers(values: np.ndarray, tolerance: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return np.empty(0, dtype=float)
    sorted_values = np.sort(values)
    groups: List[List[float]] = [[float(sorted_values[0])]]
    for value in sorted_values[1:]:
        current = float(value)
        if abs(current - groups[-1][-1]) <= tolerance:
            groups[-1].append(current)
        else:
            groups.append([current])
    return np.asarray([float(np.mean(group)) for group in groups], dtype=float)


def _infer_irregular_axis_widths(local_values: np.ndarray, parent_extent: float, tolerance: float) -> np.ndarray:
    parent_extent = float(parent_extent)
    local_values = np.asarray(local_values, dtype=float)
    if len(local_values) == 0:
        return np.empty(0, dtype=float)
    if parent_extent <= 0 or not math.isfinite(parent_extent):
        raise ValueError("Parent block extent must be a finite positive value for irregular width inference.")

    hierarchy_widths: List[float] = []
    width = parent_extent
    min_width = max(float(tolerance) * 4.0, abs(parent_extent) * 1e-9, 1e-12)
    for _ in range(24):
        hierarchy_widths.append(float(width))
        width *= 0.5
        if width < min_width:
            break

    inferred = np.empty(len(local_values), dtype=float)
    for row_index, center_value in enumerate(local_values):
        center = float(center_value)
        selected_width = None
        for candidate_width in hierarchy_widths:
            axis_tolerance = max(float(tolerance), abs(candidate_width) * 1e-6, 1e-9)
            lower = center - 0.5 * candidate_width
            upper = center + 0.5 * candidate_width
            if lower < -axis_tolerance or upper > parent_extent + axis_tolerance:
                continue
            lattice_position = (center / candidate_width) - 0.5
            if abs(lattice_position - round(lattice_position)) <= axis_tolerance / max(candidate_width, 1e-12):
                selected_width = candidate_width
                break
        if selected_width is None:
            raise ValueError(
                "Could not infer a hierarchical block width from centroid-only coordinates. "
                "Valid inferred widths must be parent_size, parent_size/2, parent_size/4, ... along each axis."
            )
        inferred[row_index] = selected_width
    return inferred


def _hierarchy_level(value: float, quantum: float, tolerance: float) -> int | None:
    return _hierarchy_level_for_multiple(value, quantum, tolerance)


def _format_hierarchy_size(value: float) -> str:
    numeric = float(value)
    if not math.isfinite(numeric):
        return str(numeric)
    text = f"{numeric:.6f}".rstrip("0").rstrip(".")
    return text or "0"


def _build_irregular_hierarchy_report(
    prepared: Mapping[str, object],
    xyz_cols: Sequence[str],
    index_tolerance: float,
) -> Dict[str, object] | None:
    if not prepared.get("is_irregular") or prepared.get("irregular_widths") is None:
        return None

    df = prepared.get("df")
    if not isinstance(df, pd.DataFrame):
        return None

    widths = np.asarray(prepared["irregular_widths"], dtype=float)
    parent_size = np.asarray(prepared["cell"], dtype=float)
    coords = df[list(xyz_cols)].to_numpy(dtype=float)
    if widths.ndim != 2 or widths.shape[1] != 3 or parent_size.shape != (3,):
        return None

    axis_reports: Dict[str, object] = {}
    for axis, axis_name in enumerate(("x", "y", "z")):
        axis_widths = widths[:, axis]
        finite_widths = axis_widths[np.isfinite(axis_widths) & (axis_widths > 0)]
        if len(finite_widths) == 0:
            continue

        quantum = float(np.min(finite_widths))
        base = float(parent_size[axis])
        base_level = _hierarchy_level(base, quantum, index_tolerance)
        if base_level is None:
            base_level = max(_hierarchy_level(value, quantum, index_tolerance) or 0 for value in finite_widths)

        row_counts: Counter[int] = Counter()
        for value in finite_widths:
            level = _hierarchy_level(float(value), quantum, index_tolerance)
            if level is not None:
                row_counts[level] += 1

        spacing_counts: Counter[int] = Counter()
        unique_coords = np.unique(coords[np.isfinite(coords[:, axis]), axis])
        diffs = np.diff(np.sort(unique_coords)) if len(unique_coords) >= 2 else np.empty(0, dtype=float)
        diffs = diffs[np.isfinite(diffs) & (diffs > max(float(index_tolerance), 1e-9))]
        for value in diffs:
            level = _hierarchy_level(float(value), quantum, index_tolerance)
            if level is not None:
                spacing_counts[level] += 1

        max_level = max([int(base_level), *(row_counts.keys() or [0]), *(spacing_counts.keys() or [0])])
        rows = []
        max_row_count = max(row_counts.values(), default=0)
        max_spacing_count = max(spacing_counts.values(), default=0)
        for level in range(max_level + 1):
            size = quantum * (2.0 ** level)
            row_count = int(row_counts.get(level, 0))
            spacing_count = int(spacing_counts.get(level, 0))
            row_bar = "#" * int(round((row_count / max_row_count) * 24)) if max_row_count else ""
            spacing_bar = "#" * int(round((spacing_count / max_spacing_count) * 24)) if max_spacing_count else ""
            marker = "base" if level == int(base_level) else ("sub" if level < int(base_level) else "gap")
            rows.append({
                "level": int(level),
                "size": float(size),
                "multiple": int(2 ** level),
                "role": marker,
                "row_count": row_count,
                "spacing_count": spacing_count,
                "row_histogram": row_bar,
                "spacing_histogram": spacing_bar,
            })

        axis_reports[axis_name] = {
            "quantum": float(quantum),
            "base": float(base),
            "base_level": int(base_level),
            "levels": rows,
        }

    return {
        "source": prepared.get("irregular_width_source"),
        "axes": axis_reports,
    }


def _format_irregular_hierarchy_report(report: Mapping[str, object] | None) -> str:
    if not report:
        return ""
    axes = report.get("axes")
    if not isinstance(axes, Mapping) or not axes:
        return ""
    lines = ["Inferred irregular block hierarchy:"]
    for axis_name in ("x", "y", "z"):
        axis_report = axes.get(axis_name)
        if not isinstance(axis_report, Mapping):
            continue
        quantum = _format_hierarchy_size(float(axis_report.get("quantum", float("nan"))))
        base = _format_hierarchy_size(float(axis_report.get("base", float("nan"))))
        lines.append(f"  Axis {axis_name.upper()}: quantum={quantum}, base={base}")
        levels = axis_report.get("levels", [])
        if not isinstance(levels, list):
            continue
        for row in levels:
            if not isinstance(row, Mapping):
                continue
            size = _format_hierarchy_size(float(row.get("size", 0.0)))
            role = str(row.get("role", ""))
            row_count = int(row.get("row_count", 0))
            spacing_count = int(row.get("spacing_count", 0))
            row_bar = str(row.get("row_histogram", ""))
            spacing_bar = str(row.get("spacing_histogram", ""))
            lines.append(
                f"    {size:>10} ({role:>4}) rows={row_count:>8} {row_bar} | spacings={spacing_count:>8} {spacing_bar}"
            )
    return "\n".join(lines)


def _prepare_irregular_blocks(
    df: pd.DataFrame,
    xyz_cols: Sequence[str],
    metadata: Mapping[str, object],
    cell_size: Sequence[float] | None,
    origin: Sequence[float] | None,
    index_tolerance: float,
    geometry_size_cols: Mapping[str, str] | None = None,
    geometry_extent_cols: Mapping[str, str] | None = None,
) -> Dict[str, object]:
    clean = _to_numeric_xyz(df, xyz_cols)
    if clean.empty:
        raise ValueError("No valid XYZ rows found after numeric conversion.")

    extent_col_map = dict(geometry_extent_cols or {})
    if not extent_col_map and all(column in clean.columns for column in GEOMETRY_EXTENT_COLUMNS):
        extent_col_map = dict(zip(GEOMETRY_EXTENT_COLUMNS, GEOMETRY_EXTENT_COLUMNS))
    has_extent_columns = all(extent_col_map.get(role) in clean.columns for role in GEOMETRY_EXTENT_COLUMNS)
    if has_extent_columns:
        extent_source_cols = [extent_col_map[role] for role in GEOMETRY_EXTENT_COLUMNS]
        extent_values = clean[extent_source_cols].apply(pd.to_numeric, errors="coerce")
        valid_extents = extent_values.notna().all(axis=1)
        clean = clean.loc[valid_extents].reset_index(drop=True)
        extent_values = extent_values.loc[valid_extents].reset_index(drop=True)
        if clean.empty:
            raise ValueError("No valid irregular lower/upper extent rows found in CSV.")
        lower = extent_values[[extent_col_map["__lower_x"], extent_col_map["__lower_y"], extent_col_map["__lower_z"]]].to_numpy(dtype=float)
        upper = extent_values[[extent_col_map["__upper_x"], extent_col_map["__upper_y"], extent_col_map["__upper_z"]]].to_numpy(dtype=float)
        widths = upper - lower
        xyz = (lower + upper) * 0.5
        for axis, col in enumerate(xyz_cols):
            clean[col] = xyz[:, axis]
        if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper)) or np.any(widths <= 0):
            raise ValueError("Invalid irregular lower/upper extents in CSV.")

        parent_size_source = metadata.get("parent_block_size") or cell_size
        if parent_size_source is not None:
            parent_size = np.asarray(parent_size_source, dtype=float)
        else:
            parent_size = np.nanmax(widths, axis=0)
        if parent_size.shape != (3,) or np.any(parent_size <= 0) or np.any(~np.isfinite(parent_size)):
            raise ValueError(f"Invalid parent block size for irregular export: {parent_size_source}")
        if origin is not None:
            origin_arr = np.asarray(origin, dtype=float)
        elif metadata.get("minimum_corner") is not None:
            origin_arr = np.asarray(metadata["minimum_corner"], dtype=float)
        else:
            origin_arr = np.nanmin(lower, axis=0)
        if origin_arr.shape != (3,):
            raise ValueError(f"Invalid origin values for irregular export: {origin_arr}")
        parent_idx = np.floor((xyz - origin_arr) / parent_size + float(index_tolerance)).astype(np.int64)
        if metadata.get("size_in_parent_blocks") is not None:
            dims = np.asarray(metadata["size_in_parent_blocks"], dtype=np.int64)
        elif len(parent_idx):
            dims = np.maximum(parent_idx.max(axis=0) + 1, 1)
        else:
            dims = np.ones(3, dtype=np.int64)
        extents = np.maximum(np.nanmax(upper, axis=0) - origin_arr, dims * parent_size)
        canonical_extent_source = all(extent_col_map[role] == role for role in GEOMETRY_EXTENT_COLUMNS)
        clean = clean.drop(columns=[col for col in set(extent_source_cols) if col in clean.columns and col not in xyz_cols], errors="ignore")
        return {
            "df": clean.reset_index(drop=True),
            "idx": parent_idx,
            "grid_tuples": [tuple(map(int, row)) for row in parent_idx],
            "origin": origin_arr,
            "cell": parent_size,
            "dims": dims.astype(int),
            "extents": extents,
            "duplicates": 0,
            "max_index_error": 0.0,
            "is_irregular": True,
            "irregular_lower": lower,
            "irregular_upper": upper,
            "irregular_widths": widths,
            "irregular_width_source": "lower_upper_columns" if canonical_extent_source else "mapped_lower_upper_columns",
            "geometry_column_mapping": {"extent_cols": dict(extent_col_map)},
            "grid_index_source": "irregular_lower_upper_columns" if canonical_extent_source else "irregular_mapped_lower_upper_columns",
        }

    size_col_map = dict(geometry_size_cols or {})
    has_size_columns = all(size_col_map.get(role) in clean.columns for role in GEOMETRY_SIZE_COLUMN_ROLES)
    if has_size_columns:
        size_source_cols = [size_col_map[role] for role in GEOMETRY_SIZE_COLUMN_ROLES]
        size_values = clean[size_source_cols].apply(pd.to_numeric, errors="coerce")
        valid_sizes = size_values.notna().all(axis=1)
        clean = clean.loc[valid_sizes].reset_index(drop=True)
        size_values = size_values.loc[valid_sizes].reset_index(drop=True)
        if clean.empty:
            raise ValueError("No valid irregular dimension rows found in CSV.")
        widths = size_values[size_source_cols].to_numpy(dtype=float)
        xyz = clean[list(xyz_cols)].to_numpy(dtype=float)
        lower = xyz - 0.5 * widths
        upper = xyz + 0.5 * widths
        if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper)) or np.any(widths <= 0):
            raise ValueError("Invalid irregular dimension columns in CSV.")

        parent_size_source = metadata.get("parent_block_size") or cell_size
        if parent_size_source is not None:
            parent_size = np.asarray(parent_size_source, dtype=float)
        else:
            parent_size = np.nanmax(widths, axis=0)
        if parent_size.shape != (3,) or np.any(parent_size <= 0) or np.any(~np.isfinite(parent_size)):
            raise ValueError(f"Invalid parent block size for irregular export: {parent_size_source}")
        if origin is not None:
            origin_arr = np.asarray(origin, dtype=float)
        elif metadata.get("minimum_corner") is not None:
            origin_arr = np.asarray(metadata["minimum_corner"], dtype=float)
        else:
            origin_arr = np.nanmin(lower, axis=0)
        if origin_arr.shape != (3,):
            raise ValueError(f"Invalid origin values for irregular export: {origin_arr}")
        parent_idx = np.floor((xyz - origin_arr) / parent_size + float(index_tolerance)).astype(np.int64)
        if metadata.get("size_in_parent_blocks") is not None:
            dims = np.asarray(metadata["size_in_parent_blocks"], dtype=np.int64)
        elif len(parent_idx):
            dims = np.maximum(parent_idx.max(axis=0) + 1, 1)
        else:
            dims = np.ones(3, dtype=np.int64)
        extents = np.maximum(np.nanmax(upper, axis=0) - origin_arr, dims * parent_size)
        clean = clean.drop(columns=[col for col in set(size_source_cols) if col in clean.columns and col not in xyz_cols], errors="ignore")
        return {
            "df": clean.reset_index(drop=True),
            "idx": parent_idx,
            "grid_tuples": [tuple(map(int, row)) for row in parent_idx],
            "origin": origin_arr,
            "cell": parent_size,
            "dims": dims.astype(int),
            "extents": extents,
            "duplicates": 0,
            "max_index_error": 0.0,
            "is_irregular": True,
            "irregular_lower": lower,
            "irregular_upper": upper,
            "irregular_widths": widths,
            "irregular_width_source": "dimension_columns",
            "geometry_column_mapping": {"size_cols": dict(size_col_map)},
            "grid_index_source": "irregular_dimension_columns",
        }

    parent_size_source = metadata.get("parent_block_size") or cell_size
    if parent_size_source is None:
        raise ValueError("Irregular tbms-config-text export requires parent block size metadata or explicit Cell Size X/Y/Z.")
    parent_size = np.asarray(parent_size_source, dtype=float)
    if parent_size.shape != (3,) or np.any(parent_size <= 0):
        raise ValueError(f"Invalid parent block size for irregular export: {parent_size_source}")

    if origin is not None:
        origin_arr = np.asarray(origin, dtype=float)
    elif metadata.get("minimum_corner") is not None:
        origin_arr = np.asarray(metadata["minimum_corner"], dtype=float)
    elif metadata.get("minimum_parent_centroid") is not None:
        origin_arr = np.asarray(metadata["minimum_parent_centroid"], dtype=float) - 0.5 * parent_size
    else:
        xyz_min = clean[list(xyz_cols)].to_numpy(dtype=float).min(axis=0)
        origin_arr = np.floor(xyz_min / parent_size) * parent_size
    if origin_arr.shape != (3,):
        raise ValueError(f"Invalid origin values for irregular export: {origin_arr}")

    xyz = clean[list(xyz_cols)].to_numpy(dtype=float)
    raw_parent_idx = (xyz - origin_arr) / parent_size
    parent_idx = np.floor(raw_parent_idx + index_tolerance).astype(np.int64)
    local = xyz - (origin_arr + parent_idx.astype(float) * parent_size)
    local = np.clip(local, 0.0, parent_size)

    widths = np.empty_like(xyz, dtype=float)
    parent_keys, inverse = np.unique(parent_idx, axis=0, return_inverse=True)
    axis_tolerance = max(float(index_tolerance), 1e-9)
    for group_index in range(len(parent_keys)):
        group_positions = np.flatnonzero(inverse == group_index)
        group_local = local[group_positions]
        for axis in range(3):
            widths[group_positions, axis] = _infer_irregular_axis_widths(
                group_local[:, axis],
                float(parent_size[axis]),
                axis_tolerance,
            )

    lower = xyz - 0.5 * widths
    upper = xyz + 0.5 * widths
    if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper)) or np.any(widths <= 0):
        raise ValueError("Could not infer valid irregular block extents from CSV coordinates.")

    if metadata.get("size_in_parent_blocks") is not None:
        dims = np.asarray(metadata["size_in_parent_blocks"], dtype=np.int64)
    elif len(parent_idx):
        dims = parent_idx.max(axis=0) + 1
    else:
        dims = np.ones(3, dtype=np.int64)
    dims = np.maximum(dims.astype(np.int64), 1)
    extents = dims * parent_size
    grid_tuples = [tuple(map(int, row)) for row in parent_idx]
    width_source = "centroid_hierarchy_inference"
    grid_index_source = "irregular_xyz_metadata_parent" if _metadata_supports_tbms_irregular(metadata) else "irregular_xyz_hierarchy_inference"

    return {
        "df": clean.reset_index(drop=True),
        "idx": parent_idx,
        "grid_tuples": grid_tuples,
        "origin": origin_arr,
        "cell": parent_size,
        "dims": dims.astype(int),
        "extents": extents,
        "duplicates": 0,
        "max_index_error": 0.0,
        "is_irregular": True,
        "irregular_lower": lower,
        "irregular_upper": upper,
        "irregular_widths": widths,
        "irregular_width_source": width_source,
        "grid_index_source": grid_index_source,
    }


def _modal_nonblank_value(series: pd.Series):
    counts: Dict[object, int] = {}
    first_values: Dict[object, object] = {}
    first_order: Dict[object, int] = {}
    for position, value in enumerate(series):
        if pd.isna(value):
            continue
        if isinstance(value, str):
            key = value.strip()
            if not key:
                continue
        else:
            try:
                hash(value)
                key = value
            except TypeError:
                key = str(value)
        if key not in counts:
            counts[key] = 0
            first_values[key] = value
            first_order[key] = position
        counts[key] += 1
    if not counts:
        return np.nan
    best_key = max(counts, key=lambda key: (counts[key], -first_order[key]))
    return first_values[best_key]


def _regularize_to_base_cell_grid(
    df: pd.DataFrame,
    xyz_cols: Sequence[str],
    cell_size: Sequence[float] | None,
    origin: Sequence[float] | None,
    value_cols: Sequence[str] | None,
    column_types: Mapping[str, object] | None,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    if cell_size is None:
        raise ValueError(
            "Regularize to base block size requires explicit Cell Size X/Y/Z values. "
            "Set the BMF cell size to the intended base block dimensions, or disable regularization."
        )
    cell = np.asarray(cell_size, dtype=float)
    if cell.shape != (3,) or np.any(cell <= 0):
        raise ValueError(f"Invalid base block cell size values for regularization: {cell_size}")

    clean = _to_numeric_xyz(df, xyz_cols)
    if clean.empty:
        raise ValueError("No valid XYZ rows found after numeric conversion.")

    xyz = clean[list(xyz_cols)].to_numpy(dtype=float)
    if origin is None:
        base_origin = np.floor(xyz.min(axis=0) / cell) * cell
    else:
        base_origin = np.asarray(origin, dtype=float)
        if base_origin.shape != (3,):
            raise ValueError(f"Invalid origin values for regularization: {origin}")

    idx = np.floor((xyz - base_origin) / cell + 1e-9).astype(np.int64)
    if np.any(idx < 0):
        raise ValueError(
            "Regularized base-cell indices include negative values. Check the explicit origin and cell size."
        )

    value_columns = _classify_columns(clean, xyz_cols, value_cols)
    work = clean.copy()
    index_cols = ["__bmf_base_i", "__bmf_base_j", "__bmf_base_k"]
    for axis, col in enumerate(index_cols):
        work[col] = idx[:, axis]

    grouped = work.groupby(index_cols, sort=True, dropna=False)
    groups = grouped.size().reset_index(name="__bmf_source_rows")
    out = pd.DataFrame(index=groups.index)
    centers = groups[index_cols].to_numpy(dtype=float) * cell + base_origin + 0.5 * cell
    for axis, col in enumerate(xyz_cols):
        out[col] = centers[:, axis]

    normalized_overrides = _normalize_export_type_overrides(column_types)
    aggregation_modes: Dict[str, str] = {}

    for col in value_columns:
        forced_type = normalized_overrides.get(col)
        numeric_values = pd.to_numeric(work[col], errors="coerce")
        missing_mask = work[col].isna() | work[col].astype(str).str.strip().eq("")
        numeric_compatible = bool((numeric_values.notna() | missing_mask).all() and numeric_values.notna().any())
        should_average = forced_type == "double" or (forced_type is None and numeric_compatible)
        if should_average:
            numeric_work = work[index_cols].copy()
            numeric_work["__bmf_value"] = numeric_values
            out[col] = numeric_work.groupby(index_cols, sort=True, dropna=False)["__bmf_value"].mean().to_numpy()
            aggregation_modes[col] = "mean"
        else:
            out[col] = grouped[col].agg(_modal_nonblank_value).to_numpy()
            aggregation_modes[col] = "mode"

    summary = {
        "enabled": True,
        "input_rows": int(len(clean)),
        "output_rows": int(len(out)),
        "merged_rows": int(len(clean) - len(out)),
        "origin": [float(x) for x in base_origin],
        "cell_size": [float(x) for x in cell],
        "aggregation": aggregation_modes,
    }
    return out, summary


def _classify_columns(df: pd.DataFrame, xyz_cols: Sequence[str], value_cols: Sequence[str] | None) -> List[str]:
    geometry_cols = set(GRID_INDEX_COLUMNS) | set(GEOMETRY_EXTENT_COLUMNS)
    if value_cols:
        cols = [c for c in value_cols if c in df.columns and c not in xyz_cols and c not in geometry_cols]
        if not cols:
            raise ValueError("No valid value columns were provided.")
        return cols
    return [c for c in df.columns if c not in xyz_cols and c not in geometry_cols]


def _to_builtin(value):
    if isinstance(value, np.generic):
        return value.item()
    return value


def _align_offset(offset: int, alignment: int = 8) -> int:
    if alignment <= 1:
        return offset
    return ((offset + alignment - 1) // alignment) * alignment


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _json_bytes(payload: object) -> bytes:
    return json.dumps(payload, indent=2, sort_keys=True, default=_json_default).encode("utf-8")


def _normalize_string_records(df: pd.DataFrame, cols: Sequence[str]) -> Dict[str, List[str | None]]:
    records: Dict[str, List[str | None]] = {}
    for col in cols:
        series = df[col]
        records[col] = [None if pd.isna(v) else str(v) for v in series]
    return records


def _tbms_escape_string(value: str) -> str:
    text = str(value)
    text = text.replace("\\", "\\\\")
    text = text.replace('"', '\\"')
    text = text.replace("\r", "\\r")
    text = text.replace("\n", "\\n")
    text = text.replace("\t", "\\t")
    return f'"{text}"'


def _tbms_config_text(value: object, indent: int = 0) -> str:
    pad = " " * indent
    if isinstance(value, dict):
        if not value:
            return "{}"
        parts = ["{"]
        items = list(value.items())
        for index, (key, item) in enumerate(items):
            rendered = _tbms_config_text(item, indent + 2)
            suffix = "," if index < len(items) - 1 else ""
            parts.append(f"{pad}  {_tbms_escape_string(str(key))} = {rendered}{suffix}")
        parts.append(f"{pad}}}")
        return "\n".join(parts)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (np.bool_,)):
        return "true" if bool(value) else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError("TBMS config cannot serialize non-finite floats.")
        return format(numeric, ".15g")
    return _tbms_escape_string(str(value))


def _tbms_linear_indices(prepared: Dict[str, object]) -> np.ndarray:
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    idx = np.asarray(prepared["idx"], dtype=np.int64)
    return idx[:, 0] + dims[0] * (idx[:, 1] + dims[1] * idx[:, 2])


def _tbms_dense_field_values(
    prepared: Dict[str, object],
    series: pd.Series,
    default_value: object,
    dtype: np.dtype,
) -> Dict[str, object]:
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    row_count = int(np.prod(dims, dtype=np.int64))
    linear = _tbms_linear_indices(prepared)
    if len(linear) != len(series):
        raise ValueError("BMF field values do not align with prepared grid rows.")
    order = np.argsort(linear, kind="stable")
    return {
        "row_count": row_count,
        "linear": linear[order].astype(np.int64, copy=False),
        "source_values": np.asarray(series, dtype=dtype)[order],
    }


def _tbms_export_field_values(
    prepared: Dict[str, object],
    series: pd.Series,
    default_value: object,
    dtype: np.dtype,
) -> object:
    if prepared.get("is_irregular"):
        return np.asarray(series.reset_index(drop=True), dtype=dtype)
    return _tbms_dense_field_values(prepared, series, default_value, dtype)


def _metadata_supports_tbms_irregular(metadata: Mapping[str, object]) -> bool:
    return bool(metadata.get("parent_block_size") and metadata.get("subblock_factors"))


def _backend_requires_dense_grid(backend: str) -> bool:
    return str(backend or "").strip() in DENSE_EXPORT_BACKENDS


def _validate_dense_export_size(prepared: Dict[str, object], value_cols: Sequence[str], backend: str | None = None) -> None:
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    row_count = int(np.prod(dims, dtype=np.int64))
    field_count = max(len(value_cols), 1)
    estimated_bytes = row_count * 8 * field_count
    if estimated_bytes <= MAX_DENSE_EXPORT_BYTES:
        return

    cell = np.asarray(prepared.get("cell", [1.0, 1.0, 1.0]), dtype=float)
    valid_rows = len(prepared.get("df", []))
    backend_text = f" for backend {backend!r}" if backend else ""
    raise MemoryError(
        f"BMF export would create an oversized dense regular-grid file{backend_text}. "
        f"Inferred dimensions={tuple(int(value) for value in dims)} ({row_count:,} cells) from "
        f"{valid_rows:,} valid CSV rows using cell_size={cell.tolist()}. "
        f"Estimated uncompressed dense value payload for {field_count} field(s) is {_format_byte_size(estimated_bytes)}. "
        "This limit is about the dense output implied by the CSV geometry, not just peak RAM. The writer streams pages, "
        "but a dense Vulcan-style BMF still needs storage for every implied regular-grid cell. The CSV likely has very "
        "fine coordinate spacing, irregular/sub-block coordinates, or an incorrect XYZ/cell-size setup. "
        "For Vulcan-compatible dense BMF output, export a regular base-block grid or set Cell Size X/Y/Z to the intended "
        "base block dimensions, not sub-block spacing. For a row-indexed Anterpolator container that avoids dense grid "
        "allocation, choose the tbms-experimental backend."
    )


def _validate_selected_csv_read_size(
    path: Path,
    header_line: int,
    selected_read_cols: Sequence[str] | None,
    backend: str,
    regularize_to_base_block: bool,
    source_row_mode: bool = False,
    progress_callback=None,
) -> None:
    if not selected_read_cols:
        return
    row_count = _count_csv_data_rows(path, header_line, progress_callback, 3, 5)
    if row_count is None:
        return
    column_count = len(selected_read_cols)
    estimated_object_block_bytes = int(row_count) * int(column_count) * 8
    if source_row_mode:
        estimated_coordinate_bytes = int(row_count) * 3 * 8
        if estimated_coordinate_bytes > MAX_SOURCE_ROW_PREP_BYTES:
            raise MemoryError(
                f"The source-row tbms-config-text export is too large for the current in-memory irregular row preparation. "
                f"It would need to allocate at least {_format_byte_size(estimated_coordinate_bytes)} just for the "
                f"XYZ coordinate matrix ({row_count:,} rows x 3 float64 columns) before writing. "
                "This path does not create an implied dense grid, but this particular export could not use the streaming "
                "irregular writer. Streaming source-row BMF export requires explicit row geometry, Leapfrog parent/sub-block "
                "metadata, or an explicit parent Cell Size X/Y/Z. Increasing ANTERPOLATOR_BMF_MAX_SOURCE_ROW_PREP_BYTES "
                "only helps on machines with enough RAM."
            )
    if estimated_object_block_bytes <= MAX_SELECTED_CSV_OBJECT_BYTES:
        return

    if regularize_to_base_block:
        operation = "regularized dense BMF export"
    elif source_row_mode:
        operation = "source-row BMF export"
    elif _backend_requires_dense_grid(backend):
        operation = "dense BMF export"
    else:
        operation = "source-row BMF export"
    raise MemoryError(
        f"The selected CSV columns are too wide for {operation} with backend {backend!r}. "
        f"The exporter would need to load {row_count:,} row(s) x {column_count:,} selected column(s) before writing, "
        f"which is at least {_format_byte_size(estimated_object_block_bytes)} just for pandas object references. "
        "This source-row path avoids creating an implied dense grid, but the current standalone writer still has to "
        "materialize the selected CSV columns before writing when the streaming irregular writer cannot be used. "
        "Streaming source-row BMF export requires explicit row geometry, Leapfrog parent/sub-block metadata, or an explicit "
        "parent Cell Size X/Y/Z. Select fewer value columns and export in smaller batches, or increase "
        "ANTERPOLATOR_BMF_MAX_SELECTED_CSV_BYTES only if the machine has enough RAM."
    )


def _tbms_encode_export_field(
    prepared: Dict[str, object],
    series: pd.Series,
    null_float: float,
    forced_type: str | None = None,
) -> Dict[str, object]:
    series = series.reset_index(drop=True)
    normalized_forced_type = _normalize_export_field_type(forced_type)

    def missing_or_blank_mask(values: pd.Series) -> pd.Series:
        return values.isna() | values.astype(str).str.strip().eq("")

    def invalid_numeric_values(values: pd.Series, numeric_values: pd.Series) -> List[str]:
        invalid_mask = ~missing_or_blank_mask(values) & numeric_values.isna()
        return sorted(pd.unique(values.loc[invalid_mask].astype(str)))[:5]

    if normalized_forced_type == "boolean":
        boolean_series = series.fillna(False)
        if pd.api.types.is_bool_dtype(series):
            boolean_series = boolean_series.astype(bool)
        else:
            truthy = {"1", "true", "t", "yes", "y"}
            falsy = {"0", "false", "f", "no", "n", ""}
            normalized_text = boolean_series.astype(str).str.strip().str.lower()
            invalid_mask = ~series.isna() & ~normalized_text.isin(truthy | falsy)
            if invalid_mask.any():
                invalid_values = sorted(pd.unique(series.loc[invalid_mask].astype(str)))[:5]
                raise ValueError(
                    f"Column {series.name!r} cannot be exported as boolean. Invalid values include: {invalid_values}"
                )
            boolean_series = normalized_text.isin(truthy)
        encoded = boolean_series.astype(np.uint8)
        return {
            "field_type": "boolean",
            "default": 0,
            "dtype": np.dtype("u1"),
            "values": _tbms_export_field_values(prepared, encoded, 0, np.dtype("u1")),
        }

    numeric = pd.to_numeric(series, errors="coerce")
    missing_mask = missing_or_blank_mask(series)
    numeric_compatible = bool((numeric.notna() | missing_mask).all())
    if normalized_forced_type == "int":
        if not numeric_compatible:
            invalid_values = invalid_numeric_values(series, numeric)
            raise ValueError(
                f"Column {series.name!r} cannot be exported as int because it contains non-numeric values. "
                f"Invalid values include: {invalid_values}"
            )
        finite = numeric.dropna().to_numpy(dtype=float)
        if finite.size > 0 and not np.allclose(finite, np.rint(finite), atol=1e-9, rtol=0.0):
            raise ValueError(f"Column {series.name!r} cannot be exported as int because it contains non-integer numeric values.")
        min_value = int(np.min(finite)) if finite.size else 0
        max_value = int(np.max(finite)) if finite.size else 0
        if not (-2147483648 <= min_value and max_value <= 2147483647):
            raise ValueError(f"Column {series.name!r} exceeds the supported int32 range for int export.")
        converted = numeric.fillna(int(null_float)).astype(np.int32)
        return {
            "field_type": "int",
            "default": int(null_float),
            "dtype": np.dtype("<i4"),
            "values": _tbms_export_field_values(prepared, converted, int(null_float), np.dtype("<i4")),
        }

    if normalized_forced_type == "double":
        if not numeric_compatible:
            invalid_values = invalid_numeric_values(series, numeric)
            raise ValueError(
                f"Column {series.name!r} cannot be exported as double because it contains non-numeric values. "
                f"Invalid values include: {invalid_values}"
            )
        converted = numeric.fillna(float(null_float)).astype(np.float64)
        return {
            "field_type": "double",
            "default": float(null_float),
            "dtype": np.dtype("<f8"),
            "values": _tbms_export_field_values(prepared, converted, float(null_float), np.dtype("<f8")),
        }

    if normalized_forced_type == "string":
        labels = [str(value) for value in pd.unique(series.dropna())]
        if len(labels) > 32766:
            raise ValueError(f"Column {series.name!r} has too many categorical values for string export.")
        code_map = {label: index + 1 for index, label in enumerate(labels)}
        codes = series.map(lambda value: 0 if pd.isna(value) else code_map[str(value)]).astype(np.int16)
        field_info = {
            "field_type": "namedshort",
            "default": 0,
            "dtype": np.dtype("<i2"),
            "values": _tbms_export_field_values(prepared, codes, 0, np.dtype("<i2")),
        }
        for label, code in code_map.items():
            field_info[f"string_{code}"] = label
        field_info["string_0"] = ""
        return field_info

    if pd.api.types.is_bool_dtype(series):
        encoded = series.fillna(False).astype(bool).astype(np.uint8)
        return {
            "field_type": "boolean",
            "default": 0,
            "dtype": np.dtype("u1"),
            "values": _tbms_export_field_values(prepared, encoded, 0, np.dtype("u1")),
        }

    if numeric_compatible:
        finite = numeric.dropna().to_numpy(dtype=float)
        integer_like = finite.size > 0 and np.allclose(finite, np.rint(finite), atol=1e-9, rtol=0.0)
        if integer_like:
            min_value = int(np.min(finite)) if finite.size else 0
            max_value = int(np.max(finite)) if finite.size else 0
            if -2147483648 <= min_value and max_value <= 2147483647:
                converted = numeric.fillna(int(null_float)).astype(np.int32)
                return {
                    "field_type": "int",
                    "default": int(null_float),
                    "dtype": np.dtype("<i4"),
                    "values": _tbms_export_field_values(prepared, converted, int(null_float), np.dtype("<i4")),
                }
        converted = numeric.fillna(float(null_float)).astype(np.float64)
        return {
            "field_type": "double",
            "default": float(null_float),
            "dtype": np.dtype("<f8"),
            "values": _tbms_export_field_values(prepared, converted, float(null_float), np.dtype("<f8")),
        }

    labels = [str(value) for value in pd.unique(series.dropna())]
    if len(labels) > 32766:
        raise ValueError(f"Column {series.name!r} has too many categorical values for namedshort export.")
    code_map = {label: index + 1 for index, label in enumerate(labels)}
    codes = series.map(lambda value: 0 if pd.isna(value) else code_map[str(value)]).astype(np.int16)
    field_info = {
        "field_type": "namedshort",
        "default": 0,
        "dtype": np.dtype("<i2"),
        "values": _tbms_export_field_values(prepared, codes, 0, np.dtype("<i2")),
    }
    for label, code in code_map.items():
        field_info[f"string_{code}"] = label
    field_info["string_0"] = ""
    return field_info


def _tbms_value_page(values: np.ndarray, dtype: np.dtype, default_value: object) -> List[bytes]:
    payload_size = 2048
    capacity = payload_size // dtype.itemsize
    page_count = max(1, int(math.ceil(len(values) / capacity)))
    padded = np.full(page_count * capacity, default_value, dtype=dtype)
    padded[: len(values)] = values.astype(dtype, copy=False)

    pages: List[bytes] = []
    for page_index in range(page_count):
        start = page_index * capacity
        stop = start + capacity
        payload = padded[start:stop].astype(dtype, copy=False).tobytes(order="C")
        pages.append(struct.pack("<II", 512, 0) + payload)
    return pages


def _iter_tbms_value_pages(values: object, dtype: np.dtype, default_value: object) -> Iterable[bytes]:
    payload_size = 2048
    capacity = payload_size // dtype.itemsize
    if isinstance(values, dict):
        row_count = int(values.get("row_count", 0) or 0)
        linear = np.asarray(values.get("linear", []), dtype=np.int64)
        source_values = np.asarray(values.get("source_values", []), dtype=dtype)
        page_count = max(1, int(math.ceil(row_count / capacity)))
        position = 0
        source_count = len(linear)
        for page_index in range(page_count):
            start = page_index * capacity
            stop = min(start + capacity, row_count)
            page_values = np.full(capacity, default_value, dtype=dtype)
            while position < source_count and linear[position] < start:
                position += 1
            page_stop = position
            while page_stop < source_count and linear[page_stop] < stop:
                page_stop += 1
            if page_stop > position:
                offsets = (linear[position:page_stop] - start).astype(np.int64, copy=False)
                page_values[offsets] = source_values[position:page_stop]
            position = page_stop
            yield struct.pack("<II", 512, 0) + page_values.tobytes(order="C")
        return

    dense_values = np.asarray(values, dtype=dtype)
    page_count = max(1, int(math.ceil(len(dense_values) / capacity)))
    for page_index in range(page_count):
        start = page_index * capacity
        stop = min(start + capacity, len(dense_values))
        page_values = np.full(capacity, default_value, dtype=dtype)
        if stop > start:
            page_values[: stop - start] = dense_values[start:stop].astype(dtype, copy=False)
        yield struct.pack("<II", 512, 0) + page_values.tobytes(order="C")


def _tbms_pointer_page(child_offsets: Sequence[int]) -> bytes:
    entries = [257] + [int(offset) for offset in child_offsets[:256]]
    if len(entries) < 257:
        entries.extend([0] * (257 - len(entries)))
    return struct.pack("<257Q", *entries)


class _TbmsStreamingPageField:
    def __init__(self, dtype: np.dtype, default_value: object, allocate_page):
        self.dtype = np.dtype(dtype)
        self.default_value = default_value
        self.allocate_page = allocate_page
        self.capacity = 2048 // self.dtype.itemsize
        self.buffer = np.full(self.capacity, default_value, dtype=self.dtype)
        self.buffer_count = 0
        self.leaf_offsets: list[int] = []

    def append(self, values: Sequence[object] | np.ndarray) -> None:
        incoming = np.asarray(values, dtype=self.dtype)
        position = 0
        total = len(incoming)
        while position < total:
            available = self.capacity - self.buffer_count
            take = min(available, total - position)
            self.buffer[self.buffer_count:self.buffer_count + take] = incoming[position:position + take]
            self.buffer_count += take
            position += take
            if self.buffer_count == self.capacity:
                self._flush_full_buffer()

    def _flush_full_buffer(self) -> None:
        payload = self.buffer.astype(self.dtype, copy=False).tobytes(order="C")
        self.leaf_offsets.append(self.allocate_page(struct.pack("<II", 512, 0) + payload))
        self.buffer.fill(self.default_value)
        self.buffer_count = 0

    def finish(self) -> list[int]:
        if self.buffer_count > 0 or not self.leaf_offsets:
            payload = self.buffer.astype(self.dtype, copy=False).tobytes(order="C")
            self.leaf_offsets.append(self.allocate_page(struct.pack("<II", 512, 0) + payload))
            self.buffer.fill(self.default_value)
            self.buffer_count = 0
        return list(self.leaf_offsets)


def _tbms_page_tree_root_from_leaf_offsets(leaf_offsets: Sequence[int], allocate_page) -> int:
    current_offsets = list(leaf_offsets)
    if not current_offsets:
        current_offsets = [allocate_page(_tbms_value_page(np.empty(0, dtype="<f4"), np.dtype("<f4"), 0.0)[0])]
    while len(current_offsets) > 1:
        next_offsets: list[int] = []
        for start in range(0, len(current_offsets), 256):
            next_offsets.append(allocate_page(_tbms_pointer_page(current_offsets[start:start + 256])))
        current_offsets = next_offsets
    return int(current_offsets[0])


def _stream_missing_or_blank_mask(series: pd.Series) -> pd.Series:
    return series.isna() | series.astype(str).str.strip().eq("")


def _new_stream_field_state(name: str, forced_type: str | None) -> dict[str, object]:
    return {
        "name": name,
        "forced_type": _normalize_export_field_type(forced_type),
        "numeric_compatible": True,
        "integer_like": True,
        "has_numeric": False,
        "min": None,
        "max": None,
        "labels": {},
        "too_many_labels": False,
        "invalid_values": [],
    }


def _remember_stream_invalid_values(state: dict[str, object], values: pd.Series) -> None:
    invalid_values = state.setdefault("invalid_values", [])
    if not isinstance(invalid_values, list) or len(invalid_values) >= 5:
        return
    for value in pd.unique(values.astype(str)):
        if len(invalid_values) >= 5:
            break
        text = str(value)
        if text not in invalid_values:
            invalid_values.append(text)


def _scan_stream_field_values(state: dict[str, object], series: pd.Series) -> None:
    missing_mask = _stream_missing_or_blank_mask(series)
    label_values = series.loc[series.notna()]
    labels = state.get("labels")
    if isinstance(labels, dict) and not state.get("too_many_labels"):
        for value in pd.unique(label_values.astype(str)):
            label = str(value)
            if label not in labels:
                if len(labels) >= 32766:
                    state["too_many_labels"] = True
                    break
                labels[label] = len(labels) + 1

    numeric = pd.to_numeric(series, errors="coerce")
    invalid_numeric = ~missing_mask & numeric.isna()
    if invalid_numeric.any():
        state["numeric_compatible"] = False
        _remember_stream_invalid_values(state, series.loc[invalid_numeric])
        return

    finite = numeric.loc[numeric.notna()].to_numpy(dtype=float)
    if finite.size == 0:
        return
    state["has_numeric"] = True
    min_value = int(np.min(finite)) if np.allclose(finite, np.rint(finite), atol=1e-9, rtol=0.0) else float(np.min(finite))
    max_value = int(np.max(finite)) if np.allclose(finite, np.rint(finite), atol=1e-9, rtol=0.0) else float(np.max(finite))
    state["min"] = min_value if state.get("min") is None else min(float(state["min"]), float(min_value))
    state["max"] = max_value if state.get("max") is None else max(float(state["max"]), float(max_value))
    if not np.allclose(finite, np.rint(finite), atol=1e-9, rtol=0.0):
        state["integer_like"] = False


def _finalize_stream_field_state(state: dict[str, object], null_float: float) -> dict[str, object]:
    name = str(state.get("name"))
    forced_type = state.get("forced_type")
    numeric_compatible = bool(state.get("numeric_compatible"))
    invalid_values = state.get("invalid_values") if isinstance(state.get("invalid_values"), list) else []

    if forced_type == "boolean":
        return {"field_type": "boolean", "default": 0, "dtype": np.dtype("u1")}
    if forced_type == "int":
        if not numeric_compatible:
            raise ValueError(f"Column {name!r} cannot be exported as int because it contains non-numeric values. Invalid values include: {invalid_values}")
        if not bool(state.get("integer_like", True)):
            raise ValueError(f"Column {name!r} cannot be exported as int because it contains non-integer numeric values.")
        min_value = int(float(state.get("min") or 0))
        max_value = int(float(state.get("max") or 0))
        if not (-2147483648 <= min_value and max_value <= 2147483647):
            raise ValueError(f"Column {name!r} exceeds the supported int32 range for int export.")
        return {"field_type": "int", "default": int(null_float), "dtype": np.dtype("<i4")}
    if forced_type == "double":
        if not numeric_compatible:
            raise ValueError(f"Column {name!r} cannot be exported as double because it contains non-numeric values. Invalid values include: {invalid_values}")
        return {"field_type": "double", "default": float(null_float), "dtype": np.dtype("<f8")}
    if forced_type == "string":
        if state.get("too_many_labels"):
            raise ValueError(f"Column {name!r} has too many categorical values for string export.")
        labels = dict(state.get("labels") or {})
        return {"field_type": "namedshort", "default": 0, "dtype": np.dtype("<i2"), "labels": labels}

    if numeric_compatible:
        if bool(state.get("integer_like", True)) and bool(state.get("has_numeric", False)):
            min_value = int(float(state.get("min") or 0))
            max_value = int(float(state.get("max") or 0))
            if -2147483648 <= min_value and max_value <= 2147483647:
                return {"field_type": "int", "default": int(null_float), "dtype": np.dtype("<i4")}
        return {"field_type": "double", "default": float(null_float), "dtype": np.dtype("<f8")}

    if state.get("too_many_labels"):
        raise ValueError(f"Column {name!r} has too many categorical values for namedshort export.")
    labels = dict(state.get("labels") or {})
    return {"field_type": "namedshort", "default": 0, "dtype": np.dtype("<i2"), "labels": labels}


def _encode_stream_field_chunk(series: pd.Series, field_info: Mapping[str, object], null_float: float) -> np.ndarray:
    field_type = str(field_info.get("field_type"))
    if field_type == "boolean":
        truthy = {"1", "true", "t", "yes", "y"}
        falsy = {"0", "false", "f", "no", "n", ""}
        normalized = series.fillna(False).astype(str).str.strip().str.lower()
        invalid = ~series.isna() & ~normalized.isin(truthy | falsy)
        if invalid.any():
            invalid_values = sorted(pd.unique(series.loc[invalid].astype(str)))[:5]
            raise ValueError(f"Column {series.name!r} cannot be exported as boolean. Invalid values include: {invalid_values}")
        return normalized.isin(truthy).astype(np.uint8).to_numpy(dtype="u1", copy=False)
    if field_type == "int":
        return pd.to_numeric(series, errors="coerce").fillna(int(null_float)).astype(np.int32).to_numpy(dtype="<i4", copy=False)
    if field_type == "double":
        return pd.to_numeric(series, errors="coerce").fillna(float(null_float)).astype(np.float64).to_numpy(dtype="<f8", copy=False)
    labels = field_info.get("labels") if isinstance(field_info.get("labels"), Mapping) else {}
    return series.map(lambda value: 0 if pd.isna(value) else int(labels[str(value)])).astype(np.int16).to_numpy(dtype="<i2", copy=False)


def _classify_header_value_columns(
    header_names: Sequence[str],
    xyz_cols: Sequence[str],
    value_cols: Sequence[str] | None,
    geometry_columns: Sequence[str] | None = None,
) -> list[str]:
    geometry_cols = set(GRID_INDEX_COLUMNS) | set(GEOMETRY_EXTENT_COLUMNS) | set(geometry_columns or [])
    if value_cols:
        cols = [c for c in value_cols if c in header_names and c not in xyz_cols and c not in geometry_cols]
        if not cols:
            raise ValueError("No valid value columns were provided.")
        return cols
    return [c for c in header_names if c not in xyz_cols and c not in geometry_cols]


def _prepare_stream_irregular_chunk_geometry(
    chunk: pd.DataFrame,
    xyz_cols: Sequence[str],
    metadata: Mapping[str, object],
    cell_size: Sequence[float] | np.ndarray | None,
    origin: Sequence[float] | np.ndarray | None,
    index_tolerance: float,
    geometry_size_cols: Mapping[str, str] | None = None,
    geometry_extent_cols: Mapping[str, str] | None = None,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    clean = chunk.copy()
    for col in xyz_cols:
        clean[col] = pd.to_numeric(clean[col], errors="coerce")
    clean = clean.dropna(subset=list(xyz_cols))
    if clean.empty:
        return clean, np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=np.int64), np.empty((0, 3), dtype=float)

    extent_col_map = dict(geometry_extent_cols or {})
    if not extent_col_map and all(column in clean.columns for column in GEOMETRY_EXTENT_COLUMNS):
        extent_col_map = dict(zip(GEOMETRY_EXTENT_COLUMNS, GEOMETRY_EXTENT_COLUMNS))
    if all(extent_col_map.get(role) in clean.columns for role in GEOMETRY_EXTENT_COLUMNS):
        extent_values = clean[[extent_col_map[role] for role in GEOMETRY_EXTENT_COLUMNS]].apply(pd.to_numeric, errors="coerce")
        valid = extent_values.notna().all(axis=1)
        clean = clean.loc[valid].reset_index(drop=True)
        extent_values = extent_values.loc[valid].reset_index(drop=True)
        lower = extent_values[[extent_col_map["__lower_x"], extent_col_map["__lower_y"], extent_col_map["__lower_z"]]].to_numpy(dtype=float)
        upper = extent_values[[extent_col_map["__upper_x"], extent_col_map["__upper_y"], extent_col_map["__upper_z"]]].to_numpy(dtype=float)
        xyz = (lower + upper) * 0.5
        parent_size_source = metadata.get("parent_block_size") or cell_size
        parent_size = np.asarray(parent_size_source if parent_size_source is not None else np.nanmax(upper - lower, axis=0), dtype=float)
        origin_arr = np.asarray(origin if origin is not None else metadata.get("minimum_corner", np.nanmin(lower, axis=0)), dtype=float)
        parent_idx = np.floor((xyz - origin_arr) / parent_size + float(index_tolerance)).astype(np.int64)
        return clean, lower, upper, parent_idx, parent_size

    size_col_map = dict(geometry_size_cols or {})
    if all(size_col_map.get(role) in clean.columns for role in GEOMETRY_SIZE_COLUMN_ROLES):
        size_cols = [size_col_map[role] for role in GEOMETRY_SIZE_COLUMN_ROLES]
        size_values = clean[size_cols].apply(pd.to_numeric, errors="coerce")
        valid = size_values.notna().all(axis=1)
        clean = clean.loc[valid].reset_index(drop=True)
        size_values = size_values.loc[valid].reset_index(drop=True)
        widths = size_values[size_cols].to_numpy(dtype=float)
        xyz = clean[list(xyz_cols)].to_numpy(dtype=float)
        lower = xyz - 0.5 * widths
        upper = xyz + 0.5 * widths
        parent_size_source = metadata.get("parent_block_size") or cell_size
        parent_size = np.asarray(parent_size_source if parent_size_source is not None else np.nanmax(widths, axis=0), dtype=float)
        origin_arr = np.asarray(origin if origin is not None else metadata.get("minimum_corner", np.nanmin(lower, axis=0)), dtype=float)
        parent_idx = np.floor((xyz - origin_arr) / parent_size + float(index_tolerance)).astype(np.int64)
        return clean, lower, upper, parent_idx, parent_size

    parent_size_source = metadata.get("parent_block_size") or cell_size
    if parent_size_source is None:
        raise ValueError("Streaming irregular export requires parent block size metadata or explicit Cell Size X/Y/Z when no row extent/dimension columns are available.")
    parent_size = np.asarray(parent_size_source, dtype=float)
    if origin is not None:
        origin_arr = np.asarray(origin, dtype=float)
    elif metadata.get("minimum_corner") is not None:
        origin_arr = np.asarray(metadata["minimum_corner"], dtype=float)
    elif metadata.get("minimum_parent_centroid") is not None:
        origin_arr = np.asarray(metadata["minimum_parent_centroid"], dtype=float) - 0.5 * parent_size
    else:
        xyz_min = clean[list(xyz_cols)].to_numpy(dtype=float).min(axis=0)
        origin_arr = np.floor(xyz_min / parent_size) * parent_size
    xyz = clean[list(xyz_cols)].to_numpy(dtype=float)
    parent_idx = np.floor((xyz - origin_arr) / parent_size + float(index_tolerance)).astype(np.int64)
    local = np.clip(xyz - (origin_arr + parent_idx.astype(float) * parent_size), 0.0, parent_size)
    widths = np.empty_like(xyz, dtype=float)
    axis_tolerance = max(float(index_tolerance), 1e-9)
    for axis in range(3):
        widths[:, axis] = _infer_irregular_axis_widths(local[:, axis], float(parent_size[axis]), axis_tolerance)
    lower = xyz - 0.5 * widths
    upper = xyz + 0.5 * widths
    return clean.reset_index(drop=True), lower, upper, parent_idx, parent_size


def _scan_stream_irregular_export(
    input_csv: Path,
    delimiter: str | None,
    header_line: int,
    selected_read_cols: Sequence[str],
    xyz_cols: Sequence[str],
    value_cols: Sequence[str],
    metadata: Mapping[str, object],
    cell_size: Sequence[float] | np.ndarray | None,
    origin: Sequence[float] | np.ndarray | None,
    index_tolerance: float,
    geometry_size_cols: Mapping[str, str] | None,
    geometry_extent_cols: Mapping[str, str] | None,
    column_types: Mapping[str, object] | None,
    value_exceptions: Mapping[str, Mapping[str, object]] | None,
    null_float: float,
    total_rows: int | None,
    progress_callback=None,
) -> dict[str, object]:
    normalized_overrides = _normalize_export_type_overrides(column_types)
    field_states = {
        name: _new_stream_field_state(name, normalized_overrides.get(name))
        for name in value_cols
    }
    row_count = 0
    min_lower = np.full(3, np.inf, dtype=float)
    max_upper = np.full(3, -np.inf, dtype=float)
    max_parent_idx = np.full(3, -1, dtype=np.int64)
    parent_size = None

    _emit_progress(progress_callback, 30, 100, "Scanning CSV chunks for streaming BMF schema and geometry...")
    for chunk in _iter_csv_chunks_for_export(
        input_csv,
        delimiter=delimiter,
        header_line=header_line,
        usecols=selected_read_cols,
        progress_callback=progress_callback,
        progress_start=30,
        progress_end=45,
        progress_label="Scanning CSV rows",
        total_rows=total_rows,
    ):
        if value_exceptions:
            chunk = _apply_value_exceptions(chunk, value_exceptions)
        clean, lower, upper, parent_idx, chunk_parent_size = _prepare_stream_irregular_chunk_geometry(
            chunk,
            xyz_cols=xyz_cols,
            metadata=metadata,
            cell_size=cell_size,
            origin=origin,
            index_tolerance=index_tolerance,
            geometry_size_cols=geometry_size_cols,
            geometry_extent_cols=geometry_extent_cols,
        )
        if clean.empty:
            continue
        widths = upper - lower
        if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper)) or np.any(widths <= 0):
            raise ValueError("Invalid irregular row geometry encountered while scanning CSV chunks.")
        row_count += len(clean)
        min_lower = np.minimum(min_lower, np.nanmin(lower, axis=0))
        max_upper = np.maximum(max_upper, np.nanmax(upper, axis=0))
        if len(parent_idx):
            max_parent_idx = np.maximum(max_parent_idx, np.nanmax(parent_idx, axis=0).astype(np.int64))
        chunk_parent_size = np.asarray(chunk_parent_size, dtype=float)
        parent_size = chunk_parent_size if parent_size is None else np.maximum(np.asarray(parent_size, dtype=float), chunk_parent_size)
        for name in value_cols:
            _scan_stream_field_values(field_states[name], clean[name])

    if row_count <= 0:
        raise ValueError("No valid irregular source rows found while streaming CSV export.")
    parent_size_arr = np.asarray(parent_size if parent_size is not None else cell_size, dtype=float)
    if parent_size_arr.shape != (3,) or np.any(parent_size_arr <= 0) or np.any(~np.isfinite(parent_size_arr)):
        raise ValueError("Could not determine a valid parent block size for streaming irregular export.")
    if origin is not None:
        origin_arr = np.asarray(origin, dtype=float)
    elif metadata.get("minimum_corner") is not None:
        origin_arr = np.asarray(metadata["minimum_corner"], dtype=float)
    else:
        origin_arr = min_lower
    if origin_arr.shape != (3,) or np.any(~np.isfinite(origin_arr)):
        raise ValueError("Could not determine a valid origin for streaming irregular export.")

    if metadata.get("size_in_parent_blocks") is not None:
        dims = np.asarray(metadata["size_in_parent_blocks"], dtype=np.int64)
    else:
        span_dims = np.ceil((max_upper - origin_arr) / parent_size_arr - 1e-9).astype(np.int64)
        dims = np.maximum(span_dims, max_parent_idx + 1)
    dims = np.maximum(dims.astype(np.int64), 1)
    extents = np.maximum(max_upper - origin_arr, dims.astype(float) * parent_size_arr)
    field_infos = {
        name: _finalize_stream_field_state(state, null_float)
        for name, state in field_states.items()
    }
    return {
        "row_count": int(row_count),
        "origin": origin_arr,
        "cell": parent_size_arr,
        "dims": dims,
        "extents": extents,
        "field_infos": field_infos,
    }


def _write_streaming_irregular_tbms_config_text(
    input_csv: Path,
    prepared: Mapping[str, object],
    xyz_cols: Sequence[str],
    value_cols: Sequence[str],
    selected_read_cols: Sequence[str],
    delimiter: str | None,
    header_line: int,
    metadata: Mapping[str, object],
    cell_size: Sequence[float] | np.ndarray | None,
    origin: Sequence[float] | np.ndarray | None,
    index_tolerance: float,
    geometry_size_cols: Mapping[str, str] | None,
    geometry_extent_cols: Mapping[str, str] | None,
    value_exceptions: Mapping[str, Mapping[str, object]] | None,
    out_bmf: Path,
    null_float: float,
    progress_callback=None,
    progress_start: int = 50,
    progress_end: int = 98,
    total_rows: int | None = None,
) -> dict[str, object]:
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    origin_arr = np.asarray(prepared["origin"], dtype=float)
    cell = np.asarray(prepared["cell"], dtype=float)
    extents = np.asarray(prepared["extents"], dtype=float)
    row_count = int(prepared["row_count"])
    field_infos = dict(prepared.get("field_infos") or {})

    first_page_offset = TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET
    config_pointer_header_offsets = (24, 40)
    page_size = 2056
    field_entries: list[dict[str, object]] = []
    categorical_fields = 0
    page_count = 0

    out_bmf.parent.mkdir(parents=True, exist_ok=True)
    with out_bmf.open("wb") as fh:
        _emit_progress(progress_callback, progress_start, 100, "Writing streaming BMF header...")
        header = bytearray(TBMS_HEADER_SIZE)
        header[: len(TBMS_SIGNATURE)] = TBMS_SIGNATURE
        struct.pack_into("<I", header, 12, 1)
        struct.pack_into("<I", header, 16, first_page_offset)
        fh.write(header)
        if fh.tell() < first_page_offset:
            fh.write(b"\x00" * (first_page_offset - fh.tell()))

        def allocate_page(blob: bytes) -> int:
            nonlocal page_count
            if len(blob) != page_size:
                raise ValueError("TBMS page blobs must be exactly 2056 bytes.")
            offset = first_page_offset + page_count * page_size
            if fh.tell() < offset:
                fh.write(b"\x00" * (offset - fh.tell()))
            fh.write(blob)
            page_count += 1
            return offset

        geometry_names = ["__lower_x", "__upper_x", "__lower_y", "__upper_y", "__lower_z", "__upper_z"]
        geometry_builders = {
            name: _TbmsStreamingPageField(np.dtype("<f4"), 0.0, allocate_page)
            for name in geometry_names
        }
        _emit_progress(progress_callback, progress_start, 100, "Streaming irregular geometry pages from CSV chunks...")
        for chunk in _iter_csv_chunks_for_export(
            input_csv,
            delimiter=delimiter,
            header_line=header_line,
            usecols=selected_read_cols,
            progress_callback=progress_callback,
            progress_start=progress_start,
            progress_end=_scale_progress(progress_start, progress_end, 0.35),
            progress_label="Writing geometry rows",
            total_rows=total_rows,
        ):
            if value_exceptions:
                chunk = _apply_value_exceptions(chunk, value_exceptions)
            clean, lower, upper, _parent_idx, _parent_size = _prepare_stream_irregular_chunk_geometry(
                chunk,
                xyz_cols=xyz_cols,
                metadata=metadata,
                cell_size=cell_size,
                origin=origin,
                index_tolerance=index_tolerance,
                geometry_size_cols=geometry_size_cols,
                geometry_extent_cols=geometry_extent_cols,
            )
            if clean.empty:
                continue
            geometry_builders["__lower_x"].append(lower[:, 0].astype("<f4", copy=False))
            geometry_builders["__upper_x"].append(upper[:, 0].astype("<f4", copy=False))
            geometry_builders["__lower_y"].append(lower[:, 1].astype("<f4", copy=False))
            geometry_builders["__upper_y"].append(upper[:, 1].astype("<f4", copy=False))
            geometry_builders["__lower_z"].append(lower[:, 2].astype("<f4", copy=False))
            geometry_builders["__upper_z"].append(upper[:, 2].astype("<f4", copy=False))

        for geometry_index, name in enumerate(geometry_names):
            root_offset = _tbms_page_tree_root_from_leaf_offsets(geometry_builders[name].finish(), allocate_page)
            field_entries.append({
                "entry_key": f"var_{geometry_index}",
                "entry": {
                    "name": name,
                    "type": "float",
                    "location": int(root_offset),
                    "default": 0.0,
                    "global": 0,
                    "read_only": 1,
                    "description": "Irregular block geometry extent",
                },
            })

        value_entry_offset = len(field_entries)
        value_builders = {
            name: _TbmsStreamingPageField(np.dtype(field_infos[name]["dtype"]), field_infos[name]["default"], allocate_page)
            for name in value_cols
        }
        _emit_progress(progress_callback, _scale_progress(progress_start, progress_end, 0.36), 100, "Streaming BMF value pages from CSV chunks...")
        for chunk in _iter_csv_chunks_for_export(
            input_csv,
            delimiter=delimiter,
            header_line=header_line,
            usecols=selected_read_cols,
            progress_callback=progress_callback,
            progress_start=_scale_progress(progress_start, progress_end, 0.36),
            progress_end=_scale_progress(progress_start, progress_end, 0.82),
            progress_label="Writing value rows",
            total_rows=total_rows,
        ):
            if value_exceptions:
                chunk = _apply_value_exceptions(chunk, value_exceptions)
            clean, _lower, _upper, _parent_idx, _parent_size = _prepare_stream_irregular_chunk_geometry(
                chunk,
                xyz_cols=xyz_cols,
                metadata=metadata,
                cell_size=cell_size,
                origin=origin,
                index_tolerance=index_tolerance,
                geometry_size_cols=geometry_size_cols,
                geometry_extent_cols=geometry_extent_cols,
            )
            if clean.empty:
                continue
            for name in value_cols:
                value_builders[name].append(_encode_stream_field_chunk(clean[name], field_infos[name], null_float))

        for field_index, name in enumerate(value_cols):
            encoded = field_infos[name]
            root_offset = _tbms_page_tree_root_from_leaf_offsets(value_builders[name].finish(), allocate_page)
            entry = {
                "name": name,
                "type": encoded["field_type"],
                "location": int(root_offset),
                "default": encoded["default"],
                "global": 0,
                "read_only": 0,
                "description": f"Exported from CSV column {name}",
            }
            labels = encoded.get("labels") if isinstance(encoded.get("labels"), Mapping) else {}
            for label, code in labels.items():
                entry[f"string_{int(code)}"] = str(label)
            if encoded["field_type"] == "namedshort":
                entry["string_0"] = ""
                categorical_fields += 1
            field_entries.append({"entry_key": f"var_{value_entry_offset + field_index}", "entry": entry})

    now_text = datetime.now(timezone.utc).isoformat()
    config_object: dict[str, object] = {
        "created": now_text,
        "modified": now_text,
        "history_source": "Anterpolator streaming tbms-config-text export",
        "n_blocks": row_count,
        "n_schemas": 1,
        "is_irregular": 1,
        "geometry_encoding": "row_extents_lower_upper",
        "origin_x": 0.0,
        "origin_y": 0.0,
        "origin_z": 0.0,
        "lower_x": float(origin_arr[0]),
        "lower_y": float(origin_arr[1]),
        "lower_z": float(origin_arr[2]),
        "upper_x": float(origin_arr[0] + extents[0]),
        "upper_y": float(origin_arr[1] + extents[1]),
        "upper_z": float(origin_arr[2] + extents[2]),
        "voxel_length_x": float(cell[0]),
        "voxel_length_y": float(cell[1]),
        "voxel_length_z": float(cell[2]),
        "schema_0": {
            "dim_x": int(dims[0]),
            "dim_y": int(dims[1]),
            "dim_z": int(dims[2]),
            "lower_x": float(origin_arr[0]),
            "lower_y": float(origin_arr[1]),
            "lower_z": float(origin_arr[2]),
            "upper_x": float(origin_arr[0] + extents[0]),
            "upper_y": float(origin_arr[1] + extents[1]),
            "upper_z": float(origin_arr[2] + extents[2]),
            "max_size_x": float(cell[0]),
            "max_size_y": float(cell[1]),
            "max_size_z": float(cell[2]),
        },
    }
    for field in field_entries:
        config_object[str(field["entry_key"])] = field["entry"]

    config_bytes = _tbms_config_text(config_object).encode("latin1", errors="replace")
    config_offset = _align_offset(first_page_offset + page_count * page_size, 8)
    file_size = config_offset + len(config_bytes)
    with out_bmf.open("r+b") as fh:
        fh.seek(0, os.SEEK_END)
        if fh.tell() < config_offset:
            fh.write(b"\x00" * (config_offset - fh.tell()))
        _emit_progress(progress_callback, progress_end, 100, "Writing streaming BMF configuration...")
        fh.write(config_bytes)
        fh.seek(0)
        header = bytearray(TBMS_HEADER_SIZE)
        header[: len(TBMS_SIGNATURE)] = TBMS_SIGNATURE
        struct.pack_into("<I", header, 12, 1)
        struct.pack_into("<I", header, 16, first_page_offset)
        for header_offset in config_pointer_header_offsets:
            struct.pack_into("<Q", header, header_offset, config_offset)
        struct.pack_into("<Q", header, 48, file_size)
        fh.write(header)

    return {
        "backend": "tbms-config-text",
        "streaming": True,
        "file_size": int(file_size),
        "page_count": int(page_count),
        "field_count": len(field_entries),
        "categorical_fields": categorical_fields,
        "config_offset": int(config_offset),
        "first_page_offset": int(first_page_offset),
        "row_count": int(row_count),
        "value_columns": list(value_cols),
        "xyz_columns": list(xyz_cols),
    }


def _write_tbms_config_text(
    prepared: Dict[str, object],
    xyz_cols: Sequence[str],
    value_cols: Sequence[str],
    out_bmf: Path,
    null_float: float,
    column_types: Mapping[str, object] | None = None,
    progress_callback=None,
    progress_start: int = 50,
    progress_end: int = 98,
) -> Dict[str, object]:
    is_irregular = bool(prepared.get("is_irregular"))
    if not is_irregular and int(prepared.get("duplicates") or 0) > 0:
        raise ValueError("tbms-config-text export requires unique grid cells; duplicate XYZ rows were found.")

    df: pd.DataFrame = prepared["df"]  # type: ignore[assignment]
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    origin = np.asarray(prepared["origin"], dtype=float)
    cell = np.asarray(prepared["cell"], dtype=float)
    extents = np.asarray(prepared["extents"], dtype=float)
    row_count = int(len(df)) if is_irregular else int(np.prod(dims, dtype=np.int64))

    first_page_offset = TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET
    config_pointer_header_offsets = (24, 40)
    page_size = 2056

    field_entries: List[Dict[str, object]] = []
    categorical_fields = 0
    normalized_overrides = _normalize_export_type_overrides(column_types)
    total_fields = max(len(value_cols), 1)
    encode_end = _scale_progress(progress_start, progress_end, 0.82)
    _emit_progress(
        progress_callback,
        progress_start,
        100,
        f"Encoding BMF fields from CSV data (0/{len(value_cols)})...",
    )

    page_count = 0

    out_bmf.parent.mkdir(parents=True, exist_ok=True)
    with out_bmf.open("wb") as fh:
        _emit_progress(progress_callback, progress_start, 100, "Writing BMF header...")
        header = bytearray(TBMS_HEADER_SIZE)
        header[: len(TBMS_SIGNATURE)] = TBMS_SIGNATURE
        struct.pack_into("<I", header, 12, 1)
        struct.pack_into("<I", header, 16, first_page_offset)
        fh.write(header)
        if fh.tell() < first_page_offset:
            fh.write(b"\x00" * (first_page_offset - fh.tell()))

        def allocate_page(blob: bytes) -> int:
            nonlocal page_count
            if len(blob) != page_size:
                raise ValueError("TBMS page blobs must be exactly 2056 bytes.")
            offset = first_page_offset + page_count * page_size
            if fh.tell() < offset:
                fh.write(b"\x00" * (offset - fh.tell()))
            fh.write(blob)
            page_count += 1
            return offset

        geometry_entries: List[Tuple[str, pd.Series]] = []
        if is_irregular:
            lower = np.asarray(prepared["irregular_lower"], dtype=float)
            upper = np.asarray(prepared["irregular_upper"], dtype=float)
            geometry_entries = [
                ("__lower_x", pd.Series(lower[:, 0], name="__lower_x")),
                ("__upper_x", pd.Series(upper[:, 0], name="__upper_x")),
                ("__lower_y", pd.Series(lower[:, 1], name="__lower_y")),
                ("__upper_y", pd.Series(upper[:, 1], name="__upper_y")),
                ("__lower_z", pd.Series(lower[:, 2], name="__lower_z")),
                ("__upper_z", pd.Series(upper[:, 2], name="__upper_z")),
            ]

        for geometry_index, (name, values) in enumerate(geometry_entries):
            encoded = {
                "field_type": "float",
                "default": 0.0,
                "dtype": np.dtype("<f4"),
                "values": np.asarray(values, dtype="<f4"),
            }
            leaf_offsets = [
                allocate_page(blob)
                for blob in _iter_tbms_value_pages(encoded["values"], encoded["dtype"], encoded["default"])
            ]
            current_offsets = leaf_offsets
            while len(current_offsets) > 1:
                next_offsets: List[int] = []
                for start in range(0, len(current_offsets), 256):
                    next_offsets.append(allocate_page(_tbms_pointer_page(current_offsets[start:start + 256])))
                current_offsets = next_offsets
            field_entries.append({
                "entry_key": f"var_{geometry_index}",
                "entry": {
                    "name": name,
                    "type": encoded["field_type"],
                    "location": int(current_offsets[0]),
                    "default": encoded["default"],
                    "global": 0,
                    "read_only": 1,
                    "description": "Irregular block geometry extent",
                },
            })

        value_entry_offset = len(field_entries)
        for field_index, name in enumerate(value_cols):
            field_progress = _scale_progress(progress_start, encode_end, field_index / total_fields)
            _emit_progress(
                progress_callback,
                field_progress,
                100,
                f"Encoding BMF field {field_index + 1}/{len(value_cols)}: {name}...",
            )
            encoded = _tbms_encode_export_field(prepared, df[name], null_float, forced_type=normalized_overrides.get(name))
            leaf_offsets = [
                allocate_page(blob)
                for blob in _iter_tbms_value_pages(encoded["values"], encoded["dtype"], encoded["default"])
            ]
            current_offsets = leaf_offsets
            while len(current_offsets) > 1:
                next_offsets: List[int] = []
                for start in range(0, len(current_offsets), 256):
                    next_offsets.append(allocate_page(_tbms_pointer_page(current_offsets[start:start + 256])))
                current_offsets = next_offsets
            location = int(current_offsets[0])

            entry = {
                "name": name,
                "type": encoded["field_type"],
                "location": location,
                "default": encoded["default"],
                "global": 0,
                "read_only": 0,
                "description": f"Exported from CSV column {name}",
            }
            for key, value in encoded.items():
                if str(key).startswith("string_"):
                    entry[str(key)] = value
            if encoded["field_type"] == "namedshort":
                categorical_fields += 1
            field_entries.append({"entry_key": f"var_{value_entry_offset + field_index}", "entry": entry})

        _emit_progress(
            progress_callback,
            encode_end,
            100,
            f"Encoded {len(field_entries)} BMF field(s); preparing file write...",
        )

    config_object: Dict[str, object] = {
        "created": datetime.now(timezone.utc).isoformat(),
        "modified": datetime.now(timezone.utc).isoformat(),
        "history_source": "Anterpolator tbms-config-text export",
        "n_blocks": row_count,
        "n_schemas": 1,
        "is_irregular": 1 if is_irregular else 0,
        "origin_x": 0.0,
        "origin_y": 0.0,
        "origin_z": 0.0,
        "lower_x": float(origin[0]),
        "lower_y": float(origin[1]),
        "lower_z": float(origin[2]),
        "upper_x": float(origin[0] + extents[0]),
        "upper_y": float(origin[1] + extents[1]),
        "upper_z": float(origin[2] + extents[2]),
        "voxel_length_x": float(cell[0]),
        "voxel_length_y": float(cell[1]),
        "voxel_length_z": float(cell[2]),
        "schema_0": {
            "dim_x": int(dims[0]),
            "dim_y": int(dims[1]),
            "dim_z": int(dims[2]),
            "lower_x": float(origin[0]),
            "lower_y": float(origin[1]),
            "lower_z": float(origin[2]),
            "upper_x": float(origin[0] + extents[0]),
            "upper_y": float(origin[1] + extents[1]),
            "upper_z": float(origin[2] + extents[2]),
            "max_size_x": float(cell[0]),
            "max_size_y": float(cell[1]),
            "max_size_z": float(cell[2]),
        },
    }
    if is_irregular:
        config_object["geometry_encoding"] = "row_extents_lower_upper"
    for field in field_entries:
        config_object[str(field["entry_key"])] = field["entry"]

    config_bytes = _tbms_config_text(config_object).encode("latin1", errors="replace")
    config_offset = _align_offset(first_page_offset + page_count * page_size, 8)
    file_size = config_offset + len(config_bytes)
    with out_bmf.open("r+b") as fh:
        fh.seek(0, os.SEEK_END)
        if fh.tell() < config_offset:
            fh.write(b"\x00" * (config_offset - fh.tell()))
        _emit_progress(progress_callback, progress_end, 100, "Writing BMF configuration...")
        fh.write(config_bytes)

        fh.seek(0)
        header = bytearray(TBMS_HEADER_SIZE)
        header[: len(TBMS_SIGNATURE)] = TBMS_SIGNATURE
        struct.pack_into("<I", header, 12, 1)
        struct.pack_into("<I", header, 16, first_page_offset)
        for header_offset in config_pointer_header_offsets:
            struct.pack_into("<Q", header, header_offset, config_offset)
        struct.pack_into("<Q", header, 48, file_size)
        fh.write(header)

    return {
        "backend": "tbms-config-text",
        "file_size": int(file_size),
        "page_count": int(page_count),
        "field_count": len(field_entries),
        "categorical_fields": categorical_fields,
        "config_offset": int(config_offset),
        "first_page_offset": int(first_page_offset),
        "row_count": int(row_count),
        "value_columns": list(value_cols),
        "xyz_columns": list(xyz_cols),
    }


def _write_tbms_experimental(
    prepared: Dict[str, object],
    xyz_cols: Sequence[str],
    value_cols: Sequence[str],
    out_bmf: Path,
    null_float: float,
    column_types: Mapping[str, object] | None = None,
    progress_callback=None,
    progress_start: int = 50,
    progress_end: int = 98,
) -> Dict[str, object]:
    df: pd.DataFrame = prepared["df"]  # type: ignore[assignment]
    idx = np.asarray(prepared["idx"], dtype=np.int32)
    origin = np.asarray(prepared["origin"], dtype=float)
    cell = np.asarray(prepared["cell"], dtype=float)
    dims = np.asarray(prepared["dims"], dtype=np.int32)
    extents = np.asarray(prepared["extents"], dtype=float)

    normalized_overrides = _normalize_export_type_overrides(column_types)
    _emit_progress(progress_callback, progress_start, 100, "Preparing experimental BMF sections...")
    numeric_cols = [
        col for col in value_cols
        if normalized_overrides.get(col) in {"boolean", "int", "double"}
        or (normalized_overrides.get(col) is None and pd.api.types.is_numeric_dtype(df[col]))
    ]
    string_cols = [col for col in value_cols if col not in numeric_cols]

    numeric_matrix = np.empty((len(df), len(numeric_cols)), dtype="<f8")
    if numeric_cols:
        numeric_matrix = df[numeric_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype="<f8", copy=True)

    column_spec = []
    for col in numeric_cols:
        column_spec.append({"name": col, "storage": "float64", "null_value": float(null_float), "kind": "numeric"})
    for col in string_cols:
        column_spec.append({"name": col, "storage": "json_utf8", "kind": "string"})

    now_utc = datetime.now(timezone.utc).isoformat()
    metadata = {
        "created": now_utc,
        "modified": now_utc,
        "description": "Experimental Anterpolator TBMS2.0 container. Not confirmed Vulcan-compatible.",
        "format_family": "TBMS2.0",
        "experimental_backend": "anterpolator.tbms-experimental.v1",
        "row_count": int(len(df)),
        "duplicate_grid_rows": int(prepared["duplicates"]),
        "max_index_error": float(prepared["max_index_error"]),
        "x_col": xyz_cols[0],
        "y_col": xyz_cols[1],
        "z_col": xyz_cols[2],
        "lower_x": float(origin[0]),
        "lower_y": float(origin[1]),
        "lower_z": float(origin[2]),
        "upper_x": float(origin[0] + extents[0]),
        "upper_y": float(origin[1] + extents[1]),
        "upper_z": float(origin[2] + extents[2]),
        "origin": origin.tolist(),
        "cell_size": cell.tolist(),
        "dimensions": dims.tolist(),
        "extents": extents.tolist(),
        "value_columns": list(value_cols),
        "numeric_columns": numeric_cols,
        "string_columns": string_cols,
        "null_float": float(null_float),
    }

    sections: List[Tuple[str, int, bytes]] = [
        ("metadata", TBMS_EXPERIMENTAL_SECTION_FLAGS["json_utf8"], _json_bytes(metadata)),
        ("columns", TBMS_EXPERIMENTAL_SECTION_FLAGS["json_utf8"], _json_bytes(column_spec)),
        ("indices", TBMS_EXPERIMENTAL_SECTION_FLAGS["int32_le"], idx.astype("<i4", copy=False).tobytes(order="C")),
    ]
    if numeric_cols:
        sections.append(("numeric", TBMS_EXPERIMENTAL_SECTION_FLAGS["float64_le"], numeric_matrix.tobytes(order="C")))
    if string_cols:
        sections.append(("strings", TBMS_EXPERIMENTAL_SECTION_FLAGS["json_utf8"], _json_bytes(_normalize_string_records(df, string_cols))))

    dir_size = 16 + len(sections) * TBMS_EXPERIMENTAL_SECTION_ENTRY.size
    if TBMS_HEADER_SIZE + dir_size > TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET:
        raise ValueError("Experimental section directory does not fit before first payload offset.")

    section_infos = []
    current_offset = TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET
    for name, flags, payload in sections:
        current_offset = _align_offset(current_offset, 8)
        section_infos.append({
            "name": name,
            "flags": int(flags),
            "offset": int(current_offset),
            "size": int(len(payload)),
            "payload": payload,
        })
        current_offset += len(payload)

    file_size = current_offset
    out_bmf.parent.mkdir(parents=True, exist_ok=True)
    with out_bmf.open("wb") as fh:
        _emit_progress(progress_callback, _scale_progress(progress_start, progress_end, 0.1), 100, "Writing experimental BMF header...")
        header = bytearray(TBMS_HEADER_SIZE)
        header[: len(TBMS_SIGNATURE)] = TBMS_SIGNATURE
        struct.pack_into("<I", header, 12, 1)
        struct.pack_into("<I", header, 16, TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET)
        struct.pack_into("<Q", header, 48, file_size)
        fh.write(header)

        fh.write(TBMS_EXPERIMENTAL_DIR_MAGIC)
        fh.write(struct.pack("<II", 1, len(section_infos)))
        for info in section_infos:
            name_bytes = info["name"].encode("ascii", errors="ignore")[:16].ljust(16, b"\x00")
            fh.write(TBMS_EXPERIMENTAL_SECTION_ENTRY.pack(
                name_bytes,
                int(info["offset"]),
                int(info["size"]),
                int(info["flags"]),
                0,
            ))

        if fh.tell() > TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET:
            raise ValueError("Experimental TBMS directory overflowed reserved payload boundary.")
        fh.write(b"\x00" * (TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET - fh.tell()))

        total_sections = max(len(section_infos), 1)
        for section_index, info in enumerate(section_infos, start=1):
            if fh.tell() < info["offset"]:
                fh.write(b"\x00" * (info["offset"] - fh.tell()))
            fh.write(info["payload"])
            _emit_progress(
                progress_callback,
                _scale_progress(progress_start, progress_end, section_index / total_sections),
                100,
                f"Writing experimental BMF section {section_index}/{total_sections}: {info['name']}...",
            )

    return {
        "backend": "tbms-experimental",
        "file_size": int(file_size),
        "sections": [{k: v for k, v in info.items() if k != "payload"} for info in section_infos],
    }


def _parse_tbms_experimental_sections(path: Path) -> List[Dict[str, object]] | None:
    with path.open("rb") as fh:
        fh.seek(TBMS_HEADER_SIZE)
        prefix = fh.read(16)
        if len(prefix) < 16 or prefix[:8] != TBMS_EXPERIMENTAL_DIR_MAGIC:
            return None
        version, count = struct.unpack("<II", prefix[8:16])
        if version != 1 or count < 0 or count > 64:
            return None

        sections: List[Dict[str, object]] = []
        for _ in range(count):
            raw = fh.read(TBMS_EXPERIMENTAL_SECTION_ENTRY.size)
            if len(raw) != TBMS_EXPERIMENTAL_SECTION_ENTRY.size:
                return None
            name_bytes, offset, size, flags, _reserved = TBMS_EXPERIMENTAL_SECTION_ENTRY.unpack(raw)
            name = name_bytes.split(b"\x00", 1)[0].decode("ascii", errors="ignore")
            sections.append({
                "name": name,
                "offset": int(offset),
                "size": int(size),
                "flags": int(flags),
            })

        for section in sections:
            if section["flags"] != TBMS_EXPERIMENTAL_SECTION_FLAGS["json_utf8"]:
                continue
            if section["size"] > 1024 * 1024:
                continue
            fh.seek(int(section["offset"]))
            payload = fh.read(int(section["size"]))
            try:
                section["json_preview"] = json.loads(payload.decode("utf-8"))
            except Exception:
                section["json_preview"] = None
        return sections


def _read_section_bytes(path: Path, section: Dict[str, object], limit_bytes: int | None = None) -> bytes:
    read_size = int(section["size"])
    if limit_bytes is not None:
        read_size = min(read_size, int(limit_bytes))
    with path.open("rb") as fh:
        fh.seek(int(section["offset"]))
        return fh.read(read_size)


class _TbmsConfigParser:
    def __init__(self, text: str):
        self.text = text
        self.length = len(text)
        self.pos = 0

    def parse(self) -> Dict[str, object]:
        self._skip_ws()
        value = self._parse_object()
        self._skip_ws()
        return value

    def _skip_ws(self) -> None:
        while self.pos < self.length and ord(self.text[self.pos]) <= 32:
            self.pos += 1

    def _peek(self) -> str:
        if self.pos >= self.length:
            return ""
        return self.text[self.pos]

    def _expect(self, token: str) -> None:
        self._skip_ws()
        if not self.text.startswith(token, self.pos):
            raise ValueError(f"Expected {token!r} at position {self.pos}")
        self.pos += len(token)

    def _skip_to_object_key_start(self) -> None:
        self._skip_ws()
        while self.pos < self.length:
            ch = self._peek()
            if ch in ('"', "}"):
                return
            self.pos += 1
            self._skip_ws()

    def _parse_object(self) -> Dict[str, object]:
        self._expect("{")
        result: Dict[str, object] = {}
        while True:
            self._skip_to_object_key_start()
            if self._peek() == "}":
                self.pos += 1
                return result
            key = self._parse_string()
            self._skip_ws()
            self._expect("=")
            value = self._parse_value()
            result[key] = value
            self._skip_ws()
            if self._peek() == ",":
                self.pos += 1
                continue
            if self._peek() == "}":
                self.pos += 1
                return result
            raise ValueError(f"Unexpected character {self._peek()!r} at position {self.pos}")

    def _parse_string(self) -> str:
        self._skip_ws()
        if self._peek() != '"':
            raise ValueError(f"Expected string at position {self.pos}")
        self.pos += 1
        chars: List[str] = []
        while self.pos < self.length:
            ch = self.text[self.pos]
            self.pos += 1
            if ch == '"':
                return "".join(ch2 for ch2 in chars if ord(ch2) >= 32)
            if ch == "\\" and self.pos < self.length:
                esc = self.text[self.pos]
                self.pos += 1
                escape_map = {"n": "\n", "r": "\r", "t": "\t", '"': '"', "\\": "\\"}
                chars.append(escape_map.get(esc, esc))
            else:
                chars.append(ch)
        raise ValueError("Unterminated string in TBMS config block")

    def _parse_number(self) -> object:
        start = self.pos
        while self.pos < self.length and self.text[self.pos] in "+-0123456789.eE":
            self.pos += 1
        token = self.text[start:self.pos]
        if token == "":
            raise ValueError(f"Expected number at position {start}")
        if any(ch in token for ch in ".eE"):
            return float(token)
        return int(token)

    def _parse_value(self) -> object:
        self._skip_ws()
        ch = self._peek()
        if ch == "{":
            return self._parse_object()
        if ch == '"':
            return self._parse_string()
        if ch in "+-0123456789":
            return self._parse_number()
        if self.text.startswith("true", self.pos):
            self.pos += 4
            return True
        if self.text.startswith("false", self.pos):
            self.pos += 5
            return False
        raise ValueError(f"Unsupported value at position {self.pos}")


def _find_bytes(path: Path, needle: bytes, start: int = 0, chunk_size: int = 8 * 1024 * 1024) -> int:
    overlap = max(1, len(needle) - 1)
    position = start
    previous = b""
    with path.open("rb") as fh:
        fh.seek(start)
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                return -1
            buffer = previous + chunk
            idx = buffer.find(needle)
            if idx >= 0:
                return position - len(previous) + idx
            position += len(chunk)
            previous = buffer[-overlap:]


def _extract_balanced_brace_block(path: Path, start: int, max_scan: int) -> Dict[str, object] | None:
    with path.open("rb") as fh:
        fh.seek(start)
        data = fh.read(max_scan)

    if not data or data[0] != 0x7B:
        return None

    depth = 0
    in_string = False
    escape = False
    end = None
    for offset in range(0, min(len(data), max_scan)):
        byte = data[offset]
        if in_string:
            if escape:
                escape = False
            elif byte == 0x5C:
                escape = True
            elif byte == 0x22:
                in_string = False
            continue

        if byte == 0x22:
            in_string = True
        elif byte == 0x7B:
            depth += 1
        elif byte == 0x7D:
            depth -= 1
            if depth == 0:
                end = offset
                break

    if end is None:
        return None

    text = data[:end + 1].decode("latin1", errors="replace")
    parsed = _TbmsConfigParser(text).parse()
    if not _looks_like_tbms_config(parsed):
        return None
    return {
        "start": int(start),
        "end": int(start + end),
        "text": text,
        "parsed": parsed,
    }


def _looks_like_tbms_config(parsed: object) -> bool:
    if not isinstance(parsed, dict):
        return False
    keys = {str(key) for key in parsed.keys()}
    if "created" in keys or "n_blocks" in keys:
        return True
    return any(
        key.startswith(("schema_", "var_", "special_"))
        for key in keys
    )


def _score_tbms_config_result(result: Dict[str, object] | None) -> int:
    if result is None:
        return -1
    parsed = result.get("parsed")
    if not isinstance(parsed, dict):
        return -1
    score = 0
    if "created" in parsed:
        score += 10
    if "n_blocks" in parsed:
        score += 10
    try:
        fields = _tbms_schema_fields(parsed)
        score += len(fields) * 100
        if "location" in fields.columns:
            locations = pd.to_numeric(fields["location"], errors="coerce").fillna(0)
            score += int((locations > 0).sum()) * 10
    except Exception:
        pass
    return score


def _extract_tbms_segmented_config_block(path: Path, directory_offset: int) -> Dict[str, object] | None:
    if directory_offset <= 0 or directory_offset >= path.stat().st_size:
        return None
    table = _tbms_read_qword_table(path, directory_offset)
    if table is None or table[0] not in (131329,):
        return None
    pointers = [int(pointer) for pointer in table[1:] if pointer]
    if not pointers:
        return None

    file_size = path.stat().st_size
    if any(pointer <= 0 or pointer >= file_size for pointer in pointers):
        return None

    chunks: List[bytes] = []
    with path.open("rb") as fh:
        for pointer in pointers:
            fh.seek(pointer + 8)
            chunks.append(fh.read(2048))
    raw = b"".join(chunks)
    brace_index = raw.find(b"{")
    if brace_index < 0:
        return None
    raw = raw[brace_index:]
    nul_index = raw.find(b"\x00")
    if nul_index >= 0:
        raw = raw[:nul_index]
    text = raw.decode("latin1", errors="replace")
    parsed = _TbmsConfigParser(text).parse()
    if not _looks_like_tbms_config(parsed):
        return None
    return {
        "start": int(pointers[0] + 8 + brace_index),
        "end": int(pointers[-1] + 8 + 2048),
        "text": text,
        "parsed": parsed,
        "segments": [int(pointer) for pointer in pointers],
        "directory_offset": int(directory_offset),
    }


def _extract_tbms_first_page_assignments(
    path: Path,
    first_page_offset: int = 2056,
    max_pages: int = 8,
    min_matches: int = 5,
) -> Dict[str, object] | None:
    assignment_pattern = re.compile(r'"(?P<key>[^"\r\n]+)"\s*=\s*(?:"(?P<quoted>[^"]*)"|(?P<bare>[^,\r\n}]+))')
    assignments: List[Dict[str, object]] = []
    pages: List[Dict[str, object]] = []
    file_size = path.stat().st_size
    page_size = 2056
    if first_page_offset <= 0 or first_page_offset >= file_size:
        return None

    with path.open("rb") as fh:
        for index in range(max_pages):
            page_offset = first_page_offset + index * page_size
            if page_offset + page_size > file_size:
                break
            fh.seek(page_offset)
            raw = fh.read(page_size)
            if len(raw) != page_size:
                break
            payload = raw[8:]
            text = payload.decode("latin1", errors="replace")
            matches = list(assignment_pattern.finditer(text))
            if len(matches) < min_matches:
                if assignments:
                    break
                continue
            pages.append({
                "page_offset": int(page_offset),
                "page_header": list(struct.unpack("<2I", raw[:8])),
                "assignment_count": len(matches),
            })
            for match in matches:
                value = match.group("quoted")
                if value is None:
                    value = (match.group("bare") or "").strip()
                assignments.append({
                    "entry_type": "text_assignment",
                    "entry_key": match.group("key"),
                    "name": match.group("key"),
                    "value": value,
                    "page_offset": int(page_offset),
                })

    if not assignments:
        return None
    return {
        "first_page_offset": int(first_page_offset),
        "pages": pages,
        "assignments": assignments,
    }


def _extract_tbms_text_tree_config_block(path: Path, root_offset: int, max_depth: int = 4) -> Dict[str, object] | None:
    file_size = path.stat().st_size
    page_size = 2056
    text_leaf_magic = 134218240

    def gather_leaves(offset: int, depth: int = 0, visited: set[int] | None = None) -> List[int]:
        if visited is None:
            visited = set()
        if depth > max_depth or offset <= 0 or offset >= file_size or offset in visited:
            return []
        visited.add(offset)

        with path.open("rb") as fh:
            fh.seek(offset)
            raw = fh.read(page_size)
        if len(raw) != page_size:
            return []

        magic = struct.unpack("<I", raw[:4])[0]
        if magic == text_leaf_magic:
            return [offset]

        table = struct.unpack("<257Q", raw)
        leaves: List[int] = []
        for pointer in table[1:]:
            if pointer:
                leaves.extend(gather_leaves(int(pointer), depth + 1, visited))
        return leaves

    leaves = gather_leaves(root_offset)
    if not leaves:
        return None

    chunks: List[bytes] = []
    with path.open("rb") as fh:
        for leaf in leaves:
            fh.seek(leaf + 8)
            chunks.append(fh.read(2048))
    raw = b"".join(chunks)
    span_start = min(leaves) + 8
    span_end = max(leaves) + 2056
    brace_index = raw.find(b"{")
    if brace_index < 0:
        return None

    text = raw[brace_index:].decode("latin1", errors="replace")
    parse_error: Exception | None = None
    for candidate in (text,):
        try:
            parsed = _TbmsConfigParser(candidate).parse()
            if _looks_like_tbms_config(parsed):
                return {
                    "start": int(span_start),
                    "end": int(span_end),
                    "text": candidate,
                    "parsed": parsed,
                    "segments": [int(leaf) for leaf in leaves],
                    "directory_offset": int(root_offset),
                    "source": "text-tree",
                }
        except Exception as exc:
            parse_error = exc

    failure_pos = None
    match = re.search(r"position (\d+)", str(parse_error or ""))
    if match:
        failure_pos = int(match.group(1))

    cut_candidates: List[int] = []
    search_upto = failure_pos if failure_pos is not None else len(text)
    for marker, marker_len in (('\n  },', 5), ('\n }', 3), ('\n  }', 4), ('},', 2)):
        cut = text.rfind(marker, 0, search_upto)
        if cut >= 0:
            cut_candidates.append(cut + marker_len)
    if not cut_candidates:
        return None

    trimmed = text[:max(cut_candidates)]
    if not trimmed.rstrip().endswith('}'):
        trimmed = trimmed.rstrip(',\r\n ') + '\n}\n'
    try:
        parsed = _TbmsConfigParser(trimmed).parse()
    except Exception:
        return None
    if not _looks_like_tbms_config(parsed):
        return None
    return {
        "start": int(span_start),
        "end": int(span_end),
        "text": trimmed,
        "parsed": parsed,
        "segments": [int(leaf) for leaf in leaves],
        "directory_offset": int(root_offset),
        "source": "text-tree-trimmed",
    }


def _find_tbms_config_block(path: Path, search_start: int = 0x5000, max_scan: int = 2 * 1024 * 1024) -> Dict[str, object] | None:
    candidates: List[Dict[str, object]] = []
    with path.open("rb") as fh:
        header = fh.read(64)
    if len(header) >= 56 and header.startswith(TBMS_SIGNATURE):
        header_offsets = struct.unpack("<8Q", header)
        for pointer in (header_offsets[3], header_offsets[5]):
            if pointer <= 0 or pointer >= path.stat().st_size:
                continue
            try:
                result = _extract_tbms_text_tree_config_block(path, int(pointer))
            except Exception:
                result = None
            if result is not None:
                candidates.append(result)

            try:
                result = _extract_tbms_segmented_config_block(path, int(pointer))
            except Exception:
                result = None
            if result is not None:
                candidates.append(result)

            for delta in (0, 8):
                candidate_start = int(pointer) + delta
                try:
                    with path.open("rb") as fh:
                        fh.seek(candidate_start)
                        first = fh.read(1)
                    if first != b"{":
                        continue
                    result = _extract_balanced_brace_block(path, candidate_start, max_scan=max_scan)
                except Exception:
                    result = None
                if result is not None:
                    candidates.append(result)

        if candidates:
            return max(candidates, key=_score_tbms_config_result)

    anchors = [b"created", b"n_blocks", b"origin_x", b"lower_x", b"schema_0"]
    back_window = 256 * 1024

    for anchor in anchors:
        anchor_offset = _find_bytes(path, anchor, start=search_start)
        if anchor_offset < 0:
            continue

        with path.open("rb") as fh:
            back_start = max(search_start, anchor_offset - back_window)
            fh.seek(back_start)
            back_data = fh.read(anchor_offset - back_start + len(anchor))

        brace_positions = [
            back_start + idx
            for idx, byte in enumerate(back_data)
            if byte == 0x7B
        ]
        for candidate_start in brace_positions:
            try:
                result = _extract_balanced_brace_block(path, candidate_start, max_scan=max_scan)
            except Exception:
                result = None
            if result is not None:
                return result

    first_brace = _find_bytes(path, b"{", start=search_start)
    if first_brace < 0:
        return None
    try:
        return _extract_balanced_brace_block(path, first_brace, max_scan=max_scan)
    except Exception:
        return None


def _sanitize_tbms_key(key: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z_.-]+", "", str(key))
    return cleaned or str(key)


def _sanitize_tbms_value(value: object) -> object:
    if isinstance(value, dict):
        return {_sanitize_tbms_key(k): _sanitize_tbms_value(v) for k, v in value.items()}
    return value


def _flatten_tbms_config(data: Dict[str, object], prefix: str = "") -> Dict[str, object]:
    flat: Dict[str, object] = {}
    for key, value in data.items():
        key = _sanitize_tbms_key(key)
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            flat.update(_flatten_tbms_config(value, name))
        else:
            flat[name] = value
    return flat


def _tbms_schema_fields(parsed: Dict[str, object]) -> pd.DataFrame:
    rows = []
    for key, value in parsed.items():
        key = _sanitize_tbms_key(key)
        if not isinstance(value, dict):
            continue
        value = { _sanitize_tbms_key(sub_key): sub_value for sub_key, sub_value in value.items() }
        if key.startswith("schema_"):
            row = {"entry_key": key, "entry_type": "schema"}
            rows.append(row)
            row.update(value)
            continue
        if "name" in value and ("type" in value or "location" in value or "description" in value):
            row = {"entry_key": key, "entry_type": "field"}
            row.update(value)
            rows.append(row)
    if rows:
        frame = pd.DataFrame(rows)
        preferred = [col for col in ["entry_type", "entry_key", "name", "type", "location", "default", "global", "read_only", "description"] if col in frame.columns]
        remaining = [col for col in frame.columns if col not in preferred]
        return frame[preferred + remaining]
    flat_rows = [{"name": key, "value": value} for key, value in parsed.items()]
    return pd.DataFrame(flat_rows)


def _tbms_numeric_dtype(field_type: object) -> np.dtype | None:
    mapping = {
        "float": np.dtype("<f4"),
        "double": np.dtype("<f8"),
        "int": np.dtype("<i4"),
        "integer32": np.dtype("<i4"),
        "short": np.dtype("<i2"),
        "namedshort": np.dtype("<i2"),
        "longlong": np.dtype("<i8"),
        "integer64": np.dtype("<i8"),
        "byte": np.dtype("u1"),
        "boolean": np.dtype("u1"),
    }
    return mapping.get(str(field_type or "").strip().lower())


def _tbms_read_qword_table(path: Path, offset: int) -> Tuple[int, ...] | None:
    if offset < 0:
        return None
    with path.open("rb") as fh:
        fh.seek(offset)
        raw = fh.read(2056)
    if len(raw) != 2056:
        return None
    return struct.unpack("<257Q", raw)


def _tbms_looks_like_value_page(path: Path, offset: int) -> bool:
    if offset <= 0 or offset >= path.stat().st_size:
        return False
    with path.open("rb") as fh:
        fh.seek(offset)
        raw = fh.read(8)
    if len(raw) != 8:
        return False
    payload_words, reserved = struct.unpack("<II", raw)
    return reserved == 0 and payload_words in (256, 512, 768, 1024, 1536, 2048, 2560)


def _iter_tbms_page_tree(path: Path, offset: int, depth: int = 0, visited: set[int] | None = None) -> Iterable[int]:
    if visited is None:
        visited = set()
    if depth > 5 or offset <= 0 or offset in visited or offset >= path.stat().st_size:
        return
    visited.add(offset)
    if _tbms_looks_like_value_page(path, offset):
        yield offset
        return
    table = _tbms_read_qword_table(path, offset)
    if table is None or table[0] not in (257, 258, 259, 260):
        return
    for pointer in table[1:]:
        if pointer:
            yield from _iter_tbms_page_tree(path, int(pointer), depth + 1, visited)


def _iter_tbms_field_pages(path: Path, location: int) -> Iterable[int]:
    yield from _iter_tbms_page_tree(path, location)


def _decode_tbms_field_values(
    path: Path,
    location: int,
    field_type: object,
    row_limit: int,
) -> np.ndarray | None:
    dtype = _tbms_numeric_dtype(field_type)
    if dtype is None or row_limit <= 0 or location <= 0:
        return None

    payload_size = 2048
    page_capacity = payload_size // dtype.itemsize
    parts: List[np.ndarray] = []
    loaded = 0

    with path.open("rb") as fh:
        for page_offset in _iter_tbms_field_pages(path, location):
            fh.seek(page_offset + 8)
            raw = fh.read(payload_size)
            if len(raw) != payload_size:
                break
            page_values = np.frombuffer(raw, dtype=dtype, count=page_capacity).copy()
            if dtype == np.dtype("u1") and str(field_type or "").strip().lower() == "boolean":
                page_values = page_values.astype(bool)
            parts.append(page_values)
            loaded += len(page_values)
            if loaded >= row_limit:
                break

    if not parts:
        return None
    values = np.concatenate(parts)
    return values[:row_limit]


def _tbms_grid_spec(metadata: Dict[str, object]) -> Dict[str, np.ndarray] | None:
    dims = np.asarray([
        int(metadata.get("schema_0.dim_x", metadata.get("dim_x", metadata.get("voxel_dim_x", 0))) or 0),
        int(metadata.get("schema_0.dim_y", metadata.get("dim_y", metadata.get("voxel_dim_y", 0))) or 0),
        int(metadata.get("schema_0.dim_z", metadata.get("dim_z", metadata.get("voxel_dim_z", 0))) or 0),
    ], dtype=np.int64)

    fallback_cell = []
    for axis in ("x", "y", "z"):
        upper = float(metadata.get(f"upper_{axis}", metadata.get(f"schema_0.upper_{axis}", 0.0)) or 0.0)
        lower_axis = float(metadata.get(f"lower_{axis}", metadata.get(f"schema_0.lower_{axis}", 0.0)) or 0.0)
        dim_axis = dims[{"x": 0, "y": 1, "z": 2}[axis]]
        fallback_cell.append((upper - lower_axis) / dim_axis if dim_axis > 0 else 0.0)

    cell = np.asarray([
        float(metadata.get("voxel_length_x", metadata.get("schema_0.max_size_x", fallback_cell[0])) or fallback_cell[0]),
        float(metadata.get("voxel_length_y", metadata.get("schema_0.max_size_y", fallback_cell[1])) or fallback_cell[1]),
        float(metadata.get("voxel_length_z", metadata.get("schema_0.max_size_z", fallback_cell[2])) or fallback_cell[2]),
    ], dtype=float)
    lower = np.asarray([
        float(metadata.get("schema_0.lower_x", metadata.get("lower_x", metadata.get("voxel_lower_x", 0.0))) or 0.0),
        float(metadata.get("schema_0.lower_y", metadata.get("lower_y", metadata.get("voxel_lower_y", 0.0))) or 0.0),
        float(metadata.get("schema_0.lower_z", metadata.get("lower_z", metadata.get("voxel_lower_z", 0.0))) or 0.0),
    ], dtype=float)
    origin = np.asarray([
        float(metadata.get("origin_x", 0.0) or 0.0),
        float(metadata.get("origin_y", 0.0) or 0.0),
        float(metadata.get("origin_z", 0.0) or 0.0),
    ], dtype=float)
    if np.any(dims <= 0) or np.any(cell <= 0):
        return None
    return {
        "dims": dims,
        "cell": cell,
        "lower": lower,
        "origin": origin,
    }


def _build_tbms_grid_frame(metadata: Dict[str, object], rows_to_load: int) -> pd.DataFrame:
    spec = _tbms_grid_spec(metadata)
    if spec is None or rows_to_load <= 0:
        return pd.DataFrame()

    dims = spec["dims"]
    cell = spec["cell"]
    lower = spec["lower"]
    origin = spec["origin"]

    idx = np.arange(rows_to_load, dtype=np.int64)
    plane = int(dims[0] * dims[1])
    grid_i = idx % dims[0]
    grid_j = (idx // dims[0]) % dims[1]
    grid_k = idx // plane

    local_lower_x = lower[0] + grid_i.astype(float) * cell[0]
    local_lower_y = lower[1] + grid_j.astype(float) * cell[1]
    local_lower_z = lower[2] + grid_k.astype(float) * cell[2]

    frame = pd.DataFrame({
        "grid_i": grid_i,
        "grid_j": grid_j,
        "grid_k": grid_k,
        "x": origin[0] + local_lower_x + 0.5 * cell[0],
        "y": origin[1] + local_lower_y + 0.5 * cell[1],
        "z": origin[2] + local_lower_z + 0.5 * cell[2],
    })
    return frame


def _probe_tbms_first_page(path: Path, first_offset: int | None = None) -> Dict[str, object]:
    page_offset = int(first_offset or TBMS_EXPERIMENTAL_FIRST_SECTION_OFFSET)
    with path.open("rb") as fh:
        fh.seek(page_offset)
        raw = fh.read(2056)
    if len(raw) < 8:
        return {
            "page_offset": page_offset,
            "available": False,
        }

    header_words = struct.unpack("<2I", raw[:8])
    return {
        "page_offset": page_offset,
        "available": True,
        "payload_words": int(header_words[0]),
        "reserved_word": int(header_words[1]),
        "sample_u32": list(struct.unpack("<8I", raw[:32])),
    }


def _read_tbms_value_page(path: Path, offset: int, dtype: np.dtype, row_limit: int) -> np.ndarray | None:
    if offset <= 0 or row_limit <= 0:
        return None
    payload_size = 2048
    count = payload_size // dtype.itemsize
    with path.open("rb") as fh:
        fh.seek(offset + 8)
        raw = fh.read(payload_size)
    if len(raw) != payload_size:
        return None
    return np.frombuffer(raw, dtype=dtype, count=count).copy()[:row_limit]


def _scan_tbms_flat_roots(
    path: Path,
    page_start: int = 2056,
    scan_bytes: int = 8 * 1024 * 1024,
    max_roots: int = 256,
) -> List[Dict[str, object]]:
    page_size = 2056
    roots: List[Dict[str, object]] = []
    file_size = path.stat().st_size
    if page_start <= 0 or page_start >= file_size:
        page_start = page_size
    scan_end = min(file_size, page_start + max(scan_bytes, page_size))
    with path.open("rb") as fh:
        for offset in range(page_start, scan_end - page_size + 1, page_size):
            fh.seek(offset)
            raw = fh.read(page_size)
            if len(raw) != page_size:
                break
            table = struct.unpack("<257Q", raw)
            if table[0] != 257:
                continue
            pointers = [int(value) for value in table[1:] if value]
            if not pointers:
                continue
            if any(
                pointer < page_start or pointer >= file_size or (pointer - page_start) % page_size != 0
                for pointer in pointers
            ):
                continue
            with path.open("rb") as check_fh:
                valid_pages = 0
                sample_header = None
                for pointer in pointers[:8]:
                    check_fh.seek(pointer)
                    page_header = check_fh.read(8)
                    if len(page_header) != 8:
                        continue
                    payload_words, reserved = struct.unpack("<II", page_header)
                    if reserved == 0 and payload_words in (256, 512, 768, 1024, 1536, 2048, 2560):
                        valid_pages += 1
                        sample_header = (payload_words, reserved)
                if valid_pages == 0:
                    continue
            roots.append({
                "root_offset": int(offset),
                "page_offsets": pointers,
                "sample_page_header": sample_header,
            })
            if len(roots) >= max_roots:
                break
    return roots


def _decode_tbms_flat_root_values(
    path: Path,
    page_offsets: Sequence[int],
    dtype: np.dtype,
    row_limit: int,
) -> np.ndarray | None:
    parts: List[np.ndarray] = []
    loaded = 0
    payload_size = 2048
    page_capacity = payload_size // dtype.itemsize
    with path.open("rb") as fh:
        for page_offset in page_offsets:
            fh.seek(int(page_offset) + 8)
            raw = fh.read(payload_size)
            if len(raw) != payload_size:
                break
            values = np.frombuffer(raw, dtype=dtype, count=page_capacity).copy()
            parts.append(values)
            loaded += len(values)
            if loaded >= row_limit:
                break
    if not parts:
        return None
    return np.concatenate(parts)[:row_limit]


def _infer_flat_root_dtype_and_name(index: int, float_values: np.ndarray | None, int_values: np.ndarray | None) -> Tuple[str, np.dtype]:
    standard = [
        ("__lower_x", np.dtype("<f4")),
        ("__upper_x", np.dtype("<f4")),
        ("__lower_y", np.dtype("<f4")),
        ("__upper_y", np.dtype("<f4")),
        ("__lower_z", np.dtype("<f4")),
        ("__upper_z", np.dtype("<f4")),
    ]
    if index < len(standard):
        return standard[index]

    if float_values is not None and int_values is not None:
        finite = float_values[np.isfinite(float_values)]
        abs_finite = np.abs(finite) if finite.size else np.asarray([], dtype=float)
        integer_like = np.asarray(int_values, dtype=np.int64)
        plausible_ints = integer_like[(integer_like >= -1_000_000_000) & (integer_like <= 1_000_000_000)]
        tiny_float_ratio = float(np.mean(abs_finite < 1e-20)) if abs_finite.size else 1.0
        plausible_float_ratio = float(np.mean((abs_finite > 1e-6) & (abs_finite < 1e9))) if abs_finite.size else 0.0
        plausible_int_ratio = float(len(plausible_ints) / len(integer_like)) if len(integer_like) else 0.0
        if tiny_float_ratio > 0.8 and plausible_int_ratio > 0.8:
            return f"field_{index:03d}", np.dtype("<i4")
        if plausible_float_ratio > 0.25:
            return f"field_{index:03d}", np.dtype("<f4")
    return f"field_{index:03d}", np.dtype("<f4")


def _load_tbms_headerless_flat_roots(
    path: Path,
    report: Dict[str, object],
    recognition: Dict[str, object],
    row_limit: int | None,
) -> Dict[str, object]:
    roots = recognition.get("flat_roots") or []
    rows_to_load = int(row_limit) if row_limit is not None and row_limit > 0 else 1000
    frame = pd.DataFrame()
    fields: List[Dict[str, object]] = []
    decoded_names: List[str] = []

    for index, root in enumerate(roots):
        page_offsets = list(root.get("page_offsets") or [])
        if not page_offsets:
            continue
        sample_count = min(rows_to_load, 64)
        sample_floats = _decode_tbms_flat_root_values(path, page_offsets, np.dtype("<f4"), sample_count)
        sample_ints = _decode_tbms_flat_root_values(path, page_offsets, np.dtype("<i4"), sample_count)
        name, dtype = _infer_flat_root_dtype_and_name(index, sample_floats, sample_ints)
        values = _decode_tbms_flat_root_values(path, page_offsets, dtype, rows_to_load)
        if values is None:
            continue
        frame[name] = values
        decoded_names.append(name)
        fields.append({
            "entry_type": "field",
            "entry_key": f"flat_root_{index:03d}",
            "name": name,
            "type": "float" if dtype == np.dtype("<f4") else "int",
            "location": int(root.get("root_offset") or 0),
            "page_count": len(page_offsets),
            "encoding": "tbms-flat-root",
            "inferred_name": bool(index >= 6),
        })

    metadata = {
        "file": str(path),
        "size_bytes": report.get("size_bytes"),
        "starts_with_tbms2": report.get("starts_with_tbms2"),
        "header_fields": report.get("header_fields"),
        "tbms_variant_family": "tbms-binary-flat-roots",
        "tbms_variant_evidence": recognition.get("evidence", []),
        "flat_root_count": len(roots),
        "decoded_fields": decoded_names,
        "row_limit_defaulted": row_limit is None or row_limit <= 0,
    }
    report["tbms_variant_family"] = "tbms-binary-flat-roots"
    report["tbms_variant_evidence"] = recognition.get("evidence", [])
    report["flat_roots"] = [
        {
            "root_offset": root.get("root_offset"),
            "page_count": len(root.get("page_offsets") or []),
            "first_pages": list(root.get("page_offsets") or [])[:8],
        }
        for root in roots[:50]
    ]
    return {
        "report": report,
        "metadata": metadata,
        "fields": fields,
        "dataframe": frame,
        "reader_mode": "tbms-binary-flat-roots",
        "row_count": 0,
        "rows_loaded": int(len(frame)),
    }


def _recognize_tbms_variant(path: Path, report: Dict[str, object]) -> Dict[str, object]:
    recognition = {
        "family": "unknown",
        "evidence": [],
    }

    if report.get("experimental_sections"):
        recognition["family"] = "tbms-experimental"
        recognition["evidence"].append("experimental section directory present")
        return recognition

    config_block = _find_tbms_config_block(path)
    if config_block is not None:
        recognition["family"] = "tbms-config-text"
        recognition["evidence"].append("balanced text config block parsed")
        recognition["config_block"] = config_block
        return recognition

    header_fields = report.get("header_fields") or {}
    first_offset = header_fields.get("u32_16")
    page_probe = _probe_tbms_first_page(path, first_offset if isinstance(first_offset, int) else None)
    recognition["page_probe"] = page_probe
    if report.get("starts_with_tbms2") and page_probe.get("available"):
        payload_words = int(page_probe.get("payload_words") or 0)
        reserved_word = int(page_probe.get("reserved_word") or 0)
        if payload_words in (256, 512, 768, 1024, 1536, 2048) and reserved_word == 0:
            header_fields = report.get("header_fields") or {}
            flat_roots = _scan_tbms_flat_roots(path, int(header_fields.get("first_section_offset") or 2056))
            if flat_roots:
                recognition["family"] = "tbms-binary-flat-roots"
                recognition["evidence"].append("TBMS header present with flat root pages and no text config")
                recognition["flat_roots"] = flat_roots
                return recognition
            recognition["family"] = "tbms-binary-pages"
            recognition["evidence"].append("TBMS header present with numeric first page and no text config")
            return recognition

    if report.get("starts_with_tbms2"):
        recognition["family"] = "tbms-unknown"
        recognition["evidence"].append("TBMS header present but no recognized metadata encoding")
    return recognition


def _load_tbms_config_variant(
    path: Path,
    report: Dict[str, object],
    config_block: Dict[str, object],
    row_limit: int | None,
) -> Dict[str, object]:
    parsed = config_block["parsed"]
    metadata = {
        "file": str(path),
        "size_bytes": report.get("size_bytes"),
        "starts_with_tbms2": report.get("starts_with_tbms2"),
        "header_fields": report.get("header_fields"),
    }
    if "directory_offset" in config_block:
        metadata["config_directory_offset"] = config_block.get("directory_offset")
        metadata["config_segments"] = config_block.get("segments", [])
    flat_metadata = _flatten_tbms_config(parsed)
    metadata.update(flat_metadata)
    schema_frame = _tbms_schema_fields(parsed)
    row_count = int(parsed.get("n_blocks", 0) or 0)
    if row_count <= 0:
        try:
            dim_values = [
                int(flat_metadata.get(key, flat_metadata.get(f"schema_0.{key}", 0)) or 0)
                for key in ("dim_x", "dim_y", "dim_z")
            ]
            if all(value > 0 for value in dim_values):
                row_count = int(dim_values[0] * dim_values[1] * dim_values[2])
        except Exception:
            row_count = 0
    if row_limit is not None and row_limit > 0:
        rows_to_load = min(int(row_limit), row_count)
    else:
        rows_to_load = row_count

    is_irregular = bool(int(parsed.get("is_irregular", 0) or 0))
    if is_irregular:
        frame = pd.DataFrame(index=pd.RangeIndex(rows_to_load))
    else:
        frame = _build_tbms_grid_frame(metadata, rows_to_load)
    decoded_fields: List[str] = []
    skipped_fields: List[str] = []
    decoded_columns: Dict[str, object] = {}
    if rows_to_load > 0 and not schema_frame.empty and "entry_type" in schema_frame.columns:
        field_rows = schema_frame[schema_frame["entry_type"] == "field"]
        for field in field_rows.to_dict(orient="records"):
            name = str(field.get("name") or "").strip()
            field_type = field.get("type")
            location = field.get("location")
            if not name or location in (None, ""):
                continue
            try:
                loc = int(float(location))
            except Exception:
                skipped_fields.append(name)
                continue

            values = _decode_tbms_field_values(
                path=path,
                location=loc,
                field_type=field_type,
                row_limit=rows_to_load,
            )
            if values is None or len(values) < rows_to_load:
                skipped_fields.append(name)
                continue
            if str(field_type or "").strip().lower() == "namedshort":
                labels = {
                    int(str(key).split("_", 1)[1]): value
                    for key, value in field.items()
                    if str(key).startswith("string_") and value is not None and value == value
                }
                if labels:
                    values = np.asarray([labels.get(int(value), value) for value in values], dtype=object)
            decoded_columns[name] = values
            decoded_fields.append(name)

    if decoded_columns:
        frame = pd.concat([frame, pd.DataFrame(decoded_columns, index=frame.index)], axis=1)

    if is_irregular and all(column in frame.columns for column in ("__lower_x", "__upper_x", "__lower_y", "__upper_y", "__lower_z", "__upper_z")):
        frame.insert(0, "x", (pd.to_numeric(frame["__lower_x"], errors="coerce") + pd.to_numeric(frame["__upper_x"], errors="coerce")) * 0.5)
        frame.insert(1, "y", (pd.to_numeric(frame["__lower_y"], errors="coerce") + pd.to_numeric(frame["__upper_y"], errors="coerce")) * 0.5)
        frame.insert(2, "z", (pd.to_numeric(frame["__lower_z"], errors="coerce") + pd.to_numeric(frame["__upper_z"], errors="coerce")) * 0.5)

    report["config_block"] = {
        "start": config_block["start"],
        "end": config_block["end"],
        "length": int(config_block["end"] - config_block["start"] + 1),
    }
    if "directory_offset" in config_block:
        report["config_block"]["directory_offset"] = config_block.get("directory_offset")
        report["config_block"]["segments"] = config_block.get("segments", [])
    report["decoded_fields"] = decoded_fields
    report["skipped_fields"] = skipped_fields
    return {
        "report": report,
        "metadata": metadata,
        "fields": schema_frame.to_dict(orient="records"),
        "dataframe": frame,
        "reader_mode": "tbms-config",
        "row_count": row_count,
        "rows_loaded": int(len(frame)),
    }


def _load_tbms_binary_variant(
    path: Path,
    report: Dict[str, object],
    recognition: Dict[str, object],
    row_limit: int | None = None,
) -> Dict[str, object]:
    if recognition.get("family") == "tbms-binary-flat-roots":
        return _load_tbms_headerless_flat_roots(path, report, recognition, row_limit)

    metadata = {
        "file": str(path),
        "size_bytes": report.get("size_bytes"),
        "starts_with_tbms2": report.get("starts_with_tbms2"),
        "header_fields": report.get("header_fields"),
        "tbms_variant_family": recognition.get("family", "tbms-binary-pages"),
        "tbms_variant_evidence": recognition.get("evidence", []),
        "page_probe": recognition.get("page_probe", {}),
    }
    report["tbms_variant_family"] = recognition.get("family", "tbms-binary-pages")
    report["tbms_variant_evidence"] = recognition.get("evidence", [])
    if recognition.get("page_probe"):
        report["page_probe"] = recognition["page_probe"]

    return {
        "report": report,
        "metadata": metadata,
        "fields": [],
        "dataframe": pd.DataFrame(),
        "reader_mode": "tbms-binary-pages",
        "row_count": 0,
        "rows_loaded": 0,
    }


def export_bmf(
    input_csv: str | Path,
    output_bmf: str | Path,
    backend: str = "tbms-config-text",
    x_col: str = "x",
    y_col: str = "y",
    z_col: str = "z",
    delimiter: str | None = None,
    header_line: int = 1,
    cell_size: Sequence[float] | None = None,
    origin: Sequence[float] | None = None,
    value_cols: Sequence[str] | None = None,
    size_cols: Sequence[str] | Mapping[str, str] | str | None = None,
    extent_cols: Sequence[str] | Mapping[str, str] | str | None = None,
    column_types: Mapping[str, object] | None = None,
    value_exceptions: Mapping[str, Mapping[str, object]] | None = None,
    null_float: float = -99.0,
    index_tolerance: float = 1e-3,
    regularize_to_base_block: bool = False,
    dry_run: bool = False,
    summary_json: str | Path | None = None,
    progress_callback=None,
) -> Dict[str, object]:
    in_csv = Path(input_csv)
    out_bmf = Path(output_bmf)
    if not in_csv.exists():
        raise FileNotFoundError(f"Input file not found: {in_csv}")

    user_cell_size_provided = cell_size is not None
    user_origin_provided = origin is not None
    _emit_progress(progress_callback, 0, 100, "Preparing BMF export...")
    csv_metadata = parse_leapfrog_block_metadata(in_csv)
    log_leapfrog_metadata_summary(in_csv, csv_metadata, context="BMF export")
    effective_header_line = resolve_effective_csv_header_line(in_csv, header_line)
    metadata_geometry = infer_bmf_export_geometry_from_metadata(
        csv_metadata,
        regularize_to_base_block=regularize_to_base_block,
        dense_regular_grid=_backend_requires_dense_grid(backend),
    )
    if cell_size is None and metadata_geometry.get("cell_size"):
        cell_size = metadata_geometry["cell_size"]  # type: ignore[assignment]
        _emit_progress(
            progress_callback,
            1,
            100,
            f"Using CSV metadata cell size from {metadata_geometry.get('cell_size_source')}...",
        )
    if origin is None and metadata_geometry.get("origin"):
        origin = metadata_geometry["origin"]  # type: ignore[assignment]
        _emit_progress(
            progress_callback,
            2,
            100,
            f"Using CSV metadata origin from {metadata_geometry.get('origin_source')}...",
        )
    xyz_cols = [x_col, y_col, z_col]
    csv_delimiter = delimiter or detect_csv_delimiter(in_csv)
    csv_header_names = parse_header_line(in_csv, csv_delimiter, effective_header_line)
    available_grid_index_cols = [col for col in GRID_INDEX_COLUMNS if col in csv_header_names]
    mapped_extent_cols = _resolve_geometry_column_map(csv_header_names, extent_cols, GEOMETRY_EXTENT_COLUMNS, "lower/upper extent columns")
    mapped_size_cols = _resolve_geometry_column_map(csv_header_names, size_cols, GEOMETRY_SIZE_COLUMN_ROLES, "dimension columns")
    if not mapped_extent_cols:
        mapped_extent_cols = _detect_geometry_extent_columns(csv_header_names)
    if not mapped_size_cols:
        mapped_size_cols = _detect_geometry_size_columns(csv_header_names)
    explicit_geometry_cols = list(dict.fromkeys([*mapped_extent_cols.values(), *mapped_size_cols.values()]))
    has_explicit_row_geometry = bool(mapped_extent_cols or mapped_size_cols)
    use_irregular_tbms_config = (
        backend == "tbms-config-text"
        and not regularize_to_base_block
        and (has_explicit_row_geometry or _metadata_supports_tbms_irregular(csv_metadata))
    )
    use_grid_index_columns = (
        len(available_grid_index_cols) == len(GRID_INDEX_COLUMNS)
        and not regularize_to_base_block
        and not use_irregular_tbms_config
    )
    source_row_mode = bool(backend == "tbms-config-text" and not regularize_to_base_block)
    header_value_cols = _classify_header_value_columns(
        csv_header_names,
        xyz_cols=xyz_cols,
        value_cols=value_cols,
        geometry_columns=explicit_geometry_cols,
    ) if source_row_mode else []
    selected_read_cols = (
        list(dict.fromkeys([*xyz_cols, *available_grid_index_cols, *(explicit_geometry_cols if use_irregular_tbms_config else []), *(value_cols or [])]))
        if value_cols else None
    )
    if source_row_mode and selected_read_cols is None:
        selected_read_cols = list(dict.fromkeys([*xyz_cols, *(explicit_geometry_cols if use_irregular_tbms_config else []), *header_value_cols]))
    selected_message = "selected columns" if selected_read_cols else "all columns"
    counted_source_rows = _count_csv_data_rows(in_csv, effective_header_line, progress_callback, 3, 5) if selected_read_cols else None
    selected_column_count = len(selected_read_cols or [])
    estimated_selected_bytes = int(counted_source_rows or 0) * selected_column_count * 8
    estimated_coordinate_bytes = int(counted_source_rows or 0) * 3 * 8
    can_stream_irregular_source_rows = bool(
        source_row_mode
        and selected_read_cols
        and not use_grid_index_columns
        and (use_irregular_tbms_config or cell_size is not None or csv_metadata.get("parent_block_size"))
    )
    should_stream_irregular_source_rows = bool(
        can_stream_irregular_source_rows
        and counted_source_rows is not None
        and (
            estimated_coordinate_bytes > MAX_SOURCE_ROW_PREP_BYTES
            or estimated_selected_bytes > MAX_SELECTED_CSV_OBJECT_BYTES
        )
    )

    if should_stream_irregular_source_rows:
        chosen_value_cols = header_value_cols
        if not chosen_value_cols:
            raise ValueError("No valid value columns were available for streaming source-row BMF export.")
        normalized_value_exceptions = _normalize_value_exceptions(value_exceptions)
        _emit_progress(
            progress_callback,
            6,
            100,
            "Using streaming source-row BMF writer; CSV rows will not be materialized as one table.",
        )
        parent_cell_size = cell_size
        stream_scan = _scan_stream_irregular_export(
            input_csv=in_csv,
            delimiter=delimiter,
            header_line=effective_header_line,
            selected_read_cols=selected_read_cols,
            xyz_cols=xyz_cols,
            value_cols=chosen_value_cols,
            metadata=csv_metadata,
            cell_size=parent_cell_size,
            origin=origin,
            index_tolerance=index_tolerance,
            geometry_size_cols=mapped_size_cols,
            geometry_extent_cols=mapped_extent_cols,
            column_types=column_types,
            value_exceptions=normalized_value_exceptions,
            null_float=null_float,
            total_rows=counted_source_rows,
            progress_callback=progress_callback,
        )
        prepared_summary = {
            "df": pd.DataFrame(index=range(0)),
            "origin": stream_scan["origin"],
            "cell": stream_scan["cell"],
            "dims": stream_scan["dims"],
            "extents": stream_scan["extents"],
            "duplicates": 0,
            "max_index_error": 0.0,
            "is_irregular": True,
            "irregular_width_source": (
                "mapped_lower_upper_columns" if mapped_extent_cols and any(mapped_extent_cols.get(role) != role for role in GEOMETRY_EXTENT_COLUMNS)
                else "lower_upper_columns" if mapped_extent_cols
                else "dimension_columns" if mapped_size_cols
                else "centroid_hierarchy_inference"
            ),
            "grid_index_source": (
                "irregular_mapped_lower_upper_columns" if mapped_extent_cols and any(mapped_extent_cols.get(role) != role for role in GEOMETRY_EXTENT_COLUMNS)
                else "irregular_lower_upper_columns" if mapped_extent_cols
                else "irregular_dimension_columns" if mapped_size_cols
                else "irregular_xyz_metadata_parent"
            ),
        }
        if mapped_extent_cols:
            prepared_summary["geometry_column_mapping"] = {"extent_cols": dict(mapped_extent_cols)}
        elif mapped_size_cols:
            prepared_summary["geometry_column_mapping"] = {"size_cols": dict(mapped_size_cols)}
        block_size_determination = _build_block_size_determination_report(
            prepared=prepared_summary,
            metadata_geometry=metadata_geometry,
            statistical_base_geometry={},
            user_cell_size_provided=user_cell_size_provided,
            user_origin_provided=user_origin_provided,
            use_grid_index_columns=False,
            regularize_to_base_block=regularize_to_base_block,
        )
        block_size_message = (
            f"Block size determination: {block_size_determination['message']} "
            f"cell_size={_format_numeric_vector(block_size_determination['cell_size'])}, "
            f"origin={_format_numeric_vector(block_size_determination['origin'])}."
        )
        print(block_size_message)
        summary = {
            "input_csv": str(in_csv),
            "output_bmf": str(out_bmf),
            "backend": backend,
            "streaming": True,
            "csv_header_line": int(effective_header_line),
            "csv_configured_header_line": int(max(int(header_line or 1), 1)),
            "csv_metadata": csv_metadata,
            "csv_metadata_geometry": metadata_geometry,
            "rows": int(stream_scan["row_count"]),
            "value_columns": chosen_value_cols,
            "value_exceptions": normalized_value_exceptions,
            "csv_statistical_base_block": {},
            "block_size_determination": block_size_determination,
            "regularization": {"enabled": False},
            "grid": {
                "origin": [float(x) for x in np.asarray(stream_scan["origin"], dtype=float)],
                "cell_size": [float(x) for x in np.asarray(stream_scan["cell"], dtype=float)],
                "dimensions": [int(x) for x in np.asarray(stream_scan["dims"], dtype=int)],
                "extents": [float(x) for x in np.asarray(stream_scan["extents"], dtype=float)],
                "duplicate_grid_rows": 0,
                "max_index_error": 0.0,
                "index_source": str(prepared_summary.get("grid_index_source")),
                "is_irregular": True,
                "irregular_width_source": prepared_summary.get("irregular_width_source"),
                "irregular_hierarchy": None,
                "irregular_hierarchy_text": "",
            },
        }
        if summary_json:
            summary_path = Path(summary_json)
            summary_path.parent.mkdir(parents=True, exist_ok=True)
            summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
            _emit_progress(progress_callback, 49, 100, "BMF export summary written.")
        backend_summary = None
        if not dry_run:
            backend_summary = _write_streaming_irregular_tbms_config_text(
                input_csv=in_csv,
                prepared=stream_scan,
                xyz_cols=xyz_cols,
                value_cols=chosen_value_cols,
                selected_read_cols=selected_read_cols,
                delimiter=delimiter,
                header_line=effective_header_line,
                metadata=csv_metadata,
                cell_size=parent_cell_size,
                origin=origin,
                index_tolerance=index_tolerance,
                geometry_size_cols=mapped_size_cols,
                geometry_extent_cols=mapped_extent_cols,
                value_exceptions=normalized_value_exceptions,
                out_bmf=out_bmf,
                null_float=null_float,
                progress_callback=progress_callback,
                progress_start=50,
                progress_end=98,
                total_rows=counted_source_rows,
            )
        else:
            _emit_progress(progress_callback, 98, 100, "Dry-run streaming BMF export validation complete.")
        _emit_progress(progress_callback, 100, 100, "BMF export complete.")
        return {"summary": summary, "backend_summary": backend_summary, "dry_run": bool(dry_run)}

    _validate_selected_csv_read_size(
        in_csv,
        effective_header_line,
        selected_read_cols,
        backend,
        regularize_to_base_block,
        source_row_mode=source_row_mode,
        progress_callback=None if counted_source_rows is not None else progress_callback,
    )
    _emit_progress(progress_callback, 5, 100, f"Reading CSV {selected_message}...")
    df = _auto_read_csv(
        in_csv,
        delimiter=delimiter,
        header_line=effective_header_line,
        usecols=selected_read_cols,
        progress_callback=progress_callback,
        progress_start=5,
        progress_end=30,
    )
    _emit_progress(
        progress_callback,
        30,
        100,
        f"CSV read complete: {len(df):,} rows, {len(df.columns):,} columns.",
    )
    missing = [col for col in xyz_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing coordinate columns: {missing}")

    statistical_base_geometry = {}
    if cell_size is None:
        statistical_base_geometry = infer_base_block_size_from_centroid_statistics(
            df,
            xyz_cols,
            index_tolerance=index_tolerance,
        )
    statistical_cell_size = statistical_base_geometry.get("cell_size") if isinstance(statistical_base_geometry, Mapping) else None
    statistical_origin = statistical_base_geometry.get("origin") if isinstance(statistical_base_geometry, Mapping) else None
    use_provided_parent_irregular_tbms_config = bool(
        backend == "tbms-config-text"
        and not regularize_to_base_block
        and not use_irregular_tbms_config
        and not use_grid_index_columns
        and cell_size is not None
    )
    use_statistical_irregular_tbms_config = bool(
        backend == "tbms-config-text"
        and not regularize_to_base_block
        and not use_irregular_tbms_config
        and not use_grid_index_columns
        and not use_provided_parent_irregular_tbms_config
        and statistical_cell_size
    )

    normalized_value_exceptions = _normalize_value_exceptions(value_exceptions)
    regularization_summary = {"enabled": False}
    if regularize_to_base_block:
        if cell_size is None and statistical_base_geometry.get("cell_size"):
            cell_size = statistical_base_geometry["cell_size"]  # type: ignore[assignment]
            if origin is None and statistical_base_geometry.get("origin"):
                origin = statistical_base_geometry["origin"]  # type: ignore[assignment]
            _emit_progress(
                progress_callback,
                31,
                100,
                "Using statistically inferred base block size for regularization...",
            )
        pre_regularization_exceptions = _filter_value_exceptions_for_regularization(
            normalized_value_exceptions,
            include_in_regularization=True,
        )
        if pre_regularization_exceptions:
            _emit_progress(progress_callback, 32, 100, "Applying BMF value exception rules before regularization...")
            df = _apply_value_exceptions(df, pre_regularization_exceptions)

        _emit_progress(progress_callback, 34, 100, "Regularizing CSV rows to base block cells...")
        df, regularization_summary = _regularize_to_base_cell_grid(
            df=df,
            xyz_cols=xyz_cols,
            cell_size=cell_size,
            origin=origin,
            value_cols=value_cols,
            column_types=column_types,
        )

        post_regularization_exceptions = _filter_value_exceptions_for_regularization(
            normalized_value_exceptions,
            include_in_regularization=False,
        )
        if post_regularization_exceptions:
            _emit_progress(progress_callback, 35, 100, "Applying BMF value exception rules after regularization...")
            df = _apply_value_exceptions(df, post_regularization_exceptions)
    elif normalized_value_exceptions:
        _emit_progress(progress_callback, 32, 100, "Applying BMF value exception rules...")
        df = _apply_value_exceptions(df, normalized_value_exceptions)

    if use_irregular_tbms_config or use_provided_parent_irregular_tbms_config or use_statistical_irregular_tbms_config:
        _emit_progress(progress_callback, 36, 100, "Preparing irregular BMF blocks from source CSV row extents...")
        parent_cell_size = cell_size or statistical_cell_size
        try:
            prepared = _prepare_irregular_blocks(
                df=df,
                xyz_cols=xyz_cols,
                metadata=csv_metadata,
                cell_size=parent_cell_size,
                origin=origin or statistical_origin,
                index_tolerance=index_tolerance,
                geometry_size_cols=mapped_size_cols,
                geometry_extent_cols=mapped_extent_cols,
            )
        except ValueError as irregular_error:
            if use_provided_parent_irregular_tbms_config and parent_cell_size is not None:
                raise ValueError(
                    _build_centroid_hierarchy_error_message(irregular_error, df, xyz_cols, parent_cell_size)
                ) from irregular_error
            raise
        if use_statistical_irregular_tbms_config:
            prepared["statistical_base_block_source"] = "centroid_spacing_statistics"
    elif use_grid_index_columns and all(col in df.columns for col in GRID_INDEX_COLUMNS):
        _emit_progress(progress_callback, 36, 100, "Preparing regular BMF grid from CSV grid_i/grid_j/grid_k columns...")
        prepared = _prepare_grid_from_index_columns(
            df=df,
            xyz_cols=xyz_cols,
            cell_size=cell_size,
            origin=origin,
        )
    else:
        _emit_progress(progress_callback, 36, 100, "Preparing regular BMF grid from XYZ columns...")
        try:
            prepared = _prepare_grid(
                df=df,
                xyz_cols=xyz_cols,
                cell_size=cell_size,
                origin=origin,
                index_tolerance=index_tolerance,
            )
        except ValueError as regular_grid_error:
            inferred_cell_size = statistical_base_geometry.get("cell_size") if isinstance(statistical_base_geometry, Mapping) else None
            inferred_origin = statistical_base_geometry.get("origin") if isinstance(statistical_base_geometry, Mapping) else None
            can_use_statistical_base = bool(inferred_cell_size and backend == "tbms-config-text" and not regularize_to_base_block)
            if backend != "tbms-config-text" or regularize_to_base_block or (cell_size is None and not can_use_statistical_base):
                raise
            _emit_progress(
                progress_callback,
                37,
                100,
                "Regular grid alignment failed; inferring irregular row extents from centroid hierarchy...",
            )
            try:
                prepared = _prepare_irregular_blocks(
                    df=df,
                    xyz_cols=xyz_cols,
                    metadata=csv_metadata,
                    cell_size=cell_size or inferred_cell_size,
                    origin=origin or inferred_origin,
                    index_tolerance=index_tolerance,
                    geometry_size_cols=mapped_size_cols,
                    geometry_extent_cols=mapped_extent_cols,
                )
                if cell_size is None and can_use_statistical_base:
                    prepared["statistical_base_block_source"] = "centroid_spacing_statistics"
            except ValueError as irregular_error:
                raise regular_grid_error from irregular_error
    dims_text = " x ".join(str(int(value)) for value in np.asarray(prepared["dims"], dtype=int))
    geometry_text = "Irregular blocks prepared" if prepared.get("is_irregular") else "Grid prepared"
    _emit_progress(
        progress_callback,
        45,
        100,
        f"{geometry_text}: {dims_text} parent cells from {len(prepared['df']):,} valid CSV rows.",
    )

    chosen_value_cols = _classify_columns(prepared["df"], xyz_cols, value_cols)
    _emit_progress(progress_callback, 47, 100, f"Selected {len(chosen_value_cols)} BMF value field(s).")
    if backend not in SUPPORTED_EXPORT_BACKENDS:
        raise ValueError(f"Unsupported backend: {backend}")
    if _backend_requires_dense_grid(backend) and not prepared.get("is_irregular"):
        _validate_dense_export_size(prepared, chosen_value_cols, backend=backend)
    irregular_hierarchy_report = _build_irregular_hierarchy_report(
        prepared,
        xyz_cols=xyz_cols,
        index_tolerance=index_tolerance,
    )
    irregular_hierarchy_text = _format_irregular_hierarchy_report(irregular_hierarchy_report)
    if irregular_hierarchy_text:
        print(irregular_hierarchy_text)
    block_size_determination = _build_block_size_determination_report(
        prepared=prepared,
        metadata_geometry=metadata_geometry,
        statistical_base_geometry=statistical_base_geometry,
        user_cell_size_provided=user_cell_size_provided,
        user_origin_provided=user_origin_provided,
        use_grid_index_columns=use_grid_index_columns,
        regularize_to_base_block=regularize_to_base_block,
    )
    block_size_message = (
        f"Block size determination: {block_size_determination['message']} "
        f"cell_size={_format_numeric_vector(block_size_determination['cell_size'])}, "
        f"origin={_format_numeric_vector(block_size_determination['origin'])}."
    )
    print(block_size_message)
    _emit_progress(progress_callback, 46, 100, block_size_message)
    summary = {
        "input_csv": str(in_csv),
        "output_bmf": str(out_bmf),
        "backend": backend,
        "csv_header_line": int(effective_header_line),
        "csv_configured_header_line": int(max(int(header_line or 1), 1)),
        "csv_metadata": csv_metadata,
        "csv_metadata_geometry": metadata_geometry,
        "rows": int(len(prepared["df"])),
        "value_columns": chosen_value_cols,
        "value_exceptions": normalized_value_exceptions,
        "csv_statistical_base_block": statistical_base_geometry,
        "block_size_determination": block_size_determination,
        "regularization": regularization_summary,
        "grid": {
            "origin": [float(x) for x in np.asarray(prepared["origin"], dtype=float)],
            "cell_size": [float(x) for x in np.asarray(prepared["cell"], dtype=float)],
            "dimensions": [int(x) for x in np.asarray(prepared["dims"], dtype=int)],
            "extents": [float(x) for x in np.asarray(prepared["extents"], dtype=float)],
            "duplicate_grid_rows": int(prepared["duplicates"]),
            "max_index_error": float(prepared["max_index_error"]),
            "index_source": str(prepared.get("grid_index_source") or "xyz"),
            "is_irregular": bool(prepared.get("is_irregular", False)),
            "irregular_width_source": prepared.get("irregular_width_source"),
            "irregular_hierarchy": irregular_hierarchy_report,
            "irregular_hierarchy_text": irregular_hierarchy_text,
        },
    }

    if summary_json:
        summary_path = Path(summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        _emit_progress(progress_callback, 49, 100, "BMF export summary written.")

    backend_summary = None
    if not dry_run:
        out_bmf.parent.mkdir(parents=True, exist_ok=True)
        if backend == "vulcan":
            _emit_progress(progress_callback, 50, 100, "Writing BMF through Vulcan backend...")
            _export_with_vulcan(
                prepared=prepared,
                xyz_cols=xyz_cols,
                value_cols=chosen_value_cols,
                out_bmf=out_bmf,
                null_float=null_float,
                column_types=column_types,
            )
            backend_summary = {"backend": "vulcan"}
            _emit_progress(progress_callback, 98, 100, "Vulcan BMF write complete.")
        elif backend == "tbms-config-text":
            backend_summary = _write_tbms_config_text(
                prepared=prepared,
                xyz_cols=xyz_cols,
                value_cols=chosen_value_cols,
                out_bmf=out_bmf,
                null_float=null_float,
                column_types=column_types,
                progress_callback=progress_callback,
                progress_start=50,
                progress_end=98,
            )
        elif backend == "tbms-experimental":
            backend_summary = _write_tbms_experimental(
                prepared=prepared,
                xyz_cols=xyz_cols,
                value_cols=chosen_value_cols,
                out_bmf=out_bmf,
                null_float=null_float,
                column_types=column_types,
                progress_callback=progress_callback,
                progress_start=50,
                progress_end=98,
            )
        else:
            raise ValueError(f"Unsupported backend: {backend}")
    else:
        _emit_progress(progress_callback, 98, 100, "Dry-run BMF export validation complete.")

    _emit_progress(progress_callback, 100, 100, "BMF export complete.")

    return {
        "summary": summary,
        "backend_summary": backend_summary,
        "dry_run": bool(dry_run),
    }


def inspect_bmf_file(
    input_bmf: str | Path,
    scan_bytes: int = 2 * 1024 * 1024,
    min_string_len: int = 4,
    max_strings: int = 30,
    top_chunks: int = 20,
) -> Dict[str, object]:
    path = Path(input_bmf)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    with path.open("rb") as fh:
        head = fh.read(scan_bytes)
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - 128), os.SEEK_SET)
        tail = fh.read(128)

    header_ascii = head[:16].decode("ascii", errors="replace")
    starts_tbms = head.startswith(TBMS_SIGNATURE)

    ints = []
    for off in range(8, min(132, len(head) - 3), 4):
        ints.append({"offset": off, "int32_le": int(struct.unpack_from("<i", head, off)[0])})

    strings = _extract_strings(head, min_len=min_string_len)
    strings = strings[: max_strings]

    chunks = [head[i : i + 4] for i in range(0, max(0, len(head) - 3))]
    common_4b = [
        {"hex": chunk.hex(), "count": count}
        for chunk, count in Counter(chunks).most_common(top_chunks)
    ]

    report = {
        "file": str(path),
        "size_bytes": int(size),
        "header_hex": head[:32].hex(),
        "header_ascii": header_ascii,
        "starts_with_tbms2": bool(starts_tbms),
        "header_fields": {
            "u32_12": int(struct.unpack_from("<I", head, 12)[0]) if len(head) >= 16 else None,
            "u32_16": int(struct.unpack_from("<I", head, 16)[0]) if len(head) >= 20 else None,
            "u64_48": int(struct.unpack_from("<Q", head, 48)[0]) if len(head) >= 56 else None,
        },
        "int32_candidates": ints,
        "top_4byte_chunks": common_4b,
        "strings": [{"offset": offset, "value": value} for offset, value in strings],
        "tail_hex": tail.hex(),
    }

    experimental_sections = _parse_tbms_experimental_sections(path)
    if experimental_sections is not None:
        report["experimental_sections"] = experimental_sections

    return report


def load_bmf_table(input_bmf: str | Path, row_limit: int | None = None) -> Dict[str, object]:
    path = Path(input_bmf)
    report = inspect_bmf_file(path, scan_bytes=256 * 1024)
    result: Dict[str, object] = {
        "report": report,
        "metadata": {},
        "fields": [],
        "dataframe": pd.DataFrame(),
        "reader_mode": "basic-inspect",
        "row_count": 0,
        "rows_loaded": 0,
    }

    sections = report.get("experimental_sections")
    if not sections:
        recognition = _recognize_tbms_variant(path, report)
        report["tbms_variant_family"] = recognition.get("family", "unknown")
        report["tbms_variant_evidence"] = recognition.get("evidence", [])

        if recognition.get("family") == "tbms-config-text" and recognition.get("config_block") is not None:
            return _load_tbms_config_variant(path, report, recognition["config_block"], row_limit)

        if recognition.get("family") in ("tbms-binary-pages", "tbms-binary-flat-roots"):
            return _load_tbms_binary_variant(path, report, recognition, row_limit)

        metadata = {
            "file": str(path),
            "size_bytes": report.get("size_bytes"),
            "starts_with_tbms2": report.get("starts_with_tbms2"),
            "header_fields": report.get("header_fields"),
            "tbms_variant_family": recognition.get("family", "unknown"),
            "tbms_variant_evidence": recognition.get("evidence", []),
        }
        result["metadata"] = metadata
        return result

    section_map = {section["name"]: section for section in sections if isinstance(section, dict) and section.get("name")}
    metadata = dict(section_map.get("metadata", {}).get("json_preview") or {})
    field_specs = list(section_map.get("columns", {}).get("json_preview") or [])
    row_count = int(metadata.get("row_count", 0))
    if row_limit is not None and row_limit > 0:
        rows_to_load = min(int(row_limit), row_count)
    else:
        rows_to_load = row_count

    numeric_columns = [field["name"] for field in field_specs if field.get("kind") == "numeric"]
    string_columns = [field["name"] for field in field_specs if field.get("kind") == "string"]

    x_col = str(metadata.get("x_col", "x"))
    y_col = str(metadata.get("y_col", "y"))
    z_col = str(metadata.get("z_col", "z"))

    frame = pd.DataFrame()
    if rows_to_load > 0 and "indices" in section_map:
        indices_blob = _read_section_bytes(path, section_map["indices"], limit_bytes=rows_to_load * 3 * 4)
        indices = np.frombuffer(indices_blob, dtype="<i4").reshape((-1, 3))
        origin = np.asarray(metadata.get("origin", [0.0, 0.0, 0.0]), dtype=float)
        cell_size = np.asarray(metadata.get("cell_size", [1.0, 1.0, 1.0]), dtype=float)
        centers = origin + (indices.astype(float) + 0.5) * cell_size

        frame["grid_i"] = indices[:, 0]
        frame["grid_j"] = indices[:, 1]
        frame["grid_k"] = indices[:, 2]
        frame[x_col] = centers[:, 0]
        frame[y_col] = centers[:, 1]
        frame[z_col] = centers[:, 2]

        if numeric_columns and "numeric" in section_map:
            numeric_blob = _read_section_bytes(
                path,
                section_map["numeric"],
                limit_bytes=rows_to_load * len(numeric_columns) * 8,
            )
            numeric = np.frombuffer(numeric_blob, dtype="<f8").reshape((-1, len(numeric_columns)))
            numeric_df = pd.DataFrame(numeric, columns=numeric_columns)
            frame = pd.concat([frame.reset_index(drop=True), numeric_df.reset_index(drop=True)], axis=1)

        if string_columns and "strings" in section_map:
            strings_payload = _read_section_bytes(path, section_map["strings"])
            string_records = json.loads(strings_payload.decode("utf-8"))
            for column_name in string_columns:
                values = list(string_records.get(column_name, []))[:rows_to_load]
                if len(values) < rows_to_load:
                    values.extend([None] * (rows_to_load - len(values)))
                frame[column_name] = values

    result.update({
        "metadata": metadata,
        "fields": field_specs,
        "dataframe": frame,
        "reader_mode": "tbms-experimental",
        "row_count": row_count,
        "rows_loaded": int(len(frame)),
    })
    return result


def _export_with_vulcan(
    prepared: Dict[str, object],
    xyz_cols: Sequence[str],
    value_cols: Sequence[str],
    out_bmf: Path,
    null_float: float,
    column_types: Mapping[str, object] | None = None,
) -> None:
    try:
        import vulcan  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "Vulcan Python module is not available. BMF export requires Maptek Vulcan runtime."
        ) from exc

    df: pd.DataFrame = prepared["df"]  # type: ignore[assignment]
    origin = np.asarray(prepared["origin"], dtype=float)
    dims = np.asarray(prepared["dims"], dtype=int)
    extents = np.asarray(prepared["extents"], dtype=float)
    normalized_overrides = _normalize_export_type_overrides(column_types)

    bm = vulcan.block_model()
    bm.create_regular(
        str(out_bmf),
        0.0,
        0.0,
        0.0,
        float(extents[0]),
        float(extents[1]),
        float(extents[2]),
        int(dims[0]),
        int(dims[1]),
        int(dims[2]),
    )
    bm.set_model_origin(float(origin[0]), float(origin[1]), float(origin[2]))

    dtype_map: Dict[str, str] = {}
    for col in value_cols:
        col_series = df[col]
        forced_type = normalized_overrides.get(col)
        if forced_type in {"boolean", "int", "double"} or (forced_type is None and pd.api.types.is_numeric_dtype(col_series)):
            bm.add_variable(col, "float", str(null_float), "")
            dtype_map[col] = "float"
        else:
            bm.add_variable(col, "name", "n", "")
            dtype_map[col] = "name"

    bm.write()
    bm.index_model()
    bm.write()

    coords = df[list(xyz_cols)].to_numpy(dtype=float)
    for i, world in enumerate(coords):
        if bm.find_world_xyz(float(world[0]), float(world[1]), float(world[2])):
            continue
        row = df.iloc[i]
        for col in value_cols:
            mode = dtype_map[col]
            if mode == "float":
                v = pd.to_numeric(pd.Series([row[col]]), errors="coerce").iloc[0]
                if pd.isna(v):
                    bm.put(col, float(null_float))
                else:
                    bm.put(col, float(v))
            else:
                text = "" if pd.isna(row[col]) else str(row[col])
                bm.put_string(col, text)

    bm.write()


def cmd_export(args: argparse.Namespace) -> int:
    in_csv = Path(args.input_csv)
    out_bmf = Path(args.output_bmf)
    if not in_csv.exists():
        print(f"Input file not found: {in_csv}", file=sys.stderr)
        return 2
    try:
        result = export_bmf(
            input_csv=in_csv,
            output_bmf=out_bmf,
            backend=args.backend,
            x_col=args.x_col,
            y_col=args.y_col,
            z_col=args.z_col,
            cell_size=args.cell_size,
            origin=args.origin,
            value_cols=args.value_cols,
            size_cols=args.size_cols,
            extent_cols=args.extent_cols,
            null_float=args.null_float,
            index_tolerance=args.index_tolerance,
            regularize_to_base_block=getattr(args, "regularize_to_base_block", False),
            dry_run=args.dry_run,
            summary_json=args.summary_json,
        )
    except Exception as exc:
        print(f"BMF export failed: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(result["summary"], indent=2))
    if result.get("backend_summary"):
        print(json.dumps(result["backend_summary"], indent=2))
    if args.dry_run:
        print("Dry-run mode: no BMF file written.")
        return 0

    print(f"BMF export completed: {out_bmf}")
    return 0


def _extract_strings(blob: bytes, min_len: int = 4) -> List[Tuple[int, str]]:
    strings: List[Tuple[int, str]] = []
    start = None
    for i, b in enumerate(blob):
        if 32 <= b <= 126:
            if start is None:
                start = i
        else:
            if start is not None and (i - start) >= min_len:
                strings.append((start, blob[start:i].decode("ascii", errors="ignore")))
            start = None
    if start is not None and (len(blob) - start) >= min_len:
        strings.append((start, blob[start:].decode("ascii", errors="ignore")))
    return strings


def cmd_inspect(args: argparse.Namespace) -> int:
    path = Path(args.input_bmf)
    if not path.exists():
        print(f"Input file not found: {path}", file=sys.stderr)
        return 2

    with path.open("rb") as fh:
        head = fh.read(args.scan_bytes)
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - 128), os.SEEK_SET)
        tail = fh.read(128)

    header_ascii = head[:16].decode("ascii", errors="replace")
    starts_tbms = head.startswith(b"TBMS2.0\x00")

    ints = []
    for off in range(8, min(132, len(head) - 3), 4):
        ints.append({"offset": off, "int32_le": int(struct.unpack_from("<i", head, off)[0])})

    strings = _extract_strings(head, min_len=args.min_string_len)
    strings = strings[: args.max_strings]

    chunks = [head[i : i + 4] for i in range(0, max(0, len(head) - 3))]
    common_4b = [
        {"hex": b.hex(), "count": c}
        for b, c in Counter(chunks).most_common(args.top_chunks)
    ]

    report = {
        "file": str(path),
        "size_bytes": int(size),
        "header_hex": head[:32].hex(),
        "header_ascii": header_ascii,
        "starts_with_tbms2": bool(starts_tbms),
        "header_fields": {
            "u32_12": int(struct.unpack_from("<I", head, 12)[0]) if len(head) >= 16 else None,
            "u32_16": int(struct.unpack_from("<I", head, 16)[0]) if len(head) >= 20 else None,
            "u64_48": int(struct.unpack_from("<Q", head, 48)[0]) if len(head) >= 56 else None,
        },
        "int32_candidates": ints,
        "top_4byte_chunks": common_4b,
        "strings": [{"offset": o, "value": s} for o, s in strings],
        "tail_hex": tail.hex(),
    }

    experimental_sections = _parse_tbms_experimental_sections(path)
    if experimental_sections is not None:
        report["experimental_sections"] = experimental_sections

    print(json.dumps(report, indent=2))
    if args.output_json:
        out = Path(args.output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Standalone Anterpolator grid-to-BMF tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    export_p = sub.add_parser("export", help="Export regular-grid CSV to BMF")
    export_p.add_argument("input_csv", help="Input grid CSV path")
    export_p.add_argument("output_bmf", help="Output BMF path")
    export_p.add_argument(
        "--backend",
        choices=["vulcan", "tbms-config-text", "tbms-experimental"],
        default="tbms-config-text",
        help="Export backend. tbms-config-text writes a reverse-engineered TBMS config/page layout that this tool can reopen; tbms-experimental writes the older standalone container.",
    )
    export_p.add_argument("--x-col", default="x", help="X coordinate column")
    export_p.add_argument("--y-col", default="y", help="Y coordinate column")
    export_p.add_argument("--z-col", default="z", help="Z coordinate column")
    export_p.add_argument(
        "--cell-size",
        nargs=3,
        type=float,
        metavar=("DX", "DY", "DZ"),
        default=None,
        help="Explicit cell size; if omitted, inferred from coordinates",
    )
    export_p.add_argument(
        "--origin",
        nargs=3,
        type=float,
        metavar=("OX", "OY", "OZ"),
        default=None,
        help="Explicit model origin; if omitted, inferred from min XYZ",
    )
    export_p.add_argument(
        "--value-cols",
        nargs="+",
        default=None,
        help="Columns to export as variables; default is all non-XYZ columns",
    )
    export_p.add_argument(
        "--size-cols",
        nargs=3,
        metavar=("DX_COL", "DY_COL", "DZ_COL"),
        default=None,
        help="CSV dimension columns for row sizes; converted to lower/upper irregular extents",
    )
    export_p.add_argument(
        "--extent-cols",
        nargs=6,
        metavar=("LOWER_X", "UPPER_X", "LOWER_Y", "UPPER_Y", "LOWER_Z", "UPPER_Z"),
        default=None,
        help="CSV lower/upper extent columns, in X/Y/Z lower/upper order",
    )
    export_p.add_argument(
        "--null-float",
        type=float,
        default=-99.0,
        help="Null value for numeric variables",
    )
    export_p.add_argument(
        "--index-tolerance",
        type=float,
        default=1e-3,
        help="Tolerance for coordinate-to-grid index alignment",
    )
    export_p.add_argument(
        "--regularize-to-base-block",
        action="store_true",
        help="Aggregate rows into explicit cell-size base blocks before writing dense BMF output",
    )
    export_p.add_argument("--summary-json", default=None, help="Write summary JSON")
    export_p.add_argument("--dry-run", action="store_true", help="Validate and summarize only")
    export_p.set_defaults(func=cmd_export)

    inspect_p = sub.add_parser("inspect", help="Inspect BMF binary signatures")
    inspect_p.add_argument("input_bmf", help="Input BMF path")
    inspect_p.add_argument("--scan-bytes", type=int, default=2 * 1024 * 1024, help="Bytes to scan from start")
    inspect_p.add_argument("--min-string-len", type=int, default=4, help="Minimum printable string length")
    inspect_p.add_argument("--max-strings", type=int, default=30, help="Maximum strings in report")
    inspect_p.add_argument("--top-chunks", type=int, default=20, help="Top repeated 4-byte chunks")
    inspect_p.add_argument("--output-json", default=None, help="Optional JSON report path")
    inspect_p.set_defaults(func=cmd_inspect)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
