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
        return pd.concat(chunks, ignore_index=True)
    return pd.DataFrame()


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
    header = max(int(header_line or 1), 1) - 1
    read_kwargs = {
        "header": header,
        "usecols": list(usecols) if usecols else None,
        "dtype": str,
    }
    total_rows = None
    if progress_callback is not None:
        total_rows = _count_csv_data_rows(path, header_line, progress_callback, progress_start, min(progress_start + 3, progress_end))
        read_progress_start = min(progress_start + 3, progress_end)
    else:
        read_progress_start = progress_start
    if delimiter:
        try:
            if progress_callback is not None:
                return _read_csv_with_chunk_progress(
                    path,
                    total_rows,
                    progress_callback=progress_callback,
                    progress_start=read_progress_start,
                    progress_end=progress_end,
                    sep=delimiter,
                    low_memory=True,
                    memory_map=True,
                    **read_kwargs,
                )
            return pd.read_csv(path, sep=delimiter, low_memory=True, memory_map=True, **read_kwargs)
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
                    sep=delimiter,
                    engine="python",
                    **read_kwargs,
                )
            return pd.read_csv(path, sep=delimiter, engine="python", **read_kwargs)
    try:
        if progress_callback is not None:
            return _read_csv_with_chunk_progress(
                path,
                total_rows,
                progress_callback=progress_callback,
                progress_start=read_progress_start,
                progress_end=progress_end,
                sep=None,
                engine="python",
                **read_kwargs,
            )
        return pd.read_csv(path, sep=None, engine="python", **read_kwargs)
    except Exception:
        if progress_callback is not None:
            return _read_csv_with_chunk_progress(
                path,
                total_rows,
                progress_callback=progress_callback,
                progress_start=read_progress_start,
                progress_end=progress_end,
                low_memory=True,
                memory_map=True,
                **read_kwargs,
            )
        return pd.read_csv(path, low_memory=True, memory_map=True, **read_kwargs)


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


def _format_numeric_vector(values: Sequence[float] | np.ndarray, precision: int = 6) -> str:
    formatted = []
    for value in np.asarray(values, dtype=float):
        if not math.isfinite(float(value)):
            formatted.append(str(float(value)))
            continue
        text = f"{float(value):.{precision}f}".rstrip("0").rstrip(".")
        formatted.append(text or "0")
    return "[" + ", ".join(formatted) + "]"


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
    if value_cols:
        cols = [c for c in value_cols if c in df.columns and c not in xyz_cols]
        if not cols:
            raise ValueError("No valid value columns were provided.")
        return cols
    return [c for c in df.columns if c not in xyz_cols]


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
) -> np.ndarray:
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    row_count = int(np.prod(dims, dtype=np.int64))
    required_bytes = row_count * np.dtype(dtype).itemsize
    if required_bytes > MAX_DENSE_EXPORT_BYTES:
        cell = np.asarray(prepared.get("cell", [1.0, 1.0, 1.0]), dtype=float)
        raise MemoryError(
            "BMF export would require a dense grid that is too large for this writer: "
            f"dimensions={tuple(int(value) for value in dims)}, cells={row_count:,}, "
            f"field_size={_format_byte_size(required_bytes)}, cell_size={cell.tolist()}. "
            "This usually means the inferred cell size is too small for the CSV coordinates, or the CSV contains "
            "irregular/sub-block coordinates. Export a coarser regular-grid CSV or specify an explicit cell size."
        )
    linear = _tbms_linear_indices(prepared)
    dense = np.full(row_count, default_value, dtype=dtype)
    dense[linear] = np.asarray(series, dtype=dtype)
    return dense


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
        f"BMF export would create an oversized dense grid{backend_text} before writing any values. "
        f"Inferred dimensions={tuple(int(value) for value in dims)} ({row_count:,} cells) from "
        f"{valid_rows:,} valid CSV rows using cell_size={cell.tolist()}. "
        f"Estimated dense value storage for {field_count} field(s) is {_format_byte_size(estimated_bytes)}. "
        "The CSV likely has very fine coordinate spacing, irregular/sub-block coordinates, or an incorrect XYZ/cell-size setup. "
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
    progress_callback=None,
) -> None:
    if not selected_read_cols or not _backend_requires_dense_grid(backend):
        return
    row_count = _count_csv_data_rows(path, header_line, progress_callback, 3, 5)
    if row_count is None:
        return
    column_count = len(selected_read_cols)
    estimated_object_block_bytes = int(row_count) * int(column_count) * 8
    if estimated_object_block_bytes <= MAX_SELECTED_CSV_OBJECT_BYTES:
        return

    operation = "regularized dense BMF export" if regularize_to_base_block else "dense BMF export"
    raise MemoryError(
        f"The selected CSV columns are too wide for {operation} with backend {backend!r}. "
        f"The exporter would need to load {row_count:,} row(s) x {column_count:,} selected column(s) before writing, "
        f"which is at least {_format_byte_size(estimated_object_block_bytes)} just for pandas object references. "
        "Select fewer value columns and export in smaller batches, or use the tbms-experimental backend if you need a "
        "row-indexed Anterpolator container rather than a dense Vulcan-style BMF."
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
            "values": _tbms_dense_field_values(prepared, encoded, 0, np.dtype("u1")),
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
            "values": _tbms_dense_field_values(prepared, converted, int(null_float), np.dtype("<i4")),
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
            "values": _tbms_dense_field_values(prepared, converted, float(null_float), np.dtype("<f8")),
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
            "values": _tbms_dense_field_values(prepared, codes, 0, np.dtype("<i2")),
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
            "values": _tbms_dense_field_values(prepared, encoded, 0, np.dtype("u1")),
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
                    "values": _tbms_dense_field_values(prepared, converted, int(null_float), np.dtype("<i4")),
                }
        converted = numeric.fillna(float(null_float)).astype(np.float64)
        return {
            "field_type": "double",
            "default": float(null_float),
            "dtype": np.dtype("<f8"),
            "values": _tbms_dense_field_values(prepared, converted, float(null_float), np.dtype("<f8")),
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
        "values": _tbms_dense_field_values(prepared, codes, 0, np.dtype("<i2")),
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


def _iter_tbms_value_pages(values: np.ndarray, dtype: np.dtype, default_value: object) -> Iterable[bytes]:
    payload_size = 2048
    capacity = payload_size // dtype.itemsize
    page_count = max(1, int(math.ceil(len(values) / capacity)))
    for page_index in range(page_count):
        start = page_index * capacity
        stop = min(start + capacity, len(values))
        page_values = np.full(capacity, default_value, dtype=dtype)
        if stop > start:
            page_values[: stop - start] = values[start:stop].astype(dtype, copy=False)
        yield struct.pack("<II", 512, 0) + page_values.tobytes(order="C")


def _tbms_pointer_page(child_offsets: Sequence[int]) -> bytes:
    entries = [257] + [int(offset) for offset in child_offsets[:256]]
    if len(entries) < 257:
        entries.extend([0] * (257 - len(entries)))
    return struct.pack("<257Q", *entries)


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
    if int(prepared.get("duplicates") or 0) > 0:
        raise ValueError("tbms-config-text export requires unique grid cells; duplicate XYZ rows were found.")

    df: pd.DataFrame = prepared["df"]  # type: ignore[assignment]
    dims = np.asarray(prepared["dims"], dtype=np.int64)
    origin = np.asarray(prepared["origin"], dtype=float)
    cell = np.asarray(prepared["cell"], dtype=float)
    extents = np.asarray(prepared["extents"], dtype=float)
    row_count = int(np.prod(dims, dtype=np.int64))

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
            field_entries.append({"entry_key": f"var_{field_index}", "entry": entry})

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
        "is_irregular": 0,
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

    _emit_progress(progress_callback, 0, 100, "Preparing BMF export...")
    xyz_cols = [x_col, y_col, z_col]
    selected_read_cols = list(dict.fromkeys([*xyz_cols, *(value_cols or [])])) if value_cols else None
    selected_message = "selected columns" if selected_read_cols else "all columns"
    _validate_selected_csv_read_size(
        in_csv,
        header_line,
        selected_read_cols,
        backend,
        regularize_to_base_block,
        progress_callback=progress_callback,
    )
    _emit_progress(progress_callback, 5, 100, f"Reading CSV {selected_message}...")
    df = _auto_read_csv(
        in_csv,
        delimiter=delimiter,
        header_line=header_line,
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

    normalized_value_exceptions = _normalize_value_exceptions(value_exceptions)
    regularization_summary = {"enabled": False}
    if regularize_to_base_block:
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

    _emit_progress(progress_callback, 36, 100, "Preparing regular BMF grid from XYZ columns...")
    prepared = _prepare_grid(
        df=df,
        xyz_cols=xyz_cols,
        cell_size=cell_size,
        origin=origin,
        index_tolerance=index_tolerance,
    )
    dims_text = " x ".join(str(int(value)) for value in np.asarray(prepared["dims"], dtype=int))
    _emit_progress(
        progress_callback,
        45,
        100,
        f"Grid prepared: {dims_text} cells from {len(prepared['df']):,} valid CSV rows.",
    )

    chosen_value_cols = _classify_columns(prepared["df"], xyz_cols, value_cols)
    _emit_progress(progress_callback, 47, 100, f"Selected {len(chosen_value_cols)} BMF value field(s).")
    if backend not in SUPPORTED_EXPORT_BACKENDS:
        raise ValueError(f"Unsupported backend: {backend}")
    if _backend_requires_dense_grid(backend):
        _validate_dense_export_size(prepared, chosen_value_cols, backend=backend)
    summary = {
        "input_csv": str(in_csv),
        "output_bmf": str(out_bmf),
        "backend": backend,
        "rows": int(len(prepared["df"])),
        "value_columns": chosen_value_cols,
        "value_exceptions": normalized_value_exceptions,
        "regularization": regularization_summary,
        "grid": {
            "origin": [float(x) for x in np.asarray(prepared["origin"], dtype=float)],
            "cell_size": [float(x) for x in np.asarray(prepared["cell"], dtype=float)],
            "dimensions": [int(x) for x in np.asarray(prepared["dims"], dtype=int)],
            "extents": [float(x) for x in np.asarray(prepared["extents"], dtype=float)],
            "duplicate_grid_rows": int(prepared["duplicates"]),
            "max_index_error": float(prepared["max_index_error"]),
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

    df = _auto_read_csv(in_csv)
    xyz_cols = [args.x_col, args.y_col, args.z_col]
    missing = [c for c in xyz_cols if c not in df.columns]
    if missing:
        print(f"Missing coordinate columns: {missing}", file=sys.stderr)
        return 2

    regularization_summary = {"enabled": False}
    if getattr(args, "regularize_to_base_block", False):
        df, regularization_summary = _regularize_to_base_cell_grid(
            df=df,
            xyz_cols=xyz_cols,
            cell_size=args.cell_size,
            origin=args.origin,
            value_cols=args.value_cols,
            column_types=None,
        )

    prepared = _prepare_grid(
        df=df,
        xyz_cols=xyz_cols,
        cell_size=args.cell_size,
        origin=args.origin,
        index_tolerance=args.index_tolerance,
    )

    value_cols = _classify_columns(prepared["df"], xyz_cols, args.value_cols)

    summary = {
        "input_csv": str(in_csv),
        "output_bmf": str(out_bmf),
        "backend": args.backend,
        "rows": int(len(prepared["df"])),
        "value_columns": value_cols,
        "regularization": regularization_summary,
        "grid": {
            "origin": [float(x) for x in np.asarray(prepared["origin"], dtype=float)],
            "cell_size": [float(x) for x in np.asarray(prepared["cell"], dtype=float)],
            "dimensions": [int(x) for x in np.asarray(prepared["dims"], dtype=int)],
            "extents": [float(x) for x in np.asarray(prepared["extents"], dtype=float)],
            "duplicate_grid_rows": int(prepared["duplicates"]),
            "max_index_error": float(prepared["max_index_error"]),
        },
    }

    if args.summary_json:
        summary_path = Path(args.summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))

    if args.dry_run:
        print("Dry-run mode: no BMF file written.")
        return 0

    try:
        out_bmf.parent.mkdir(parents=True, exist_ok=True)
        if args.backend == "vulcan":
            _export_with_vulcan(
                prepared=prepared,
                xyz_cols=xyz_cols,
                value_cols=value_cols,
                out_bmf=out_bmf,
                null_float=args.null_float,
            )
        elif args.backend == "tbms-config-text":
            config_summary = _write_tbms_config_text(
                prepared=prepared,
                xyz_cols=xyz_cols,
                value_cols=value_cols,
                out_bmf=out_bmf,
                null_float=args.null_float,
            )
            print(json.dumps(config_summary, indent=2))
        else:
            experimental_summary = _write_tbms_experimental(
                prepared=prepared,
                xyz_cols=xyz_cols,
                value_cols=value_cols,
                out_bmf=out_bmf,
                null_float=args.null_float,
            )
            print(json.dumps(experimental_summary, indent=2))
    except Exception as exc:
        print(f"BMF export failed: {exc}", file=sys.stderr)
        return 1

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
