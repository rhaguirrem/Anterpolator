import colorsys
import os
import threading
import time
import traceback
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import taichi as ti


TI_INITIALIZED = False


def ensure_taichi_initialized():
    global TI_INITIALIZED
    if TI_INITIALIZED:
        return
    try:
        ti.init(arch=ti.gpu)
    except Exception:
        ti.init(arch=ti.cpu)
    TI_INITIALIZED = True


def detect_csv_delimiter(path):
    if not os.path.isfile(path):
        return ','
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = []
            for _ in range(10):
                try:
                    line = next(f)
                except StopIteration:
                    break
                if line.strip() and not line.lstrip().startswith(('#', '//')):
                    lines.append(line)
                if len(lines) >= 3:
                    break
            sample = ''.join(lines)
    except StopIteration:
        sample = ''
    counts = {d: sample.count(d) for d in [',', ';', '\t', '|']}
    delim = max(counts, key=counts.get)
    return delim if counts[delim] > 0 else ','


def parse_header_line(path, delimiter, line_number):
    if not os.path.isfile(path):
        raise ValueError(f"File not found: {path}")
    if line_number < 1:
        raise ValueError("Header line number must be >= 1")
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        for current, line in enumerate(f, start=1):
            if current == line_number:
                raw = line.strip('\n\r')
                tokens = [t.strip() for t in raw.split(delimiter)]
                tokens = [t for t in tokens if t != '']
                if not tokens:
                    raise ValueError(
                        f"Parsed header line {line_number} in '{os.path.basename(path)}' produced no tokens."
                    )
                return tokens
    raise ValueError(f"Header line {line_number} exceeds total lines in file '{os.path.basename(path)}'.")


def build_unique_column_names(headers):
    name_counts = {}
    final_names = []
    for header in headers:
        key = header if header else 'Unnamed'
        count = name_counts.get(key, 0)
        final_names.append(f"{key}_{count}" if count > 0 else key)
        name_counts[key] = count + 1
    return final_names


def read_csv_with_selected_header(path, delimiter, header_line, expected_min_cols=1):
    headers = parse_header_line(path, delimiter, header_line)
    final_names = build_unique_column_names(headers)
    df = pd.read_csv(
        path,
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=header_line - 1,
        comment='#',
    )
    if df.shape[0] and all(
        str(df.iloc[0, i]).strip() == final_names[i]
        for i in range(min(len(final_names), df.shape[1]))
    ):
        df = df.iloc[1:].reset_index(drop=True)
    empty_cols = [
        c for c in df.columns
        if df[c].isna().all() or (df[c].astype(str).str.strip() == '').all()
    ]
    if empty_cols:
        df = df.drop(columns=empty_cols)
    if df.shape[1] < expected_min_cols:
        raise ValueError(f"File '{path}' has fewer than {expected_min_cols} non-empty columns after cleanup.")
    df._detected_delimiter = delimiter
    return df, final_names


def read_autodetect_csv(path, min_cols=1, forced_delimiter=None):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")
    delim = forced_delimiter if forced_delimiter else detect_csv_delimiter(path)
    try:
        df = pd.read_csv(path, delimiter=delim, comment='#')
    except Exception:
        alt = ';' if delim == ',' else (',' if delim == ';' else None)
        if alt:
            df = pd.read_csv(path, delimiter=alt, comment='#')
            delim = alt
        else:
            raise
    if df.shape[1] == 1:
        text_sample = str(df.iloc[0, 0]) if len(df) else ''
        candidates = ['; ', ',', '\t', '|']
        for cand in candidates:
            if cand.strip() in text_sample and cand != delim:
                df_try = pd.read_csv(path, delimiter=cand, comment='#')
                if df_try.shape[1] > 1:
                    df = df_try
                    delim = cand
                    break
    empty_cols = [
        c for c in df.columns
        if df[c].isna().all() or (df[c].astype(str).str.strip() == '').all()
    ]
    if empty_cols:
        df = df.drop(columns=empty_cols)
    if df.shape[1] < min_cols:
        raise ValueError(f"File '{path}' has fewer than {min_cols} non-empty columns after cleanup.")
    df._detected_delimiter = delim
    return df


def load_samples_dataframe(samples_file, samples_delimiter=None, samples_header_line=1,
                           sample_x_col=None, sample_y_col=None, sample_z_col=None, sample_value_col=None):
    explicit_mapping = all([sample_x_col, sample_y_col, sample_z_col, sample_value_col])
    if explicit_mapping:
        delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
        df, _ = read_csv_with_selected_header(
            samples_file,
            delimiter,
            samples_header_line or 1,
            expected_min_cols=4,
        )
        rename_map = {
            sample_x_col: 'x',
            sample_y_col: 'y',
            sample_z_col: 'z',
            sample_value_col: 'Value',
        }
        df = df.rename(columns=rename_map)
        return df, rename_map

    if samples_header_line and samples_header_line != 1 and samples_delimiter:
        df, _ = read_csv_with_selected_header(
            samples_file,
            samples_delimiter,
            samples_header_line,
            expected_min_cols=4,
        )
    else:
        df = read_autodetect_csv(samples_file, forced_delimiter=samples_delimiter)
    return df, None


def load_lfc_colormap(lfc_file):
    if not lfc_file or not os.path.isfile(lfc_file):
        return [], [], []
    try:
        tree = ET.parse(lfc_file)
        root = tree.getroot()
        ranges = root.find('ranges')
        if ranges is None:
            return [], [], []
        entries = ranges.findall('entry')
        colormap = []
        thresholds = []
        boundaries = []
        labels = []
        for entry in entries:
            colour_elem = entry.find('colour')
            if colour_elem is None:
                continue
            colour_text = colour_elem.text.strip()
            colour = tuple(map(float, colour_text.split()))
            colormap.append(colour)
            value_elem = entry.find('value')
            min_elem = entry.find('min_value')
            max_elem = entry.find('max_value')
            label_elem = entry.find('label')
            lbl = label_elem.text.strip() if label_elem is not None and label_elem.text else None
            labels.append(lbl)
            if value_elem is not None:
                try:
                    thresholds.append(float(value_elem.text.strip()))
                except Exception:
                    pass
            else:
                try:
                    vmin = float(min_elem.text.strip()) if min_elem is not None else None
                except Exception:
                    vmin = None
                try:
                    vmax = float(max_elem.text.strip()) if max_elem is not None else None
                except Exception:
                    vmax = None
                boundaries.append((vmin, vmax))
        end_colour_elem = root.find('end_colour')
        if end_colour_elem is not None:
            end_colour_text = end_colour_elem.find('colour').text.strip()
            end_colour = tuple(map(float, end_colour_text.split()))
            colormap.append(end_colour)
        flat = [c for tpl in colormap for c in tpl]
        max_val = max(flat) if flat else 1.0
        scale = 255.0 if max_val > 1.0 else 1.0
        normalized = []
        for rgba in colormap:
            if len(rgba) == 3:
                rgba = (*rgba, 1.0)
            elif len(rgba) > 4:
                rgba = rgba[:4]
            r, g, b, a = rgba
            normalized.append((r / scale, g / scale, b / scale, a if scale == 1.0 else a / scale if a > 1 else a))
        if thresholds:
            numeric_edges = sorted(set(thresholds))
        else:
            numeric_edges = []
            for vmin, vmax in boundaries:
                if vmin is not None:
                    if not numeric_edges:
                        numeric_edges.append(vmin)
                if vmax is not None:
                    numeric_edges.append(vmax)
            numeric_edges = sorted(set(numeric_edges))
        final_labels = []
        for i, lbl in enumerate(labels):
            if lbl:
                final_labels.append(lbl)
            else:
                if i < len(boundaries):
                    vmin, vmax = boundaries[i]
                    final_labels.append(f"{vmin} - {vmax}")
                else:
                    final_labels.append(f"Class {i}")
        return normalized, numeric_edges, final_labels
    except Exception:
        return [], [], []


def build_lfc_tick_labels(colors, bins):
    if not colors or not bins:
        return []
    bins_np = np.array(sorted(bins), dtype=np.float32)
    n_colors = len(colors)
    threshold_style = (n_colors == len(bins_np) + 1)
    tick_labels = []
    if threshold_style:
        for i in range(n_colors):
            if i == 0:
                rng = f"< {bins_np[0]}"
            elif i < len(bins_np):
                rng = f"{bins_np[i-1]} - {bins_np[i]}"
            else:
                rng = f">= {bins_np[-1]}"
            tick_labels.append(rng)
    elif len(bins_np) == n_colors + 1:
        for i in range(n_colors):
            tick_labels.append(f"{bins_np[i]} - {bins_np[i+1]}")
    else:
        for i in range(n_colors):
            tick_labels.append(f"Class {i}")
    return tick_labels


def map_values_to_colors(values, colors, bins):
    if values is None or len(values) == 0:
        return np.zeros((0, 3), dtype=np.float32)
    if not colors:
        return np.full((len(values), 3), 0.6, dtype=np.float32)
    colors_np = np.array(colors, dtype=np.float32)
    if colors_np.shape[1] == 4:
        colors_np = colors_np[:, :3]

    values = np.asarray(values, dtype=np.float32)
    vmin, vmax = float(np.min(values)), float(np.max(values))
    if bins:
        bins_np = np.array(bins, dtype=np.float32)
        n_colors = len(colors_np)
        if n_colors == len(bins_np) + 1:
            idx = np.searchsorted(bins_np, values, side='right')
        elif len(bins_np) == n_colors + 1:
            idx = np.digitize(values, bins_np, right=False) - 1
        else:
            idx = np.digitize(values, bins_np, right=False) - 1
        idx = np.clip(idx, 0, n_colors - 1)
        return colors_np[idx]
    if vmax <= vmin:
        return np.tile(colors_np[0], (len(values), 1))
    scaled = (values - vmin) / (vmax - vmin)
    idx = np.clip((scaled * (len(colors_np) - 1)).astype(int), 0, len(colors_np) - 1)
    return colors_np[idx]


def prepare_points(df, value_col='Value'):
    for c in ['x', 'y', 'z']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df = df.dropna(subset=['x', 'y', 'z'])
    values = None
    if value_col in df.columns:
        df[value_col] = pd.to_numeric(df[value_col], errors='coerce')
        df = df.dropna(subset=[value_col])
        values = df[value_col].values.astype(np.float32)
    points = df[['x', 'y', 'z']].values.astype(np.float32)
    return points, values, df


def build_domain_color_map(domains):
    unique_domains = [str(domain).strip() for domain in domains if str(domain).strip()]
    unique_domains = sorted(set(unique_domains))
    color_map = {}
    golden_ratio = 0.61803398875
    hue = 0.05
    for domain in unique_domains:
        hue = (hue + golden_ratio) % 1.0
        sat = 0.65
        val = 0.95
        color_map[domain] = np.array(colorsys.hsv_to_rgb(hue, sat, val), dtype=np.float32)
    if 'Undomained' not in color_map:
        color_map['Undomained'] = np.array([0.5, 0.5, 0.5], dtype=np.float32)
    return color_map


class TaichiInterpolationViewer:
    def __init__(
        self,
        sample_points,
        sample_values,
        blocks_data,
        block_size,
        value_filter=0.0,
        lfc_colors=None,
        lfc_bins=None,
        lfc_tick_labels=None,
        sample_domains=None,
        config=None,
        window_title='Anterpolator Taichi Viewer',
    ):
        ensure_taichi_initialized()
        self.sample_points = np.asarray(sample_points, dtype=np.float32)
        self.sample_values = None if sample_values is None else np.asarray(sample_values, dtype=np.float32)
        self.sample_domains = None if sample_domains is None else np.asarray(sample_domains, dtype=object)
        self.blocks_data = blocks_data
        self.block_size = np.asarray(block_size, dtype=np.float32)
        self.block_info = getattr(blocks_data, '_block_info', {})
        self.min_bounds = np.asarray(self.block_info.get('min_bounds', np.zeros(3)), dtype=np.float32)
        self.positions_are_world = bool(self.block_info.get('positions_are_world', False))
        self.value_filter = float(value_filter)
        self.config = config or {}
        self.window_title = window_title
        self.lfc_colors = lfc_colors or []
        self.lfc_bins = lfc_bins or []
        self.lfc_tick_labels = lfc_tick_labels or build_lfc_tick_labels(self.lfc_colors, self.lfc_bins)
        self.domain_data_available = self._has_domain_data()
        self.color_mode = 'domain' if self._should_default_to_domain_mode() else 'value'
        self.show_samples = True
        self.show_blocks = True
        self.show_filled_block_faces = True
        self.show_legend_overlay = True
        self._worker = None
        self._pending_result = None
        self._result_lock = threading.Lock()
        self._raw_snapshot = None
        self._last_legend_mode = None
        self._last_status = ''
        self._render_center = np.zeros(3, dtype=np.float32)
        self._render_scale = 1.0
        self._world_span = 1.0
        self._sample_count = 0
        self._block_count = 0
        self._use_flat_block_lighting = True

        self.sample_fields = (None, None)
        self.block_fields = (None, None)
        self.block_mesh_fields = (None, None, None)
        self.block_line_fields = (None, None)

        self._domain_color_map = build_domain_color_map(self._collect_known_domains()) if self.domain_data_available else {}
        self._refresh_snapshot()

    def _extract_interpolators(self):
        blocks = self.blocks_data
        interpolators = []
        if hasattr(blocks, '_interpolators') and blocks._interpolators:
            for entry in blocks._interpolators.values():
                if isinstance(entry, list):
                    if entry:
                        interpolators.append(entry[0])
                else:
                    interpolators.append(entry)
        elif hasattr(blocks, '_ant_colony'):
            interpolators = [blocks._ant_colony]
        return interpolators

    def _has_domain_data(self):
        if self.sample_domains is not None and len(self.sample_domains) > 0:
            return True
        if self.block_info.get('domain_mapping'):
            return True
        for interp in self._extract_interpolators():
            if getattr(interp, 'domain_mapping', None):
                return True
        return False

    def _should_default_to_domain_mode(self):
        for interp in self._extract_interpolators():
            if str(getattr(interp, 'interpolation_target', 'value')).strip().lower() == 'domain':
                return True
        return False

    def _collect_known_domains(self):
        domains = []
        if self.sample_domains is not None:
            domains.extend([str(v).strip() for v in self.sample_domains if str(v).strip()])
        if self.block_info.get('domain_mapping'):
            domains.extend(str(v).strip() for v in self.block_info['domain_mapping'].values() if str(v).strip())
        for interp in self._extract_interpolators():
            domain_mapping = getattr(interp, 'domain_mapping', None)
            if domain_mapping:
                domains.extend(str(v).strip() for v in domain_mapping.values() if str(v).strip())
        return domains

    def _grid_to_world(self, pos):
        if self.positions_are_world:
            return np.asarray(pos, dtype=np.float32)
        return self.min_bounds + np.asarray(pos, dtype=np.float32) * self.block_size

    def _collect_block_snapshot(self):
        sample_blocks = dict(getattr(self.blocks_data, '_sample_blocks', {}) or {})
        merged_values = dict(sample_blocks)
        merged_domains = dict(self.block_info.get('domain_mapping', {}) or {})
        for interp in self._extract_interpolators():
            try:
                merged_values.update(interp.get_interpolated_values())
            except Exception:
                blocks = getattr(interp, 'blocks', {}) or {}
                if isinstance(blocks, dict):
                    for pos, block in blocks.items():
                        if hasattr(block, 'value'):
                            merged_values[pos] = block.value
                        elif isinstance(block, dict) and 'value' in block:
                            merged_values[pos] = block['value']
            domain_mapping = getattr(interp, 'domain_mapping', None)
            if domain_mapping:
                merged_domains.update(domain_mapping)

        positions = []
        values = []
        domains = []
        is_sample = []
        for pos, value in merged_values.items():
            positions.append(self._grid_to_world(pos))
            values.append(float(value))
            domains.append(str(merged_domains.get(pos, 'Undomained')))
            is_sample.append(pos in sample_blocks)

        return {
            'positions': np.asarray(positions, dtype=np.float32),
            'values': np.asarray(values, dtype=np.float32),
            'domains': np.asarray(domains, dtype=object),
            'is_sample': np.asarray(is_sample, dtype=bool),
        }

    def _sample_colors(self):
        if self.color_mode == 'domain' and self.sample_domains is not None and len(self.sample_domains) == len(self.sample_points):
            colors = [self._domain_color_map.get(str(domain).strip() or 'Undomained', self._domain_color_map['Undomained']) for domain in self.sample_domains]
            return np.asarray(colors, dtype=np.float32)
        if self.sample_values is None:
            return np.full((len(self.sample_points), 3), 0.6, dtype=np.float32)
        return map_values_to_colors(self.sample_values, self.lfc_colors, self.lfc_bins)

    def _block_display_arrays(self):
        snapshot = self._raw_snapshot
        if snapshot is None or len(snapshot['positions']) == 0:
            return np.zeros((0, 3), dtype=np.float32), np.zeros((0, 3), dtype=np.float32)

        if self.color_mode == 'domain' and self.domain_data_available:
            mask = np.ones(len(snapshot['positions']), dtype=bool)
            colors = [self._domain_color_map.get(str(domain).strip() or 'Undomained', self._domain_color_map['Undomained']) for domain in snapshot['domains']]
            return snapshot['positions'][mask], np.asarray(colors, dtype=np.float32)[mask]

        values = snapshot['values']
        mask = snapshot['is_sample'] | (values >= self.value_filter)
        colors = map_values_to_colors(values, self.lfc_colors, self.lfc_bins)
        return snapshot['positions'][mask], colors[mask]

    def _build_fields(self, points, colors):
        if points is None or len(points) == 0:
            return None, None
        pos_field = ti.Vector.field(3, dtype=ti.f32, shape=len(points))
        col_field = ti.Vector.field(3, dtype=ti.f32, shape=len(points))
        pos_field.from_numpy(np.asarray(points, dtype=np.float32))
        col_field.from_numpy(np.asarray(colors, dtype=np.float32))
        return pos_field, col_field

    def _build_block_line_fields(self, block_centers, block_colors):
        if block_centers is None or len(block_centers) == 0:
            return None, None

        half_size = (self.block_size.astype(np.float32) * self._render_scale) / 2.0
        corners_template = np.array([
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ], dtype=np.float32)
        edge_pairs = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]

        vertex_count = len(block_centers) * len(edge_pairs) * 2
        vertices = np.empty((vertex_count, 3), dtype=np.float32)
        colors = np.empty((vertex_count, 3), dtype=np.float32)

        cursor = 0
        scaled_corners = corners_template * half_size
        for center, color in zip(block_centers, block_colors):
            block_corners = scaled_corners + center
            for start_idx, end_idx in edge_pairs:
                vertices[cursor] = block_corners[start_idx]
                vertices[cursor + 1] = block_corners[end_idx]
                colors[cursor] = color
                colors[cursor + 1] = color
                cursor += 2

        return self._build_fields(vertices, colors)

    def _build_block_mesh_fields(self, block_centers, block_colors):
        if block_centers is None or len(block_centers) == 0:
            return None, None, None

        half_size = (self.block_size.astype(np.float32) * self._render_scale) / 2.0
        corners_template = np.array([
            [-1.0, -1.0, -1.0],
            [1.0, -1.0, -1.0],
            [1.0, 1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 1.0, 1.0],
        ], dtype=np.float32)
        triangle_indices = np.array([
            0, 1, 2, 0, 2, 3,
            4, 6, 5, 4, 7, 6,
            0, 4, 5, 0, 5, 1,
            1, 5, 6, 1, 6, 2,
            2, 6, 7, 2, 7, 3,
            3, 7, 4, 3, 4, 0,
        ], dtype=np.int32)

        vertex_count = len(block_centers) * 8
        index_count = len(block_centers) * len(triangle_indices)
        vertices = np.empty((vertex_count, 3), dtype=np.float32)
        colors = np.empty((vertex_count, 3), dtype=np.float32)
        indices = np.empty(index_count, dtype=np.int32)

        scaled_corners = corners_template * half_size
        for block_idx, (center, color) in enumerate(zip(block_centers, block_colors)):
            vertex_offset = block_idx * 8
            index_offset = block_idx * len(triangle_indices)
            block_corners = scaled_corners + center
            vertices[vertex_offset:vertex_offset + 8] = block_corners
            # Taichi scene.mesh does not expose alpha; lighten faces to approximate translucency.
            face_color = 0.55 * np.asarray(color, dtype=np.float32) + 0.45
            colors[vertex_offset:vertex_offset + 8] = np.clip(face_color, 0.0, 1.0)
            indices[index_offset:index_offset + len(triangle_indices)] = triangle_indices + vertex_offset

        vertex_field = ti.Vector.field(3, dtype=ti.f32, shape=vertex_count)
        color_field = ti.Vector.field(3, dtype=ti.f32, shape=vertex_count)
        index_field = ti.field(dtype=ti.i32, shape=index_count)
        vertex_field.from_numpy(vertices)
        color_field.from_numpy(colors)
        index_field.from_numpy(indices)
        return vertex_field, color_field, index_field

    def _update_render_transform(self, all_points):
        if all_points is None or len(all_points) == 0:
            self._render_center = np.zeros(3, dtype=np.float32)
            self._render_scale = 1.0
            self._world_span = 1.0
            return
        mins = np.min(all_points, axis=0).astype(np.float32)
        maxs = np.max(all_points, axis=0).astype(np.float32)
        self._render_center = ((mins + maxs) / 2.0).astype(np.float32)
        self._world_span = float(np.max(maxs - mins)) if np.any(maxs > mins) else 1.0
        self._render_scale = 1.6 / max(self._world_span, 1e-6)

    def _to_render_space(self, points):
        if points is None or len(points) == 0:
            return np.zeros((0, 3), dtype=np.float32)
        points = np.asarray(points, dtype=np.float32)
        return (points - self._render_center) * self._render_scale

    def _refresh_render_data(self):
        sample_colors = self._sample_colors()
        block_points, block_colors = self._block_display_arrays()
        if len(self.sample_points) and len(block_points):
            all_points = np.vstack([self.sample_points, block_points])
        elif len(self.sample_points):
            all_points = self.sample_points
        else:
            all_points = block_points
        self._update_render_transform(all_points)

        render_sample_points = self._to_render_space(self.sample_points)
        render_block_points = self._to_render_space(block_points)
        self.sample_fields = self._build_fields(render_sample_points, sample_colors)
        self.block_fields = self._build_fields(render_block_points, block_colors)
        self.block_mesh_fields = self._build_block_mesh_fields(render_block_points, block_colors)
        self.block_line_fields = self._build_block_line_fields(render_block_points, block_colors)
        self._sample_count = len(render_sample_points)
        self._block_count = len(render_block_points)

    def _refresh_snapshot(self):
        self._raw_snapshot = self._collect_block_snapshot()
        self._refresh_render_data()
        self._print_legend_if_needed(force=True)

    def _run_interpolators_once(self):
        dims = tuple(self.block_info.get('dims', (0, 0, 0)))
        changes_made = False
        for interp in self._extract_interpolators():
            if not hasattr(interp, 'next_block_id'):
                existing_ids = []
                blocks = getattr(interp, 'blocks', {}) or {}
                if isinstance(blocks, dict):
                    for block in blocks.values():
                        if hasattr(block, 'block_id'):
                            existing_ids.append(int(block.block_id))
                        elif isinstance(block, dict) and 'Block_ID' in block:
                            existing_ids.append(int(block['Block_ID']))
                interp.next_block_id = (max(existing_ids) + 1) if existing_ids else 1
            try:
                previous_blocks = len(interp.get_interpolated_values())
            except Exception:
                previous_blocks = len(getattr(interp, 'blocks', {}) or {})
            if interp.run_iteration(dims):
                changes_made = True
            else:
                try:
                    current_blocks = len(interp.get_interpolated_values())
                except Exception:
                    current_blocks = len(getattr(interp, 'blocks', {}) or {})
                if current_blocks > previous_blocks:
                    changes_made = True
            if self.config.get('fill_unvisited_domainwise', False):
                try:
                    if hasattr(interp, 'fill_unvisited_blocks_domainwise'):
                        created, assigned = interp.fill_unvisited_blocks_domainwise(dims)
                        if created or assigned:
                            changes_made = True
                except Exception:
                    pass
        return changes_made

    def request_update(self):
        if self._worker is not None and self._worker.is_alive():
            self._set_status('Update already running...')
            return

        def _worker():
            started = time.perf_counter()
            try:
                changed = self._run_interpolators_once()
                snapshot = self._collect_block_snapshot()
                result = {
                    'changed': changed,
                    'snapshot': snapshot,
                    'elapsed': time.perf_counter() - started,
                    'error': None,
                }
            except Exception:
                result = {
                    'changed': False,
                    'snapshot': None,
                    'elapsed': time.perf_counter() - started,
                    'error': traceback.format_exc(),
                }
            with self._result_lock:
                self._pending_result = result

        self._set_status('Running interpolation update...')
        self._worker = threading.Thread(target=_worker, daemon=True)
        self._worker.start()

    def _poll_worker_result(self):
        with self._result_lock:
            result = self._pending_result
            self._pending_result = None
        if result is None:
            return
        if result['error']:
            print(result['error'])
            self._set_status('Update failed. See console.')
            return
        if result['snapshot'] is not None:
            self._raw_snapshot = result['snapshot']
            self._refresh_render_data()
        if result['changed']:
            self._set_status(
                f"Updated in {result['elapsed']:.2f}s. Blocks now: {len(self._raw_snapshot['positions'])}."
            )
        else:
            self._set_status(f"No blocks to update ({result['elapsed']:.2f}s).")

    def _set_status(self, message):
        if message != self._last_status:
            self._last_status = message
            print(message)

    def _print_legend_if_needed(self, force=False):
        if not force and self._last_legend_mode == self.color_mode:
            return
        self._last_legend_mode = self.color_mode
        if self.color_mode == 'domain' and self.domain_data_available:
            print('\nDomain Legend:')
            for domain in sorted(self._domain_color_map.keys()):
                color = self._domain_color_map[domain]
                print(f"{domain}: RGB({color[0]:.3f},{color[1]:.3f},{color[2]:.3f})")
            return
        if self.lfc_colors and self.lfc_bins and self.lfc_tick_labels:
            print('\nLFC Discrete Legend:')
            colors = np.array(self.lfc_colors, dtype=np.float32)
            if colors.shape[1] == 4:
                colors = colors[:, :3]
            for idx, label in enumerate(self.lfc_tick_labels):
                if idx >= len(colors):
                    break
                color = colors[idx]
                print(f"Class {idx}: {label} | RGB({color[0]:.3f},{color[1]:.3f},{color[2]:.3f})")
        elif self.sample_values is not None and len(self.sample_values):
            print(
                f"\nContinuous Value Range: min={float(np.min(self.sample_values)):.3f} "
                f"max={float(np.max(self.sample_values)):.3f}"
            )

    def _get_overlay_entries(self, max_items=12):
        if self.color_mode == 'domain' and self.domain_data_available:
            items = sorted(self._domain_color_map.items(), key=lambda item: item[0])
            entries = [
                (domain, tuple(np.asarray(color, dtype=np.float32)[:3].tolist()))
                for domain, color in items[:max_items]
            ]
            extra_count = max(0, len(items) - max_items)
            return 'Domain Legend', entries, extra_count

        if self.lfc_colors and self.lfc_tick_labels:
            entries = []
            for idx, label in enumerate(self.lfc_tick_labels[:max_items]):
                if idx >= len(self.lfc_colors):
                    break
                color = tuple(np.asarray(self.lfc_colors[idx], dtype=np.float32)[:3].tolist())
                entries.append((label, color))
            extra_count = max(0, len(self.lfc_tick_labels) - len(entries))
            return 'Value Legend', entries, extra_count

        if self.sample_values is not None and len(self.sample_values):
            min_val = float(np.min(self.sample_values))
            max_val = float(np.max(self.sample_values))
            return 'Value Legend', [(f"Range: {min_val:.3f} to {max_val:.3f}", (0.4, 0.4, 0.4))], 0

        return 'Legend', [('No legend data available', (0.4, 0.4, 0.4))], 0

    def _overlay_text_color(self, preferred=None, min_luminance=0.72):
        if preferred is None:
            return (0.92, 0.92, 0.92)
        color = np.asarray(preferred, dtype=np.float32)[:3]
        luminance = float(0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2])
        if luminance >= min_luminance:
            return tuple(np.clip(color, 0.0, 1.0).tolist())
        boosted = color + (min_luminance - luminance)
        boosted = np.clip(boosted, min_luminance, 1.0)
        return tuple(boosted.tolist())

    def _draw_overlay(self, gui):
        if not self.show_legend_overlay:
            return

        legend_title, legend_entries, extra_count = self._get_overlay_entries()
        if self.color_mode == 'domain' and self.domain_data_available:
            base_height = 0.34
        else:
            base_height = 0.30

        with gui.sub_window('Viewer Legend', 0.015, 0.015, 0.32, base_height):
            gui.text('Anterpolator Taichi', color=(0.96, 0.96, 0.96))
            gui.text(f"Mode: {self.color_mode}", color=(0.90, 0.90, 0.90))
            gui.text(f"Samples: {'on' if self.show_samples else 'off'} | Blocks: {'on' if self.show_blocks else 'off'}", color=(0.84, 0.84, 0.84))
            gui.text(f"Filled faces: {'on' if self.show_filled_block_faces else 'off'}", color=(0.84, 0.84, 0.84))
            gui.text(f"Flat faces: {'on' if self._use_flat_block_lighting else 'off'}", color=(0.84, 0.84, 0.84))
            gui.text(f"Visible points: samples={self._sample_count} blocks={self._block_count}", color=(0.84, 0.84, 0.84))
            gui.text(f"Filter: {self.value_filter:.2f}", color=(0.84, 0.84, 0.84))
            if self._last_status:
                gui.text(self._last_status[:72], color=(0.78, 0.90, 0.98))

            gui.text('')
            gui.text(legend_title, color=(0.95, 0.95, 0.95))
            for label, color in legend_entries:
                gui.text(f"■ {label}", color=self._overlay_text_color(color))
            if extra_count:
                gui.text(f"... and {extra_count} more", color=(0.74, 0.74, 0.74))

            gui.text('')
            gui.text('Keys: b blocks, v samples, f faces, c mode, l legend, i update', color=(0.74, 0.74, 0.74))

    def toggle_color_mode(self):
        if not self.domain_data_available:
            self._set_status('No domain data available for domain coloring.')
            return
        self.color_mode = 'domain' if self.color_mode == 'value' else 'value'
        self._refresh_render_data()
        self._print_legend_if_needed(force=True)
        self._set_status(f"Color mode: {self.color_mode}")

    def run(self):
        window = ti.ui.Window(self.window_title, (1280, 720), vsync=True)
        canvas = window.get_canvas()
        canvas.set_background_color((1.0, 1.0, 1.0))
        camera = ti.ui.Camera()

        camera.position(0.0, -2.2, 1.4)
        camera.lookat(0.0, 0.0, 0.0)
        camera.up(0.0, 0.0, 1.0)

        min_block_dim_render = float(max(float(np.min(self.block_size)) * self._render_scale, 1e-6))
        sample_radius = float(np.clip(min_block_dim_render * 0.25, 0.0015, 0.03))
        self._set_status("Controls: RMB drag camera, b toggle blocks, v toggle samples, f toggle faces, i update, c toggle color mode, l toggle legend, [ ] value filter, q quit")
        self._print_legend_if_needed(force=True)

        while window.running:
            self._poll_worker_result()
            scene = window.get_scene()
            min_block_dim_render = float(max(float(np.min(self.block_size)) * self._render_scale, 1e-6))
            sample_radius = float(np.clip(min_block_dim_render * 0.25, 0.0015, 0.03))

            for event in window.get_events(ti.ui.PRESS):
                if event.key == 'b':
                    self.show_blocks = not self.show_blocks
                elif event.key == 'v':
                    self.show_samples = not self.show_samples
                elif event.key == 'f':
                    self.show_filled_block_faces = not self.show_filled_block_faces
                    self._set_status(f"Filled block faces: {'on' if self.show_filled_block_faces else 'off'}")
                elif event.key == 'i':
                    self.request_update()
                elif event.key == 'c':
                    self.toggle_color_mode()
                elif event.key == 'l':
                    self.show_legend_overlay = not self.show_legend_overlay
                    self._set_status(f"Legend overlay: {'on' if self.show_legend_overlay else 'off'}")
                elif event.key == '[':
                    self.value_filter -= 1.0
                    self._refresh_render_data()
                    self._set_status(f"value_filter={self.value_filter}")
                elif event.key == ']':
                    self.value_filter += 1.0
                    self._refresh_render_data()
                    self._set_status(f"value_filter={self.value_filter}")
                elif event.key in ('q', ti.ui.ESCAPE):
                    window.running = False

            camera.track_user_inputs(window, movement_speed=0.03, hold_key=ti.ui.RMB)
            scene.set_camera(camera)
            if self._use_flat_block_lighting:
                scene.ambient_light((1.0, 1.0, 1.0))
            else:
                scene.ambient_light((0.55, 0.55, 0.55))
                scene.point_light((0.0, 0.0, 2.5), (1.0, 1.0, 1.0))

            sample_points_field, sample_colors_field = self.sample_fields
            if self.show_samples and sample_points_field is not None:
                scene.particles(sample_points_field, per_vertex_color=sample_colors_field, radius=sample_radius)

            block_mesh_vertices, block_mesh_colors, block_mesh_indices = self.block_mesh_fields
            if self.show_blocks and self.show_filled_block_faces and block_mesh_vertices is not None:
                scene.mesh(
                    block_mesh_vertices,
                    indices=block_mesh_indices,
                    per_vertex_color=block_mesh_colors,
                    two_sided=True,
                )

            block_line_points_field, block_line_colors_field = self.block_line_fields
            if self.show_blocks and block_line_points_field is not None:
                scene.lines(block_line_points_field, width=1.25, per_vertex_color=block_line_colors_field)

            canvas.scene(scene)
            self._draw_overlay(window.GUI)
            window.show()


def create_viewer_from_files(
    samples_file,
    color_file=None,
    interpolation_file=None,
    blocks_file=None,
    samples_delimiter=None,
    blocks_delimiter=None,
    samples_header_line=1,
    blocks_header_line=1,
    sample_x_col=None,
    sample_y_col=None,
    sample_z_col=None,
    sample_value_col=None,
    block_size=(10.0, 10.0, 10.0),
    value_filter=0.0,
):
    colors, bins, _ = load_lfc_colormap(color_file)
    tick_labels = build_lfc_tick_labels(colors, bins)
    df_samples, sample_map = load_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        sample_x_col=sample_x_col,
        sample_y_col=sample_y_col,
        sample_z_col=sample_z_col,
        sample_value_col=sample_value_col,
    )
    if sample_map is None and all([sample_x_col, sample_y_col, sample_z_col, sample_value_col]):
        df_samples = df_samples.rename(columns={
            sample_x_col: 'x',
            sample_y_col: 'y',
            sample_z_col: 'z',
            sample_value_col: 'Value',
        })
    sample_domains = df_samples['Domain'].values if 'Domain' in df_samples.columns else None
    sample_points, sample_values, _ = prepare_points(df_samples, value_col='Value')

    block_positions = {}
    block_domains = {}
    if interpolation_file and os.path.isfile(interpolation_file):
        df_blocks = read_autodetect_csv(interpolation_file, forced_delimiter=blocks_delimiter)
        block_points, block_values, df_blocks = prepare_points(df_blocks, value_col='Value')
        for point, value in zip(block_points, block_values):
            block_positions[tuple(point)] = value
        if 'Domain' in df_blocks.columns:
            for point, domain in zip(block_points, df_blocks['Domain'].fillna('Undomained').astype(str)):
                block_domains[tuple(point)] = domain
    elif blocks_file and os.path.isfile(blocks_file):
        df_blocks = read_autodetect_csv(blocks_file, forced_delimiter=blocks_delimiter)
        block_points, block_values, df_blocks = prepare_points(df_blocks, value_col='Value')
        for point, value in zip(block_points, block_values if block_values is not None else np.zeros(len(block_points), dtype=np.float32)):
            block_positions[tuple(point)] = value
        if 'Domain' in df_blocks.columns:
            for point, domain in zip(block_points, df_blocks['Domain'].fillna('Undomained').astype(str)):
                block_domains[tuple(point)] = domain

    class SimpleBlocks:
        pass

    blocks_data = SimpleBlocks()
    blocks_data._sample_blocks = {
        tuple(point.astype(np.float32)): float(value)
        for point, value in zip(sample_points, sample_values)
    }
    blocks_data._block_info = {
        'min_bounds': np.min(sample_points, axis=0) if len(sample_points) else np.zeros(3, dtype=np.float32),
        'block_size': list(block_size),
        'domain_mapping': block_domains,
        'positions_are_world': True,
    }

    class StaticInterpolator:
        def __init__(self, values, domains):
            self._values = values
            self.domain_mapping = domains
            self.interpolation_target = 'value'

        def get_interpolated_values(self):
            return self._values

        def run_iteration(self, dims):
            return False

    blocks_data._ant_colony = StaticInterpolator(block_positions, block_domains)

    return TaichiInterpolationViewer(
        sample_points=sample_points,
        sample_values=sample_values,
        sample_domains=sample_domains,
        blocks_data=blocks_data,
        block_size=block_size,
        value_filter=value_filter,
        lfc_colors=colors,
        lfc_bins=bins,
        lfc_tick_labels=tick_labels,
        config={},
    )
