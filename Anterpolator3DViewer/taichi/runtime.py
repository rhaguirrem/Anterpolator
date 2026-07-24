import colorsys
import importlib
import os
import sys
import threading
import time
import traceback
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

from provenance_utils import (
    copy_interpolator_provenance,
    finalize_phase_provenance,
    seed_original_sample_provenance,
    snapshot_interpolator_state,
)

try:
    from pynput import mouse as pynput_mouse
except Exception:
    pynput_mouse = None


TI_INITIALIZED = False
TI_MODULE = None
VOLUME_RENDERER_MODULE = None
VOLUME_RENDERER_CANVAS_MODULE = None
ti = None


def get_taichi_module():
    global TI_MODULE, ti
    if TI_MODULE is not None:
        ti = TI_MODULE
        return TI_MODULE

    existing_module = sys.modules.get('taichi')
    if existing_module is not None and hasattr(existing_module, 'init'):
        TI_MODULE = existing_module
        ti = TI_MODULE
        return TI_MODULE

    runtime_dir = os.path.normcase(os.path.abspath(os.path.dirname(__file__)))
    viewer_dir = os.path.normcase(os.path.abspath(os.path.dirname(runtime_dir)))
    original_sys_path = list(sys.path)

    try:
        sys.modules.pop('taichi', None)
        sys.path = [
            entry for entry in original_sys_path
            if os.path.normcase(os.path.abspath(entry or os.curdir)) not in {runtime_dir, viewer_dir}
        ]
        module = importlib.import_module('taichi')
    except ModuleNotFoundError as exc:
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        version_hint = ''
        if sys.version_info >= (3, 14):
            version_hint = (
                f" The active interpreter is Python {python_version}; Taichi wheels may not be available for that "
                "version yet, so use Python 3.13 or earlier for the Taichi viewer."
            )
        raise ImportError(
            "The external 'taichi' package is not installed in the active Python environment. "
            "This viewer also has a local 'taichi' folder, so a plain 'import taichi' resolves to the local "
            "namespace package instead of the third-party library. Install the 'taichi' package in the "
            f"active environment before launching the Taichi viewer.{version_hint}"
        ) from exc
    finally:
        sys.path = original_sys_path

    if not hasattr(module, 'init'):
        raise ImportError(
            "Resolved module 'taichi' does not expose 'init'. The active environment is still loading the wrong "
            "module instead of the third-party Taichi package."
        )

    sys.modules['taichi'] = module
    TI_MODULE = module
    ti = TI_MODULE
    return TI_MODULE


def ensure_taichi_initialized():
    global TI_INITIALIZED
    if TI_INITIALIZED:
        return
    ti = get_taichi_module()
    try:
        ti.init(arch=ti.gpu)
    except Exception:
        ti.init(arch=ti.cpu)
    TI_INITIALIZED = True


def import_volume_renderer():
    global VOLUME_RENDERER_MODULE, VOLUME_RENDERER_CANVAS_MODULE
    if VOLUME_RENDERER_MODULE is not None and VOLUME_RENDERER_CANVAS_MODULE is not None:
        return VOLUME_RENDERER_MODULE, VOLUME_RENDERER_CANVAS_MODULE

    get_taichi_module()

    import taichi_volume_renderer as volume_renderer
    from taichi_volume_renderer import canvas as volume_canvas

    VOLUME_RENDERER_MODULE = volume_renderer
    VOLUME_RENDERER_CANVAS_MODULE = volume_canvas
    return VOLUME_RENDERER_MODULE, VOLUME_RENDERER_CANVAS_MODULE


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
        external_state_callback=None,
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
        self._overlay_mouse_capture = False
        self._overlay_ini_bounds = {}
        self._overlay_ini_last_check = 0.0
        self._overlay_ini_mtime = None
        self._overlay_ini_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'imgui.ini')
        self._worker = None
        self._pending_result = None
        self._result_lock = threading.Lock()
        self._result_condition = threading.Condition(self._result_lock)
        self._raw_snapshot = None
        self._last_legend_mode = None
        self._last_status = ''
        self._render_center = np.zeros(3, dtype=np.float32)
        self._render_scale = 1.0
        self._world_span = 1.0
        self._render_sample_points = np.zeros((0, 3), dtype=np.float32)
        self._render_block_points = np.zeros((0, 3), dtype=np.float32)
        self._sample_count = 0
        self._block_count = 0
        self._use_flat_block_lighting = True
        self._external_state_callback = external_state_callback
        self._external_state_last_check = 0.0
        self._external_state_poll_interval = 0.75
        self._external_wheel_listener = None
        self._external_wheel_lock = threading.Lock()
        self._external_wheel_delta = 0.0
        self.block_render_mode = self._normalize_block_render_mode(
            self.config.get('taichi_block_render_mode', 'mesh')
        )
        self._interpolation_groups = []
        self._interpolation_domain_index = 0
        self._interpolation_pass_index = 0
        self._interpolation_pass_iteration = 0
        self._interpolation_total_phases = 0
        self._interpolation_completed_phases = 0
        self._interpolation_complete = False
        self._phase_provenance_snapshot = None
        self._phase_provenance_interp = None

        self.sample_fields = (None, None)
        self.block_fields = (None, None)
        self.block_mesh_fields = (None, None, None)
        self.block_line_fields = (None, None)

        self._domain_color_map = build_domain_color_map(self._collect_known_domains()) if self.domain_data_available else {}
        self._reset_interpolation_progress()
        self._refresh_snapshot()

    def apply_viewer_state(self, state):
        previous_blocks_data = self.blocks_data
        previous_color_mode = self.color_mode
        self.sample_points = np.asarray(state['sample_points'], dtype=np.float32)
        sample_values = state.get('sample_values')
        sample_domains = state.get('sample_domains')
        self.sample_values = None if sample_values is None else np.asarray(sample_values, dtype=np.float32)
        self.sample_domains = None if sample_domains is None else np.asarray(sample_domains, dtype=object)
        self.blocks_data = state['blocks_data']
        self.block_size = np.asarray(state['block_size'], dtype=np.float32)
        self.value_filter = float(state.get('value_filter', 0.0))
        self.config = state.get('config', {}) or {}
        self.window_title = state.get('window_title', self.window_title)
        self.lfc_colors = state.get('lfc_colors', []) or []
        self.lfc_bins = state.get('lfc_bins', []) or []
        self.lfc_tick_labels = state.get('lfc_tick_labels', []) or build_lfc_tick_labels(self.lfc_colors, self.lfc_bins)
        self.block_info = getattr(self.blocks_data, '_block_info', {})
        self.min_bounds = np.asarray(self.block_info.get('min_bounds', np.zeros(3)), dtype=np.float32)
        self.positions_are_world = bool(self.block_info.get('positions_are_world', False))
        self.domain_data_available = self._has_domain_data()
        self._domain_color_map = build_domain_color_map(self._collect_known_domains()) if self.domain_data_available else {}
        if previous_color_mode == 'domain' and self.domain_data_available:
            self.color_mode = 'domain'
        elif previous_color_mode == 'value':
            self.color_mode = 'value'
        else:
            self.color_mode = 'domain' if self._should_default_to_domain_mode() else 'value'
        self.block_render_mode = self._normalize_block_render_mode(
            self.config.get('taichi_block_render_mode', self.block_render_mode)
        )
        self._last_legend_mode = None
        if previous_blocks_data is not self.blocks_data:
            self._reset_interpolation_progress()
        self._refresh_snapshot()

    def _maybe_reload_external_state(self):
        if self._external_state_callback is None:
            return None
        now = time.perf_counter()
        if (now - self._external_state_last_check) < self._external_state_poll_interval:
            return None
        self._external_state_last_check = now
        try:
            previous_mode = self.block_render_mode
            state = self._external_state_callback()
            if state is None:
                return None
            self.apply_viewer_state(state)
            if self.block_render_mode != previous_mode:
                return 'restart'
            self._set_status('Viewer reloaded from config.')
            return 'reloaded'
        except Exception:
            print('Failed to reload external viewer state:')
            traceback.print_exc()
            return None

    def _normalize_block_render_mode(self, mode):
        return 'mesh'

    def _volume_transparency_percent(self):
        try:
            return float(np.clip(self.config.get('taichi_transparency_percent', 50.0), 0.0, 100.0))
        except Exception:
            return 50.0

    def _volume_smoke_density_factor(self, metadata):
        if 'taichi_transparency_percent' in self.config:
            opacity = np.clip(1.0 - (self._volume_transparency_percent() / 100.0), 0.0, 1.0)
            if opacity <= 0.0:
                return 0.0

            # Convert the requested block opacity into a per-step extinction value
            # over one block thickness in the voxelized volume.
            block_density = max(float(self.config.get('taichi_volume_block_density', 1.0)), 1e-6)
            volume_shape = metadata.get('volume_shape') or (1, 1, 1)
            max_volume_dim = max(float(np.max(volume_shape)), 1.0)
            block_voxel_thickness = max(int(metadata.get('upsample', 1)), 1)
            max_supported_opacity = 0.995
            target_opacity = min(float(opacity), max_supported_opacity)
            target_transmittance = max(1.0 - target_opacity, 1e-6)
            per_step_transmittance = target_transmittance ** (1.0 / float(block_voxel_thickness))
            per_step_extinction = np.clip(1.0 - per_step_transmittance, 0.0, 0.999999)
            return float((per_step_extinction * max_volume_dim) / block_density)
        return float(
            self.config.get(
                'taichi_volume_smoke_density_factor',
                max(12.0, metadata['upsample'] * 8.0),
            )
        )

    def _volume_resolution(self):
        value = self.config.get('taichi_volume_resolution', (1280, 720))
        if isinstance(value, (list, tuple)) and len(value) == 2:
            width = max(int(value[0]), 320)
            height = max(int(value[1]), 240)
            return width, height
        return 1280, 720

    def _volume_background(self):
        value = self.config.get('taichi_volume_background', (0.97, 0.97, 0.99))
        if isinstance(value, (list, tuple)) and len(value) >= 3:
            return tuple(float(np.clip(channel, 0.0, 1.0)) for channel in value[:3])
        return (0.97, 0.97, 0.99)

    def _create_volume_display_window(self, volume_renderer, volume_canvas, resolution=None):
        smoke_density, smoke_color, metadata = self._build_volume_fields(volume_canvas)
        resolution = resolution or self._volume_resolution()
        window = volume_renderer.DisplayWindow(
            smoke_density=smoke_density,
            smoke_color=smoke_color,
            point_lights_pos=np.asarray([[0.0, 0.0, 5.0]], dtype=np.float32),
            point_lights_intensity=np.asarray([[50.0, 50.0, 50.0]], dtype=np.float32),
            lighting=bool(self.config.get('taichi_volume_lighting', True)),
            resolution=resolution,
            background=self._volume_background(),
            init_taichi=False,
            smoke_density_factor=self._volume_smoke_density_factor(metadata),
            ray_tracing_step_size_factor=float(self.config.get('taichi_volume_step_size_factor', 1.0)),
            light_ray_tracing_step_size_factor=float(self.config.get('taichi_volume_light_step_size_factor', 3.0)),
            ray_tracing_stop_threshold=float(self.config.get('taichi_volume_stop_threshold', 0.01)),
            ray_tracing_max_steps=int(self.config.get('taichi_volume_max_steps', 10000)),
        )
        window.scene.set_camera_phi(float(self.config.get('taichi_volume_camera_phi', -90.0)))
        window.scene.set_camera_theta(float(self.config.get('taichi_volume_camera_theta', 32.5)))
        window.scene.camera_distance = float(self.config.get('taichi_volume_camera_distance', 3.0))
        return window, metadata

    def _build_volume_fields(self, volume_canvas):
        block_points, block_colors = self._block_display_arrays()
        if block_points is None or len(block_points) == 0:
            raise ValueError('No visible blocks are available for volume rendering.')

        block_points = np.asarray(block_points, dtype=np.float32)
        block_colors = np.asarray(block_colors, dtype=np.float32)
        block_size = np.maximum(np.asarray(self.block_size, dtype=np.float32), 1e-6)
        origin = np.min(block_points, axis=0).astype(np.float32)
        block_indices = np.rint((block_points - origin) / block_size).astype(np.int32)
        block_indices -= np.min(block_indices, axis=0)

        padding = max(int(self.config.get('taichi_volume_padding_blocks', 1)), 0)
        scale = max(int(self.config.get('taichi_volume_upsample', 4)), 1)
        max_memory_mb = float(self.config.get('taichi_volume_max_memory_mb', 384.0))
        max_memory_bytes = max_memory_mb * 1024.0 * 1024.0
        grid_shape = np.max(block_indices, axis=0) + 1
        volume_shape = None
        estimated_bytes = 0.0

        while scale >= 1:
            scaled_shape = ((grid_shape + padding * 2) * scale).astype(np.int32)
            volume_shape = tuple(int(v) for v in scaled_shape.tolist())
            estimated_bytes = float(np.prod(np.asarray(volume_shape, dtype=np.int64))) * 16.0
            if estimated_bytes <= max_memory_bytes or scale == 1:
                break
            scale -= 1

        if volume_shape is None:
            raise ValueError('Failed to compute a valid volume shape.')
        if estimated_bytes > max_memory_bytes:
            raise ValueError(
                f"Volume would require about {estimated_bytes / (1024.0 * 1024.0):.1f} MB; "
                f"increase taichi_volume_max_memory_mb or reduce taichi_volume_upsample."
            )

        smoke_density = ti.field(dtype=ti.f32, shape=volume_shape)
        smoke_color = ti.Vector.field(3, dtype=ti.f32, shape=volume_shape)
        volume_canvas.clean(smoke_density, smoke_color)

        block_density = float(self.config.get('taichi_volume_block_density', 1.0))
        cube_scale = np.full(3, scale, dtype=np.float32)
        for block_index, color in zip(block_indices, block_colors):
            start = np.asarray((block_index + padding) * scale, dtype=np.float32)
            color_vec = np.asarray(color[:3], dtype=np.float32)
            volume_canvas.fill_rectangle(
                smoke_density,
                smoke_color,
                start,
                cube_scale,
                block_density,
                color_vec,
            )

        sample_count = 0
        if bool(self.config.get('taichi_volume_show_samples', True)) and len(self.sample_points):
            sample_positions = (
                ((self.sample_points.astype(np.float32) - origin) / block_size) + padding + 0.5
            ) * np.float32(scale)
            sample_positions = np.asarray(sample_positions, dtype=np.float32)
            sample_shape = np.asarray(volume_shape, dtype=np.float32)
            inside = np.all(sample_positions >= 0.0, axis=1) & np.all(sample_positions < sample_shape, axis=1)
            if np.any(inside):
                sample_colors = np.asarray(self._sample_colors(), dtype=np.float32)[inside]
                volume_canvas.draw_particles(
                    smoke_density,
                    np.asarray(sample_positions[inside], dtype=np.float32),
                    smoke_color,
                    densities=float(self.config.get('taichi_volume_sample_density', block_density * 0.75)),
                    colors=np.asarray(sample_colors, dtype=np.float32),
                )
                sample_count = int(np.count_nonzero(inside))

        return smoke_density, smoke_color, {
            'volume_shape': volume_shape,
            'upsample': int(scale),
            'estimated_memory_mb': estimated_bytes / (1024.0 * 1024.0),
            'block_count': int(len(block_indices)),
            'sample_count': int(sample_count),
        }

    def _run_volume_renderer(self):
        try:
            volume_renderer, volume_canvas = import_volume_renderer()
        except Exception as exc:
            self._set_status(
                f"Transparent block mode needs taichi-volume-renderer ({exc}). Falling back to mesh blocks."
            )
            self.block_render_mode = 'mesh'
            self._run_native_scene()
            return

        try:
            window, metadata = self._create_volume_display_window(volume_renderer, volume_canvas)
        except Exception as exc:
            self._set_status(f"Transparent block mode unavailable: {exc}. Falling back to mesh blocks.")
            self.block_render_mode = 'mesh'
            return 'restart'

        print(
            'Transparent block volume renderer: '
            f"shape={metadata['volume_shape']} upsample={metadata['upsample']} "
            f"~{metadata['estimated_memory_mb']:.1f}MB blocks={metadata['block_count']} "
            f"samples={metadata['sample_count']} transparency={self._volume_transparency_percent():.0f}%"
        )

        class _RestartVolumeWindow(Exception):
            pass

        def _volume_callback(_iteration, _scene):
            reload_status = self._maybe_reload_external_state()
            if reload_status == 'restart':
                raise _RestartVolumeWindow()
            if reload_status == 'reloaded':
                refreshed_window, refreshed_meta = self._create_volume_display_window(
                    volume_renderer,
                    volume_canvas,
                    resolution=window.resolution,
                )
                window.scene = refreshed_window.scene
                print(
                    'Transparent block volume renderer: '
                    f"shape={refreshed_meta['volume_shape']} upsample={refreshed_meta['upsample']} "
                    f"~{refreshed_meta['estimated_memory_mb']:.1f}MB blocks={refreshed_meta['block_count']} "
                    f"samples={refreshed_meta['sample_count']} transparency={self._volume_transparency_percent():.0f}%"
                )

        try:
            window.show(title=f'{self.window_title} (Transparent Blocks)', callback=_volume_callback)
        except _RestartVolumeWindow:
            return 'restart'
        return None

    def _run_native_scene(self):
        window = ti.ui.Window(self.window_title, (1280, 720), vsync=True)
        canvas = window.get_canvas()
        canvas.set_background_color((1.0, 1.0, 1.0))
        camera = ti.ui.Camera()
        movement_speed = 0.01
        min_movement_speed = 0.00075
        max_movement_speed = 0.25
        speed_step_factor = 1.6
        wheel_zoom_factor = 6.0
        yaw_speed = 2.0
        pitch_speed = 2.0
        click_drag_threshold = 0.0035
        initial_camera_position = np.array([0.0, -1.35, 0.85], dtype=np.float32)
        initial_camera_lookat = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        initial_camera_up = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        orbit_pivot = initial_camera_lookat.copy()
        lmb_last_mouse = None
        rmb_last_mouse = None
        mmb_press_pos = None
        mmb_last_mouse = None
        mmb_dragging = False
        lmb_press_pos = None
        lmb_dragging = False
        previous_lmb_pressed = False
        previous_mmb_pressed = False
        previous_rmb_pressed = False
        previous_speed_up_pressed = False
        previous_speed_down_pressed = False
        last_time_ns = None

        camera.position(*initial_camera_position)
        camera.lookat(*initial_camera_lookat)
        camera.up(*initial_camera_up)
        camera.fov(float(self.config.get('taichi_mesh_camera_fov', 55.0)))

        sample_radius = self._sample_render_radius()
        self._start_external_wheel_listener()
        self._set_status("Controls: RMB orbit pivot, WASD/E/Q move along current view, Wheel/Z/X zoom, +/- speed, MMB click pivot, MMB drag pan, LMB drag look/click info, I update, Shift+I all, H reset")
        self._print_legend_if_needed(force=True)

        try:
            while window.running:
                reload_status = self._maybe_reload_external_state()
                if reload_status == 'restart':
                    return 'restart'
                self._poll_worker_result()
                scene = window.get_scene()
                sample_radius = self._sample_render_radius()
                current_time_ns = time.perf_counter_ns()
                if last_time_ns is None:
                    last_time_ns = current_time_ns
                time_elapsed = (current_time_ns - last_time_ns) * 1e-9
                last_time_ns = current_time_ns

                for event in window.get_events(ti.ui.PRESS):
                    event_key = getattr(event, 'key', None)
                    event_key_normalized = str(event_key).strip().upper() if event_key is not None else ''
                    event_key_lower = str(event_key).strip().lower() if event_key is not None else ''
                    shift_pressed = self._event_has_shift(event, window)

                    if event_key == 'b':
                        self.show_blocks = not self.show_blocks
                    elif event_key == 'v':
                        self.show_samples = not self.show_samples
                    elif event_key == 'f':
                        self.show_filled_block_faces = not self.show_filled_block_faces
                        self._set_status(f"Filled block faces: {'on' if self.show_filled_block_faces else 'off'}")
                    elif event_key_lower == 'i':
                        if shift_pressed:
                            self.request_update_all()
                        else:
                            self.request_update()
                    elif event_key == 'c':
                        self.toggle_color_mode()
                    elif event_key == 'l':
                        self.show_legend_overlay = not self.show_legend_overlay
                        self._overlay_mouse_capture = False
                        self._set_status(f"Legend overlay: {'on' if self.show_legend_overlay else 'off'}")
                    elif event_key == 'h':
                        camera.position(*initial_camera_position)
                        camera.lookat(*initial_camera_lookat)
                        camera.up(*initial_camera_up)
                        orbit_pivot = initial_camera_lookat.copy()
                        lmb_last_mouse = None
                        rmb_last_mouse = None
                        mmb_press_pos = None
                        mmb_last_mouse = None
                        mmb_dragging = False
                        lmb_press_pos = None
                        lmb_dragging = False
                        previous_lmb_pressed = False
                        previous_mmb_pressed = False
                        previous_rmb_pressed = False
                        self._clear_external_wheel_delta()
                        self._set_status('View reset to initial camera position.')
                    elif event_key == '[':
                        self.value_filter -= 1.0
                        self._refresh_render_data()
                        self._set_status(f"value_filter={self.value_filter}")
                    elif event_key == ']':
                        self.value_filter += 1.0
                        self._refresh_render_data()
                        self._set_status(f"value_filter={self.value_filter}")
                    elif event_key_lower == 'z':
                        self._apply_zoom_delta(camera, orbit_pivot, 1.0, movement_speed, wheel_zoom_factor)
                        self._set_status('Zoom in')
                    elif event_key_lower == 'x':
                        self._apply_zoom_delta(camera, orbit_pivot, -1.0, movement_speed, wheel_zoom_factor)
                        self._set_status('Zoom out')
                    elif event_key in ('=', '+') or event_key_normalized in {'KP_ADD', 'ADD', 'PLUS', 'EQUAL', 'EQUALS', 'NUMPADADD'}:
                        movement_speed = min(max_movement_speed, movement_speed * speed_step_factor)
                        self._set_status(f"movement_speed={movement_speed:.4f}")
                    elif event_key in ('-', '_') or event_key_normalized in {'KP_SUBTRACT', 'SUBTRACT', 'MINUS', 'UNDERSCORE', 'NUMPADSUBTRACT'}:
                        movement_speed = max(min_movement_speed, movement_speed / speed_step_factor)
                        self._set_status(f"movement_speed={movement_speed:.4f}")
                    elif event_key == ti.ui.ESCAPE:
                        window.running = False

                window_shape = window.get_window_shape()
                curr_mouse_x, curr_mouse_y = window.get_cursor_pos()
                self._sync_overlay_bounds_from_ini()
                overlay_hovered = self._cursor_over_overlay(curr_mouse_x, curr_mouse_y, window_shape)

                for event in window.get_events():
                    event_key = getattr(event, 'key', None)
                    event_key_normalized = str(event_key).strip().upper() if event_key is not None else ''
                    delta = getattr(event, 'delta', None)
                    if 'WHEEL' not in event_key_normalized and 'SCROLL' not in event_key_normalized:
                        continue
                    if isinstance(delta, (tuple, list, np.ndarray)):
                        wheel_delta = float(delta[1] if len(delta) > 1 else delta[0])
                    else:
                        wheel_delta = float(delta if delta is not None else 0.0)
                    if abs(wheel_delta) <= 1e-9:
                        continue
                    if overlay_hovered or self._overlay_mouse_capture:
                        continue
                    self._apply_zoom_delta(camera, orbit_pivot, wheel_delta, movement_speed, wheel_zoom_factor)

                current_cursor = np.array([curr_mouse_x, curr_mouse_y], dtype=np.float32)
                current_lmb_pressed = window.is_pressed(ti.ui.LMB)
                current_mmb_pressed = window.is_pressed(ti.ui.MMB)
                current_rmb_pressed = window.is_pressed(ti.ui.RMB)
                mouse_press_started = (
                    (current_lmb_pressed and not previous_lmb_pressed)
                    or (current_mmb_pressed and not previous_mmb_pressed)
                    or (current_rmb_pressed and not previous_rmb_pressed)
                )
                if self._overlay_mouse_capture:
                    if not (current_lmb_pressed or current_mmb_pressed or current_rmb_pressed):
                        self._overlay_mouse_capture = False
                elif overlay_hovered and mouse_press_started:
                    self._overlay_mouse_capture = True

                if self._overlay_mouse_capture:
                    lmb_press_pos = None
                    lmb_last_mouse = None
                    lmb_dragging = False
                    mmb_press_pos = None
                    mmb_last_mouse = None
                    mmb_dragging = False
                    rmb_last_mouse = None

                external_wheel_delta = self._consume_external_wheel_delta()
                if abs(external_wheel_delta) > 1e-9:
                    if 0.0 <= curr_mouse_x <= 1.0 and 0.0 <= curr_mouse_y <= 1.0 and not overlay_hovered and not self._overlay_mouse_capture:
                        self._apply_zoom_delta(camera, orbit_pivot, external_wheel_delta, movement_speed, wheel_zoom_factor)

                current_speed_up_pressed = any(
                    window.is_pressed(key)
                    for key in ['=', '+', 'KP_ADD', 'Add', 'PLUS']
                )
                current_speed_down_pressed = any(
                    window.is_pressed(key)
                    for key in ['-', '_', 'KP_SUBTRACT', 'Subtract', 'MINUS']
                )
                if current_speed_up_pressed and not previous_speed_up_pressed:
                    movement_speed = min(max_movement_speed, movement_speed * speed_step_factor)
                    self._set_status(f"movement_speed={movement_speed:.4f}")
                if current_speed_down_pressed and not previous_speed_down_pressed:
                    movement_speed = max(min_movement_speed, movement_speed / speed_step_factor)
                    self._set_status(f"movement_speed={movement_speed:.4f}")

                if self._overlay_mouse_capture:
                    pass
                elif current_mmb_pressed:
                    lmb_press_pos = None
                    lmb_last_mouse = None
                    lmb_dragging = False

                if self._overlay_mouse_capture:
                    pass
                elif not current_mmb_pressed and current_lmb_pressed and not previous_lmb_pressed:
                    lmb_press_pos = current_cursor.copy()
                    lmb_last_mouse = current_cursor.copy()
                    lmb_dragging = False
                elif not current_mmb_pressed and current_lmb_pressed and lmb_last_mouse is not None:
                    dx = float(current_cursor[0] - lmb_last_mouse[0])
                    dy = float(current_cursor[1] - lmb_last_mouse[1])
                    if abs(dx) > 1e-6 or abs(dy) > 1e-6:
                        if lmb_press_pos is not None:
                            total_drag = float(np.linalg.norm(current_cursor - lmb_press_pos))
                            if total_drag > click_drag_threshold:
                                lmb_dragging = True
                        position, lookat, front, right, up = self._camera_basis(camera)
                        lookat_offset = lookat - position
                        rotation_scale = float(time_elapsed * 60.0)
                        yaw_angle = float(-dx * yaw_speed * rotation_scale)
                        pitch_angle = float(dy * pitch_speed * rotation_scale)
                        rotated_lookat_offset = self._rotate_vector(lookat_offset, up, yaw_angle)
                        rotated_up = up.copy()
                        rotated_right = self._rotate_vector(right, up, yaw_angle)
                        rotated_right_norm = float(np.linalg.norm(rotated_right))
                        if rotated_right_norm > 1e-6:
                            rotated_right = (rotated_right / rotated_right_norm).astype(np.float32)
                        else:
                            rotated_right = right
                        candidate_lookat_offset = self._rotate_vector(rotated_lookat_offset, rotated_right, pitch_angle)
                        candidate_up = self._rotate_vector(rotated_up, rotated_right, pitch_angle)
                        candidate_front_norm = float(np.linalg.norm(candidate_lookat_offset))
                        if candidate_front_norm > 1e-6:
                            candidate_front = (candidate_lookat_offset / candidate_front_norm).astype(np.float32)
                        else:
                            candidate_front = front
                        vertical_component = float(np.clip(candidate_front[1], -1.0, 1.0))
                        if abs(vertical_component) < 0.995:
                            new_lookat_offset = candidate_lookat_offset
                            new_up = candidate_up
                        else:
                            new_lookat_offset = rotated_lookat_offset
                            new_up = rotated_up
                        camera.lookat(*(position + new_lookat_offset))
                        camera.up(*new_up)
                    lmb_last_mouse = current_cursor.copy()
                elif previous_lmb_pressed and not current_lmb_pressed:
                    if not lmb_dragging and lmb_press_pos is not None:
                        picked = self._pick_visible_object(tuple(lmb_press_pos.tolist()), camera, window.get_window_shape())
                        if picked is not None:
                            obj_type = 'sample' if picked['is_sample'] else 'block'
                            pos = picked['world_position']
                            message = (
                                f"{obj_type} x={pos[0]:.2f} y={pos[1]:.2f} z={pos[2]:.2f} "
                                f"value={picked['value']:.3f} domain={picked['domain']}"
                            )
                            self._set_status(message)
                            print(
                                f"Clicked {obj_type}: world=({pos[0]:.6g}, {pos[1]:.6g}, {pos[2]:.6g}) "
                                f"value={picked['value']:.6g} domain={picked['domain']}"
                            )
                        else:
                            self._set_status('No visible object under cursor.')
                    lmb_press_pos = None
                    lmb_last_mouse = None
                    lmb_dragging = False

                if self._overlay_mouse_capture:
                    pass
                elif current_mmb_pressed and not previous_mmb_pressed:
                    mmb_press_pos = current_cursor.copy()
                    mmb_last_mouse = current_cursor.copy()
                    mmb_dragging = False
                elif current_mmb_pressed and mmb_last_mouse is not None:
                    dx = float(current_cursor[0] - mmb_last_mouse[0])
                    dy = float(current_cursor[1] - mmb_last_mouse[1])
                    if abs(dx) > 1e-6 or abs(dy) > 1e-6:
                        if mmb_press_pos is not None:
                            total_drag = float(np.linalg.norm(current_cursor - mmb_press_pos))
                            if total_drag > click_drag_threshold:
                                mmb_dragging = True
                        if mmb_dragging:
                            position, lookat, front, right, up = self._camera_basis(camera)
                            aspect = float(window.get_window_shape()[0]) / float(window.get_window_shape()[1])
                            look_distance = max(float(np.linalg.norm(lookat - position)), 0.05)
                            fov_rad = np.deg2rad(float(self.config.get('taichi_mesh_camera_fov', 55.0)))
                            pan_scale_y = look_distance * np.tan(fov_rad * 0.5) * 2.0
                            pan_scale_x = pan_scale_y * aspect
                            offset = (-right * dx * pan_scale_x) - (up * dy * pan_scale_y)
                            new_position = position + offset
                            camera.position(*new_position)
                            camera.lookat(*(new_position + front * look_distance))
                            camera.up(*up)
                    mmb_last_mouse = current_cursor.copy()
                elif previous_mmb_pressed and not current_mmb_pressed:
                    if not mmb_dragging and mmb_press_pos is not None:
                        picked = self._pick_visible_object(tuple(mmb_press_pos.tolist()), camera, window.get_window_shape())
                        if picked is not None:
                            orbit_pivot = np.asarray(picked.get('pivot_render_position', picked['render_position']), dtype=np.float32)
                            camera.lookat(*orbit_pivot)
                            self._set_status('Pivot set from picked object and centered without changing zoom.')
                        else:
                            self._set_status('No visible object under cursor for pivot selection.')
                    mmb_press_pos = None
                    mmb_last_mouse = None
                    mmb_dragging = False

                if self._overlay_mouse_capture:
                    pass
                elif current_rmb_pressed and not previous_rmb_pressed:
                    rmb_last_mouse = current_cursor.copy()
                elif current_rmb_pressed and rmb_last_mouse is not None:
                    dx = float(current_cursor[0] - rmb_last_mouse[0])
                    dy = float(current_cursor[1] - rmb_last_mouse[1])
                    if abs(dx) > 1e-6 or abs(dy) > 1e-6:
                        position, _lookat, front, right, up = self._camera_basis(camera)
                        lookat = np.asarray(camera.curr_lookat, dtype=np.float32)
                        position_offset = position - orbit_pivot
                        lookat_offset = lookat - orbit_pivot
                        rotation_scale = float(time_elapsed * 60.0)
                        yaw_angle = float(-dx * yaw_speed * rotation_scale)
                        pitch_angle = float(dy * pitch_speed * rotation_scale)
                        rotated_position_offset = self._rotate_vector(position_offset, up, yaw_angle)
                        rotated_lookat_offset = self._rotate_vector(lookat_offset, up, yaw_angle)
                        rotated_up = up.copy()
                        rotated_right = self._rotate_vector(right, up, yaw_angle)
                        rotated_right_norm = float(np.linalg.norm(rotated_right))
                        if rotated_right_norm > 1e-6:
                            rotated_right = (rotated_right / rotated_right_norm).astype(np.float32)
                        else:
                            rotated_right = right
                        candidate_position_offset = self._rotate_vector(rotated_position_offset, rotated_right, pitch_angle)
                        candidate_lookat_offset = self._rotate_vector(rotated_lookat_offset, rotated_right, pitch_angle)
                        candidate_up = self._rotate_vector(rotated_up, rotated_right, pitch_angle)
                        candidate_front = candidate_lookat_offset - candidate_position_offset
                        candidate_front_norm = float(np.linalg.norm(candidate_front))
                        if candidate_front_norm > 1e-6:
                            candidate_front = (candidate_front / candidate_front_norm).astype(np.float32)
                        else:
                            candidate_front = front
                        vertical_component = float(np.clip(candidate_front[1], -1.0, 1.0))
                        if abs(vertical_component) < 0.995:
                            new_position_offset = candidate_position_offset
                            new_lookat_offset = candidate_lookat_offset
                            new_up = candidate_up
                        else:
                            new_position_offset = rotated_position_offset
                            new_lookat_offset = rotated_lookat_offset
                            new_up = rotated_up
                        new_position = orbit_pivot + new_position_offset
                        new_lookat = orbit_pivot + new_lookat_offset
                        camera.position(*new_position)
                        camera.lookat(*new_lookat)
                        camera.up(*new_up)
                    rmb_last_mouse = current_cursor.copy()
                else:
                    rmb_last_mouse = None

                movement_keys_pressed = any(window.is_pressed(key) for key in ['w', 'a', 's', 'd', 'e', 'q'])
                if movement_keys_pressed and not current_rmb_pressed and not current_mmb_pressed and not self._overlay_mouse_capture:
                    position, lookat, front, right, up = self._camera_basis(camera)
                    scaled_speed = float(movement_speed * time_elapsed * 60.0)
                    movement_delta = np.zeros(3, dtype=np.float32)
                    if window.is_pressed('w'):
                        movement_delta += front * scaled_speed
                    if window.is_pressed('s'):
                        movement_delta -= front * scaled_speed
                    if window.is_pressed('a'):
                        movement_delta -= right * scaled_speed
                    if window.is_pressed('d'):
                        movement_delta += right * scaled_speed
                    if window.is_pressed('e'):
                        movement_delta += up * scaled_speed
                    if window.is_pressed('q'):
                        movement_delta -= up * scaled_speed
                    if float(np.linalg.norm(movement_delta)) > 0.0:
                        camera.position(*(position + movement_delta))
                        camera.lookat(*(lookat + movement_delta))

                previous_lmb_pressed = current_lmb_pressed
                previous_mmb_pressed = current_mmb_pressed
                previous_rmb_pressed = current_rmb_pressed
                previous_speed_up_pressed = current_speed_up_pressed
                previous_speed_down_pressed = current_speed_down_pressed

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
        finally:
            self._stop_external_wheel_listener()
        return None

    def _extract_interpolators(self):
        interpolators = []
        for group in self._extract_interpolator_groups():
            interpolators.extend(group['passes'])
        return interpolators

    def _extract_interpolator_groups(self):
        blocks = self.blocks_data
        groups = []
        if hasattr(blocks, '_interpolators') and blocks._interpolators:
            for domain, entry in blocks._interpolators.items():
                if isinstance(entry, list):
                    passes = [interp for interp in entry if interp is not None]
                elif entry is not None:
                    passes = [entry]
                else:
                    passes = []
                if passes:
                    groups.append({'domain': str(domain), 'passes': passes})
        elif hasattr(blocks, '_ant_colony') and getattr(blocks, '_ant_colony', None) is not None:
            groups.append({'domain': 'Global', 'passes': [blocks._ant_colony]})
        return groups

    def _reset_interpolation_progress(self):
        self._interpolation_groups = self._extract_interpolator_groups()
        self._interpolation_domain_index = 0
        self._interpolation_pass_index = 0
        self._interpolation_pass_iteration = 0
        self._interpolation_total_phases = sum(len(group['passes']) for group in self._interpolation_groups)
        self._interpolation_completed_phases = 0
        self._interpolation_complete = self._interpolation_total_phases == 0
        self._phase_provenance_snapshot = None
        self._phase_provenance_interp = None

    def _current_interpolation_group(self):
        while self._interpolation_domain_index < len(self._interpolation_groups):
            group = self._interpolation_groups[self._interpolation_domain_index]
            if self._interpolation_pass_index < len(group['passes']):
                return group
            self._interpolation_domain_index += 1
            self._interpolation_pass_index = 0
            self._interpolation_pass_iteration = 0
        self._interpolation_complete = True
        return None

    def _current_interpolator(self):
        group = self._current_interpolation_group()
        if group is None:
            return None, None
        return group, group['passes'][self._interpolation_pass_index]

    def _format_phase_label(self, domain, pass_index, total_passes, interpolator):
        algo_name = getattr(interpolator, 'get_algorithm_name', lambda: interpolator.__class__.__name__)()
        if total_passes > 1:
            return f"{domain} - Pass {pass_index} ({algo_name})"
        if domain == 'Global':
            return str(algo_name)
        return f"{domain} ({algo_name})"

    def _ensure_interpolator_next_block_id(self, interp):
        if hasattr(interp, 'next_block_id'):
            return
        existing_ids = []
        blocks = getattr(interp, 'blocks', {}) or {}
        if isinstance(blocks, dict):
            for block in blocks.values():
                if hasattr(block, 'block_id'):
                    existing_ids.append(int(block.block_id))
                elif isinstance(block, dict) and 'Block_ID' in block:
                    existing_ids.append(int(block['Block_ID']))
        interp.next_block_id = (max(existing_ids) + 1) if existing_ids else 1

    def _run_single_interpolator_iteration(self, interp, dims):
        self._ensure_interpolator_next_block_id(interp)
        try:
            previous_blocks = len(interp.get_interpolated_values())
        except Exception:
            previous_blocks = len(getattr(interp, 'blocks', {}) or {})

        should_continue = bool(interp.run_iteration(dims))
        changed = should_continue

        if not should_continue:
            try:
                current_blocks = len(interp.get_interpolated_values())
            except Exception:
                current_blocks = len(getattr(interp, 'blocks', {}) or {})
            if current_blocks > previous_blocks:
                changed = True

        try:
            converged = bool(interp.is_converged())
        except Exception:
            converged = False

        return changed, should_continue, converged

    def _get_domain_post_process_mode(self, domain):
        if not domain or str(domain).strip().lower() == 'global':
            return 'skip'
        domain_overrides = self.config.get('domain_algorithm_overrides') or {}
        domain_cfg = domain_overrides.get(domain) or {}
        if domain_cfg.get('skip', False):
            return 'skip'
        mode = str(domain_cfg.get('post_process', 'skip')).strip().lower()
        if mode in ('fill_with_average', 'fill average', 'fill_average', 'fill'):
            return 'fill_with_average'
        return 'skip'

    def _run_domain_post_process(self, domain, interp, dims):
        if self._get_domain_post_process_mode(domain) != 'fill_with_average':
            return 0, 0
        if not hasattr(interp, 'fill_unvisited_blocks_domainwise'):
            return 0, 0
        try:
            created, assigned = interp.fill_unvisited_blocks_domainwise(dims)
        except Exception:
            return 0, 0
        return int(created or 0), int(assigned or 0)

    def _prepare_next_pass(self, completed_interp, next_interp):
        dims = tuple(self.block_info.get('dims', (0, 0, 0)))
        min_bounds = self.block_info.get('min_bounds', np.zeros(3))
        block_size = self.block_info.get('block_size', self.block_size)
        pass_values = completed_interp.get_interpolated_values()
        use_mapping = bool(getattr(next_interp, 'allowed_grid_override', None) is not None)
        next_interp.initialize_blocks(pass_values, dims, min_bounds, block_size, use_domain_mapping=use_mapping)
        copy_interpolator_provenance(completed_interp, next_interp)
        if hasattr(next_interp, 'create_ants'):
            next_interp.create_ants()

    def _advance_interpolation_phase(self, completed_interp):
        group = self._current_interpolation_group()
        if group is None:
            self._interpolation_complete = True
            return None

        self._interpolation_completed_phases = min(
            self._interpolation_total_phases,
            self._interpolation_completed_phases + 1,
        )
        next_pass_index = self._interpolation_pass_index + 1
        self._interpolation_pass_iteration = 0

        if next_pass_index < len(group['passes']):
            next_interp = group['passes'][next_pass_index]
            self._prepare_next_pass(completed_interp, next_interp)
            self._interpolation_pass_index = next_pass_index
            return self._format_phase_label(
                group['domain'],
                self._interpolation_pass_index + 1,
                len(group['passes']),
                next_interp,
            )

        self._interpolation_domain_index += 1
        self._interpolation_pass_index = 0
        next_group = self._current_interpolation_group()
        if next_group is None:
            self._interpolation_complete = True
            return None

        next_interp = next_group['passes'][0]
        return self._format_phase_label(
            next_group['domain'],
            1,
            len(next_group['passes']),
            next_interp,
        )

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

    def _sample_render_radius(self):
        diameter_world = float(self.config.get('taichi_sample_diameter', 1.0))
        diameter_world = max(diameter_world, 1e-6)
        return max((diameter_world * self._render_scale) / 2.0, 1e-6)

    def _camera_basis(self, camera):
        position = np.asarray(camera.curr_position, dtype=np.float32)
        lookat = np.asarray(camera.curr_lookat, dtype=np.float32)
        front = lookat - position
        front_norm = float(np.linalg.norm(front))
        if front_norm <= 1e-6:
            front = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        else:
            front = (front / front_norm).astype(np.float32)

        up = np.asarray(camera.curr_up, dtype=np.float32)
        up_norm = float(np.linalg.norm(up))
        if up_norm <= 1e-6:
            up = np.array([0.0, 1.0, 0.0], dtype=np.float32)
        else:
            up = (up / up_norm).astype(np.float32)

        right = np.cross(front, up).astype(np.float32)
        right_norm = float(np.linalg.norm(right))
        if right_norm <= 1e-6:
            right = np.array([1.0, 0.0, 0.0], dtype=np.float32)
        else:
            right = (right / right_norm).astype(np.float32)
        up = np.cross(right, front).astype(np.float32)
        up_norm = float(np.linalg.norm(up))
        if up_norm > 1e-6:
            up = (up / up_norm).astype(np.float32)

        return position, lookat, front, right, up

    def _queue_external_wheel_delta(self, delta):
        if abs(float(delta)) <= 1e-9:
            return
        with self._external_wheel_lock:
            self._external_wheel_delta += float(delta)

    def _worker_running(self):
        return self._worker is not None and self._worker.is_alive()

    def _queue_worker_result(self, result):
        with self._result_condition:
            while self._pending_result is not None:
                self._result_condition.wait(timeout=0.05)
            self._pending_result = result
            self._result_condition.notify_all()

    def _event_has_shift(self, event, window):
        modifier = getattr(event, 'modifier', None)
        if modifier is not None:
            if isinstance(modifier, (list, tuple, set)):
                if any('SHIFT' in str(item).strip().upper() for item in modifier):
                    return True
            elif 'SHIFT' in str(modifier).strip().upper():
                return True
        try:
            if window.is_pressed(ti.ui.SHIFT):
                return True
        except Exception:
            pass
        event_key = getattr(event, 'key', None)
        return isinstance(event_key, str) and len(event_key) == 1 and event_key.isalpha() and event_key.isupper()

    def _consume_external_wheel_delta(self):
        with self._external_wheel_lock:
            delta = self._external_wheel_delta
            self._external_wheel_delta = 0.0
        return float(delta)

    def _clear_external_wheel_delta(self):
        with self._external_wheel_lock:
            self._external_wheel_delta = 0.0

    def _start_external_wheel_listener(self):
        if self._external_wheel_listener is not None or pynput_mouse is None:
            return False
        self._clear_external_wheel_delta()
        try:
            def _on_scroll(_x, _y, _dx, dy):
                self._queue_external_wheel_delta(dy)
                return True

            listener = pynput_mouse.Listener(on_scroll=_on_scroll)
            if hasattr(listener, 'daemon'):
                listener.daemon = True
            listener.start()
            self._external_wheel_listener = listener
            return True
        except Exception as exc:
            self._external_wheel_listener = None
            print(f"External wheel listener unavailable: {exc}")
            return False

    def _stop_external_wheel_listener(self):
        listener = self._external_wheel_listener
        self._external_wheel_listener = None
        self._clear_external_wheel_delta()
        if listener is None:
            return
        try:
            listener.stop()
        except Exception:
            pass
        try:
            listener.join(0.5)
        except Exception:
            pass

    def _apply_zoom_delta(self, camera, orbit_pivot, zoom_steps, movement_speed, wheel_zoom_factor):
        zoom_steps = float(zoom_steps)
        if abs(zoom_steps) <= 1e-9:
            return
        position, lookat, front, _right, _up = self._camera_basis(camera)
        dolly_distance = zoom_steps * float(movement_speed) * float(wheel_zoom_factor)
        if float(np.linalg.norm(lookat - orbit_pivot)) <= 1e-4:
            pivot_distance = float(np.linalg.norm(orbit_pivot - position))
            if dolly_distance > 0.0:
                dolly_distance = min(dolly_distance, max(pivot_distance - 0.05, 0.0))
            camera.position(*(position + front * dolly_distance))
            camera.lookat(*orbit_pivot)
        else:
            offset = front * dolly_distance
            camera.position(*(position + offset))
            camera.lookat(*(lookat + offset))

    def _rotate_vector(self, vector, axis, angle):
        vector = np.asarray(vector, dtype=np.float32)
        axis = np.asarray(axis, dtype=np.float32)
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm <= 1e-6 or abs(angle) <= 1e-9:
            return vector.copy()
        axis = (axis / axis_norm).astype(np.float32)
        cos_angle = float(np.cos(angle))
        sin_angle = float(np.sin(angle))
        rotated = (
            vector * cos_angle
            + np.cross(axis, vector).astype(np.float32) * sin_angle
            + axis * np.dot(axis, vector) * (1.0 - cos_angle)
        )
        return np.asarray(rotated, dtype=np.float32)

    def _front_to_yaw_pitch(self, front):
        front = np.asarray(front, dtype=np.float32)
        norm = float(np.linalg.norm(front))
        if norm <= 1e-6:
            front = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        else:
            front = front / norm
        pitch = float(np.arcsin(np.clip(front[1], -1.0, 1.0)))
        yaw = float(np.arctan2(front[0], front[2]))
        return yaw, pitch

    def _yaw_pitch_to_front(self, yaw, pitch):
        cos_pitch = float(np.cos(pitch))
        return np.array([
            np.sin(yaw) * cos_pitch,
            np.sin(pitch),
            np.cos(yaw) * cos_pitch,
        ], dtype=np.float32)

    def _pick_visible_object(self, cursor_pos, camera, resolution):
        snapshot = self._raw_snapshot
        if snapshot is None or len(snapshot['positions']) == 0:
            return None

        mask = np.zeros(len(snapshot['positions']), dtype=bool)
        if self.show_blocks:
            if self.color_mode == 'domain' and self.domain_data_available:
                mask |= np.ones(len(snapshot['positions']), dtype=bool)
            else:
                mask |= snapshot['is_sample'] | (snapshot['values'] >= self.value_filter)
        if self.show_samples and not self.show_blocks:
            mask |= snapshot['is_sample']
        if not np.any(mask):
            return None

        points_world = np.asarray(snapshot['positions'][mask], dtype=np.float32)
        points_render = self._to_render_space(points_world)
        values = np.asarray(snapshot['values'][mask], dtype=np.float32)
        domains = np.asarray(snapshot['domains'][mask], dtype=object)
        is_sample = np.asarray(snapshot['is_sample'][mask], dtype=bool)

        width, height = resolution
        if width <= 0 or height <= 0:
            return None

        position, _lookat, front, right, up = self._camera_basis(camera)
        relative = points_render - position[None, :]
        depth = relative @ front
        valid = depth > 1e-5
        if not np.any(valid):
            return None

        fov_deg = float(self.config.get('taichi_mesh_camera_fov', 55.0))
        tan_half_fov = np.tan(np.deg2rad(fov_deg) * 0.5)
        aspect = float(width) / float(height)
        x_camera = relative @ right
        y_camera = relative @ up
        x_ndc = x_camera / (depth * tan_half_fov * aspect)
        y_ndc = y_camera / (depth * tan_half_fov)
        screen_x = 0.5 + 0.5 * x_ndc
        screen_y = 0.5 + 0.5 * y_ndc
        screen_dist = np.sqrt((screen_x - float(cursor_pos[0])) ** 2 + (screen_y - float(cursor_pos[1])) ** 2)
        threshold = 28.0 / float(min(width, height))
        candidate_mask = valid & (np.abs(x_ndc) <= 1.1) & (np.abs(y_ndc) <= 1.1) & (screen_dist <= threshold)
        if not np.any(candidate_mask):
            return None

        candidate_indices = np.where(candidate_mask)[0]
        order = np.lexsort((depth[candidate_indices], screen_dist[candidate_indices]))
        picked_index = candidate_indices[order[0]]
        cursor_x_ndc = (float(cursor_pos[0]) - 0.5) * 2.0
        cursor_y_ndc = (float(cursor_pos[1]) - 0.5) * 2.0
        ray_direction = (
            front
            + right * (cursor_x_ndc * tan_half_fov * aspect)
            + up * (cursor_y_ndc * tan_half_fov)
        ).astype(np.float32)
        ray_norm = float(np.linalg.norm(ray_direction))
        if ray_norm > 1e-6:
            ray_direction = (ray_direction / ray_norm).astype(np.float32)
        depth_along_ray = float(np.dot(ray_direction, front))
        if depth_along_ray > 1e-6:
            pivot_render_position = position + ray_direction * (float(depth[picked_index]) / depth_along_ray)
        else:
            pivot_render_position = np.asarray(points_render[picked_index], dtype=np.float32)
        return {
            'render_position': np.asarray(points_render[picked_index], dtype=np.float32),
            'pivot_render_position': np.asarray(pivot_render_position, dtype=np.float32),
            'world_position': np.asarray(points_world[picked_index], dtype=np.float32),
            'value': float(values[picked_index]),
            'domain': str(domains[picked_index]),
            'is_sample': bool(is_sample[picked_index]),
        }

    def _pick_mesh_pivot(self, cursor_pos, camera, resolution):
        picked = self._pick_visible_object(cursor_pos, camera, resolution)
        if picked is None:
            return None
        return np.asarray(picked['render_position'], dtype=np.float32)

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
        self._render_sample_points = render_sample_points
        self._render_block_points = render_block_points
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
        target_iterations = max(int(self.config.get('iterations', 500)), 1)
        group, interp = self._current_interpolator()
        if group is None or interp is None:
            self._interpolation_complete = True
            return {
                'changed': False,
                'all_finished': True,
                'phase_finished': True,
                'phase_label': 'Interpolation complete',
                'phase_iteration': 0,
                'phase_iterations': target_iterations,
                'next_phase_label': None,
                'completed_phases': self._interpolation_completed_phases,
                'total_phases': self._interpolation_total_phases,
            }

        phase_label = self._format_phase_label(
            group['domain'],
            self._interpolation_pass_index + 1,
            len(group['passes']),
            interp,
        )

        if self._interpolation_pass_iteration == 0 or self._phase_provenance_interp is not interp:
            if self._interpolation_pass_index == 0:
                seed_original_sample_provenance(interp)
            self._phase_provenance_snapshot = snapshot_interpolator_state(interp)
            self._phase_provenance_interp = interp

        changed, should_continue, converged = self._run_single_interpolator_iteration(interp, dims)
        self._interpolation_pass_iteration += 1
        phase_iteration = self._interpolation_pass_iteration
        phase_finished = (not should_continue) or converged or (phase_iteration >= target_iterations)
        next_phase_label = None
        if phase_finished:
            phase_source = 'First Pass' if self._interpolation_pass_index == 0 else 'Second Pass'
            finalize_phase_provenance(interp, phase_source, interp.get_algorithm_name(), self._phase_provenance_snapshot)
            is_last_pass_in_domain = (self._interpolation_pass_index + 1) >= len(group['passes'])
            if is_last_pass_in_domain:
                post_process_snapshot = snapshot_interpolator_state(interp)
                created, assigned = self._run_domain_post_process(group['domain'], interp, dims)
                if created or assigned:
                    finalize_phase_provenance(interp, 'Post-process', 'Fill with Average', post_process_snapshot)
                    changed = True
            next_phase_label = self._advance_interpolation_phase(interp)
            self._phase_provenance_snapshot = None
            self._phase_provenance_interp = None

        return {
            'changed': changed,
            'all_finished': self._interpolation_complete,
            'phase_finished': phase_finished,
            'phase_label': phase_label,
            'phase_iteration': phase_iteration,
            'phase_iterations': target_iterations,
            'next_phase_label': next_phase_label,
            'completed_phases': self._interpolation_completed_phases,
            'total_phases': self._interpolation_total_phases,
        }

    def request_update(self):
        if self._worker_running():
            self._set_status('Update already running...')
            return

        def _worker():
            started = time.perf_counter()
            try:
                step = None
                snapshot = None
                worker_steps = 0
                while True:
                    step = self._run_interpolators_once()
                    snapshot = self._collect_block_snapshot()
                    worker_steps += 1
                    if step['all_finished'] or step['phase_finished']:
                        break

                if step is None:
                    step = {
                        'changed': False,
                        'all_finished': True,
                        'phase_finished': True,
                        'phase_label': 'Interpolation complete',
                        'phase_iteration': 0,
                        'phase_iterations': max(int(self.config.get('iterations', 500)), 1),
                        'next_phase_label': None,
                        'completed_phases': self._interpolation_completed_phases,
                        'total_phases': self._interpolation_total_phases,
                    }
                    snapshot = self._collect_block_snapshot()
                self._queue_worker_result({
                    'mode': 'single',
                    **step,
                    'snapshot': snapshot,
                    'step': worker_steps,
                    'elapsed': time.perf_counter() - started,
                    'error': None,
                })
            except Exception:
                self._queue_worker_result({
                    'mode': 'single',
                    'all_finished': False,
                    'phase_finished': False,
                    'phase_label': None,
                    'phase_iteration': 0,
                    'phase_iterations': max(int(self.config.get('iterations', 500)), 1),
                    'next_phase_label': None,
                    'completed_phases': self._interpolation_completed_phases,
                    'total_phases': self._interpolation_total_phases,
                    'changed': False,
                    'snapshot': None,
                    'step': 0,
                    'elapsed': time.perf_counter() - started,
                    'error': traceback.format_exc(),
                })

        self._set_status('Running next interpolation step...')
        self._worker = threading.Thread(target=_worker, daemon=True)
        self._worker.start()

    def request_update_all(self):
        if self._worker_running():
            self._set_status('Update already running...')
            return

        target_iterations = max(int(self.config.get('iterations', 500)), 1)

        def _worker():
            started = time.perf_counter()
            completed_steps = 0
            any_changed = False
            try:
                while True:
                    step = self._run_interpolators_once()
                    snapshot = self._collect_block_snapshot()
                    completed_steps += 1
                    any_changed = any_changed or step['changed']
                    self._queue_worker_result({
                        'mode': 'batch-progress',
                        **step,
                        'snapshot': snapshot,
                        'elapsed': time.perf_counter() - started,
                        'step': completed_steps,
                        'error': None,
                    })
                    if step['all_finished']:
                        break

                self._queue_worker_result({
                    'mode': 'batch-complete',
                    'changed': any_changed,
                    'snapshot': None,
                    'elapsed': time.perf_counter() - started,
                    'step': completed_steps,
                    'completed_phases': self._interpolation_completed_phases,
                    'total_phases': self._interpolation_total_phases,
                    'error': None,
                })
            except Exception:
                self._queue_worker_result({
                    'mode': 'batch-complete',
                    'changed': any_changed,
                    'snapshot': None,
                    'elapsed': time.perf_counter() - started,
                    'step': completed_steps,
                    'completed_phases': self._interpolation_completed_phases,
                    'total_phases': self._interpolation_total_phases,
                    'error': traceback.format_exc(),
                })

        self._set_status('Running all remaining interpolation steps...')
        self._worker = threading.Thread(target=_worker, daemon=True)
        self._worker.start()

    def _poll_worker_result(self):
        with self._result_condition:
            result = self._pending_result
            self._pending_result = None
            self._result_condition.notify_all()
        if result is None:
            return
        if result['error']:
            print(result['error'])
            self._set_status('Update failed. See console.')
            return
        if result['snapshot'] is not None:
            self._raw_snapshot = result['snapshot']
            self._refresh_render_data()
        if result['mode'] == 'single' and result['changed']:
            if result.get('all_finished') and not result.get('next_phase_label'):
                self._set_status(
                    f"Completed {result['phase_label']} in {result['elapsed']:.2f}s over {result.get('step', 0)} iterations. Interpolation sequence complete."
                )
            elif result.get('phase_finished'):
                next_suffix = f" Next: {result['next_phase_label']}." if result.get('next_phase_label') else ''
                self._set_status(
                    f"Completed {result['phase_label']} at {result['phase_iteration']}/{result['phase_iterations']} in {result['elapsed']:.2f}s over {result.get('step', 0)} iterations.{next_suffix}"
                )
            else:
                self._set_status(
                    f"Advanced {result['phase_label']} to {result['phase_iteration']}/{result['phase_iterations']} in {result['elapsed']:.2f}s."
                )
        elif result['mode'] == 'single':
            if result.get('all_finished'):
                self._set_status(f"Interpolation sequence already complete ({result['elapsed']:.2f}s).")
            elif result.get('phase_finished') and result.get('next_phase_label'):
                self._set_status(
                    f"Completed {result['phase_label']} at {result['phase_iteration']}/{result['phase_iterations']} in {result['elapsed']:.2f}s over {result.get('step', 0)} iterations. Next: {result['next_phase_label']}."
                )
            elif result.get('phase_finished'):
                self._set_status(
                    f"Completed {result['phase_label']} at {result['phase_iteration']}/{result['phase_iterations']} in {result['elapsed']:.2f}s over {result.get('step', 0)} iterations."
                )
            else:
                self._set_status(
                    f"{result['phase_label']} {result['phase_iteration']}/{result['phase_iterations']} made no visible changes ({result['elapsed']:.2f}s)."
                )
        elif result['mode'] == 'batch-progress':
            block_count = 0 if self._raw_snapshot is None else len(self._raw_snapshot['positions'])
            self._set_status(
                f"{result['phase_label']} {result['phase_iteration']}/{result['phase_iterations']} refreshed in {result['elapsed']:.2f}s. Blocks: {block_count}."
            )
        elif result['mode'] == 'batch-complete':
            self._set_status(
                f"All-iterations update finished in {result['elapsed']:.2f}s. Completed {result.get('completed_phases', 0)}/{result.get('total_phases', 0)} phases over {result.get('step', 0)} steps."
            )

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

    def _overlay_window_height(self, line_count, min_height=0.12, max_height=0.42):
        estimated = 0.045 + (0.021 * max(1, int(line_count)))
        return max(min_height, min(max_height, estimated))

    def _overlay_layout(self):
        panel_x = 0.015
        panel_y = 0.015
        panel_width = 0.32
        panel_gap = 0.012
        info_line_count = 7 + (1 if self._last_status else 0)
        legend_title, legend_entries, extra_count = self._get_overlay_entries()
        legend_line_count = 1 + len(legend_entries) + (1 if extra_count else 0)
        info_height = self._overlay_window_height(info_line_count, min_height=0.20, max_height=0.26)
        legend_height = self._overlay_window_height(legend_line_count, min_height=0.14, max_height=0.34)
        return {
            'Viewer Panel': (panel_x, panel_y, panel_width, info_height),
            'Viewer Legend': (panel_x, panel_y + info_height + panel_gap, panel_width, legend_height),
        }

    def _sync_overlay_bounds_from_ini(self):
        now = time.monotonic()
        if now - self._overlay_ini_last_check < 0.25:
            return
        self._overlay_ini_last_check = now
        if not self._overlay_ini_path or not os.path.isfile(self._overlay_ini_path):
            return

        try:
            current_mtime = os.path.getmtime(self._overlay_ini_path)
        except OSError:
            return

        if self._overlay_ini_mtime is not None and current_mtime <= self._overlay_ini_mtime:
            return

        parsed_bounds = {}
        current_window = None
        try:
            with open(self._overlay_ini_path, 'r', encoding='utf-8') as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if line.startswith('[Window][') and line.endswith(']'):
                        current_window = line[len('[Window]['):-1]
                        continue
                    if current_window not in {'Viewer Panel', 'Viewer Legend'}:
                        continue
                    if line.startswith('Pos='):
                        try:
                            x_str, y_str = line[4:].split(',', 1)
                            left = float(x_str)
                            top = float(y_str)
                        except ValueError:
                            continue
                        parsed_bounds.setdefault(current_window, {})['left'] = left
                        parsed_bounds.setdefault(current_window, {})['top'] = top
                    elif line.startswith('Size='):
                        try:
                            width_str, height_str = line[5:].split(',', 1)
                            width = float(width_str)
                            height = float(height_str)
                        except ValueError:
                            continue
                        parsed_bounds.setdefault(current_window, {})['width'] = width
                        parsed_bounds.setdefault(current_window, {})['height'] = height
        except OSError:
            return

        normalized_bounds = {}
        for name, values in parsed_bounds.items():
            if {'left', 'top', 'width', 'height'} <= values.keys():
                normalized_bounds[name] = (
                    values['left'],
                    values['top'],
                    values['width'],
                    values['height'],
                )

        if normalized_bounds:
            self._overlay_ini_bounds = normalized_bounds
            self._overlay_ini_mtime = current_mtime

    def _overlay_input_bounds(self, window_shape):
        width = max(float(window_shape[0]), 1.0)
        height = max(float(window_shape[1]), 1.0)
        layout = self._overlay_layout()
        bounds = {
            name: (left * width, top * height, panel_width * width, panel_height * height)
            for name, (left, top, panel_width, panel_height) in layout.items()
        }
        if self._overlay_ini_bounds:
            bounds.update(self._overlay_ini_bounds)
        return bounds

    def _cursor_over_overlay(self, cursor_x, cursor_y, window_shape):
        if not self.show_legend_overlay:
            return False
        cursor_px_x = float(cursor_x) * float(window_shape[0])
        cursor_px_y = (1.0 - float(cursor_y)) * float(window_shape[1])
        for left, top, width, height in self._overlay_input_bounds(window_shape).values():
            if left <= cursor_px_x <= left + width and top <= cursor_px_y <= top + height:
                return True
        return False

    def _draw_overlay(self, gui):
        if not self.show_legend_overlay:
            return

        legend_title, legend_entries, extra_count = self._get_overlay_entries()
        layout = self._overlay_layout()
        panel_x, panel_y, panel_width, info_height = layout['Viewer Panel']
        _legend_x, legend_y, _legend_width, legend_height = layout['Viewer Legend']

        with gui.sub_window('Viewer Panel', panel_x, panel_y, panel_width, info_height):
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
            gui.text('Keys: b blocks, v samples, f faces, c mode, l legend, i update, Shift+I all', color=(0.74, 0.74, 0.74))

        with gui.sub_window('Viewer Legend', panel_x, legend_y, panel_width, legend_height):
            gui.text(legend_title, color=(0.95, 0.95, 0.95))
            for label, color in legend_entries:
                gui.text(f"■ {label}", color=self._overlay_text_color(color))
            if extra_count:
                gui.text(f"... and {extra_count} more", color=(0.74, 0.74, 0.74))

    def toggle_color_mode(self):
        if not self.domain_data_available:
            self._set_status('No domain data available for domain coloring.')
            return
        self.color_mode = 'domain' if self.color_mode == 'value' else 'value'
        self._refresh_render_data()
        self._print_legend_if_needed(force=True)
        self._set_status(f"Color mode: {self.color_mode}")

    def run(self):
        while True:
            result = self._run_native_scene()
            if result != 'restart':
                return


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
