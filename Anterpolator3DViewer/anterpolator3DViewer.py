import pandas as pd
import numpy as np
import csv
import ast
import atexit
import signal
from collections import Counter
import importlib.util
import io
import threading
import time
import tempfile
import subprocess
import traceback
import warnings
import math
import re
import hashlib
from decimal import Decimal, InvalidOperation
from tqdm import tqdm
import sys
import os
from datetime import datetime
import xml.etree.ElementTree as ET
from matplotlib.colors import ListedColormap
import json
from PyQt5 import QtWidgets, QtCore
sys.path.append("C:/Projects/Anterpolator")

from provenance_utils import (
    copy_interpolator_provenance,
    finalize_phase_provenance,
    get_export_provenance,
    normalize_provenance_algorithm_name,
    seed_original_sample_provenance,
    snapshot_interpolator_state,
)

pv = None
taichi_runtime_module = None
bmf_tools_module = None
LARGE_BLOCK_FILE_THRESHOLD = 512 * 1024 * 1024
INITIAL_BLOCK_RENDER_THRESHOLD = 5000
SAMPLE_BLOCK_CACHE_VERSION = 1
INVALID_FILENAME_CHARS = str.maketrans({ch: '_' for ch in '<>:"/\\|?*'})
_LOGGED_LEAPFROG_METADATA_SIGNATURES = set()


class LightweightBlocks(list):
    pass


def _require_pyvista():
    global pv
    if pv is None:
        import pyvista as _pyvista

        pv = _pyvista
    return pv


def _load_taichi_runtime_module():
    global taichi_runtime_module
    if taichi_runtime_module is not None:
        return taichi_runtime_module

    module_path = os.path.join(os.path.dirname(__file__), 'taichi', 'runtime.py')
    spec = importlib.util.spec_from_file_location('anterpolator_taichi_runtime', module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load Taichi runtime from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    taichi_runtime_module = module
    return taichi_runtime_module


def _get_bmf_tools_module():
    global bmf_tools_module
    if bmf_tools_module is None:
        import bmf_standalone_exporter as _bmf_tools

        bmf_tools_module = _bmf_tools
    return bmf_tools_module


def is_bmf_file(path):
    return str(path or '').strip().lower().endswith('.bmf')


def _load_bmf_dataframe(path, row_limit=None, progress_label=None, progress_callback=None):
    label = str(progress_label or 'Reading BMF file').strip()
    if progress_callback is not None:
        progress_callback(0, 100, f'{label}...')
    result = _get_bmf_tools_module().load_bmf_table(path, row_limit=row_limit)
    df = result.get('dataframe', pd.DataFrame()).copy()
    df._detected_delimiter = 'bmf'
    df.attrs['bmf_result'] = result
    if progress_callback is not None:
        progress_callback(100, 100, f'{label} complete.')
    return df, result


def _normalize_viewer_backend(viewer_backend):
    backend = str(viewer_backend or 'taichi').strip().lower()
    if backend == 'pyvista':
        print("PyVista viewer is deprecated. Falling back to Taichi viewer.")
        return 'taichi'
    return 'taichi'


def _normalize_taichi_block_render_mode(block_render_mode):
    return 'mesh'


def _write_json_atomic(path, data):
    directory = os.path.dirname(path) or '.'
    with tempfile.NamedTemporaryFile('w', suffix='.json', prefix='anterpolator_cfg_', dir=directory, delete=False, encoding='utf-8') as handle:
        json.dump(data, handle, indent=4)
        temp_path = handle.name
    os.replace(temp_path, path)


class BackgroundOperationWorker(QtCore.QObject):
    progress = QtCore.pyqtSignal(int, int, str)
    finished = QtCore.pyqtSignal(object)
    failed = QtCore.pyqtSignal(str, str)

    def __init__(self, operation, kwargs):
        super().__init__()
        self._operation = operation
        self._kwargs = dict(kwargs)

    @QtCore.pyqtSlot()
    def run(self):
        try:
            result = self._operation(progress_callback=self._emit_progress, **self._kwargs)
        except Exception as exc:
            self.failed.emit(str(exc), traceback.format_exc())
            return
        self.finished.emit(result)

    def _emit_progress(self, value, maximum=100, message=''):
        self.progress.emit(int(value or 0), max(int(maximum or 0), 1), str(message or 'Working...'))


class TrimmedDisplayDoubleSpinBox(QtWidgets.QDoubleSpinBox):
    def textFromValue(self, value):
        decimals = max(int(self.decimals()), 0)
        text = f'{float(value):.{decimals}f}'
        if '.' in text:
            text = text.rstrip('0').rstrip('.')
        if text in {'-0', '-0.0', ''}:
            text = '0'
        decimal_point = str(self.locale().decimalPoint())
        if decimal_point != '.':
            text = text.replace('.', decimal_point)
        return text

# --- Interpolator Factory ---
def create_interpolator(config, domain=None, current_algorithm=None):
    """
    Factory function to create appropriate interpolator based on config and domain.
    
    Parameters:
    -----------
    config : dict
        Configuration dictionary
    domain : str, optional
        Domain name for domain-specific algorithm selection
    current_algorithm : str, optional
        Override algorithm selection from UI
        
    Returns:
    --------
    InterpolatorBase : Instance of appropriate interpolator
    """
    # Determine which algorithm to use
    algo_type = current_algorithm if current_algorithm else config.get('algorithm', 'ant_colony')
    
    # Check for domain-specific override (ONLY if not explicitly overridden by current_algorithm)
    if not current_algorithm and domain and 'domain_algorithm_overrides' in config:
        domain_config = config['domain_algorithm_overrides'].get(domain, {})
        if 'algorithm' in domain_config:
            algo_type = domain_config['algorithm']
            # Merge domain-specific parameters
            if algo_type == 'molecular_clock' and 'molecular_clock_params' in config:
                mc_params = {**config['molecular_clock_params'], **domain_config}
                config = {**config, 'molecular_clock_params': mc_params}
    
    if algo_type == 'ant_colony':
        from ant_colony import AntColonyInterpolator

        ant_target = config.get('ant_colony_interpolate_target', 'value')
        return AntColonyInterpolator(
            range_size=config.get('range_size', 10),
            max_pheromone=config.get('max_pheromone', 150),
            ants_per_sample=config.get('ants_per_sample', 3),
            verbose=config.get('verbose', False),
            background_value=config.get('background_value', 0.0),
            background_distance=config.get('background_distance', None),
            average_with_blocks=config.get('average_with_blocks', False),
            avoid_visited_threshold_enabled=config.get('avoid_visited_threshold_enabled', False),
            avoid_visited_threshold=config.get('avoid_visited_threshold', 100),
            ants_sampling_percentage=config.get('ants_sampling_percentage', 100.0),
            pheromone_decay_rate=config.get('pheromone_decay_rate', 1),
            interpolation_target=ant_target,
        )
    
    elif algo_type == 'molecular_clock':
        from molecular_clock_interpolator import MolecularClockInterpolator

        mc_params = config.get('molecular_clock_params', {})
        return MolecularClockInterpolator(
            spatial_weight=mc_params.get('spatial_weight', 1.0),
            attr_weight=mc_params.get('attr_weight', 1.0),
            detect_multiple_events=mc_params.get('detect_multiple_events', True),
            branch_threshold=mc_params.get('branch_threshold', 2.0),
            min_samples_per_event=mc_params.get('min_samples_per_event', 3),
            max_samples_per_event=mc_params.get('max_samples_per_event', 1000),
            interpolation_method=mc_params.get('interpolation_method', 'linear'),
            fill_background=mc_params.get('fill_background', False),
            background_value=mc_params.get('background_value', 0.0),
            ancestor_depth_offset=mc_params.get('ancestor_depth_offset', 1.0),
            verbose=config.get('verbose', False)
        )

    elif algo_type == 'gaussian_kernel':
        from gaussian_kernel_interpolator import GaussianKernelInterpolator

        gk_params = config.get('gaussian_kernel_params', {})
        return GaussianKernelInterpolator(
            bandwidth=gk_params.get('bandwidth', 3.0),
            cutoff_sigma=gk_params.get('cutoff_sigma', 3.0),
            use_nearest_fallback=gk_params.get('use_nearest_fallback', True),
            fill_background=gk_params.get('fill_background', False),
            background_value=gk_params.get('background_value', 0.0),
            verbose=config.get('verbose', False),
        )

    elif algo_type == 'adaptive_octree':
        from octree_domain_interpolator import OctreeDomainInterpolator

        octree_params = config.get('adaptive_octree_params', {})
        return OctreeDomainInterpolator(
            output_mode=octree_params.get('output_mode', 'dense_blocks_cover'),
            max_levels=octree_params.get('max_levels', 0),
            support_density_alpha=octree_params.get('support_density_alpha', 0.0),
            include_dense_provenance=octree_params.get('include_dense_provenance', True),
            verbose=config.get('verbose', False),
        )

    
    elif algo_type == 'string_theory' or algo_type == 'net_connector':
        from string_theory_interpolator import StringTheoryInterpolator

        st_params = config.get('string_theory_params', {})
        # Fallback to root config if not in sub-dict (for backward compatibility or direct config usage)
        dist_thresh = st_params.get('distance_threshold', config.get('distance_threshold', 10.0))
        grade_diff = st_params.get('grade_difference', config.get('grade_difference', 1.0))
        connect_all = st_params.get('connect_to_all', config.get('connect_to_all', True))
        max_connections = st_params.get('max_connections', config.get('max_connections', 1))
        min_connections = st_params.get('min_connections', config.get('min_connections', 1))
        collision_policy = st_params.get('collision_policy', config.get('collision_policy', 'average'))
        processing_order = st_params.get('processing_order', config.get('processing_order', 'ascending'))
        
        filter_by_frequency = st_params.get('filter_by_frequency', config.get('filter_by_frequency', False))
        min_azimuth_freq = st_params.get('min_azimuth_freq', config.get('min_azimuth_freq', 10.0))
        min_dip_freq = st_params.get('min_dip_freq', config.get('min_dip_freq', 10.0))

        interpolate_target = st_params.get('interpolate_target', 'value')
        
        return StringTheoryInterpolator(
            distance_threshold=dist_thresh,
            grade_difference=grade_diff,
            connect_to_all=connect_all,
            max_connections=max_connections,
            min_connections=min_connections,
            collision_policy=collision_policy,
            processing_order=processing_order,
            filter_by_frequency=filter_by_frequency,
            min_azimuth_freq=min_azimuth_freq,
            min_dip_freq=min_dip_freq,
            interpolation_target=interpolate_target,
            verbose=config.get('verbose', False)
        )
    
    else:
        raise ValueError(f"Unknown algorithm type: {algo_type}")

# --- New helper utilities for configurable headers/mappings ---
def parse_header_line(path, delimiter, line_number):
    """Return list of header tokens from the specified 1-based line number.
    If line_number is 1 we still parse it manually for consistency.
    Whitespace trimmed; empty tokens filtered out.
    Raises ValueError if the line does not exist or no tokens found."""
    if not os.path.isfile(path):
        raise ValueError(f"File not found: {path}")
    if line_number < 1:
        raise ValueError("Header line number must be >= 1")
    if is_bmf_file(path):
        df, _ = _load_bmf_dataframe(path, row_limit=1)
        tokens = [str(column).strip() for column in df.columns if str(column).strip()]
        if not tokens:
            raise ValueError(f"BMF file '{os.path.basename(path)}' produced no columns.")
        return tokens
    return _get_bmf_tools_module().parse_header_line(path, delimiter, line_number)


def resolve_effective_csv_header_line(path, configured_line=1):
    if is_bmf_file(path):
        return 1
    return _get_bmf_tools_module().resolve_effective_csv_header_line(path, configured_line)


def parse_effective_header_line(path, delimiter, line_number):
    if is_bmf_file(path):
        return parse_header_line(path, delimiter, 1)
    return _get_bmf_tools_module().parse_effective_header_line(path, delimiter, line_number)


def sync_csv_header_line_widget(spin_box, path, configured_line=None):
    current_line = int(configured_line if configured_line is not None else spin_box.value())
    effective_line = int(resolve_effective_csv_header_line(path, current_line))
    if effective_line != current_line:
        blocker = QtCore.QSignalBlocker(spin_box)
        spin_box.setValue(effective_line)
        del blocker
    return effective_line


def _parse_metadata_numeric_values(text, numeric_type=float, stop_at_equals=False):
    value_text = str(text or '')
    if stop_at_equals:
        value_text = value_text.split('=', 1)[0]
    values = re.findall(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', value_text)
    return [numeric_type(value) for value in values]


def parse_leapfrog_block_metadata(path, max_lines=100):
    if not path or is_bmf_file(path):
        return {}
    return _get_bmf_tools_module().parse_leapfrog_block_metadata(path, max_lines=max_lines)


def _format_metadata_vector(metadata, key):
    values = metadata.get(key)
    if not values:
        return None
    return '(' + ', '.join(f'{float(value):g}' for value in values) + ')'


def log_leapfrog_metadata_summary(path, metadata, context=''):
    if not metadata:
        return

    recognized_fields = [
        key for key in (
            'rotation_type',
            'azimuth_degrees',
            'dip_degrees',
            'pitch_degrees',
            'parent_block_size',
            'size_in_parent_blocks',
            'minimum_parent_centroid',
            'maximum_parent_centroid',
            'minimum_corner',
            'maximum_corner',
            'subblock_scheme',
            'subblock_factors',
        )
        if key in metadata
    ]
    if not recognized_fields:
        return

    signature = (os.path.abspath(path), str(context or 'default'), tuple(recognized_fields))
    if signature in _LOGGED_LEAPFROG_METADATA_SIGNATURES:
        return
    _LOGGED_LEAPFROG_METADATA_SIGNATURES.add(signature)

    summary_parts = []
    parent_size = _format_metadata_vector(metadata, 'parent_block_size')
    if parent_size:
        summary_parts.append(f'parent block size={parent_size}')
    grid_size = metadata.get('size_in_parent_blocks')
    if grid_size:
        summary_parts.append(f'parent grid={tuple(int(value) for value in grid_size)}')
    min_corner = _format_metadata_vector(metadata, 'minimum_corner')
    if min_corner:
        summary_parts.append(f'minimum corner={min_corner}')
    min_centroid = _format_metadata_vector(metadata, 'minimum_parent_centroid')
    if min_centroid:
        summary_parts.append(f'minimum parent centroid={min_centroid}')
    if any(key in metadata for key in ('azimuth_degrees', 'dip_degrees', 'pitch_degrees')):
        summary_parts.append(
            'rotation=(azimuth {azimuth:g}, dip {dip:g}, pitch {pitch:g})'.format(
                azimuth=float(metadata.get('azimuth_degrees', 0.0) or 0.0),
                dip=float(metadata.get('dip_degrees', 0.0) or 0.0),
                pitch=float(metadata.get('pitch_degrees', 0.0) or 0.0),
            )
        )
    subblock_factors = metadata.get('subblock_factors')
    if subblock_factors:
        scheme = metadata.get('subblock_scheme') or 'sub-blocks'
        summary_parts.append(f'{scheme}={tuple(int(value) for value in subblock_factors)}')

    label = f' ({context})' if context else ''
    print(f"Recognized Leapfrog block metadata{label} in '{os.path.basename(path)}': " + '; '.join(summary_parts))


def warn_leapfrog_metadata_mismatch(path, message):
    warning_text = f"Leapfrog metadata mismatch in '{os.path.basename(path)}': {message}"
    print(f"WARNING: {warning_text}")
    warnings.warn(warning_text, RuntimeWarning, stacklevel=2)


def validate_leapfrog_metadata_against_block_data(path, metadata, observed_min_bounds, observed_max_bounds,
                                                   unified_dims, min_grid_index=None, max_grid_index=None):
    if not metadata:
        return

    tolerance = 1e-6
    observed_min_bounds = np.asarray(observed_min_bounds, dtype=float)
    observed_max_bounds = np.asarray(observed_max_bounds, dtype=float)
    unified_dims = np.asarray(unified_dims, dtype=float)

    minimum_corner = np.asarray(metadata.get('minimum_corner', []), dtype=float)
    maximum_corner = np.asarray(metadata.get('maximum_corner', []), dtype=float)
    size_in_parent_blocks = np.asarray(metadata.get('size_in_parent_blocks', []), dtype=int)
    minimum_parent_centroid = np.asarray(metadata.get('minimum_parent_centroid', []), dtype=float)

    if minimum_corner.shape == (3,) and np.any(observed_min_bounds < minimum_corner - tolerance):
        warn_leapfrog_metadata_mismatch(
            path,
            f"block coordinates extend below metadata minimum corner; data min={observed_min_bounds.tolist()}, "
            f"metadata minimum corner={minimum_corner.tolist()}.",
        )
    if maximum_corner.shape == (3,) and np.any(observed_max_bounds > maximum_corner + tolerance):
        warn_leapfrog_metadata_mismatch(
            path,
            f"block coordinates extend above metadata maximum corner; data max={observed_max_bounds.tolist()}, "
            f"metadata maximum corner={maximum_corner.tolist()}.",
        )
    if minimum_corner.shape == (3,) and maximum_corner.shape == (3,) and size_in_parent_blocks.shape == (3,):
        expected_maximum_corner = minimum_corner + size_in_parent_blocks.astype(float) * unified_dims
        if not np.allclose(expected_maximum_corner, maximum_corner, atol=1e-6, rtol=0.0):
            warn_leapfrog_metadata_mismatch(
                path,
                f"metadata maximum corner does not match minimum corner + parent grid size * block size; "
                f"expected {expected_maximum_corner.tolist()}, metadata has {maximum_corner.tolist()}.",
            )
    if minimum_corner.shape == (3,) and minimum_parent_centroid.shape == (3,):
        expected_minimum_centroid = minimum_corner + unified_dims / 2.0
        if not np.allclose(expected_minimum_centroid, minimum_parent_centroid, atol=1e-6, rtol=0.0):
            warn_leapfrog_metadata_mismatch(
                path,
                f"metadata minimum parent centroid does not match minimum corner + half block size; "
                f"expected {expected_minimum_centroid.tolist()}, metadata has {minimum_parent_centroid.tolist()}.",
            )
    if size_in_parent_blocks.shape == (3,) and min_grid_index is not None and max_grid_index is not None:
        min_grid_index = np.asarray(min_grid_index, dtype=int)
        max_grid_index = np.asarray(max_grid_index, dtype=int)
        if np.any(min_grid_index < 0) or np.any(max_grid_index >= size_in_parent_blocks):
            warn_leapfrog_metadata_mismatch(
                path,
                f"block coordinates map outside metadata parent grid; observed index range "
                f"{min_grid_index.tolist()} to {max_grid_index.tolist()}, metadata grid size={size_in_parent_blocks.tolist()}.",
            )

def prepare_csv_read_kwargs(source, **read_csv_kwargs):
    prepared = dict(read_csv_kwargs)
    delimiter = prepared.get('delimiter', prepared.get('sep', ','))
    if 'engine' not in prepared:
        if delimiter is None or not isinstance(delimiter, str) or len(delimiter) != 1:
            prepared['engine'] = 'python'
    if isinstance(source, str) and prepared.get('engine', 'c') != 'python':
        prepared.setdefault('memory_map', True)
    return prepared


def _emit_progress(progress_callback, value, maximum=100, message=''):
    if progress_callback is None:
        return
    max_value = max(int(maximum or 0), 1)
    bounded_value = max(0, min(int(value or 0), max_value))
    progress_callback(bounded_value, max_value, str(message or ''))


def _format_progress_message(message, value, maximum):
    text = str(message or '').strip() or 'Working...'
    max_value = max(int(maximum or 0), 1)
    current_value = max(0, min(int(value or 0), max_value))
    if '%' in text:
        return text
    return f'{text} ({(current_value / max_value) * 100.0:.1f}%)'


def _format_row_progress_label(progress_label, processed_rows):
    row_label = 'row' if int(processed_rows) == 1 else 'rows'
    return f'{progress_label} ({int(processed_rows):,} {row_label})'


def _make_scaled_progress_callback(progress_callback, start, end, default_message=''):
    start_value = int(start)
    end_value = int(end)
    span = max(end_value - start_value, 0)

    def _callback(current, total, message=''):
        total_value = max(int(total or 0), 1)
        current_value = max(0, min(int(current or 0), total_value))
        if span == 0:
            mapped_value = start_value
        else:
            mapped_value = start_value + int(round((current_value / total_value) * span))
        _emit_progress(progress_callback, mapped_value, 100, message or default_message)

    return _callback

class ProgressTextReader:
    def __init__(self, path, label, progress_callback=None):
        self.path = path
        self.label = label
        self.progress_callback = progress_callback
        self.total_bytes = max(os.path.getsize(path), 1)
        self._raw = open(path, 'rb')
        self._text = io.TextIOWrapper(self._raw, encoding='utf-8', errors='ignore', newline='')
        self._displayed_bytes = 0
        self._started_at = time.perf_counter()
        self._last_postfix_update = self._started_at - 1.0
        self._refresh_bytes = max(512 * 1024, min(16 * 1024 * 1024, max(self.total_bytes // 2000, 1)))
        self._monitor_interval = 0.25
        self._lock = threading.Lock()
        self._stop_event = threading.Event()
        self._pbar = tqdm(
            total=self.total_bytes,
            desc=label,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        )
        self._update_postfix(force=True)
        self._monitor_thread = threading.Thread(target=self._monitor_progress, daemon=True)
        self._monitor_thread.start()

    def _sync_progress(self, force_postfix=False):
        with self._lock:
            current_bytes = self._raw.tell()
            byte_delta = current_bytes - self._displayed_bytes
            if byte_delta > 0:
                self._pbar.update(byte_delta)
                self._displayed_bytes = current_bytes
            if force_postfix or byte_delta >= self._refresh_bytes:
                self._update_postfix(force=force_postfix)
            if self.progress_callback and (byte_delta > 0 or force_postfix):
                try:
                    self.progress_callback(self._displayed_bytes, self.total_bytes, self.label)
                except Exception:
                    pass

    def _monitor_progress(self):
        while not self._stop_event.wait(self._monitor_interval):
            self._sync_progress()

    def _update_postfix(self, force=False):
        now = time.perf_counter()
        if not force and (now - self._last_postfix_update) < 1.0:
            return
        elapsed_seconds = max(int(now - self._started_at), 0)
        self._pbar.set_postfix_str(f"elapsed~{elapsed_seconds}s")
        self._last_postfix_update = now

    def _track_text(self, text):
        return text

    @property
    def handle(self):
        return self._text

    def read(self, *args, **kwargs):
        return self._track_text(self._text.read(*args, **kwargs))

    def readline(self, *args, **kwargs):
        return self._track_text(self._text.readline(*args, **kwargs))

    def __iter__(self):
        return self

    def __next__(self):
        line = self.readline()
        if line == '':
            raise StopIteration
        return line

    def __getattr__(self, name):
        return getattr(self._text, name)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def close(self):
        try:
            self._stop_event.set()
            if self._monitor_thread.is_alive():
                self._monitor_thread.join(timeout=max(self._monitor_interval * 2, 0.5))
            self._sync_progress(force_postfix=True)
        finally:
            self._pbar.close()
            self._text.close()

def read_csv_with_progress(path, progress_label, progress_callback=None, **read_csv_kwargs):
    read_csv_kwargs = prepare_csv_read_kwargs(path, **read_csv_kwargs)
    read_csv_kwargs.pop('memory_map', None)
    with ProgressTextReader(path, progress_label, progress_callback=progress_callback) as reader:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', pd.errors.DtypeWarning)
            df = pd.read_csv(reader.handle, **read_csv_kwargs)
        reader._sync_progress(force_postfix=True)
        print(f"{progress_label}: read {reader._displayed_bytes:,} bytes from {os.path.basename(path)}")
        df._approx_bytes_read = reader._displayed_bytes
    return df

def iterate_csv_with_progress(path, progress_label, progress_callback=None, **read_csv_kwargs):
    read_csv_kwargs = prepare_csv_read_kwargs(path, **read_csv_kwargs)
    read_csv_kwargs.pop('memory_map', None)

    def _generator():
        with ProgressTextReader(path, progress_label, progress_callback=progress_callback) as reader:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', pd.errors.DtypeWarning)
                for chunk in pd.read_csv(reader.handle, **read_csv_kwargs):
                    yield chunk
            reader._sync_progress(force_postfix=True)
            print(f"{progress_label}: read {reader._displayed_bytes:,} bytes from {os.path.basename(path)}")

    return _generator()


def _estimate_csv_chunk_progress(path, header_line=1, sample_bytes=8 * 1024 * 1024):
    total_bytes = max(os.path.getsize(path), 1)
    skipped_lines = max(int(header_line or 1), 0)

    with open(path, 'rb') as handle:
        for _ in range(skipped_lines):
            if handle.readline() == b'':
                return total_bytes, total_bytes, None
        data_start_offset = handle.tell()
        sample = handle.read(sample_bytes)

    row_count = sample.count(b'\n')
    if row_count <= 0:
        return total_bytes, data_start_offset, None
    return total_bytes, data_start_offset, len(sample) / row_count


def iterate_csv_path_chunks_with_progress(path, progress_label, progress_callback=None, header_line=1, **read_csv_kwargs):
    read_csv_kwargs = prepare_csv_read_kwargs(path, **read_csv_kwargs)
    total_bytes, data_start_offset, bytes_per_row = _estimate_csv_chunk_progress(
        path,
        header_line=header_line,
    )

    def _generator():
        processed_rows = 0
        displayed_bytes = 0
        pbar = tqdm(
            total=total_bytes,
            desc=progress_label,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        )
        if progress_callback:
            try:
                progress_callback(
                    min(data_start_offset, total_bytes),
                    total_bytes,
                    _format_row_progress_label(progress_label, processed_rows),
                )
            except Exception:
                pass
        try:
            initial_bytes = min(data_start_offset, total_bytes)
            if initial_bytes > 0:
                pbar.update(initial_bytes)
                displayed_bytes = initial_bytes
            pbar.set_postfix_str(f'rows~{processed_rows:,}')

            with warnings.catch_warnings():
                warnings.simplefilter('ignore', pd.errors.DtypeWarning)
                for chunk in pd.read_csv(path, **read_csv_kwargs):
                    processed_rows += len(chunk)
                    message = _format_row_progress_label(progress_label, processed_rows)
                    estimated_bytes = None
                    if bytes_per_row is not None:
                        estimated_bytes = min(total_bytes, int(data_start_offset + processed_rows * bytes_per_row))
                        byte_delta = estimated_bytes - displayed_bytes
                        if byte_delta > 0:
                            pbar.update(byte_delta)
                            displayed_bytes = estimated_bytes
                    pbar.set_postfix_str(f'rows~{processed_rows:,}')
                    if progress_callback and estimated_bytes is not None:
                        try:
                            progress_callback(estimated_bytes, total_bytes, message)
                        except Exception:
                            pass
                    yield chunk

            if displayed_bytes < total_bytes:
                pbar.update(total_bytes - displayed_bytes)
                displayed_bytes = total_bytes
            pbar.set_postfix_str(f'rows~{processed_rows:,}')
            if progress_callback:
                try:
                    progress_callback(total_bytes, total_bytes, _format_row_progress_label(progress_label, processed_rows))
                except Exception:
                    pass
        finally:
            pbar.close()
        print(f"{progress_label}: processed {processed_rows:,} rows from {os.path.basename(path)}")

    return _generator()

def build_unique_column_names(headers):
    name_counts = {}
    final_names = []
    for header in headers:
        key = header if header else 'Unnamed'
        count = name_counts.get(key, 0)
        final_names.append(f"{key}_{count}" if count > 0 else key)
        name_counts[key] = count + 1
    return final_names

def read_selected_columns_with_header(path, delimiter, header_line, selected_columns, progress_label=None,
                                      progress_callback=None):
    effective_header_line = resolve_effective_csv_header_line(path, header_line)
    headers = parse_header_line(path, delimiter, effective_header_line)
    final_names = build_unique_column_names(headers)
    missing = [col for col in selected_columns if col not in final_names]
    if missing:
        raise ValueError(f"Selected columns not found in file '{os.path.basename(path)}': {missing}")
    if is_bmf_file(path):
        df, _ = _load_bmf_dataframe(path, progress_label=progress_label, progress_callback=progress_callback)
        df = df[selected_columns].copy()
        df._detected_delimiter = 'bmf'
        return df, final_names
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=effective_header_line,
        comment='#',
        usecols=selected_columns,
    )
    if progress_label:
        df = read_csv_with_progress(path, progress_label, progress_callback=progress_callback, **read_kwargs)
    else:
        df = pd.read_csv(path, **prepare_csv_read_kwargs(path, **read_kwargs))
    df = strip_leading_non_data_rows(df)
    df._detected_delimiter = delimiter
    return df, final_names

def load_samples_dataframe(samples_file, samples_delimiter=None, samples_header_line=1,
                          sample_x_col=None, sample_y_col=None, sample_z_col=None, sample_value_col=None,
                          sample_domain_col=None, sample_filters=None, progress_label=None,
                          extra_columns=None):
    explicit_mapping = all([sample_x_col, sample_y_col, sample_z_col, sample_value_col])
    if explicit_mapping:
        delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
        selected_columns = [sample_x_col, sample_y_col, sample_z_col, sample_value_col]
        selected_domain_col = str(sample_domain_col or '').strip()
        if selected_domain_col and selected_domain_col != '(None)':
            selected_columns.append(selected_domain_col)
        for column_name in extra_columns or []:
            normalized = str(column_name or '').strip()
            if normalized and normalized != '(None)':
                selected_columns.append(normalized)
        selected_columns.extend(collect_filter_fields(sample_filters))
        selected_columns = list(dict.fromkeys(selected_columns))
        df, parsed_cols = read_selected_columns_with_header(
            samples_file,
            delimiter,
            samples_header_line or 1,
            selected_columns,
            progress_label=progress_label,
        )
        if sample_filters:
            df, _ = apply_sample_filters(df, sample_filters=sample_filters)
        rename_map = {
            sample_x_col: 'x',
            sample_y_col: 'y',
            sample_z_col: 'z',
            sample_value_col: 'Value',
        }
        df = df.rename(columns=rename_map)
        return df, parsed_cols, rename_map

    if samples_header_line and samples_header_line != 1 and samples_delimiter:
        df, parsed_cols = read_csv_with_selected_header(
            samples_file,
            samples_delimiter,
            samples_header_line,
            expected_min_cols=4,
            progress_label=progress_label,
        )
    else:
        df = read_autodetect_csv(samples_file, forced_delimiter=samples_delimiter, progress_label=progress_label)
        parsed_cols = None
    if sample_filters:
        df, _ = apply_sample_filters(df, sample_filters=sample_filters)
    return df, parsed_cols, None


def normalize_selected_sample_domain_column(df, sample_domain_col=None):
    selected_domain_col = str(sample_domain_col or '').strip()
    if selected_domain_col and selected_domain_col != '(None)':
        if selected_domain_col not in df.columns:
            raise ValueError(f"Selected sample domain column not found in samples file: {selected_domain_col}")
        if selected_domain_col != 'Domain':
            df = df.copy()
            df['Domain'] = df[selected_domain_col]
        return df

    domain_like = None
    for column_name in df.columns:
        if str(column_name).strip().lower() == 'domain':
            domain_like = column_name
            break
    if domain_like and domain_like != 'Domain':
        df = df.rename(columns={domain_like: 'Domain'})
    return df


def normalize_selected_sample_weight_column(df, sample_weight_col=None, sample_value_col=None):
    selected_weight_col = str(sample_weight_col or '').strip()
    if not selected_weight_col or selected_weight_col == '(None)':
        return df

    selected_value_col = str(sample_value_col or '').strip()
    if selected_weight_col == 'Weight' and 'Weight' in df.columns:
        return df

    weight_source = selected_weight_col
    if weight_source not in df.columns:
        if selected_value_col and selected_weight_col == selected_value_col and 'Value' in df.columns:
            weight_source = 'Value'
        else:
            raise ValueError(f"Selected sample weight column not found in samples file: {selected_weight_col}")

    if weight_source == 'Weight':
        return df

    df = df.copy()
    df['Weight'] = df[weight_source]
    return df


def infer_sample_domains_from_blocks(sample_coords, blocks_file, block_size,
                                     blocks_delimiter=None, blocks_header_line=1,
                                     block_x_col=None, block_y_col=None, block_z_col=None,
                                     block_domain_col=None, block_filters=None, progress_callback=None):
    coords = np.asarray(sample_coords, dtype=float)
    if len(coords) == 0:
        return np.empty(0, dtype=object)
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('A valid blocks file is required to infer blank sample domains from blocks.')
    if block_size is None:
        raise ValueError('Block size must be specified to infer blank sample domains from blocks.')

    domain_column_name = str(block_domain_col or '').strip()
    if not domain_column_name or domain_column_name == '(None)':
        raise ValueError('A blocks domain column must be selected to infer blank sample domains from blocks.')

    delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
    block_metadata = load_large_blocks_metadata(
        blocks_file,
        delimiter,
        blocks_header_line or 1,
        block_size,
        None,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=domain_column_name,
        block_filters=block_filters,
        config=None,
        progress_callback=progress_callback,
    )

    sample_block_indices = _compute_sample_block_indices_from_metadata(coords, block_metadata)
    domain_mapping = block_metadata['domain_mapping']

    return np.array(
        [domain_mapping.get((int(idx[0]), int(idx[1]), int(idx[2])), '') for idx in sample_block_indices],
        dtype=object,
    )


def _compute_sample_block_indices_from_metadata(coords, block_metadata):
    coords_array = np.asarray(coords, dtype=float)
    if len(coords_array) == 0:
        return np.empty((0, 3), dtype=int)

    coords_for_mapping = coords_array
    if block_metadata.get('is_rotated'):
        rotation_center = block_metadata['rotation_center']
        rotation_matrix = block_metadata['rotation_matrix']
        coords_for_mapping = (coords_for_mapping - rotation_center) @ rotation_matrix.T

    all_min_bounds = np.asarray(block_metadata.get('grid_index_origin', block_metadata['all_min_bounds']), dtype=float)
    unified_dims = np.asarray(block_metadata['unified_dims'], dtype=float)
    return np.floor((coords_for_mapping - all_min_bounds) / unified_dims + 1e-6).astype(int)


def _normalize_block_transfer_columns(block_value_cols, block_x_col=None, block_y_col=None, block_z_col=None,
                                      block_domain_col=None):
    if isinstance(block_value_cols, str):
        raw_columns = [part.strip() for part in block_value_cols.split(',')]
    else:
        raw_columns = [str(part).strip() for part in (block_value_cols or [])]

    normalized_columns = []
    seen = set()
    reserved = {
        str(block_x_col or '').strip(),
        str(block_y_col or '').strip(),
        str(block_z_col or '').strip(),
        str(block_domain_col or '').strip(),
    }
    reserved = {value for value in reserved if value and value != '(None)'}

    for column_name in raw_columns:
        if not column_name or column_name == '(None)' or column_name in seen:
            continue
        if column_name in reserved:
            raise ValueError(
                f"Block transfer column '{column_name}' conflicts with the configured coordinate/domain mapping columns."
            )
        normalized_columns.append(column_name)
        seen.add(column_name)

    if not normalized_columns:
        raise ValueError('Please select at least one block column to transfer to samples.')

    return normalized_columns


def _detect_block_transfer_column_modes(blocks_file, delimiter, header_line, block_value_cols,
                                        block_x_col=None, block_y_col=None, block_z_col=None,
                                        block_filters=None, progress_callback=None):
    effective_header_line = resolve_effective_csv_header_line(blocks_file, header_line)
    headers = parse_header_line(blocks_file, delimiter, effective_header_line)
    final_names = build_unique_column_names(headers)
    filter_fields = collect_filter_fields(block_filters or [])
    selected_columns, rename_map, domain_copy_source, mapping_mode = plan_block_file_columns(
        final_names,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=None,
        extra_columns=list(block_value_cols) + filter_fields,
    )
    base_read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=effective_header_line,
        comment='#',
        usecols=selected_columns,
        chunksize=250_000,
    )
    if mapping_mode == 'explicit':
        print(f"Applied user block column mapping for transfer scan: {rename_map}")
    else:
        print(f"Applied generic positional mapping for transfer scan: {rename_map}")

    numeric_candidates = {column_name: True for column_name in block_value_cols}
    saw_values = {column_name: False for column_name in block_value_cols}
    for chunk in iterate_csv_with_progress(
        blocks_file,
        'Scanning block transfer columns',
        progress_callback=progress_callback,
        **base_read_kwargs,
    ):
        if block_filters:
            chunk, _ = apply_dataframe_filters(
                chunk,
                filters=block_filters,
                filter_subject='block',
                source_label='blocks file chunk',
                emit_logs=False,
            )
        chunk, _ = normalize_block_chunk(chunk, rename_map, domain_copy_source, extra_keep_columns=block_value_cols)
        if len(chunk) == 0:
            continue
        for column_name in block_value_cols:
            if column_name not in chunk.columns or not numeric_candidates[column_name]:
                continue
            values = chunk[column_name].fillna('').astype(str).str.strip()
            nonblank = values[(values != '') & (values.str.lower() != 'nan')]
            if len(nonblank) == 0:
                continue
            saw_values[column_name] = True
            if pd.to_numeric(nonblank, errors='coerce').isna().any():
                numeric_candidates[column_name] = False

    return {
        column_name: 'numeric' if numeric_candidates[column_name] and saw_values[column_name] else 'categorical'
        for column_name in block_value_cols
    }


def load_block_value_mappings(blocks_file, delimiter, header_line, block_size, block_value_cols,
                              block_x_col=None, block_y_col=None, block_z_col=None,
                              block_filters=None, block_metadata=None, progress_callback=None):
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')

    normalized_columns = _normalize_block_transfer_columns(
        block_value_cols,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
    )
    delimiter = delimiter or detect_csv_delimiter(blocks_file)
    header_line = header_line or 1

    block_metadata = block_metadata or load_large_blocks_metadata(
        blocks_file,
        delimiter,
        header_line,
        block_size,
        None,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=None,
        block_filters=block_filters,
        config=None,
        progress_callback=progress_callback,
    )

    column_modes = _detect_block_transfer_column_modes(
        blocks_file,
        delimiter,
        header_line,
        normalized_columns,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_filters=block_filters,
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 100, 'Scanning block transfer columns...'),
    )

    effective_header_line = resolve_effective_csv_header_line(blocks_file, header_line)
    headers = parse_header_line(blocks_file, delimiter, effective_header_line)
    final_names = build_unique_column_names(headers)
    filter_fields = collect_filter_fields(block_filters or [])
    selected_columns, rename_map, domain_copy_source, mapping_mode = plan_block_file_columns(
        final_names,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=None,
        extra_columns=normalized_columns + filter_fields,
    )
    base_read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=effective_header_line,
        comment='#',
        usecols=selected_columns,
        chunksize=250_000,
    )
    if mapping_mode == 'explicit':
        print(f"Applied user block column mapping for transfer aggregation: {rename_map}")
    else:
        print(f"Applied generic positional mapping for transfer aggregation: {rename_map}")

    numeric_summaries = {column_name: {} for column_name, mode in column_modes.items() if mode == 'numeric'}
    categorical_counts = {column_name: {} for column_name, mode in column_modes.items() if mode != 'numeric'}

    for chunk in iterate_csv_with_progress(
        blocks_file,
        'Reading grid file (block transfer mapping)',
        progress_callback=progress_callback,
        **base_read_kwargs,
    ):
        if block_filters:
            chunk, _ = apply_dataframe_filters(
                chunk,
                filters=block_filters,
                filter_subject='block',
                source_label='blocks file chunk',
                emit_logs=False,
            )
        chunk, _ = normalize_block_chunk(chunk, rename_map, domain_copy_source, extra_keep_columns=normalized_columns)
        if len(chunk) == 0:
            continue

        block_indices = _compute_sample_block_indices_from_metadata(
            chunk[['x', 'y', 'z']].to_numpy(copy=False),
            block_metadata,
        )

        for column_name, summaries in numeric_summaries.items():
            if column_name not in chunk.columns:
                continue
            numeric_values = pd.to_numeric(chunk[column_name], errors='coerce')
            valid_mask = numeric_values.notna().to_numpy()
            if not valid_mask.any():
                continue
            for idx, value in zip(block_indices[valid_mask], numeric_values.loc[valid_mask].to_numpy(copy=False)):
                block_idx = (int(idx[0]), int(idx[1]), int(idx[2]))
                current_sum, current_count = summaries.get(block_idx, (0.0, 0))
                summaries[block_idx] = (current_sum + float(value), current_count + 1)

        for column_name, counts_by_block in categorical_counts.items():
            if column_name not in chunk.columns:
                continue
            values = chunk[column_name].fillna('').astype(str).str.strip()
            valid_mask = ((values != '') & (values.str.lower() != 'nan')).to_numpy()
            if not valid_mask.any():
                continue
            for idx, value in zip(block_indices[valid_mask], values.loc[valid_mask].to_numpy(copy=False)):
                block_idx = (int(idx[0]), int(idx[1]), int(idx[2]))
                counts = counts_by_block.setdefault(block_idx, Counter())
                counts[str(value)] += 1

    value_mappings = {}
    for column_name, summaries in numeric_summaries.items():
        for block_idx, (total_value, count) in summaries.items():
            if count <= 0:
                continue
            value_mappings.setdefault(block_idx, {})[column_name] = total_value / count

    for column_name, counts_by_block in categorical_counts.items():
        for block_idx, counts in counts_by_block.items():
            if not counts:
                continue
            winner = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
            value_mappings.setdefault(block_idx, {})[column_name] = winner

    return {
        'value_mappings': value_mappings,
        'column_modes': column_modes,
        'block_metadata': block_metadata,
        'columns': normalized_columns,
    }


def apply_blank_sample_domain_behavior(df, blank_domain_behavior='skip', domain_col='Domain',
                                       x_col='x', y_col='y', z_col='z',
                                       blocks_file=None, blocks_delimiter=None, blocks_header_line=1,
                                       block_x_col=None, block_y_col=None, block_z_col=None,
                                       block_domain_col=None, block_size=None, block_filters=None,
                                       progress_callback=None):
    if domain_col not in df.columns:
        return df, {'blank_domains': 0, 'inferred_domains': 0, 'remaining_blank_domains': 0}

    updated_df = df.copy()
    stripped_domains = updated_df[domain_col].fillna('').astype(str).str.strip()
    blank_mask = updated_df[domain_col].isna() | (stripped_domains == '') | (stripped_domains.str.lower() == 'nan')
    blank_count = int(blank_mask.sum())
    updated_df[domain_col] = stripped_domains
    if blank_count == 0:
        return updated_df, {'blank_domains': 0, 'inferred_domains': 0, 'remaining_blank_domains': 0}

    inferred_count = 0
    behavior = str(blank_domain_behavior or 'skip').strip().lower()
    if behavior == 'infer_from_blocks':
        requirements_ready = bool(
            blocks_file and str(blocks_file).strip() and block_size is not None
            and str(block_domain_col or '').strip() and str(block_domain_col or '').strip() != '(None)'
        )
        if not requirements_ready:
            print(
                'Blank sample domain behavior is set to infer from blocks, but block metadata settings are incomplete. '
                'Blank-domain samples will be skipped.'
            )
        else:
            coord_frame = updated_df.loc[blank_mask, [x_col, y_col, z_col]].apply(pd.to_numeric, errors='coerce')
            valid_coord_mask = coord_frame.notna().all(axis=1)
            inferable_index = coord_frame.index[valid_coord_mask]
            skipped_invalid = int((~valid_coord_mask).sum())
            if skipped_invalid:
                print(
                    f'Skipping {skipped_invalid:,} blank-domain samples with invalid coordinates during domain inference.'
                )
            if len(inferable_index) > 0:
                print(f'Attempting to infer domains for {len(inferable_index):,} blank-domain samples from blocks...')
                inferred_domains = infer_sample_domains_from_blocks(
                    coord_frame.loc[inferable_index].to_numpy(copy=False),
                    blocks_file,
                    block_size,
                    blocks_delimiter=blocks_delimiter,
                    blocks_header_line=blocks_header_line,
                    block_x_col=block_x_col,
                    block_y_col=block_y_col,
                    block_z_col=block_z_col,
                    block_domain_col=block_domain_col,
                    block_filters=block_filters,
                    progress_callback=progress_callback,
                )
                inferred_series = pd.Series(inferred_domains, index=inferable_index, dtype=object).fillna('').astype(str).str.strip()
                inferred_mask = (inferred_series != '') & (inferred_series.str.lower() != 'nan')
                if inferred_mask.any():
                    updated_df.loc[inferred_series.index[inferred_mask], domain_col] = inferred_series.loc[inferred_mask]
                    inferred_count = int(inferred_mask.sum())
                    print(f'Inferred domains for {inferred_count:,} blank-domain samples from blocks.')
                unresolved_count = int((~inferred_mask).sum())
                if unresolved_count:
                    print(f'Could not infer domains for {unresolved_count:,} blank-domain samples; they will be skipped.')

    final_domains = updated_df[domain_col].fillna('').astype(str).str.strip()
    updated_df[domain_col] = final_domains
    remaining_blank_mask = (final_domains == '') | (final_domains.str.lower() == 'nan')
    return updated_df, {
        'blank_domains': int(blank_count),
        'inferred_domains': int(inferred_count),
        'remaining_blank_domains': int(remaining_blank_mask.sum()),
    }


def ensure_sample_domains_for_domain_operations(df, sample_domain_col=None, blank_domain_behavior='skip',
                                                x_col='x', y_col='y', z_col='z',
                                                blocks_file=None, blocks_delimiter=None, blocks_header_line=1,
                                                block_x_col=None, block_y_col=None, block_z_col=None,
                                                block_domain_col=None, block_size=None, block_filters=None,
                                                progress_callback=None):
    updated_df = normalize_selected_sample_domain_column(df, sample_domain_col=sample_domain_col)
    if 'Domain' not in updated_df.columns:
        print('No sample domain column available. Inferring domains for all samples from blocks...')
        updated_df = updated_df.copy()
        updated_df['Domain'] = ''
        coord_frame = updated_df[[x_col, y_col, z_col]].apply(pd.to_numeric, errors='coerce')
        valid_coord_mask = coord_frame.notna().all(axis=1)
        inferable_index = coord_frame.index[valid_coord_mask]
        skipped_invalid = int((~valid_coord_mask).sum())
        if skipped_invalid:
            print(f'Skipping {skipped_invalid:,} samples with invalid coordinates during full domain inference.')
        if len(inferable_index) > 0:
            inferred_domains = infer_sample_domains_from_blocks(
                coord_frame.loc[inferable_index].to_numpy(copy=False),
                blocks_file,
                block_size,
                blocks_delimiter=blocks_delimiter,
                blocks_header_line=blocks_header_line,
                block_x_col=block_x_col,
                block_y_col=block_y_col,
                block_z_col=block_z_col,
                block_domain_col=block_domain_col,
                block_filters=block_filters,
                progress_callback=progress_callback,
            )
            inferred_series = pd.Series(inferred_domains, index=inferable_index, dtype=object).fillna('').astype(str).str.strip()
            updated_df.loc[inferred_series.index, 'Domain'] = inferred_series
            inferred_count = int(((inferred_series != '') & (inferred_series.str.lower() != 'nan')).sum())
            print(f'Inferred domains for {inferred_count:,} samples from blocks.')

    return apply_blank_sample_domain_behavior(
        updated_df,
        blank_domain_behavior=blank_domain_behavior,
        domain_col='Domain',
        x_col=x_col,
        y_col=y_col,
        z_col=z_col,
        blocks_file=blocks_file,
        blocks_delimiter=blocks_delimiter,
        blocks_header_line=blocks_header_line,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=block_domain_col,
        block_size=block_size,
        block_filters=block_filters,
        progress_callback=progress_callback,
    )


def should_resolve_sample_domains_for_interpolation(wants_domain_any, blocks_file=None, block_domain_col=None):
    return bool(
        wants_domain_any
        or (
            blocks_file
            and str(block_domain_col or '').strip()
            and str(block_domain_col or '').strip() != '(None)'
        )
    )


def compute_domain_sensitive_assignment_mask(block_indices, allowed_grid, domain_mapping=None, sample_domains=None):
    indices_array = np.asarray(block_indices, dtype=int)
    allowed_mask = np.array([tuple(idx) in allowed_grid for idx in indices_array], dtype=bool)
    if sample_domains is None or domain_mapping is None:
        return allowed_mask

    domains_array = np.asarray(sample_domains, dtype=object)
    if len(domains_array) != len(indices_array):
        raise ValueError('Sample domains length must match block indices length for domain-sensitive assignment.')

    domain_match_mask = np.zeros(len(indices_array), dtype=bool)
    for index, idx in enumerate(indices_array):
        if not allowed_mask[index]:
            continue
        sample_domain = str(domains_array[index]).strip() if domains_array[index] is not None else ''
        if not sample_domain or sample_domain.lower() == 'nan':
            continue
        block_domain = str(domain_mapping.get(tuple(idx), '')).strip()
        domain_match_mask[index] = sample_domain == block_domain
    return allowed_mask & domain_match_mask


def _compute_sample_block_assignment_state(points_array, grid_origin, block_dims,
                                           allowed_grid=None, domain_mapping=None, sample_domains=None):
    block_indices = np.floor((points_array - grid_origin) / block_dims + 1e-6).astype(int)
    if allowed_grid is None:
        assigned_mask = np.ones(len(points_array), dtype=bool)
        domain_mismatch_count = 0
    else:
        allowed_mask = np.array([tuple(idx) in allowed_grid for idx in block_indices], dtype=bool)
        assigned_mask = compute_domain_sensitive_assignment_mask(
            block_indices,
            allowed_grid,
            domain_mapping=domain_mapping,
            sample_domains=sample_domains,
        )
        domain_mismatch_count = int(np.count_nonzero(allowed_mask & ~assigned_mask))
    return block_indices, assigned_mask, domain_mismatch_count


def _get_sample_block_cache_file_signature(path):
    normalized_path = os.path.abspath(path) if path else ''
    if not normalized_path:
        return {'path': '', 'exists': False}
    if not os.path.isfile(normalized_path):
        return {'path': normalized_path, 'exists': False}
    stat_result = os.stat(normalized_path)
    mtime_ns = getattr(stat_result, 'st_mtime_ns', int(stat_result.st_mtime * 1_000_000_000))
    return {
        'path': normalized_path,
        'exists': True,
        'size': int(stat_result.st_size),
        'mtime_ns': int(mtime_ns),
    }


def _normalize_sample_block_cache_block_size(block_size):
    if isinstance(block_size, np.ndarray):
        return [float(v) for v in block_size.tolist()]
    if isinstance(block_size, (list, tuple)):
        return [float(v) for v in block_size]
    return [float(block_size)] * 3


def _build_sample_block_cache_identity(config, mode='blocks_file'):
    cfg = config or {}
    return {
        'version': SAMPLE_BLOCK_CACHE_VERSION,
        'mode': str(mode or 'blocks_file'),
        'samples_file': os.path.abspath(str(cfg.get('samples_file') or '')) if cfg.get('samples_file') else '',
        'blocks_file': os.path.abspath(str(cfg.get('blocks_file') or '')) if cfg.get('blocks_file') else '',
        'samples_delimiter': str(cfg.get('samples_delimiter') or ''),
        'blocks_delimiter': str(cfg.get('blocks_delimiter') or ''),
        'samples_header_line': int(cfg.get('samples_header_line', 1) or 1),
        'blocks_header_line': int(cfg.get('blocks_header_line', 1) or 1),
        'sample_x_col': str(cfg.get('sample_x_col') or ''),
        'sample_y_col': str(cfg.get('sample_y_col') or ''),
        'sample_z_col': str(cfg.get('sample_z_col') or ''),
        'sample_value_col': str(cfg.get('sample_value_col') or ''),
        'sample_domain_col': str(cfg.get('sample_domain_col') or ''),
        'sample_weight_col': str(cfg.get('sample_weight_col') or ''),
        'block_x_col': str(cfg.get('block_x_col') or ''),
        'block_y_col': str(cfg.get('block_y_col') or ''),
        'block_z_col': str(cfg.get('block_z_col') or ''),
        'block_domain_col': str(cfg.get('block_domain_col') or ''),
        'blank_sample_domain_behavior': str(cfg.get('blank_sample_domain_behavior') or ''),
        'block_size': _normalize_sample_block_cache_block_size(cfg.get('block_size', ())),
        'sample_filters': [dict(entry) for entry in get_configured_sample_filters(cfg)],
        'block_filters': [dict(entry) for entry in get_configured_block_filters(cfg)],
        'domain_algorithm_overrides': {
            str(domain): dict(settings)
            for domain, settings in (cfg.get('domain_algorithm_overrides') or {}).items()
        },
        'subblock_domain_policy': str(cfg.get('subblock_domain_policy', 'majority') or 'majority'),
    }


def _build_sample_block_cache_manifest(config, points_count, sample_domains, sample_weights, mode='blocks_file'):
    identity = _build_sample_block_cache_identity(config, mode=mode)
    cfg = config or {}
    manifest = {
        'cache_version': SAMPLE_BLOCK_CACHE_VERSION,
        'identity': identity,
        'samples_file_signature': _get_sample_block_cache_file_signature(cfg.get('samples_file')),
        'blocks_file_signature': _get_sample_block_cache_file_signature(cfg.get('blocks_file')),
        'points_count': int(points_count),
        'uses_sample_domains': bool(sample_domains is not None),
        'uses_sample_weights': bool(sample_weights is not None),
    }
    return manifest


def _resolve_sample_block_cache_paths(config, mode='blocks_file'):
    cfg = config or {}
    samples_file = str(cfg.get('samples_file') or '').strip()
    if not samples_file:
        return None

    sample_dir = os.path.dirname(os.path.abspath(samples_file)) or os.getcwd()
    cache_dir = os.path.join(sample_dir, 'AnterpolatorCache')
    base_name = os.path.splitext(os.path.basename(samples_file))[0] or 'samples'
    identity = _build_sample_block_cache_identity(cfg, mode=mode)
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode('utf-8')).hexdigest()[:16]
    stem = f"{base_name}_sample_blocks_cache_{digest}"
    return {
        'cache_dir': cache_dir,
        'csv_path': os.path.join(cache_dir, f'{stem}.csv'),
        'manifest_path': os.path.join(cache_dir, f'{stem}.json'),
    }


def _load_sample_block_cache(config, points_count, sample_domains, sample_weights, mode='blocks_file'):
    paths = _resolve_sample_block_cache_paths(config, mode=mode)
    if not paths:
        return None
    if not os.path.isfile(paths['csv_path']) or not os.path.isfile(paths['manifest_path']):
        return None

    expected_manifest = _build_sample_block_cache_manifest(
        config,
        points_count,
        sample_domains,
        sample_weights,
        mode=mode,
    )
    try:
        with open(paths['manifest_path'], 'r', encoding='utf-8') as handle:
            stored_manifest = json.load(handle)
    except Exception as exc:
        print(f"Sample-block cache manifest could not be read; regenerating sample blocks. ({exc})")
        return None

    if stored_manifest != expected_manifest:
        return None

    try:
        cache_df = pd.read_csv(paths['csv_path'])
    except Exception as exc:
        print(f"Sample-block cache could not be read; regenerating sample blocks. ({exc})")
        return None

    required_columns = {'block_ix', 'block_iy', 'block_iz', 'value', 'sample_count', 'sample_weight_sum'}
    if not required_columns.issubset(cache_df.columns):
        print("Sample-block cache is missing required columns; regenerating sample blocks.")
        return None

    sample_block_values = {}
    sample_block_counts = {}
    sample_block_weight_sums = {}
    for row in cache_df.itertuples(index=False):
        block_idx = (int(row.block_ix), int(row.block_iy), int(row.block_iz))
        sample_block_values[block_idx] = float(row.value)
        sample_block_counts[block_idx] = int(row.sample_count)
        sample_block_weight_sums[block_idx] = float(row.sample_weight_sum)

    return {
        'paths': paths,
        'sample_block_values': sample_block_values,
        'sample_block_counts': sample_block_counts,
        'sample_block_weight_sums': sample_block_weight_sums,
    }


def _write_sample_block_cache(config, sample_block_values, sample_block_counts, sample_block_weight_sums,
                              points_count, sample_domains, sample_weights, mode='blocks_file'):
    paths = _resolve_sample_block_cache_paths(config, mode=mode)
    if not paths:
        return None

    os.makedirs(paths['cache_dir'], exist_ok=True)
    cache_rows = []
    for block_idx in sorted(sample_block_values):
        cache_rows.append({
            'block_ix': int(block_idx[0]),
            'block_iy': int(block_idx[1]),
            'block_iz': int(block_idx[2]),
            'value': float(sample_block_values[block_idx]),
            'sample_count': int(sample_block_counts.get(block_idx, 0)),
            'sample_weight_sum': float(sample_block_weight_sums.get(block_idx, 0.0)),
        })

    cache_df = pd.DataFrame(cache_rows)
    cache_df.to_csv(paths['csv_path'], index=False)
    _write_json_atomic(
        paths['manifest_path'],
        _build_sample_block_cache_manifest(
            config,
            points_count,
            sample_domains,
            sample_weights,
            mode=mode,
        ),
    )
    return paths


def _should_force_rebuild_sample_blocks(config):
    return bool((config or {}).get('force_rebuild_sample_blocks', False))


def aggregate_samples_into_blocks(points, values, grid_index_origin, unified_dims,
                                  allowed_grid=None, domain_mapping=None, sample_domains=None,
                                  sample_ids=None, sample_weights=None,
                                  progress_label='Assigning points to blocks',
                                  block_indices=None, assigned_mask=None, domain_mismatch_count=None):
    points_array = np.asarray(points, dtype=float)
    values_array = np.asarray(values, dtype=float)
    grid_origin = np.asarray(grid_index_origin, dtype=float)
    block_dims = np.asarray(unified_dims, dtype=float)

    if len(points_array) != len(values_array):
        raise ValueError('Points and values must have the same length for sample-block aggregation.')
    if sample_ids is not None and len(sample_ids) != len(points_array):
        raise ValueError('Sample IDs must have the same length as points for sample-block aggregation.')
    if sample_weights is not None and len(sample_weights) != len(points_array):
        raise ValueError('Sample weights must have the same length as points for sample-block aggregation.')

    weights_array = None
    if sample_weights is not None:
        weights_array = np.asarray(sample_weights, dtype=float)
        if np.any(~np.isfinite(weights_array)):
            raise ValueError('Sample weights must be finite numeric values for sample-block aggregation.')
        if np.any(weights_array <= 0.0):
            raise ValueError('Sample weights must be greater than zero for sample-block aggregation.')

    if block_indices is None or assigned_mask is None or domain_mismatch_count is None:
        block_indices, assigned_mask, domain_mismatch_count = _compute_sample_block_assignment_state(
            points_array,
            grid_origin,
            block_dims,
            allowed_grid=allowed_grid,
            domain_mapping=domain_mapping,
            sample_domains=sample_domains,
        )

    sample_block_sums = {}
    sample_block_counts = {}
    sample_block_weight_sums = {}
    sample_block_domain_counts = {}
    sample_block_ids = {}
    for index in tqdm(range(len(points_array)), desc=progress_label):
        if not assigned_mask[index]:
            continue

        block_idx = tuple(int(v) for v in block_indices[index])
        sample_weight = 1.0 if weights_array is None else float(weights_array[index])
        sample_block_sums[block_idx] = sample_block_sums.get(block_idx, 0.0) + float(values_array[index]) * sample_weight
        sample_block_counts[block_idx] = sample_block_counts.get(block_idx, 0) + 1
        sample_block_weight_sums[block_idx] = sample_block_weight_sums.get(block_idx, 0.0) + sample_weight

        if sample_domains is not None:
            sample_domain = str(sample_domains[index]).strip() if sample_domains[index] is not None else ''
            if sample_domain and sample_domain.lower() != 'nan':
                domain_counts = sample_block_domain_counts.setdefault(block_idx, {})
                domain_counts[sample_domain] = domain_counts.get(sample_domain, 0) + 1

        if sample_ids is not None:
            sample_id = '' if sample_ids[index] is None else str(sample_ids[index]).strip()
            if sample_id:
                sample_block_ids.setdefault(block_idx, []).append(sample_id)

    sample_block_values = {
        block_idx: sample_block_sums[block_idx] / float(sample_block_weight_sums[block_idx])
        for block_idx in sample_block_counts
        if sample_block_weight_sums.get(block_idx, 0.0) > 0.0
    }

    return {
        'block_indices': block_indices,
        'assigned_mask': assigned_mask,
        'domain_mismatch_count': domain_mismatch_count,
        'sample_block_values': sample_block_values,
        'sample_block_counts': sample_block_counts,
        'sample_block_weight_sums': sample_block_weight_sums,
        'sample_block_domain_counts': sample_block_domain_counts,
        'sample_block_ids': sample_block_ids,
    }


def _select_majority_sample_block_domains(sample_block_domain_counts):
    resolved_domains = {}
    for block_idx, domain_counts in (sample_block_domain_counts or {}).items():
        if not domain_counts:
            continue
        max_count = max(domain_counts.values())
        tied_domains = sorted(domain for domain, count in domain_counts.items() if count == max_count)
        if tied_domains:
            resolved_domains[block_idx] = tied_domains[0]
    return resolved_domains


def _build_sample_block_rows(sample_block_values, sample_block_counts, reference_origin, block_size,
                             rotation_matrix=None, rotation_center=None, domain_mapping=None,
                             value_column_name='Value', domain_column_name='Domain',
                             sample_count_column_name='Sample_Count', sample_ids_by_block=None,
                             sample_ids_column_name='SampleIDs'):
    rows = []
    origin = np.asarray(reference_origin, dtype=float)
    dims = np.asarray(block_size, dtype=float)
    for block_idx in sorted(sample_block_values):
        point = origin + (np.asarray(block_idx, dtype=float) + 0.5) * dims
        if rotation_matrix is not None and rotation_center is not None:
            point = point @ rotation_matrix + rotation_center

        row = {
            'x': float(point[0]),
            'y': float(point[1]),
            'z': float(point[2]),
            value_column_name: float(sample_block_values[block_idx]),
            sample_count_column_name: int(sample_block_counts.get(block_idx, 0)),
        }
        if sample_ids_by_block is not None:
            ordered_unique_ids = list(dict.fromkeys(sample_ids_by_block.get(block_idx, [])))
            row[sample_ids_column_name] = ' ; '.join(ordered_unique_ids)
        if domain_mapping is not None:
            domain_value = str(domain_mapping.get(block_idx, '')).strip()
            if domain_value:
                row[domain_column_name] = domain_value
        rows.append(row)
    return rows


def _load_block_assignment_metadata(blocks_file, block_size, blocks_delimiter=None, blocks_header_line=1,
                                    block_x_col=None, block_y_col=None, block_z_col=None,
                                    block_domain_col=None, block_filters=None, config=None,
                                    progress_callback=None):
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')
    if block_size is None:
        raise ValueError('Block size must be specified for sample-block export.')

    if not is_bmf_file(blocks_file):
        delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
        metadata = load_large_blocks_metadata(
            blocks_file,
            delimiter,
            blocks_header_line or 1,
            block_size,
            None,
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
            block_domain_col=block_domain_col,
            block_filters=block_filters,
            config=config,
            progress_callback=progress_callback,
        )
        metadata['allowed_grid'] = set(metadata['domain_mapping'].keys())
        return metadata

    df_blocks, _ = load_full_blocks_dataframe(
        blocks_file,
        blocks_delimiter=blocks_delimiter,
        blocks_header_line=blocks_header_line,
        block_filters=block_filters,
        progress_label='Reading blocks file',
        progress_callback=progress_callback,
    )
    block_x_col, block_y_col, block_z_col = resolve_block_coordinate_columns(
        list(df_blocks.columns),
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
    )
    coord_frame = df_blocks[[block_x_col, block_y_col, block_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1)
    if not valid_mask.any():
        raise ValueError('No valid block coordinates remain after numeric coercion.')

    aligned_coords = coord_frame.loc[valid_mask].to_numpy(copy=False)
    rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(aligned_coords, block_size_hint=block_size)
    if is_rotated:
        aligned_coords = (aligned_coords - rotation_center) @ rotation_matrix.T

    if isinstance(block_size, (list, tuple, np.ndarray)):
        unified_dims = np.asarray(block_size, dtype=float)
    else:
        unified_dims = np.array([block_size, block_size, block_size], dtype=float)

    all_min_bounds = np.min(aligned_coords, axis=0)
    all_max_bounds = np.max(aligned_coords, axis=0)
    grid_indices = np.floor((aligned_coords - all_min_bounds) / unified_dims + 1e-6).astype(int)
    grid_index_origin = all_min_bounds - (unified_dims / 2.0)

    if block_domain_col and str(block_domain_col).strip() and str(block_domain_col).strip() != '(None)':
        if block_domain_col not in df_blocks.columns:
            raise ValueError(f'Selected domain column not found in blocks file: {block_domain_col}')
        domains = df_blocks.loc[valid_mask, block_domain_col].fillna('Undomained').astype(str).str.strip()
        domains = domains.replace('', 'Undomained')
    else:
        domains = pd.Series(['Undomained'] * int(valid_mask.sum()))

    domain_mapping, subblock_counts, mixed_domain_blocks = resolve_base_block_domains(
        grid_indices,
        domains,
        policy='majority',
    )
    return {
        'all_min_bounds': np.array(all_min_bounds, copy=True),
        'all_max_bounds': np.array(all_max_bounds, copy=True),
        'unified_dims': np.array(unified_dims, copy=True),
        'domain_mapping': domain_mapping,
        'subblock_counts': subblock_counts,
        'mixed_domain_blocks': mixed_domain_blocks,
        'rotation_matrix': rotation_matrix,
        'rotation_center': rotation_center,
        'is_rotated': is_rotated,
        'grid_reference': np.array(grid_index_origin, copy=True),
        'min_quantized_idx': np.min(grid_indices, axis=0),
        'grid_index_origin': np.array(grid_index_origin, copy=True),
        'leapfrog_metadata': {},
        'source_blocks_header_line': blocks_header_line,
        'allowed_grid': set(domain_mapping.keys()),
    }

def build_domain_catalog_cache_signature(blocks_file, delimiter, header_line, domain_col, block_filters=None):
    stats = os.stat(blocks_file)
    return {
        'path': os.path.abspath(blocks_file),
        'size': int(stats.st_size),
        'mtime_ns': int(getattr(stats, 'st_mtime_ns', int(stats.st_mtime * 1_000_000_000))),
        'delimiter': delimiter,
        'header_line': int(header_line or 1),
        'domain_col': str(domain_col or ''),
        'block_filters': json.dumps(block_filters or [], sort_keys=True),
    }

def load_block_domain_catalog(blocks_file, delimiter, header_line, domain_col,
                              block_filters=None, chunksize=250_000, progress_callback=None):
    effective_header_line = resolve_effective_csv_header_line(blocks_file, header_line)
    headers = parse_header_line(blocks_file, delimiter, effective_header_line)
    final_names = build_unique_column_names(headers)
    if domain_col not in final_names:
        raise ValueError(f'Domain column "{domain_col}" not found in blocks file.')

    domains = set()
    filter_fields = collect_filter_fields(block_filters)
    selected_columns = list(dict.fromkeys([domain_col] + filter_fields))

    use_direct_read = (not block_filters) or (os.path.getsize(blocks_file) < LARGE_BLOCK_FILE_THRESHOLD)

    if use_direct_read:
        _emit_progress(progress_callback, 0, 100, 'Reading domain column...')
        df, _ = read_selected_columns_with_header(
            blocks_file,
            delimiter,
            effective_header_line,
            selected_columns,
        )
        if block_filters:
            _emit_progress(progress_callback, 80, 100, 'Applying block filters...')
            df, _ = apply_dataframe_filters(
                df,
                filters=block_filters,
                filter_subject='block',
                source_label='blocks file',
                emit_logs=False,
            )
        if domain_col not in df.columns or len(df) == 0:
            _emit_progress(progress_callback, 100, 100, 'Reading domain column complete.')
            return []
        values = df[domain_col].dropna().astype(str).str.strip()
        values = values[(values != '') & (values.str.lower() != 'nan')]
        _emit_progress(progress_callback, 100, 100, 'Reading domain column complete.')
        return sorted(set(values.tolist()))

    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=effective_header_line,
        comment='#',
        usecols=selected_columns,
        chunksize=chunksize,
    )

    for chunk in iterate_csv_path_chunks_with_progress(
        blocks_file,
        'Reading domain column',
        progress_callback=progress_callback,
        header_line=effective_header_line,
        **read_kwargs,
    ):
        if block_filters:
            chunk, _ = apply_dataframe_filters(
                chunk,
                filters=block_filters,
                filter_subject='block',
                source_label='blocks file chunk',
                emit_logs=False,
            )
        if domain_col not in chunk.columns or len(chunk) == 0:
            continue
        values = chunk[domain_col].dropna().astype(str).str.strip()
        values = values[(values != '') & (values.str.lower() != 'nan')]
        domains.update(values.tolist())

    return sorted(domains)

def quantize_grid_indices(coords, reference_origin, block_size):
    scaled = (coords - reference_origin) / block_size
    rounded = np.rint(scaled)
    return rounded.astype(np.int64)


def _series_is_all_empty(series):
    return series.isna().all() or (series.astype(str).str.strip() == '').all()


def _list_all_empty_columns(df):
    return [column_name for column_name in df.columns if _series_is_all_empty(df[column_name])]

def read_csv_with_selected_header(path, delimiter, header_line, expected_min_cols=1,
                                  progress_label=None, progress_callback=None,
                                  preserve_empty_columns=False):
    """Read CSV using a specific header line (1-based). Returns DataFrame.
    Uses manual header parsing to build names and then pandas read_csv with skiprows.
    """
    if is_bmf_file(path):
        df, _ = _load_bmf_dataframe(path, progress_label=progress_label, progress_callback=progress_callback)
        if df.shape[1] < expected_min_cols:
            raise ValueError(f"File '{path}' has fewer than {expected_min_cols} non-empty columns after cleanup.")
        df._detected_delimiter = 'bmf'
        return df, list(df.columns)

    effective_header_line = resolve_effective_csv_header_line(path, header_line)
    headers = parse_header_line(path, delimiter, effective_header_line)
    final_names = build_unique_column_names(headers)
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=effective_header_line - 1,
        comment='#'
    )
    if progress_label:
        df = read_csv_with_progress(path, progress_label, progress_callback=progress_callback, **read_kwargs)
    else:
        df = pd.read_csv(path, **prepare_csv_read_kwargs(path, **read_kwargs))
    if df.shape[0] and all(str(df.iloc[0, i]).strip() == final_names[i] for i in range(min(len(final_names), df.shape[1]))):
        df = df.iloc[1:].reset_index(drop=True)
    df = strip_leading_non_data_rows(df)
    empty_cols = _list_all_empty_columns(df)
    non_empty_column_count = int(df.shape[1] - len(empty_cols))
    if empty_cols and not preserve_empty_columns:
        df = df.drop(columns=empty_cols)
    if non_empty_column_count < expected_min_cols:
        raise ValueError(f"File '{path}' has fewer than {expected_min_cols} non-empty columns after cleanup.")
    df._detected_delimiter = delimiter
    return df, final_names

def detect_csv_delimiter(path):
    if is_bmf_file(path):
        return ','
    if not os.path.isfile(path):
        return ','  # default
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = []
            for _ in range(10):
                try:
                    line = next(f)
                except StopIteration:
                    break
                if line.strip() and not line.lstrip().startswith(('#','//')):
                    lines.append(line)
                if len(lines) >= 3:
                    break
            sample = ''.join(lines)
    except StopIteration:
        sample = ''
    counts = {d: sample.count(d) for d in [',',';','\t','|']}
    delim = max(counts, key=counts.get)
    return delim if counts[delim] > 0 else ','

def read_autodetect_csv(path, min_cols=1, forced_delimiter=None, progress_label=None,
            progress_callback=None, preserve_empty_columns=False):
    """Read a CSV file detecting delimiter (comma, semicolon, tab, pipe) unless forced.
    Drops empty columns. Returns DataFrame."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")
    if is_bmf_file(path):
        df, _ = _load_bmf_dataframe(path, progress_label=progress_label, progress_callback=progress_callback)
        if df.shape[1] < min_cols:
            raise ValueError(f"File '{path}' has fewer than {min_cols} non-empty columns after cleanup.")
        df._detected_delimiter = 'bmf'
        return df
    delim = forced_delimiter if forced_delimiter else detect_csv_delimiter(path)

    def base_read(delimiter):
        read_kwargs = dict(delimiter=delimiter, comment='#')
        if progress_label:
            return read_csv_with_progress(path, progress_label, progress_callback=progress_callback, **read_kwargs)
        return pd.read_csv(path, **prepare_csv_read_kwargs(path, **read_kwargs))

    try:
        df = base_read(delim)
    except Exception:
        alt = ';' if delim == ',' else (',' if delim == ';' else None)
        if alt:
            try:
                df = base_read(alt)
                delim = alt
            except Exception as e2:
                raise ValueError(f"Failed to read CSV with delimiter '{delim}' and fallback '{alt}': {e2}")
        else:
            raise

    if df.shape[1] == 1:
        text_sample = str(df.iloc[0, 0]) if len(df) else ''
        candidates = ['; ', ',', '\t', '|']
        for cand in candidates:
            if cand.strip() in text_sample and cand != delim:
                try:
                    df_try = base_read(cand)
                    if df_try.shape[1] > 1:
                        df = df_try
                        delim = cand
                        print(f"Reparsed '{os.path.basename(path)}' with delimiter '{cand}' -> {df.shape[1]} columns")
                        break
                except Exception:
                    pass
    df = strip_leading_non_data_rows(df)
    empty_cols = _list_all_empty_columns(df)
    non_empty_column_count = int(df.shape[1] - len(empty_cols))
    if empty_cols and not preserve_empty_columns:
        df = df.drop(columns=empty_cols)
    if non_empty_column_count < min_cols:
        raise ValueError(f"File '{path}' has fewer than {min_cols} non-empty columns after cleanup.")
    df._detected_delimiter = delim
    return df

def format_point_info(point, value):
    return f"Position: ({point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f})\nValue: {value:.2f}"

def strip_leading_non_data_rows(df):
    if df.empty:
        return df

    metadata_prefixes = ('variable descriptions', 'variable types', 'variable defaults')

    def _is_empty_value(value):
        if pd.isna(value):
            return True
        text = str(value).strip()
        return text == '' or text.lower() == 'nan'

    rows_to_drop = 0
    for row_index in range(len(df)):
        row = df.iloc[row_index]
        values = row.tolist()
        if all(_is_empty_value(value) for value in values):
            rows_to_drop += 1
            continue

        first_text = ''
        for value in values:
            if _is_empty_value(value):
                continue
            first_text = str(value).strip().lower().rstrip(':')
            break

        if any(first_text.startswith(prefix) for prefix in metadata_prefixes):
            rows_to_drop += 1
            continue

        break

    if rows_to_drop:
        df = df.iloc[rows_to_drop:].reset_index(drop=True)
    return df

def _collect_same_level_vectors(points, max_points=4096, max_points_per_level=256, neighbors_per_point=8):
    if len(points) == 0:
        return np.empty((0, 3), dtype=float)

    points = np.asarray(points, dtype=float)

    z_keys = np.round(points[:, 2], decimals=6)
    unique_levels, inverse = np.unique(z_keys, return_inverse=True)
    vectors = []

    max_levels = max(1, max_points // max_points_per_level)
    if len(unique_levels) <= max_levels:
        selected_level_indices = np.arange(len(unique_levels), dtype=int)
    else:
        selected_level_indices = np.linspace(0, len(unique_levels) - 1, num=max_levels, dtype=int)
        selected_level_indices = np.unique(selected_level_indices)

    for level_idx in selected_level_indices:
        level_points = points[inverse == level_idx]
        if len(level_points) < 2:
            continue

        order = np.lexsort((level_points[:, 1], level_points[:, 0]))
        level_points = level_points[order]

        if len(level_points) > max_points_per_level:
            level_points = level_points[:max_points_per_level]

        xy = level_points[:, :2]
        deltas = xy[:, None, :] - xy[None, :, :]
        distance_sq = np.einsum('ijk,ijk->ij', deltas, deltas)
        np.fill_diagonal(distance_sq, np.inf)

        neighbor_count = min(neighbors_per_point, len(level_points) - 1)
        if neighbor_count <= 0:
            continue

        nearest_indices = np.argsort(distance_sq, axis=1)[:, :neighbor_count]
        for point_idx, neighbor_indices in enumerate(nearest_indices):
            base_point = level_points[point_idx]
            for neighbor_idx in neighbor_indices:
                vector = level_points[neighbor_idx] - base_point
                if vector[0] != 0 or vector[1] != 0:
                    vectors.append(vector)

    if not vectors:
        return np.empty((0, 3), dtype=float)
    return np.asarray(vectors, dtype=float)

def middle_click_callback(plotter, event_id):
    if hasattr(plotter, 'picked_point'):
        point = plotter.picked_point
        if point is not None:
            plotter.camera_position = [
                plotter.camera.position,
                point,
                plotter.camera.up
            ]
            plotter.renderer.GetActiveCamera().SetFocalPoint(point)
            plotter.render()

def detect_grid_rotation(points, block_size_hint=None, sample_size=50000):
    """
    Detect rotation around Z-axis by analyzing vectors between adjacent blocks.
    Returns (rotation_matrix, center, is_rotated).
    """
    if len(points) < 10:
        return np.eye(3), np.zeros(3), False

    print(f"Detecting rotation on {len(points)} points...")

    vectors = _collect_same_level_vectors(points, max_points=min(len(points), sample_size, 4096))

    if len(vectors) == 0:
        print("No vectors found between neighbors on same Z level.")
        return np.eye(3), np.zeros(3), False
    norms = np.linalg.norm(vectors, axis=1)

    # 4. Determine target block size
    target_length = None
    
    # Try to use hint if valid
    if block_size_hint is not None:
        if isinstance(block_size_hint, (list, tuple, np.ndarray)):
            bs = float(block_size_hint[0])
        else:
            bs = float(block_size_hint)
        if bs > 0:
            target_length = bs
            print(f"Using provided block size hint: {target_length}")

    # If no hint or hint didn't yield vectors, try auto-detection
    if target_length is None:
        hist, bin_edges = np.histogram(norms, bins=50)
        peak_idx = np.argmax(hist)
        target_length = (bin_edges[peak_idx] + bin_edges[peak_idx+1]) / 2
        print(f"Auto-detected common block length: {target_length}")
    
    # Filter vectors matching target length
    # Allow 10% tolerance
    valid_mask = np.abs(norms - target_length) < (target_length * 0.1)
    grid_vectors = vectors[valid_mask]
    
    # If hint failed (e.g. sub-blocked model where neighbors are smaller), fallback to auto
    if len(grid_vectors) < 10 and block_size_hint is not None:
        print("Hint yielded few vectors, falling back to auto-detected length...")
        hist, bin_edges = np.histogram(norms, bins=50)
        peak_idx = np.argmax(hist)
        auto_length = (bin_edges[peak_idx] + bin_edges[peak_idx+1]) / 2
        if target_length and auto_length > (target_length * 1.5):
            print(
                "Auto-detected fallback length is much larger than the provided block size hint; "
                "treating rotation estimate as unreliable and skipping rotation correction."
            )
            return np.eye(3), np.zeros(3), False
        valid_mask = np.abs(norms - auto_length) < (auto_length * 0.1)
        grid_vectors = vectors[valid_mask]
        print(f"Fallback auto-detected length: {auto_length}")

    print(f"Found {len(grid_vectors)} grid vectors matching length.")
    
    if len(grid_vectors) == 0:
        print("No grid vectors found.")
        return np.eye(3), np.zeros(3), False

    # 5. Calculate angles of these vectors relative to East (0 degrees)
    angles = np.arctan2(grid_vectors[:, 1], grid_vectors[:, 0]) # -pi to pi
    
    # Convert to degrees for easier clustering
    angles_deg = np.degrees(angles) % 360
    
    # Fold into 0-90 range
    angles_folded = angles_deg % 90
    
    # Find the peak in 0-90
    hist_ang, bin_ang = np.histogram(angles_folded, bins=90, range=(0, 90))
    peak_ang_idx = np.argmax(hist_ang)
    grid_offset = (bin_ang[peak_ang_idx] + bin_ang[peak_ang_idx+1]) / 2
    
    print(f"Detected grid offset (mod 90): {grid_offset:.2f}")
    
    # Determine actual Azimuth
    candidates = [grid_offset, grid_offset+90, grid_offset+180, grid_offset+270]
    axis_match_tolerance_deg = 5.0
    
    best_y_angle = None
    
    for cand in candidates:
        cand_norm = cand % 360
        # Count vectors close to this axis in either direction.
        diff = np.abs(angles_deg - cand_norm)
        diff = np.minimum(diff, 360 - diff)
        opposite_angle = (cand_norm + 180) % 360
        diff_opposite = np.abs(angles_deg - opposite_angle)
        diff_opposite = np.minimum(diff_opposite, 360 - diff_opposite)
        count = np.sum((diff < axis_match_tolerance_deg) | (diff_opposite < axis_match_tolerance_deg))
        
        dist_to_north = abs(cand_norm - 90)
        if dist_to_north > 180: dist_to_north = 360 - dist_to_north
        
        # We prefer the axis that actually has vectors
        if count > len(grid_vectors) * 0.05: 
             # If we haven't picked one, or this one is closer to North (90)
             if best_y_angle is None:
                 best_y_angle = cand_norm
             else:
                 current_dist = abs(best_y_angle - 90)
                 if current_dist > 180: current_dist = 360 - current_dist
                 
                 if dist_to_north < current_dist:
                     best_y_angle = cand_norm
    
    if best_y_angle is None:
        print("No dominant axis found aligned with grid.")
        return np.eye(3), np.zeros(3), False

    print(f"Best Y-axis angle: {best_y_angle:.2f}")
    
    azimuth = (90 - best_y_angle) % 360
    print(f"Calculated Azimuth: {azimuth:.2f}")

    theta_deg = best_y_angle - 90
    theta_rad = np.radians(theta_deg)
    
    # Construct rotation matrix
    y_axis = np.array([np.cos(np.radians(best_y_angle)), np.sin(np.radians(best_y_angle)), 0])
    x_axis = np.array([y_axis[1], -y_axis[0], 0])
    z_axis = np.array([0, 0, 1])
    
    rotation_matrix = np.vstack([x_axis, y_axis, z_axis])

    aligned_diff = np.abs(angles_deg - best_y_angle)
    aligned_diff = np.minimum(aligned_diff, 360 - aligned_diff)
    opposite_aligned_diff = np.abs(angles_deg - ((best_y_angle + 180) % 360))
    opposite_aligned_diff = np.minimum(opposite_aligned_diff, 360 - opposite_aligned_diff)
    aligned_count = np.sum((aligned_diff < axis_match_tolerance_deg) | (opposite_aligned_diff < axis_match_tolerance_deg))
    alignment_ratio = aligned_count / float(len(grid_vectors))

    normalized_theta = ((theta_deg + 180.0) % 360.0) - 180.0
    min_significant_rotation_deg = 3.0
    min_axis_alignment_ratio = 0.10
    
    # Check significance
    if abs(normalized_theta) < min_significant_rotation_deg:
        return np.eye(3), np.zeros(3), False
    if alignment_ratio < min_axis_alignment_ratio:
        return np.eye(3), np.zeros(3), False
        
    center = np.mean(points, axis=0)
    return rotation_matrix, center, True

def resolve_base_block_domains(grid_indices, domains, policy='majority'):
    grouped_domains = {}
    for idx, domain in zip(grid_indices, domains):
        base_idx = tuple(int(v) for v in idx)
        grouped_domains.setdefault(base_idx, []).append(str(domain))

    return resolve_base_block_domains_from_counts(grouped_domains, policy=policy)

def resolve_base_block_domains_from_counts(grouped_domains, policy='majority'):
    domain_mapping = {}
    subblock_counts = {}
    mixed_domain_blocks = {}

    for base_idx, block_domains in grouped_domains.items():
        counts = block_domains if isinstance(block_domains, Counter) else Counter(block_domains)
        subblock_counts[base_idx] = int(sum(counts.values()))

        if len(counts) == 1:
            domain_mapping[base_idx] = next(iter(counts))
            continue

        sorted_counts = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        mixed_domain_blocks[base_idx] = dict(sorted_counts)

        if policy == 'strict':
            continue
        if policy != 'majority':
            raise ValueError(f"Unsupported sub-block domain policy: {policy}")

        domain_mapping[base_idx] = sorted_counts[0][0]

    if policy == 'strict' and mixed_domain_blocks:
        examples = list(mixed_domain_blocks.items())[:5]
        detail = '; '.join(f"{idx}: {counts}" for idx, counts in examples)
        raise ValueError(
            f"Found {len(mixed_domain_blocks)} base blocks with mixed sub-block domains under strict policy. "
            f"Examples: {detail}"
        )

    return domain_mapping, subblock_counts, mixed_domain_blocks

def plan_block_file_columns(header_names, block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                            extra_columns=None):
    available_cols = [c for c in header_names if str(c).strip() != '']
    if block_x_col and block_y_col and block_z_col:
        rename_map = {}
        for src, tgt in [(block_x_col, 'x'), (block_y_col, 'y'), (block_z_col, 'z')]:
            if src not in header_names:
                raise ValueError(f"Selected blocks column '{src}' not found in file.")
            rename_map[src] = tgt
        selected_columns = [block_x_col, block_y_col, block_z_col]
        domain_copy_source = None
        if block_domain_col:
            if block_domain_col not in header_names:
                raise ValueError(f"Selected domain column '{block_domain_col}' not found in blocks file.")
            if block_domain_col in selected_columns:
                domain_copy_source = block_domain_col
            else:
                selected_columns.append(block_domain_col)
                rename_map[block_domain_col] = 'Domain'
        mapping_mode = 'explicit'
    else:
        if len(available_cols) < 3:
            raise ValueError(f"Blocks file must provide at least 3 non-empty columns. Parsed headers: {header_names}")
        selected_columns = available_cols[:4]
        rename_map = {src: tgt for src, tgt in zip(selected_columns, ['x', 'y', 'z', 'Domain'])}
        domain_copy_source = None
        mapping_mode = 'positional'
    for column_name in extra_columns or []:
        if column_name not in header_names:
            raise ValueError(f"Selected block filter column '{column_name}' not found in blocks file.")
        if column_name not in selected_columns:
            selected_columns.append(column_name)
    return selected_columns, rename_map, domain_copy_source, mapping_mode

def normalize_block_chunk(chunk, rename_map, domain_copy_source=None, extra_keep_columns=None):
    if domain_copy_source and 'Domain' not in chunk.columns and domain_copy_source in chunk.columns:
        chunk['Domain'] = chunk[domain_copy_source]
    if rename_map:
        chunk = chunk.rename(columns=rename_map)
    keep_columns = [c for c in ['x', 'y', 'z', 'Domain'] if c in chunk.columns]
    for column_name in extra_keep_columns or []:
        if column_name in chunk.columns and column_name not in keep_columns:
            keep_columns.append(column_name)
    chunk = chunk.loc[:, keep_columns].copy()
    rows_before = len(chunk)
    for column in ['x', 'y', 'z']:
        if column not in chunk.columns:
            raise ValueError(f"Blocks file missing required coordinate column '{column}' after mapping.")
    coords = chunk[['x', 'y', 'z']].apply(pd.to_numeric, errors='coerce')
    valid_mask = coords.notna().all(axis=1)
    chunk = chunk.loc[valid_mask].copy()
    if len(chunk) > 0:
        chunk.loc[:, ['x', 'y', 'z']] = coords.loc[valid_mask].to_numpy(dtype=float, copy=False)
    chunk = chunk.astype({'x': float, 'y': float, 'z': float}, copy=False)
    return chunk, rows_before - len(chunk)

def load_large_blocks_metadata(blocks_file, delimiter, header_line, block_size, points,
                               block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                               block_filters=None, config=None, progress_callback=None):
    leapfrog_metadata = parse_leapfrog_block_metadata(blocks_file)
    if leapfrog_metadata:
        recognized_keys = sorted(key for key in leapfrog_metadata if key != 'raw_lines')
        print(f"Detected Leapfrog block metadata: {recognized_keys}")
        log_leapfrog_metadata_summary(blocks_file, leapfrog_metadata, context='streaming loader')

    if leapfrog_metadata.get('parent_block_size'):
        metadata_block_size = np.asarray(leapfrog_metadata['parent_block_size'], dtype=float)
        if block_size is not None:
            configured_block_size = np.asarray(
                block_size if isinstance(block_size, (list, tuple, np.ndarray)) else [block_size, block_size, block_size],
                dtype=float,
            )
            if configured_block_size.shape == (3,) and not np.allclose(configured_block_size, metadata_block_size):
                warn_leapfrog_metadata_mismatch(
                    blocks_file,
                    f"configured block size {configured_block_size.tolist()} differs from metadata parent block size "
                    f"{metadata_block_size.tolist()}; using metadata value.",
                )
        block_size = metadata_block_size.tolist()

    effective_header_line = resolve_effective_csv_header_line(blocks_file, header_line)
    headers = parse_header_line(blocks_file, delimiter, effective_header_line)
    final_names = build_unique_column_names(headers)
    effective_block_filters = [dict(entry) for entry in (block_filters or get_configured_block_filters(config))]
    filter_fields = collect_filter_fields(effective_block_filters)
    selected_columns, rename_map, domain_copy_source, mapping_mode = plan_block_file_columns(
        final_names,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=block_domain_col,
        extra_columns=filter_fields,
    )
    chunksize = 250_000
    skipped_domains = set()
    if config and 'domain_algorithm_overrides' in config:
        skipped_domains = {
            domain for domain, cfg in config['domain_algorithm_overrides'].items()
            if cfg.get('skip', False)
        }

    print(f"Streaming blocks file ({os.path.getsize(blocks_file) / 1024**3:.2f} GiB) with chunks of {chunksize:,} rows.")
    if mapping_mode == 'explicit':
        print(f"Applied user block column mapping: {rename_map}")
    else:
        print(f"Applied generic positional mapping to first columns: {rename_map}")
    print(f"Blocks selected headers: {selected_columns}")

    base_read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=effective_header_line,
        comment='#',
        usecols=selected_columns,
        chunksize=chunksize,
    )

    metadata_angles = [
        float(leapfrog_metadata.get(name, 0.0) or 0.0)
        for name in ('azimuth_degrees', 'dip_degrees', 'pitch_degrees')
    ]
    has_rotation_metadata = any(name in leapfrog_metadata for name in ('azimuth_degrees', 'dip_degrees', 'pitch_degrees'))
    if has_rotation_metadata and all(abs(angle) < 1e-9 for angle in metadata_angles):
        print("Leapfrog metadata reports zero azimuth/dip/pitch; skipping rotation sampling.")
        rotation_matrix = np.eye(3)
        rotation_center = np.zeros(3)
        is_rotated = False
        _emit_progress(progress_callback, 22, 100, 'Using Leapfrog rotation metadata...')
    else:
        rotation_sample_target = 10_000
        rotation_samples = []
        print("Reading grid file (rotation sample)")
        rotation_progress = _make_scaled_progress_callback(progress_callback, 0, 20, 'Reading grid file (rotation sample)')
        for chunk in iterate_csv_with_progress(
            blocks_file,
            "Reading grid file (rotation sample)",
            progress_callback=rotation_progress,
            **base_read_kwargs,
        ):
            if effective_block_filters:
                chunk, _ = apply_dataframe_filters(
                    chunk,
                    filters=effective_block_filters,
                    filter_subject='block',
                    source_label='blocks file chunk',
                    emit_logs=False,
                )
            chunk, dropped = normalize_block_chunk(chunk, rename_map, domain_copy_source, extra_keep_columns=filter_fields)
            if len(chunk) == 0:
                continue
            rotation_samples.append(chunk[['x', 'y', 'z']].to_numpy())
            if sum(len(array) for array in rotation_samples) >= rotation_sample_target:
                break

        if not rotation_samples:
            raise ValueError("Blocks file had no valid coordinate rows to sample for rotation detection.")

        sample_points = np.concatenate(rotation_samples, axis=0)
        _emit_progress(progress_callback, 22, 100, 'Detecting grid rotation...')
        rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(sample_points, block_size_hint=block_size)

    if block_size is not None:
        if isinstance(block_size, (list, tuple, np.ndarray)):
            unified_dims = np.array(block_size, dtype=float)
        else:
            unified_dims = np.array([block_size, block_size, block_size], dtype=float)
        print(f"Using configured block size: {unified_dims}")
    else:
        raise ValueError("Block size must be specified when streaming large blocks files.")

    metadata_grid_origin = None
    if leapfrog_metadata.get('minimum_corner'):
        metadata_grid_origin = np.asarray(leapfrog_metadata['minimum_corner'], dtype=float)
        if metadata_grid_origin.shape == (3,):
            print(f"Using Leapfrog minimum corner as parent-grid origin: {metadata_grid_origin.tolist()}")
        else:
            metadata_grid_origin = None

    print("Building bounds and domain mapping...")
    all_min_bounds = None
    all_max_bounds = None
    min_quantized_idx = None
    max_quantized_idx = None

    full_read_kwargs = dict(base_read_kwargs)

    grouped_domain_counts = {}
    skipped_count_total = 0
    resolved_rows = 0
    total_block_rows = 0
    filtered_block_rows = 0

    mapping_progress = _make_scaled_progress_callback(progress_callback, 25, 90, 'Reading grid file (bounds + domain mapping)')
    for chunk in iterate_csv_with_progress(
        blocks_file,
        "Reading grid file (bounds + domain mapping)",
        progress_callback=mapping_progress,
        **full_read_kwargs,
    ):
        total_block_rows += int(len(chunk))
        if effective_block_filters:
            chunk, _ = apply_dataframe_filters(
                chunk,
                filters=effective_block_filters,
                filter_subject='block',
                source_label='blocks file chunk',
                emit_logs=False,
            )
            filtered_block_rows += int(len(chunk))
        chunk, dropped = normalize_block_chunk(chunk, rename_map, domain_copy_source, extra_keep_columns=filter_fields)
        if len(chunk) == 0:
            continue

        coords = chunk[['x', 'y', 'z']].to_numpy(copy=False)
        if is_rotated:
            coords = (coords - rotation_center) @ rotation_matrix.T

        chunk_min = coords.min(axis=0)
        chunk_max = coords.max(axis=0)
        if all_min_bounds is None:
            all_min_bounds = chunk_min
            all_max_bounds = chunk_max
        else:
            all_min_bounds = np.minimum(all_min_bounds, chunk_min)
            all_max_bounds = np.maximum(all_max_bounds, chunk_max)

        if metadata_grid_origin is not None and not is_rotated:
            absolute_grid_indices = np.floor((coords - metadata_grid_origin) / unified_dims + 1e-6).astype(int)
        else:
            absolute_grid_indices = np.floor(coords / unified_dims + 1e-6).astype(int)
        chunk_min_idx = absolute_grid_indices.min(axis=0)
        chunk_max_idx = absolute_grid_indices.max(axis=0)
        if min_quantized_idx is None:
            min_quantized_idx = chunk_min_idx
            max_quantized_idx = chunk_max_idx
        else:
            min_quantized_idx = np.minimum(min_quantized_idx, chunk_min_idx)
            max_quantized_idx = np.maximum(max_quantized_idx, chunk_max_idx)

        if 'Domain' in chunk.columns:
            domains = chunk['Domain'].fillna("Undomained").astype(str).str.strip().replace('', "Undomained")
        else:
            domains = pd.Series(["Undomained"] * len(chunk))

        if skipped_domains:
            keep_mask = ~domains.isin(skipped_domains)
            skipped_count_total += int((~keep_mask).sum())
            absolute_grid_indices = absolute_grid_indices[keep_mask.to_numpy()]
            domains = domains[keep_mask].reset_index(drop=True)

        if len(domains) == 0:
            continue

        resolved_rows += len(domains)
        grouped = pd.DataFrame(
            {
                'ix': absolute_grid_indices[:, 0],
                'iy': absolute_grid_indices[:, 1],
                'iz': absolute_grid_indices[:, 2],
                'Domain': domains.to_numpy(copy=False),
            }
        ).groupby(['ix', 'iy', 'iz', 'Domain'], sort=False).size()

        for (ix, iy, iz, domain), count in grouped.items():
            base_idx = (int(ix), int(iy), int(iz))
            domain_counts = grouped_domain_counts.setdefault(base_idx, Counter())
            domain_counts[str(domain)] += int(count)

    if all_min_bounds is None or all_max_bounds is None or min_quantized_idx is None or max_quantized_idx is None:
        raise ValueError('Could not determine grid bounds from blocks file.')

    _emit_progress(progress_callback, 92, 100, 'Calculating grid dimensions...')
    observed_min_bounds = np.array(all_min_bounds, copy=True)
    observed_max_bounds = np.array(all_max_bounds, copy=True)
    validate_leapfrog_metadata_against_block_data(
        blocks_file,
        leapfrog_metadata,
        observed_min_bounds,
        observed_max_bounds,
        unified_dims,
        min_grid_index=min_quantized_idx if metadata_grid_origin is not None and not is_rotated else None,
        max_grid_index=max_quantized_idx if metadata_grid_origin is not None and not is_rotated else None,
    )
    if metadata_grid_origin is not None and not is_rotated:
        if leapfrog_metadata.get('size_in_parent_blocks'):
            dims_grid = np.asarray(leapfrog_metadata['size_in_parent_blocks'], dtype=int)
        else:
            dims_grid = (max_quantized_idx - min_quantized_idx + 1).astype(int)
        all_min_bounds = metadata_grid_origin.astype(float)
        min_quantized_idx = np.zeros(3, dtype=int)
    else:
        dims_grid = (max_quantized_idx - min_quantized_idx + 1).astype(int)
        all_min_bounds = min_quantized_idx.astype(float) * unified_dims
    all_max_bounds = all_min_bounds + dims_grid * unified_dims
    print('Calculated grid dimensions:', dims_grid)
    print("Finalizing domain mapping...")

    if np.any(min_quantized_idx != 0):
        shifted_grouped_domain_counts = {}
        min_quantized_idx_tuple = tuple(int(v) for v in min_quantized_idx)
        for base_idx, domain_counts in grouped_domain_counts.items():
            shifted_idx = tuple(int(base_idx[axis] - min_quantized_idx_tuple[axis]) for axis in range(3))
            shifted_grouped_domain_counts[shifted_idx] = domain_counts
        grouped_domain_counts = shifted_grouped_domain_counts

    if skipped_domains:
        print(f"Skipping domains: {skipped_domains}")
        if skipped_count_total > 0:
            print(f"Skipped {skipped_count_total:,} sub-block rows due to domain overrides.")

    if effective_block_filters:
        print(
            f"Block filtering complete: {filtered_block_rows:,} of {total_block_rows:,} rows satisfy all block filters."
        )

    subblock_domain_policy = 'majority'
    if config:
        subblock_domain_policy = str(config.get('subblock_domain_policy', 'majority')).strip().lower() or 'majority'

    domain_mapping, subblock_counts, mixed_domain_blocks = resolve_base_block_domains_from_counts(
        grouped_domain_counts,
        policy=subblock_domain_policy,
    )

    print(
        f"Resolved {resolved_rows:,} sub-block rows into {len(domain_mapping):,} base blocks "
        f"using '{subblock_domain_policy}' policy."
    )

    return {
        'all_min_bounds': all_min_bounds,
        'all_max_bounds': all_max_bounds,
        'unified_dims': unified_dims,
        'domain_mapping': domain_mapping,
        'subblock_counts': subblock_counts,
        'mixed_domain_blocks': mixed_domain_blocks,
        'rotation_matrix': rotation_matrix,
        'rotation_center': rotation_center,
        'is_rotated': is_rotated,
        'grid_reference': np.array(all_min_bounds, copy=True),
        'min_quantized_idx': np.array(min_quantized_idx, copy=True),
        'grid_index_origin': np.array(metadata_grid_origin if metadata_grid_origin is not None and not is_rotated else all_min_bounds, copy=True),
        'leapfrog_metadata': leapfrog_metadata,
        'source_blocks_header_line': effective_header_line,
    }

def create_blocks(points, values, block_size=10, verbose=False, range_size=10, max_pheromone=150,
                  ants_per_sample=3, blocks_file=None, background_value=0.0, background_distance=None, average_with_blocks=False,
                  blocks_delimiter=None,
                  avoid_visited_threshold_enabled=False,
                  avoid_visited_threshold=100,
                  blocks_header_line=1,
                  block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                  config=None,
                  sample_domains=None,
                  sample_weights=None,
                  build_visual_blocks=True):
    pv = _require_pyvista() if build_visual_blocks else None
    original_points_array = np.array(points, copy=True)
    original_values_array = np.array(values, copy=True)
    original_domains_array = np.array(sample_domains, copy=True) if sample_domains is not None else None
    original_weights_array = np.array(sample_weights, copy=True) if sample_weights is not None else None
    # Domain interpolation in String Theory must ignore blocks_file as input.
    st_target = None
    if config and config.get('algorithm') in ('string_theory', 'net_connector'):
        st_target = config.get('string_theory_params', {}).get('interpolate_target', 'value')
    if isinstance(st_target, str) and st_target.strip().lower() == 'domain':
        blocks_file = None

    if blocks_file and str(blocks_file).strip():
        print("Loading predefined cells from blocks_file...")
        load_pbar = tqdm(total=7, desc="Blocks loading phases", unit="phase")
        load_pbar.set_postfix_str("reading blocks file")
        block_filters = get_configured_block_filters(config)
        use_streaming_blocks = (not is_bmf_file(blocks_file)) and os.path.getsize(blocks_file) >= LARGE_BLOCK_FILE_THRESHOLD
        if use_streaming_blocks:
            if not blocks_delimiter:
                blocks_delimiter = detect_csv_delimiter(blocks_file)
            metadata = load_large_blocks_metadata(
                blocks_file,
                blocks_delimiter,
                blocks_header_line or 1,
                block_size,
                points,
                block_x_col=block_x_col,
                block_y_col=block_y_col,
                block_z_col=block_z_col,
                block_domain_col=block_domain_col,
                config=config,
            )
            all_min_bounds = metadata['all_min_bounds']
            all_max_bounds = metadata['all_max_bounds']
            unified_dims = metadata['unified_dims']
            domain_mapping = metadata['domain_mapping']
            subblock_counts = metadata['subblock_counts']
            mixed_domain_blocks = metadata['mixed_domain_blocks']
            rotation_matrix = metadata['rotation_matrix']
            rotation_center = metadata['rotation_center']
            is_rotated = metadata['is_rotated']
            grid_index_origin = metadata.get('grid_index_origin', all_min_bounds)
            source_blocks_header_line = metadata.get('source_blocks_header_line', blocks_header_line)
            if is_rotated:
                points = (np.asarray(points, dtype=float) - rotation_center) @ rotation_matrix.T
                print("Applied rotation to sample points.")
            load_pbar.update(4)
        elif (not is_bmf_file(blocks_file)) and blocks_header_line and blocks_header_line != 1 and blocks_delimiter:
            # Use custom header line parsing
            df_blocks, parsed_cols = read_csv_with_selected_header(
                blocks_file,
                blocks_delimiter,
                blocks_header_line,
                expected_min_cols=3,
                progress_label="Reading grid file",
            )
            print(f"Blocks file (custom header line {blocks_header_line}) parsed columns: {parsed_cols}")
        else:
            df_blocks = read_autodetect_csv(
                blocks_file,
                forced_delimiter=blocks_delimiter,
                progress_label="Reading grid file",
            )
            print(f"Blocks file delimiter used: '{df_blocks._detected_delimiter}'")
        if not use_streaming_blocks:
            # If only one column but delimiter appears inside, attempt alternate reparse
            if df_blocks.shape[1] == 1:
                sole_col = df_blocks.columns[0]
                sample_val = str(df_blocks.iloc[0, 0]) if len(df_blocks) else ''
                if (';' in sample_val and df_blocks._detected_delimiter != ';') or (',' in sample_val and df_blocks._detected_delimiter != ','):
                    alt = ';' if df_blocks._detected_delimiter == ',' else ','
                    try:
                        df_blocks = read_autodetect_csv(
                            blocks_file,
                            forced_delimiter=alt,
                            progress_label="Reading grid file",
                        )
                        print(f"Reparsed blocks file with alternate delimiter '{alt}' -> columns: {list(df_blocks.columns)}")
                        print(f"Predefined cells read from blocks_file: {len(df_blocks):,}")
                    except Exception as e:
                        print(f"Alternate delimiter parse failed: {e}")
            if block_filters:
                df_blocks, _ = apply_dataframe_filters(
                    df_blocks,
                    filters=block_filters,
                    filter_subject='block',
                    source_label='blocks file',
                )
            load_pbar.update(1)
            load_pbar.set_postfix_str("normalizing headers")
            # Robust header normalization for blocks
            # Apply explicit user mapping if provided
            if block_x_col and block_y_col and block_z_col:
                rename_map = {}
                for src, tgt in [(block_x_col, 'x'), (block_y_col, 'y'), (block_z_col, 'z')]:
                    if src not in df_blocks.columns:
                        raise ValueError(f"Selected blocks column '{src}' not found in file.")
                    rename_map[src] = tgt
                if block_domain_col:
                    if block_domain_col in [block_x_col, block_y_col, block_z_col]:
                        if block_domain_col in df_blocks.columns:
                            df_blocks['Domain'] = df_blocks[block_domain_col]
                    else:
                        if block_domain_col not in df_blocks.columns:
                            raise ValueError(f"Selected domain column '{block_domain_col}' not found in blocks file.")
                        rename_map[block_domain_col] = 'Domain'
                df_blocks = df_blocks.rename(columns=rename_map)
                print(f"Applied user block column mapping: {rename_map}")
            else:
                needed = ['x', 'y', 'z']
                if not all(n in df_blocks.columns for n in needed):
                    cols = [c for c in df_blocks.columns if str(c).strip() != '']
                    if len(cols) >= 3:
                        mapping = {}
                        for new, old in zip(['x', 'y', 'z'], cols[:3]):
                            mapping[old] = new
                        if len(cols) >= 4:
                            mapping[cols[3]] = 'Domain'
                        df_blocks = df_blocks.rename(columns=mapping)
                        print(f"Blocks header auto-mapped: {mapping}")
                cols = df_blocks.columns.tolist()
                if len(cols) >= 4 and not all(name in df_blocks.columns for name in ['x', 'y', 'z']):
                    positional_targets = ['x', 'y', 'z', 'Domain']
                    if any(cols[i] != positional_targets[i] for i in range(4)):
                        pos_map = {cols[i]: positional_targets[i] for i in range(4)}
                        df_blocks = df_blocks.rename(columns=pos_map)
                        print(f"Applied generic positional mapping to first four columns: {pos_map}")

            print(f"Blocks final headers: {list(df_blocks.columns)}")
            coord_before = len(df_blocks)
            for c in ['x', 'y', 'z']:
                if c in df_blocks.columns:
                    df_blocks[c] = pd.to_numeric(df_blocks[c], errors='coerce')
            df_blocks = df_blocks.dropna(subset=['x', 'y', 'z'])
            dropped = coord_before - len(df_blocks)
            if dropped > 0:
                print(f"Dropped {dropped} block rows with non-numeric coordinates.")
            if len(df_blocks) == 0:
                raise ValueError("All block rows have non-numeric coordinates after conversion.")
            missing_coords = [c for c in ['x', 'y', 'z'] if c not in df_blocks.columns]
            if missing_coords:
                raise ValueError(f"Blocks file missing required coordinate columns after mapping: {missing_coords}. Parsed headers: {list(df_blocks.columns)}")
            load_pbar.update(1)
            load_pbar.set_postfix_str("detecting rotation")
            centroids_all = df_blocks[['x', 'y', 'z']].values

            rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(centroids_all, block_size_hint=block_size)

            theta_rad = np.arctan2(rotation_matrix[1, 1], rotation_matrix[1, 0])
            theta_deg = np.degrees(theta_rad)
            azimuth = (90 - theta_deg) % 360

            if is_rotated:
                print(f"Detected rotated block model (Azimuth: {azimuth:.2f}°). Applying rotation correction to align with axes.")
                print(f"Rotation center: {rotation_center}")
                print(f"Rotation matrix:\n{rotation_matrix}")
                centroids_all = (centroids_all - rotation_center) @ rotation_matrix.T
                df_blocks['x'] = centroids_all[:, 0]
                df_blocks['y'] = centroids_all[:, 1]
                df_blocks['z'] = centroids_all[:, 2]
                points[:] = (points - rotation_center) @ rotation_matrix.T
                print("Applied rotation to sample points.")
            else:
                print(f"No significant rotation detected (Azimuth: {azimuth:.2f}°).")

            load_pbar.update(1)
            load_pbar.set_postfix_str("mapping domains")

            if block_size is not None:
                if isinstance(block_size, (list, tuple, np.ndarray)):
                    unified_dims = np.array(block_size, dtype=float)
                else:
                    unified_dims = np.array([block_size, block_size, block_size], dtype=float)
                print(f"Using configured block size: {unified_dims}")
            else:
                dims_all = []
                for col in ['x', 'y', 'z']:
                    uniq = np.sort(df_blocks[col].unique())
                    diffs = np.diff(uniq)
                    positive_diffs = diffs[diffs > 0]
                    dim = positive_diffs.min() if len(positive_diffs) > 0 else 10
                    dims_all.append(dim)
                unified_dims = np.array(dims_all)
                print(f"Auto-detected block dimensions: {unified_dims}")

            observed_min_centroid = np.min(centroids_all, axis=0)
            observed_max_centroid = np.max(centroids_all, axis=0)
            grid_index_origin = observed_min_centroid - (unified_dims / 2.0)
            all_min_bounds = np.array(grid_index_origin, copy=True)
            all_max_bounds = observed_max_centroid + (unified_dims / 2.0)
            dims_grid = np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int)
            print("Calculated grid dimensions:", dims_grid)
            print("Building domain mapping...")
            centroids = df_blocks[['x', 'y', 'z']].values
            grid_indices = np.floor((centroids - all_min_bounds) / unified_dims + 1e-6).astype(int)

            if 'Domain' in df_blocks.columns:
                domains = df_blocks['Domain'].fillna("Undomained").astype(str).str.strip()
                domains = domains.replace("", "Undomained")
            else:
                domains = pd.Series(["Undomained"] * len(df_blocks))

            if config and 'domain_algorithm_overrides' in config:
                skipped_domains = {domain for domain, cfg in config['domain_algorithm_overrides'].items()
                                 if cfg.get('skip', False)}
                if skipped_domains:
                    print(f"Skipping domains: {skipped_domains}")
                    keep_mask = ~domains.isin(skipped_domains)
                    skipped_count = int((~keep_mask).sum())
                    grid_indices = grid_indices[keep_mask.to_numpy()]
                    domains = domains[keep_mask].reset_index(drop=True)
                    if skipped_count > 0:
                        print(f"Skipped {skipped_count:,} sub-block rows due to domain overrides.")

            subblock_domain_policy = 'majority'
            if config:
                subblock_domain_policy = str(config.get('subblock_domain_policy', 'majority')).strip().lower() or 'majority'

            domain_mapping, subblock_counts, mixed_domain_blocks = resolve_base_block_domains(
                grid_indices,
                domains,
                policy=subblock_domain_policy,
            )
            source_blocks_header_line = blocks_header_line

            print(
                f"Resolved {len(domains):,} sub-block rows into {len(domain_mapping):,} base blocks "
                f"using '{subblock_domain_policy}' policy."
            )
        
        if mixed_domain_blocks:
            print(f"Warning: {len(mixed_domain_blocks):,} base blocks contain mixed sub-block domains.")
            preview = list(mixed_domain_blocks.items())[:5]
            for idx, counts in preview:
                chosen = domain_mapping.get(idx)
                print(f"  Base block {idx}: {counts} -> selected '{chosen}'")
            if len(mixed_domain_blocks) > len(preview):
                print(f"  ... and {len(mixed_domain_blocks) - len(preview):,} more mixed base blocks")

        allowed_grid = set(domain_mapping.keys())
        
        # Count blocks per domain for debugging
        domain_block_counts = {}
        for idx, domain in domain_mapping.items():
            domain_block_counts[domain] = domain_block_counts.get(domain, 0) + 1
        
        print(f"\nDomain block distribution (from blocks file):")
        for domain, count in sorted(domain_block_counts.items()):
            print(f"  {domain}: {count} blocks")
        print(f"  Total allowed blocks: {len(allowed_grid)}")
        load_pbar.update(1)
        load_pbar.set_postfix_str("assigning points")

        # Group sample points using the same grid origin used to index the block model.
        print("Assigning points to blocks...")
        block_indices, assigned_mask, domain_mismatch_count = _compute_sample_block_assignment_state(
            np.asarray(points, dtype=float),
            np.asarray(grid_index_origin, dtype=float),
            np.asarray(unified_dims, dtype=float),
            allowed_grid=allowed_grid,
            domain_mapping=domain_mapping,
            sample_domains=sample_domains,
        )
        sample_block_cache = None
        if not _should_force_rebuild_sample_blocks(config):
            sample_block_cache = _load_sample_block_cache(
                config,
                len(points),
                sample_domains,
                sample_weights,
                mode='blocks_file',
            )
        sample_block_data = None
        if sample_block_cache is not None:
            sample_block_data = {
                'block_indices': block_indices,
                'assigned_mask': assigned_mask,
                'domain_mismatch_count': domain_mismatch_count,
                'sample_block_values': sample_block_cache['sample_block_values'],
                'sample_block_counts': sample_block_cache['sample_block_counts'],
                'sample_block_weight_sums': sample_block_cache['sample_block_weight_sums'],
                'sample_block_domain_counts': {},
                'sample_block_ids': {},
            }
            print(
                f"Using cached sample blocks from {sample_block_cache['paths']['csv_path']} "
                f"({len(sample_block_data['sample_block_values']):,} blocks)."
            )
        else:
            sample_block_data = aggregate_samples_into_blocks(
                points,
                values,
                grid_index_origin,
                unified_dims,
                allowed_grid=allowed_grid,
                domain_mapping=domain_mapping,
                sample_domains=sample_domains,
                sample_weights=sample_weights,
                progress_label='Assigning points to blocks',
                block_indices=block_indices,
                assigned_mask=assigned_mask,
                domain_mismatch_count=domain_mismatch_count,
            )
            cache_paths = _write_sample_block_cache(
                config,
                sample_block_data['sample_block_values'],
                sample_block_data['sample_block_counts'],
                sample_block_data['sample_block_weight_sums'],
                len(points),
                sample_domains,
                sample_weights,
                mode='blocks_file',
            )
            if cache_paths is not None:
                print(
                    f"Saved sample-block cache to {cache_paths['csv_path']} "
                    f"({len(sample_block_data['sample_block_values']):,} blocks)."
                )
        sample_block_values = sample_block_data['sample_block_values']
        sample_block_counts = sample_block_data['sample_block_counts']
        if domain_mismatch_count:
            print(f"Rejected {domain_mismatch_count:,} samples whose domain does not match their target block domain.")
        
        # Count sample blocks per domain
        sample_domain_counts = {}
        for idx in sample_block_values.keys():
            domain = domain_mapping.get(idx, "Undomained")
            sample_domain_counts[domain] = sample_domain_counts.get(domain, 0) + 1
        
        print(f"\nSample block distribution (samples assigned to blocks):")
        for domain, count in sorted(sample_domain_counts.items()):
            print(f"  {domain}: {count} sample blocks")
        print(f"  Total sample blocks: {len(sample_block_values)}")
        if build_visual_blocks:
            load_pbar.update(1)
            load_pbar.set_postfix_str("creating blocks")

            block_data = []
            for idx in tqdm(sample_block_values.keys(), desc="Creating blocks"):
                corner = all_min_bounds + np.array(idx) * unified_dims
                cell = pv.Box(bounds=(
                    corner[0], corner[0] + unified_dims[0],
                    corner[1], corner[1] + unified_dims[1],
                    corner[2], corner[2] + unified_dims[2]
                ))
                avg_value = sample_block_values[idx]
                cell.cell_data['Value'] = np.full(cell.n_cells, avg_value)
                cell.cell_data['Raw_Value'] = np.full(cell.n_cells, avg_value)
                cell.cell_data['Is_Sample'] = np.full(cell.n_cells, True)
                cell.cell_data['Block_ID'] = np.full(cell.n_cells, 0)  # to be set later
                domain = domain_mapping.get(idx, "Undomained")
                cell.cell_data['Domain'] = np.full(cell.n_cells, domain)
                block_data.append(cell)
            multiblock = pv.MultiBlock(block_data)
        else:
            load_pbar.update(1)
            load_pbar.set_postfix_str("skipping visual blocks")
            print("Skipping sample block geometry materialization (interpolation-only mode).")
            multiblock = LightweightBlocks()
        block_info = {
            'min_bounds': all_min_bounds,
            'dims': np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int),
            'block_size': unified_dims.tolist(),
            'grid_index_origin': np.asarray(grid_index_origin, dtype=float).tolist(),
            'allowed_grid': list(allowed_grid),
            'rotation_matrix': rotation_matrix if is_rotated else None,
            'rotation_center': rotation_center if is_rotated else None,
            'domain_mapping': domain_mapping,
            'subblock_counts': subblock_counts,
            'mixed_domain_blocks': mixed_domain_blocks,
            'source_blocks_file': blocks_file,
            'source_blocks_delimiter': blocks_delimiter,
            'source_blocks_header_line': source_blocks_header_line,
            'source_block_x_col': block_x_col,
            'source_block_y_col': block_y_col,
            'source_block_z_col': block_z_col,
            'source_block_filters': [dict(entry) for entry in (block_filters or [])],
            'expand_interpolation_exports_to_subblocks': bool((config or {}).get('expand_interpolation_exports_to_subblocks', True)),
        }
        # Store metadata on multiblock with private-style names to avoid PyVista attribute restrictions
        multiblock._block_info = block_info
        multiblock._sample_blocks = dict(sample_block_values)
        multiblock._sample_assignment_data = {
            'points': original_points_array,
            'values': original_values_array,
            'domains': original_domains_array,
            'block_indices': block_indices,
            'assigned_mask': assigned_mask,
            'sample_block_counts': dict(sample_block_data['sample_block_counts']),
            'sample_block_weight_sums': dict(sample_block_data['sample_block_weight_sums']),
        }
        multiblock._sample_block_counts = dict(sample_block_data['sample_block_counts'])
        multiblock._sample_block_weight_sums = dict(sample_block_data['sample_block_weight_sums'])
        load_pbar.update(1)
        load_pbar.set_postfix_str("initializing interpolator")
        
        # Check if we should process domains sequentially with different algorithms
        process_sequentially = config and config.get('process_domains_sequentially', False)
        has_domain_overrides = config and 'domain_algorithm_overrides' in config and len(config['domain_algorithm_overrides']) > 0
        
        # We must split by domain if overrides exist (Scenario 3) OR if user explicitly requested sequential processing
        if has_domain_overrides or process_sequentially:
            print("\n=== Processing domains sequentially with different algorithms ===")
            # Group sample blocks by domain
            samples_by_domain = {}
            for idx, val in multiblock._sample_blocks.items():
                domain = domain_mapping.get(idx, "Undomained")
                if domain not in samples_by_domain:
                    samples_by_domain[domain] = {}
                samples_by_domain[domain][idx] = val
            
            # Store multiple interpolators (list per domain)
            multiblock._interpolators = {}
            
            for domain, domain_samples in samples_by_domain.items():
                # Determine which algorithm for this domain
                domain_config = config['domain_algorithm_overrides'].get(domain, {})
                
                # First Pass
                if domain_config.get('skip', False):
                    print(f"  Skipping domain: {domain}")
                    continue
                
                algo1 = domain_config.get('algorithm', config.get('algorithm', 'ant_colony'))
                algo2 = domain_config.get('second_pass_algorithm', 'skip')
                
                print(f"  Domain '{domain}': {len(domain_samples)} samples")
                print(f"    Pass 1: {algo1}")
                if algo2 != 'skip':
                    print(f"    Pass 2: {algo2}")
                
                domain_interpolators = []
                
                # --- Create Pass 1 ---
                domain_cfg1 = config.copy()
                domain_cfg1['algorithm'] = algo1
                interp1 = create_interpolator(domain_cfg1, domain=domain)
                
                # Filter allowed_grid to this domain only
                domain_allowed_grid = {pos for pos in allowed_grid if domain_mapping.get(pos) == domain}
                domain_sample_mapping = {pos: domain for pos in domain_samples}
                print(f"  Domain '{domain}' volume: {len(domain_allowed_grid)} blocks")
                
                # Check if we have any non-sample blocks in this domain (sparse vs dense definition)
                # If the allowed grid is just the samples, we shouldn't enforce it as a constraint.
                has_extra_blocks = len(domain_allowed_grid) > len(domain_samples) * 1.2
                print(f"  Domain '{domain}' sparse check: allowed={len(domain_allowed_grid)}, samples={len(domain_samples)}, has_extra={has_extra_blocks}")
                
                # Attach domain-specific attributes
                if has_extra_blocks:
                    interp1.allowed_grid_override = domain_allowed_grid
                    interp1.domain_mapping = {pos: domain for pos in domain_allowed_grid}
                    use_mapping = True
                else:
                    interp1.allowed_grid_override = None
                    interp1.domain_mapping = None
                    use_mapping = False
                
                # Initialize with only this domain's samples
                interp1.initialize_blocks(domain_samples, tuple(block_info['dims']),
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping,
                                             sample_domain_mapping=domain_sample_mapping,
                                             sample_block_counts={pos: multiblock._sample_block_counts.get(pos, 1) for pos in domain_samples},
                                             sample_block_weight_sums={pos: multiblock._sample_block_weight_sums.get(pos, 1.0) for pos in domain_samples})
                
                if hasattr(interp1, 'create_ants'):
                    interp1.create_ants()
                
                domain_interpolators.append(interp1)
                
                # --- Create Pass 2 (if configured) ---
                if algo2 != 'skip':
                    domain_cfg2 = config.copy()
                    domain_cfg2['algorithm'] = algo2
                    interp2 = create_interpolator(domain_cfg2, domain=domain, current_algorithm=algo2)
                    
                    if has_extra_blocks:
                        interp2.allowed_grid_override = domain_allowed_grid
                        interp2.domain_mapping = {pos: domain for pos in domain_allowed_grid}
                        use_mapping_2 = True
                    else:
                        interp2.allowed_grid_override = None
                        interp2.domain_mapping = None
                        use_mapping_2 = False
                    
                    # Initialize with original samples (will be augmented later)
                    interp2.initialize_blocks(domain_samples, tuple(block_info['dims']),
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping_2,
                                             sample_domain_mapping=domain_sample_mapping,
                                             sample_block_counts={pos: multiblock._sample_block_counts.get(pos, 1) for pos in domain_samples},
                                             sample_block_weight_sums={pos: multiblock._sample_block_weight_sums.get(pos, 1.0) for pos in domain_samples})
                    
                    if hasattr(interp2, 'create_ants'):
                        interp2.create_ants()
                    
                    domain_interpolators.append(interp2)
                
                multiblock._interpolators[domain] = domain_interpolators
            
            # Keep reference to primary interpolator (for compatibility)
            if multiblock._interpolators:
                # Just take the first pass of the first domain
                first_domain_list = list(multiblock._interpolators.values())[0]
                multiblock._ant_colony = first_domain_list[0]
            
        else:
            # Original single interpolator approach (potentially with 2 passes now)
            if config is None:
                config = {
                    'algorithm': 'ant_colony',
                    'range_size': range_size,
                    'max_pheromone': max_pheromone,
                    'ants_per_sample': ants_per_sample,
                    'verbose': verbose,
                    'background_value': background_value,
                    'background_distance': background_distance,
                    'average_with_blocks': average_with_blocks,
                    'avoid_visited_threshold_enabled': avoid_visited_threshold_enabled,
                    'avoid_visited_threshold': avoid_visited_threshold
                }
            
            # Check for second pass
            algo2 = config.get('second_pass_algorithm', 'skip')
            
            if algo2 != 'skip':
                print(f"Configuring undomained two-pass interpolation: {config.get('algorithm')} -> {algo2}")
                multiblock._interpolators = {}
                domain_interpolators = []
                
                # Check if we have any non-sample blocks (sparse vs dense definition)
                has_extra_blocks = len(allowed_grid) > len(multiblock._sample_blocks) * 1.2
                print(f"  Undomained (2-pass) sparse check: allowed={len(allowed_grid)}, samples={len(multiblock._sample_blocks)}, has_extra={has_extra_blocks}")
                
                # Pass 1
                interp1 = create_interpolator(config)
                algo_type = config.get('algorithm', 'ant_colony')
                
                if has_extra_blocks:
                    if algo_type == 'ant_colony':
                        interp1.allowed_grid_override = allowed_grid
                        interp1.domain_mapping = domain_mapping
                    elif algo_type == 'molecular_clock':
                        interp1.allowed_grid_override = allowed_grid
                        interp1.domain_mapping = domain_mapping
                    elif algo_type == 'gaussian_kernel':
                        interp1.allowed_grid_override = allowed_grid
                        interp1.domain_mapping = domain_mapping
                    elif algo_type == 'adaptive_octree':
                        interp1.allowed_grid_override = allowed_grid
                        interp1.domain_mapping = domain_mapping
                    elif algo_type == 'string_theory':
                        interp1.allowed_grid_override = allowed_grid
                        interp1.domain_mapping = domain_mapping
                    use_mapping_1 = True
                else:
                    interp1.allowed_grid_override = None
                    interp1.domain_mapping = None
                    use_mapping_1 = False
                
                interp1.initialize_blocks(multiblock._sample_blocks, tuple(block_info['dims']),
                                         all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping_1,
                                         sample_block_counts=getattr(multiblock, '_sample_block_counts', {}),
                                         sample_block_weight_sums=getattr(multiblock, '_sample_block_weight_sums', {}))
                if hasattr(interp1, 'create_ants'):
                    interp1.create_ants()
                domain_interpolators.append(interp1)
                
                # Pass 2
                config2 = config.copy()
                config2['algorithm'] = algo2
                interp2 = create_interpolator(config2)
                
                if has_extra_blocks:
                    if algo2 == 'ant_colony':
                        interp2.allowed_grid_override = allowed_grid
                        interp2.domain_mapping = domain_mapping
                    elif algo2 == 'molecular_clock':
                        interp2.allowed_grid_override = allowed_grid
                        interp2.domain_mapping = domain_mapping
                    elif algo2 == 'gaussian_kernel':
                        interp2.allowed_grid_override = allowed_grid
                        interp2.domain_mapping = domain_mapping
                    elif algo2 == 'adaptive_octree':
                        interp2.allowed_grid_override = allowed_grid
                        interp2.domain_mapping = domain_mapping
                    elif algo2 == 'string_theory':
                        interp2.allowed_grid_override = allowed_grid
                        interp2.domain_mapping = domain_mapping
                    use_mapping_2 = True
                else:
                    interp2.allowed_grid_override = None
                    interp2.domain_mapping = None
                    use_mapping_2 = False
                
                interp2.initialize_blocks(multiblock._sample_blocks, tuple(block_info['dims']),
                                         all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping_2,
                                         sample_block_counts=getattr(multiblock, '_sample_block_counts', {}),
                                         sample_block_weight_sums=getattr(multiblock, '_sample_block_weight_sums', {}))
                if hasattr(interp2, 'create_ants'):
                    interp2.create_ants()
                domain_interpolators.append(interp2)
                
                multiblock._interpolators["Undomained"] = domain_interpolators
                multiblock._ant_colony = interp1
            else:
                interpolator = create_interpolator(config)
                
                # Attach domain-specific attributes BEFORE initialize_blocks
                # For molecular clock, don't filter by allowed_grid (let it create blocks anywhere)
                # For ant colony, use allowed_grid to restrict movement
                algo_type = config.get('algorithm', 'ant_colony')
                
                # Check if we have any non-sample blocks in this domain (sparse vs dense definition)
                # If the allowed grid is just the samples, we shouldn't enforce it as a constraint.
                # In undomained mode with blocks_file=None, allowed_grid is not defined here, but we are in the blocks_file!=None branch.
                # allowed_grid comes from domain_mapping.keys()
                
                has_extra_blocks = len(allowed_grid) > len(multiblock._sample_blocks) * 1.2
                print(f"  Undomained (1-pass) sparse check: allowed={len(allowed_grid)}, samples={len(multiblock._sample_blocks)}, has_extra={has_extra_blocks}")
                
                if has_extra_blocks:
                    if algo_type == 'ant_colony':
                        interpolator.allowed_grid_override = allowed_grid
                        interpolator.domain_mapping = domain_mapping
                    elif algo_type == 'molecular_clock':
                        interpolator.allowed_grid_override = allowed_grid
                        interpolator.domain_mapping = domain_mapping
                    elif algo_type == 'gaussian_kernel':
                        interpolator.allowed_grid_override = allowed_grid
                        interpolator.domain_mapping = domain_mapping
                    elif algo_type == 'adaptive_octree':
                        interpolator.allowed_grid_override = allowed_grid
                        interpolator.domain_mapping = domain_mapping
                    elif algo_type == 'string_theory':
                        interpolator.allowed_grid_override = allowed_grid
                        interpolator.domain_mapping = domain_mapping
                    
                    use_mapping = True
                else:
                    interpolator.allowed_grid_override = None
                    interpolator.domain_mapping = None
                    use_mapping = False
                
                # Use use_domain_mapping=True to respect allowed_grid_override (geometry)
                interpolator.initialize_blocks(multiblock._sample_blocks, tuple(block_info['dims']),
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping,
                                             sample_block_counts=getattr(multiblock, '_sample_block_counts', {}),
                                             sample_block_weight_sums=getattr(multiblock, '_sample_block_weight_sums', {}))
                
                if hasattr(interpolator, 'create_ants'):
                    interpolator.create_ants()
                
                multiblock._ant_colony = interpolator
        
        load_pbar.update(1)
        load_pbar.close()
        return multiblock
    else:
        print("Assigning points to blocks...")
        min_bounds = np.min(points, axis=0)
        max_bounds = np.max(points, axis=0)
        # Ensure at least 1 cell per axis (e.g. all samples at a constant Z)
        # and include the upper boundary so points exactly on max edge still fit.
        dims = np.ceil((max_bounds - min_bounds) / np.array(block_size)).astype(int) + 1
        dims = np.maximum(dims, 1)
        blocks = {}
        block_values = {}
        block_weights = {}
        block_domains = {}
        block_info = {
            'min_bounds': min_bounds,
            'dims': dims,
            'block_size': block_size,
            'source_blocks_file': None,
            'source_blocks_delimiter': None,
            'source_blocks_header_line': 1,
            'source_block_x_col': None,
            'source_block_y_col': None,
            'source_block_z_col': None,
            'source_block_filters': [],
            'expand_interpolation_exports_to_subblocks': bool((config or {}).get('expand_interpolation_exports_to_subblocks', True)),
        }
        # Vectorized block assignment
        points_array = np.array(points)
        values_array = np.array(values)
        domains_array = None
        weights_array = None
        if sample_domains is not None:
            domains_array = np.array(sample_domains)
        if sample_weights is not None:
            weights_array = np.array(sample_weights, dtype=float)
            if len(weights_array) != len(points_array):
                raise ValueError('Sample weights must have the same length as sample points.')
            if np.any(~np.isfinite(weights_array)) or np.any(weights_array <= 0.0):
                raise ValueError('Sample weights must be finite values greater than zero.')
        block_indices = ((points_array - min_bounds) // np.array(block_size)).astype(int)
        assigned_mask = np.ones(len(points_array), dtype=bool)
        
        for i in tqdm(range(len(points_array)), desc="Assigning points to blocks"):
            block_idx = tuple(block_indices[i])
            if block_idx not in blocks:
                blocks[block_idx] = []
                block_values[block_idx] = []
                if weights_array is not None:
                    block_weights[block_idx] = []
                if domains_array is not None:
                    block_domains[block_idx] = []
            blocks[block_idx].append(points_array[i])
            block_values[block_idx].append(values_array[i])
            if weights_array is not None:
                block_weights[block_idx].append(weights_array[i])
            if domains_array is not None:
                block_domains[block_idx].append(domains_array[i])
        next_block_id = 1

        # Detect domain interpolation mode (String Theory)
        is_st_domain_interpolation = bool(
            config
            and config.get('algorithm') in ('string_theory', 'net_connector')
            and str(config.get('string_theory_params', {}).get('interpolate_target', 'value')).strip().lower() == 'domain'
        )


        # Detect domain interpolation mode (Ant Colony)
        is_ant_domain_interpolation = bool(
            config
            and config.get('algorithm') == 'ant_colony'
            and str(config.get('ant_colony_interpolate_target', 'value')).strip().lower() == 'domain'
        )

        sample_block_domain_mapping = {}
        if is_st_domain_interpolation or is_ant_domain_interpolation:
            if domains_array is None:
                raise ValueError("Domain interpolation requires a 'Domain' column in the samples file.")

            bs = np.array(block_size, dtype=float)
            for idx, pts in blocks.items():
                doms_raw = block_domains.get(idx, [])
                # Clean domains: treat blanks as missing.
                doms = []
                aligned_pts = []
                for p, d in zip(pts, doms_raw):
                    ds = str(d).strip() if d is not None else ''
                    if ds == '' or ds.lower() == 'nan':
                        continue
                    doms.append(ds)
                    aligned_pts.append(p)

                if not doms:
                    continue

                # Majority vote
                counts = {}

                for ds in doms:
                    counts[ds] = counts.get(ds, 0) + 1
                max_count = max(counts.values())
                tied = [ds for ds, c in counts.items() if c == max_count]

                if len(tied) == 1:
                    winner = tied[0]
                else:
                    # Tie-breaker: lower average distance to block centroid
                    corner = min_bounds + np.array(idx, dtype=float) * bs
                    centroid = corner + bs / 2.0
                    best = None
                    for ds in tied:
                        dists = [float(np.linalg.norm(np.array(p) - centroid)) for p, d in zip(aligned_pts, doms) if d == ds]
                        avg = float(np.mean(dists)) if dists else float('inf')
                        cand = (avg, ds)
                        if best is None or cand < best:
                            best = cand
                    winner = best[1] if best is not None else sorted(tied)[0]

                sample_block_domain_mapping[idx] = winner

        if build_visual_blocks:
            block_data = []
            for idx in tqdm(blocks.keys(), desc="Creating blocks"):
                corner = min_bounds + np.array(idx) * np.array(block_size)
                cell = pv.Box(bounds=(
                    corner[0], corner[0] + block_size[0],
                    corner[1], corner[1] + block_size[1],
                    corner[2], corner[2] + block_size[2]
                ))
                avg_value = np.average(block_values[idx], weights=block_weights[idx]) if weights_array is not None else np.mean(block_values[idx])
                cell.cell_data['Value'] = np.full(cell.n_cells, avg_value)
                cell.cell_data['Is_Sample'] = np.full(cell.n_cells, True)
                cell.cell_data['Block_ID'] = np.full(cell.n_cells, next_block_id)
                next_block_id += 1
                if is_st_domain_interpolation or is_ant_domain_interpolation:
                    dom = sample_block_domain_mapping.get(idx, "")
                    cell.cell_data['Domain'] = np.full(cell.n_cells, dom)
                block_data.append(cell)
            multiblock = pv.MultiBlock(block_data)
        else:
            print("Skipping sample block geometry materialization (interpolation-only mode).")
            multiblock = LightweightBlocks()
        sample_blocks = {
            idx: np.average(vals, weights=block_weights[idx]) if weights_array is not None else np.mean(vals)
            for idx, vals in block_values.items()
        }
        sample_block_counts = {idx: len(vals) for idx, vals in block_values.items()}
        sample_block_weight_sums = {
            idx: float(np.sum(block_weights[idx])) if weights_array is not None else float(len(vals))
            for idx, vals in block_values.items()
        }
        multiblock._block_info = block_info
        multiblock._sample_blocks = sample_blocks
        multiblock._sample_assignment_data = {
            'points': original_points_array,
            'values': original_values_array,
            'domains': original_domains_array,
            'weights': original_weights_array,
            'block_indices': block_indices,
            'assigned_mask': assigned_mask,
            'sample_block_counts': dict(sample_block_counts),
            'sample_block_weight_sums': dict(sample_block_weight_sums),
        }
        multiblock._sample_block_counts = dict(sample_block_counts)
        multiblock._sample_block_weight_sums = dict(sample_block_weight_sums)

        if is_st_domain_interpolation or is_ant_domain_interpolation:
            # Preserve sample block domains for export.
            multiblock._block_info['domain_mapping'] = dict(sample_block_domain_mapping)
        
        # Create interpolator using factory (with config if available, otherwise build from params)
        if config is None:
            config = {
                'algorithm': 'ant_colony',
                'range_size': range_size,
                'max_pheromone': max_pheromone,
                'ants_per_sample': ants_per_sample,
                'verbose': verbose,
                'background_value': background_value,
                'background_distance': background_distance,
                'average_with_blocks': average_with_blocks,
                'avoid_visited_threshold_enabled': avoid_visited_threshold_enabled,
                'avoid_visited_threshold': avoid_visited_threshold
            }
        
        # Domain interpolation (String Theory) is a special case: the domain is the output, so we don't
        # generate a fake full-grid domain mapping.
        if is_st_domain_interpolation:
            algo2 = config.get('second_pass_algorithm', 'skip')
            if algo2 != 'skip':
                print("Warning: second pass is not supported for String Theory domain interpolation; ignoring second pass.")

            interpolator = create_interpolator(config)
            # Provide sample-block domains to the interpolator
            interpolator.domain_mapping = dict(sample_block_domain_mapping)
            interpolator.initialize_blocks(
                sample_blocks,
                tuple(block_info['dims']),
                min_bounds,
                block_size,
                use_domain_mapping=False,
                sample_domain_mapping=sample_block_domain_mapping,
                sample_block_counts=sample_block_counts,
                sample_block_weight_sums=sample_block_weight_sums,
            )
            multiblock._ant_colony = interpolator
            return multiblock


        # Ant Colony domain interpolation: initialize ants with sample_block domains and allow domains to expand.
        if is_ant_domain_interpolation:
            algo2 = config.get('second_pass_algorithm', 'skip')
            if algo2 != 'skip':
                print("Warning: second pass is not supported for Ant Colony domain interpolation; ignoring second pass.")

            interpolator = create_interpolator(config)
            interpolator.initialize_blocks(
                sample_blocks,
                tuple(block_info['dims']),
                min_bounds,
                block_size,
                use_domain_mapping=False,
                sample_domain_mapping=sample_block_domain_mapping,
                sample_block_counts=sample_block_counts,
                sample_block_weight_sums=sample_block_weight_sums,
            )
            if hasattr(interpolator, 'create_ants'):
                interpolator.create_ants()
            multiblock._ant_colony = interpolator
            return multiblock

        # Generate full grid domain mapping for "Undomained" mode
        # This ensures that algorithms like String Theory (which rely on domain mapping)
        # treat the entire bounding box as a valid domain.
        print("Generating full grid domain mapping for Undomained mode...")
        full_domain_mapping = {}
        full_allowed_grid = set()
        nx, ny, nz = tuple(block_info['dims'])
        
        # Optimization: Only generate if grid is not excessively large (e.g. < 10M blocks)
        # Otherwise, we might run out of memory.
        if nx * ny * nz < 10_000_000:
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        pos = (i, j, k)
                        full_domain_mapping[pos] = "Default"
                        full_allowed_grid.add(pos)
        else:
            print("Warning: Grid too large for full domain mapping. Falling back to unconstrained mode.")
            full_domain_mapping = None
            full_allowed_grid = None

        # Check for second pass
        algo2 = config.get('second_pass_algorithm', 'skip')
        
        if algo2 != 'skip':
            print(f"Configuring undomained two-pass interpolation (auto-blocks): {config.get('algorithm')} -> {algo2}")
            multiblock._interpolators = {}
            domain_interpolators = []
            
            # Pass 1
            interp1 = create_interpolator(config)
            if full_domain_mapping:
                interp1.allowed_grid_override = full_allowed_grid
                interp1.domain_mapping = full_domain_mapping
                use_mapping = True
            else:
                use_mapping = False
                
            interp1.initialize_blocks(sample_blocks, tuple(block_info['dims']),
                                       min_bounds, block_size, use_domain_mapping=use_mapping,
                                       sample_block_counts=sample_block_counts,
                                       sample_block_weight_sums=sample_block_weight_sums)
            if hasattr(interp1, 'create_ants'):
                interp1.create_ants()
            domain_interpolators.append(interp1)
            
            # Pass 2
            config2 = config.copy()
            config2['algorithm'] = algo2
            interp2 = create_interpolator(config2)
            if full_domain_mapping:
                interp2.allowed_grid_override = full_allowed_grid
                interp2.domain_mapping = full_domain_mapping
                use_mapping = True
            else:
                use_mapping = False
                
            interp2.initialize_blocks(sample_blocks, tuple(block_info['dims']),
                                       min_bounds, block_size, use_domain_mapping=use_mapping,
                                       sample_block_counts=sample_block_counts,
                                       sample_block_weight_sums=sample_block_weight_sums)
            if hasattr(interp2, 'create_ants'):
                interp2.create_ants()
            domain_interpolators.append(interp2)
            
            multiblock._interpolators["Undomained"] = domain_interpolators
            multiblock._ant_colony = interp1
        else:
            interpolator = create_interpolator(config)
            if full_domain_mapping:
                interpolator.allowed_grid_override = full_allowed_grid
                interpolator.domain_mapping = full_domain_mapping
                use_mapping = True
            else:
                use_mapping = False
                
            interpolator.initialize_blocks(sample_blocks, tuple(block_info['dims']),
                                           min_bounds, block_size, use_domain_mapping=use_mapping,
                                           sample_block_counts=sample_block_counts,
                                           sample_block_weight_sums=sample_block_weight_sums)
            
            if hasattr(interpolator, 'create_ants'):
                interpolator.create_ants()
            
            multiblock._ant_colony = interpolator
        return multiblock

def toggle_blocks(plotter):
    if not hasattr(plotter, '_blocks_actor') or plotter._blocks_actor is None:
        if hasattr(plotter, '_blocks_data'):
            print(f"Preparing blocks display: blocks={len(plotter._blocks_data)}")
        _ensure_blocks_actor(plotter, visible=True)
        if hasattr(plotter, '_block_lookup') and plotter._block_lookup is not None:
            print(f"Blocks lookup size: {len(plotter._block_lookup)}")
        if hasattr(plotter, '_visible_blocks') and plotter._visible_blocks is not None:
            print(f"Visible blocks: {len(plotter._visible_blocks)}")
        plotter.render()
        return
    if hasattr(plotter, '_blocks_actor'):
        is_visible = not plotter._blocks_actor.GetVisibility()
        plotter._blocks_actor.SetVisibility(is_visible)
        plotter.render()

def _get_block_raw_value(block):
    if 'Raw_Value' in block.cell_data:
        return float(block.cell_data['Raw_Value'][0])
    return float(block.cell_data['Value'][0])

def _classify_block_value(plotter, raw_value):
    if getattr(plotter, '_value_is_indexed', False) and hasattr(plotter, '_lfc_bins') and isinstance(plotter._colormap, ListedColormap):
        bins_np = np.array(plotter._lfc_bins)
        n_colors = len(plotter._colormap.colors)
        threshold_style = (n_colors == len(bins_np) + 1)
        if threshold_style:
            idx_val = np.searchsorted(bins_np, [raw_value], side='right')[0]
        else:
            idx_val = np.digitize([raw_value], bins_np, right=False)[0] - 1
        return max(0, min(idx_val, n_colors - 1))
    return float(raw_value)

def _set_block_display_value(plotter, block, raw_value):
    raw_value = float(raw_value)
    display_value = _classify_block_value(plotter, raw_value)
    if 'Raw_Value' in block.cell_data:
        block.cell_data['Raw_Value'][:] = raw_value
    elif getattr(plotter, '_value_is_indexed', False):
        block.cell_data['Raw_Value'] = np.full(block.n_cells, raw_value)
    if 'Value' in block.cell_data:
        block.cell_data['Value'][:] = display_value
    else:
        block.cell_data['Value'] = np.full(block.n_cells, display_value)
    block.Modified()

def _get_block_grid_position(block, min_bounds, block_size):
    corner = np.array(block.bounds[::2], dtype=float)
    return tuple(np.floor((corner - np.array(min_bounds)) / np.array(block_size) + 1e-6).astype(int))

def _build_block_lookup(blocks, min_bounds, block_size):
    lookup = {}
    for block in blocks:
        pos = _get_block_grid_position(block, min_bounds, block_size)
        lookup[pos] = block
    return lookup

def _should_display_block(plotter, pos, block):
    if pos in plotter._blocks_data._sample_blocks:
        return True
    return _get_block_raw_value(block) >= plotter._value_filter

def _create_visible_blocks(plotter):
    pv = _require_pyvista()
    visible_blocks = pv.MultiBlock()
    visible_positions = set()
    for pos, block in plotter._block_lookup.items():
        if _should_display_block(plotter, pos, block):
            visible_blocks.append(block)
            visible_positions.add(pos)
    if len(visible_blocks) == 0 and len(plotter._block_lookup) > 0:
        for block in plotter._block_lookup.values():
            visible_blocks.append(block)
        visible_positions = set(plotter._block_lookup.keys())
    return visible_blocks, visible_positions

def _get_blocks_mesh_kwargs(plotter):
    if getattr(plotter, '_value_is_indexed', False):
        annotations = {i: plotter._lfc_tick_labels[i] for i in range(len(plotter._lfc_tick_labels))}
        return dict(
            style='surface',
            scalars='Value',
            opacity=0.5,
            show_edges=True,
            cmap=plotter._colormap,
            categories=True,
            annotations=annotations,
            clim=[-0.5, len(plotter._colormap.colors) - 0.5],
            show_scalar_bar=False,
        )
    return dict(
        style='surface',
        scalars='Value',
        opacity=0.5,
        show_edges=True,
        cmap=plotter._colormap,
        clim=list(getattr(plotter, '_scalar_range', (0.0, 1.0))),
        show_scalar_bar=False,
    )

def _prepare_blocks_for_display(plotter):
    if getattr(plotter, '_blocks_display_prepared', False):
        return
    if getattr(plotter, '_value_is_indexed', False) and hasattr(plotter, '_lfc_bins') and isinstance(plotter._colormap, ListedColormap):
        bins_np = np.array(plotter._lfc_bins)
        n_colors = len(plotter._colormap.colors)
        threshold_style = (n_colors == len(bins_np) + 1)
        for block in plotter._blocks_data:
            if 'Raw_Value' in block.cell_data:
                raw_val = block.cell_data['Raw_Value'][0]
            else:
                raw_val = block.cell_data['Value'][0]
                block.cell_data['Raw_Value'] = np.full(block.n_cells, raw_val)
            if threshold_style:
                idx_val = np.searchsorted(bins_np, [raw_val], side='right')[0]
            else:
                idx_val = np.digitize([raw_val], bins_np, right=False)[0] - 1
            idx_val = max(0, min(idx_val, n_colors - 1))
            block.cell_data['Value'][:] = idx_val
    plotter._blocks_display_prepared = True

def _ensure_blocks_actor(plotter, visible=True):
    if not hasattr(plotter, '_block_lookup') or plotter._block_lookup is None:
        plotter._block_lookup = _build_block_lookup(
            plotter._blocks_data,
            plotter._blocks_data._block_info['min_bounds'],
            plotter._blocks_data._block_info['block_size'],
        )
    if not hasattr(plotter, '_visible_blocks') or plotter._visible_blocks is None or not hasattr(plotter, '_visible_positions') or plotter._visible_positions is None:
        plotter._visible_blocks, plotter._visible_positions = _create_visible_blocks(plotter)
    _prepare_blocks_for_display(plotter)
    if not hasattr(plotter, '_blocks_actor') or plotter._blocks_actor is None:
        plotter._blocks_actor = plotter.add_mesh(plotter._visible_blocks, **_get_blocks_mesh_kwargs(plotter))
    plotter._blocks_actor.SetVisibility(visible)

def _rebuild_blocks_actor(plotter):
    _prepare_blocks_for_display(plotter)
    visible = True
    if hasattr(plotter, '_blocks_actor') and plotter._blocks_actor is not None:
        visible = bool(plotter._blocks_actor.GetVisibility())
        plotter.remove_actor(plotter._blocks_actor)
    plotter._visible_blocks, plotter._visible_positions = _create_visible_blocks(plotter)
    plotter._blocks_actor = plotter.add_mesh(plotter._visible_blocks, **_get_blocks_mesh_kwargs(plotter))
    plotter._blocks_actor.SetVisibility(visible)

def _mark_block_datasets_modified(plotter):
    if hasattr(plotter, '_blocks_data'):
        plotter._blocks_data.Modified()
    if hasattr(plotter, '_visible_blocks'):
        plotter._visible_blocks.Modified()
    if hasattr(plotter, '_blocks_actor') and hasattr(plotter._blocks_actor, 'mapper'):
        plotter._blocks_actor.mapper.Modified()

def update_interpolation(plotter):
    pv = _require_pyvista()
    print("Checking ant colony data...")
    if not hasattr(plotter, '_blocks_data'):
        print("No blocks data found")
        return

    blocks = plotter._blocks_data
    if not hasattr(blocks, '_ant_colony'):
        print("No ant colony found in blocks")
        return

    print("Found ant colony, updating interpolation...")
    interpolator = blocks._ant_colony
    dims = tuple(blocks._block_info['dims'])
    min_bounds = blocks._block_info['min_bounds']
    block_size = blocks._block_info['block_size']

    if not hasattr(plotter, '_block_lookup') or plotter._block_lookup is None:
        plotter._block_lookup = _build_block_lookup(blocks, min_bounds, block_size)
    if not hasattr(plotter, '_visible_blocks') or plotter._visible_blocks is None or not hasattr(plotter, '_visible_positions') or plotter._visible_positions is None:
        plotter._visible_blocks, plotter._visible_positions = _create_visible_blocks(plotter)

    print("Moving ants...")
    changes_made = False

    interpolators = []
    if hasattr(blocks, '_interpolators') and blocks._interpolators:
        for entry in blocks._interpolators.values():
            if isinstance(entry, list):
                if entry:
                    interpolators.append(entry[0])
            else:
                interpolators.append(entry)
    else:
        interpolators = [interpolator]

    for interp in interpolators:
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

    if changes_made:
        colony_blocks = {}
        pos_to_interpolator = {}

        for interp in interpolators:
            vals = interp.get_interpolated_values()
            colony_blocks.update(vals)
            for pos in vals:
                pos_to_interpolator[pos] = interp

        actor_rebuild_needed = False
        actor_dataset_changed = False
        new_block_count = 0
        new_visible_added = False
        scalar_min, scalar_max = getattr(plotter, '_scalar_range', (0.0, 0.0))

        for pos, value in colony_blocks.items():
            try:
                if pos in blocks._sample_blocks:
                    continue

                if pos not in plotter._block_lookup:
                    corner = min_bounds + np.array(pos) * np.array(block_size)
                    half_size = np.array(block_size) / 2
                    center = corner + half_size

                    new_block = pv.Box(bounds=(
                        center[0] - half_size[0]/2, center[0] + half_size[0]/2,
                        center[1] - half_size[1]/2, center[1] + half_size[1]/2,
                        center[2] - half_size[2]/2, center[2] + half_size[2]/2
                    ))
                    new_block.cell_data['Is_Sample'] = np.full(new_block.n_cells, False)
                    _set_block_display_value(plotter, new_block, value)

                    target_interp = pos_to_interpolator.get(pos, interpolator)
                    if not hasattr(target_interp, 'next_block_id'):
                        existing_ids = [
                            int(block.cell_data['Block_ID'][0])
                            for block in blocks
                            if 'Block_ID' in block.cell_data
                        ]
                        target_interp.next_block_id = (max(existing_ids) + 1) if existing_ids else 1
                    new_block.cell_data['Block_ID'] = np.full(new_block.n_cells, target_interp.next_block_id)
                    target_interp.next_block_id += 1

                    if hasattr(target_interp, 'domain_mapping'):
                        domain = target_interp.domain_mapping.get(pos, "Undomained")
                        new_block.cell_data['Domain'] = np.full(new_block.n_cells, domain)

                    blocks.append(new_block)
                    plotter._block_lookup[pos] = new_block
                    blocks.Modified()
                    new_block_count += 1

                    if _should_display_block(plotter, pos, new_block):
                        plotter._visible_blocks.append(new_block)
                        plotter._visible_positions.add(pos)
                        actor_dataset_changed = True
                        new_visible_added = True
                else:
                    block = plotter._block_lookup[pos]
                    base_val = _get_block_raw_value(block)
                    if abs(base_val - value) > 0.0001:
                        was_visible = pos in plotter._visible_positions
                        _set_block_display_value(plotter, block, value)
                        is_visible = _should_display_block(plotter, pos, block)
                        if was_visible and not is_visible:
                            actor_rebuild_needed = True
                        elif (not was_visible) and is_visible:
                            plotter._visible_blocks.append(block)
                            plotter._visible_positions.add(pos)
                            actor_dataset_changed = True

                scalar_min = min(scalar_min, float(value))
                scalar_max = max(scalar_max, float(value))
            except Exception as e:
                print(f"Error processing position {pos}: {str(e)}")
                continue

        if not getattr(plotter, '_value_is_indexed', False):
            current_min, current_max = getattr(plotter, '_scalar_range', (scalar_min, scalar_max))
            plotter._scalar_range = (min(current_min, scalar_min), max(current_max, scalar_max))

        if actor_rebuild_needed:
            _rebuild_blocks_actor(plotter)
        else:
            _mark_block_datasets_modified(plotter)
            if not getattr(plotter, '_value_is_indexed', False) and hasattr(plotter, '_blocks_actor') and hasattr(plotter._blocks_actor, 'mapper'):
                plotter._blocks_actor.mapper.SetScalarRange(*plotter._scalar_range)

        if new_block_count > 0 and not new_visible_added:
            print(
                "New blocks were created, but none passed the visibility filter. "
                f"Current value_filter={plotter._value_filter}."
            )

        plotter.render()
        print("Visualization updated")
    else:
        print("No blocks to update")

def resolve_block_evaluated_samples_export_path(enabled, configured_path=None, interpolation_file=None, samples_file=None):
    if not enabled:
        return None

    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if interpolation_file and str(interpolation_file).strip():
        reference_path = str(interpolation_file).strip()
    elif samples_file and str(samples_file).strip():
        reference_path = str(samples_file).strip()
    else:
        reference_path = 'interpolation.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'interpolation'
    return os.path.join(output_dir, f"{base_name}_block_evaluated_samples.csv")


def resolve_interpolation_csv_export_path(interpolation_file):
    path = str(interpolation_file or '').strip()
    if not path:
        return 'interpolation.csv'
    root, ext = os.path.splitext(path)
    if ext.lower() == '.bmf':
        return f"{root}.csv" if root else 'interpolation.csv'
    return path


def _sanitize_filename_fragment(value, fallback='Domain'):
    text = str(value or '').strip()
    if not text or text == '(None)':
        text = fallback
    return text.translate(INVALID_FILENAME_CHARS)


def resolve_sample_blocks_export_path(configured_path=None, samples_file=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if samples_file and str(samples_file).strip():
        reference_path = str(samples_file).strip()
    else:
        reference_path = 'samples.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'samples'
    return os.path.join(output_dir, f"{base_name}_SampleBlocks.csv")


def resolve_domain_samples_export_path(configured_path=None, samples_file=None, domain_col=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if samples_file and str(samples_file).strip():
        reference_path = str(samples_file).strip()
    else:
        reference_path = 'samples.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'samples'
    domain_suffix = _sanitize_filename_fragment(domain_col, fallback='Domain')
    return os.path.join(output_dir, f"{base_name}+{domain_suffix}.csv")


def resolve_block_value_transfer_export_path(configured_path=None, samples_file=None, block_value_cols=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if samples_file and str(samples_file).strip():
        reference_path = str(samples_file).strip()
    else:
        reference_path = 'samples.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'samples'
    normalized_columns = []
    seen = set()
    for value in block_value_cols or []:
        column_name = str(value or '').strip()
        if not column_name or column_name in seen:
            continue
        normalized_columns.append(column_name)
        seen.add(column_name)

    if not normalized_columns:
        suffix = 'block_values'
    else:
        preview = [_sanitize_filename_fragment(column_name, fallback='value') for column_name in normalized_columns[:3]]
        if len(normalized_columns) > 3:
            preview.append(f"plus{len(normalized_columns) - 3}")
        suffix = '+'.join(preview)

    return os.path.join(output_dir, f"{base_name}+{suffix}.csv")

def resolve_block_model_transfer_export_path(configured_path=None, target_blocks_file=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    reference_path = str(target_blocks_file or '').strip() or 'target_blocks.csv'
    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'target_blocks'
    return os.path.join(output_dir, f'{base_name}+SourceBlocks.csv')


def resolve_block_model_table_attribute_export_path(configured_path=None, block_model_file=None, table_file=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    reference_path = str(block_model_file or '').strip() or 'block_model.csv'
    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'block_model'
    table_base_name = os.path.splitext(os.path.basename(str(table_file or '').strip()))[0] or 'table'
    table_suffix = _sanitize_filename_fragment(table_base_name, fallback='table')
    return os.path.join(output_dir, f'{base_name}+{table_suffix}_attributes.csv')



def resolve_block_domain_metrics_export_path(configured_path=None, blocks_file=None, domain_col=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if blocks_file and str(blocks_file).strip():
        reference_path = str(blocks_file).strip()
    else:
        reference_path = 'blocks.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'blocks'
    domain_suffix = _sanitize_filename_fragment(domain_col, fallback='Domain')
    return os.path.join(output_dir, f"{base_name}+{domain_suffix}_sample_metrics.csv")


def resolve_block_domain_metrics_summary_export_path(metrics_file=None, blocks_file=None, domain_col=None):
    reference_path = str(metrics_file or '').strip()
    if not reference_path:
        reference_path = resolve_block_domain_metrics_export_path(
            None,
            blocks_file=blocks_file,
            domain_col=domain_col,
        )

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'blocks'
    return os.path.join(output_dir, f"{base_name}_summary.csv")


def resolve_domain_interpolation_confidence_export_path(configured_path=None, blocks_file=None, domain_col=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if blocks_file and str(blocks_file).strip():
        reference_path = str(blocks_file).strip()
    else:
        reference_path = 'blocks.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'blocks'
    domain_suffix = _sanitize_filename_fragment(domain_col, fallback='Domain')
    return os.path.join(output_dir, f"{base_name}+{domain_suffix}_interpolation_confidence.csv")


def resolve_block_volume_weighted_average_export_path(configured_path=None, blocks_file=None, value_col=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if blocks_file and str(blocks_file).strip():
        reference_path = str(blocks_file).strip()
    else:
        reference_path = 'blocks.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'blocks'
    value_suffix = _sanitize_filename_fragment(value_col, fallback='Value')
    return os.path.join(output_dir, f"{base_name}+{value_suffix}_volume_weighted.csv")


def resolve_equation_finder_export_path(configured_path=None, samples_file=None, value_col=None, domain_col=None):
    configured_path = str(configured_path or '').strip()
    if configured_path:
        return configured_path

    if samples_file and str(samples_file).strip():
        reference_path = str(samples_file).strip()
    else:
        reference_path = 'samples.csv'

    output_dir = os.path.dirname(reference_path) or '.'
    base_name = os.path.splitext(os.path.basename(reference_path))[0] or 'samples'
    value_suffix = _sanitize_filename_fragment(value_col, fallback='Value')
    domain_suffix = _sanitize_filename_fragment(domain_col, fallback='Domain')
    return os.path.join(output_dir, f"{base_name}+{domain_suffix}+{value_suffix}_equations.csv")


def load_samples_preview_dataframe(samples_file, samples_delimiter=None, samples_header_line=1, max_rows=500):
    delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
    preview_rows = max(int(max_rows or 0), 1)

    if samples_header_line and samples_header_line != 1:
        effective_header_line = resolve_effective_csv_header_line(samples_file, samples_header_line)
        headers = parse_header_line(samples_file, delimiter, effective_header_line)
        final_names = build_unique_column_names(headers)
        read_kwargs = prepare_csv_read_kwargs(
            samples_file,
            delimiter=delimiter,
            header=None,
            names=final_names,
            skiprows=max(int(effective_header_line) - 1, 0),
            comment='#',
            nrows=preview_rows + 1,
        )
        df = pd.read_csv(samples_file, **read_kwargs)
        if df.shape[0] and all(str(df.iloc[0, i]).strip() == final_names[i] for i in range(min(len(final_names), df.shape[1]))):
            df = df.iloc[1:].reset_index(drop=True)
    else:
        read_kwargs = prepare_csv_read_kwargs(
            samples_file,
            delimiter=delimiter,
            comment='#',
            nrows=preview_rows,
        )
        df = pd.read_csv(samples_file, **read_kwargs)

    df = strip_leading_non_data_rows(df)

    empty_cols = _list_all_empty_columns(df)
    if empty_cols:
        df = df.drop(columns=empty_cols)

    df._detected_delimiter = delimiter
    return df, delimiter


def infer_numeric_sample_columns(df, delimiter=None, minimum_success_ratio=0.8):
    numeric_columns = []
    for column_name in df.columns:
        series = df[column_name]
        if pd.api.types.is_bool_dtype(series):
            continue
        if pd.api.types.is_numeric_dtype(series):
            numeric_columns.append(str(column_name))
            continue

        non_blank = series.dropna()
        if len(non_blank) == 0:
            continue
        non_blank = non_blank.astype(str).str.strip()
        non_blank = non_blank[non_blank != '']
        if len(non_blank) == 0:
            continue

        converted = pd.to_numeric(non_blank, errors='coerce')
        success_ratio = float(converted.notna().mean()) if len(converted) else 0.0

        if success_ratio < float(minimum_success_ratio):
            normalized = non_blank.str.replace('\u00a0', '', regex=False).str.replace(' ', '', regex=False)
            if delimiter and delimiter != ',':
                normalized = normalized.str.replace(',', '.', regex=False)
            converted = pd.to_numeric(normalized, errors='coerce')
            success_ratio = float(converted.notna().mean()) if len(converted) else 0.0

        if success_ratio >= float(minimum_success_ratio):
            numeric_columns.append(str(column_name))
    return numeric_columns


def coerce_numeric_series(series, delimiter=None):
    numeric_series = pd.to_numeric(series, errors='coerce')
    missing_mask = numeric_series.isna()
    if not missing_mask.any():
        return numeric_series

    text_series = series.astype(str)
    normalized = text_series.str.replace('\u00a0', '', regex=False).str.replace(' ', '', regex=False)
    if delimiter and delimiter != ',':
        normalized = normalized.str.replace(',', '.', regex=False)
    alternate_numeric = pd.to_numeric(normalized, errors='coerce')
    numeric_series = numeric_series.where(~missing_mask, alternate_numeric)
    return numeric_series


def _safe_metric_r2(y_true, y_pred):
    if len(y_true) < 2:
        return np.nan
    from sklearn.metrics import r2_score

    try:
        return float(r2_score(y_true, y_pred))
    except Exception:
        return np.nan


def export_domain_symbolic_regression_equations(samples_file, output_file=None,
                                                samples_delimiter=None, samples_header_line=1,
                                                sample_value_col=None, sample_domain_col=None,
                                                predictor_cols=None, sample_filters=None,
                                                min_samples_per_domain=25,
                                                validation_fraction=0.2,
                                                max_iterations=100,
                                                timeout_seconds=60,
                                                progress_callback=None):
    if not samples_file or not os.path.isfile(samples_file):
        raise ValueError('Please select a valid samples file.')

    target_column = str(sample_value_col or '').strip()
    if not target_column or target_column == '(None)':
        raise ValueError('Please select a value column in "Samples Columns" first.')

    domain_column = str(sample_domain_col or '').strip()
    if not domain_column or domain_column == '(None)':
        raise ValueError('Please select a domain column in "Samples Columns" first.')

    predictor_columns = []
    for column_name in predictor_cols or []:
        text = str(column_name or '').strip()
        if text and text not in predictor_columns:
            predictor_columns.append(text)
    predictor_columns = [col for col in predictor_columns if col != target_column]
    if not predictor_columns:
        raise ValueError('Please select at least one predictor column.')

    try:
        pysr_module = importlib.import_module('pysr')
        PySRRegressor = pysr_module.PySRRegressor
    except Exception as exc:
        raise RuntimeError(
            'PySR is not available. Install the pysr package and ensure Julia bootstrap works on this machine before running the equation finder.'
        ) from exc

    output_file = resolve_equation_finder_export_path(
        output_file,
        samples_file=samples_file,
        value_col=target_column,
        domain_col=domain_column,
    )
    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(output_file))[0] or 'equations'
    details_dir = os.path.join(output_dir, f'{base_name}_details')
    os.makedirs(details_dir, exist_ok=True)

    selected_columns = [target_column, domain_column, *predictor_columns]
    selected_columns.extend(collect_filter_fields(sample_filters))
    selected_columns = list(dict.fromkeys(selected_columns))

    _emit_progress(progress_callback, 0, 100, 'Reading sample file...')
    delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
    df_samples, _ = read_selected_columns_with_header(
        samples_file,
        delimiter,
        samples_header_line or 1,
        selected_columns,
        progress_label='Reading sample file',
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 20, 'Reading sample file...'),
    )
    if sample_filters:
        df_samples, _ = apply_sample_filters(df_samples, sample_filters=sample_filters)

    missing_columns = [
        column_name for column_name in [target_column, domain_column, *predictor_columns]
        if column_name not in df_samples.columns
    ]
    if missing_columns:
        raise ValueError(f'Selected columns not found in samples file: {missing_columns}')

    _emit_progress(progress_callback, 22, 100, 'Validating predictors...')
    predictor_frame = pd.DataFrame(index=df_samples.index)
    for column_name in predictor_columns:
        predictor_frame[column_name] = coerce_numeric_series(df_samples[column_name], delimiter=delimiter)
    target_series = coerce_numeric_series(df_samples[target_column], delimiter=delimiter)
    domain_series = df_samples[domain_column].fillna('').astype(str).str.strip()
    domain_series = domain_series.replace('nan', '')

    invalid_predictors = []
    for column_name in predictor_columns:
        non_blank = df_samples[column_name].dropna().astype(str).str.strip()
        non_blank = non_blank[non_blank != '']
        if len(non_blank) == 0:
            invalid_predictors.append(column_name)
            continue
        converted = coerce_numeric_series(non_blank, delimiter=delimiter)
        if not converted.notna().any():
            invalid_predictors.append(column_name)

    skipped_predictors = list(invalid_predictors)
    predictor_columns = [column_name for column_name in predictor_columns if column_name not in skipped_predictors]
    if not predictor_columns:
        raise ValueError('None of the selected predictor columns contain usable numeric values after coercion.')
    if skipped_predictors:
        predictor_frame = predictor_frame[predictor_columns].copy()
        _emit_progress(
            progress_callback,
            23,
            100,
            f'Skipping {len(skipped_predictors)} predictor columns without usable numeric values...',
        )

    valid_mask = domain_series.ne('') & target_series.notna() & predictor_frame.notna().all(axis=1)
    filtered_row_count = int(valid_mask.sum())
    working_df = pd.DataFrame({
        'Domain': domain_series.loc[valid_mask].to_numpy(dtype=object, copy=False),
        'Target': target_series.loc[valid_mask].to_numpy(dtype=float, copy=False),
    })
    for column_name in predictor_columns:
        working_df[column_name] = predictor_frame.loc[valid_mask, column_name].to_numpy(dtype=float, copy=False)

    if len(working_df) == 0:
        raise ValueError('No valid rows remain after removing blank domains and non-numeric target/predictor values.')

    domain_names = sorted(str(name) for name in working_df['Domain'].dropna().unique())
    if not domain_names:
        raise ValueError('No valid domain values remain after filtering.')

    from sklearn.metrics import mean_absolute_error, mean_squared_error
    from sklearn.model_selection import train_test_split

    min_samples = max(int(min_samples_per_domain or 0), 2)
    test_fraction = float(validation_fraction or 0.0)
    if not np.isfinite(test_fraction):
        test_fraction = 0.2
    test_fraction = min(max(test_fraction, 0.0), 0.5)
    niterations = max(int(max_iterations or 0), 1)
    timeout_value = None
    if timeout_seconds is not None:
        timeout_value = max(float(timeout_seconds or 0.0), 0.0)
        if timeout_value <= 0:
            timeout_value = None

    summary_rows = []
    processed_domains = 0
    skipped_domains = 0
    for domain_index, domain_name in enumerate(domain_names, start=1):
        domain_df = working_df.loc[working_df['Domain'] == domain_name].copy()
        domain_row_count = int(len(domain_df))
        domain_slug = _sanitize_filename_fragment(domain_name, fallback=f'domain_{domain_index}')
        detail_path = os.path.join(details_dir, f'{domain_index:03d}_{domain_slug}_hall_of_fame.csv')

        base_progress = 25 + int(round(((domain_index - 1) / max(len(domain_names), 1)) * 70))
        _emit_progress(progress_callback, base_progress, 100, f'Fitting domain {domain_index}/{len(domain_names)}: {domain_name}')
        print(f'Equation finder: domain {domain_index}/{len(domain_names)} -> {domain_name} ({domain_row_count:,} valid rows)')

        if domain_row_count < min_samples:
            skipped_domains += 1
            print(f'Equation finder: skipped domain {domain_name} because it has {domain_row_count:,} valid rows; requires at least {min_samples:,}.')
            summary_rows.append({
                'Domain': domain_name,
                'Samples': domain_row_count,
                'Predictor_Count': len(predictor_columns),
                'Predictors': ', '.join(predictor_columns),
                'Skipped_Predictors': ', '.join(skipped_predictors),
                'Equation': '',
                'Complexity': np.nan,
                'Loss': np.nan,
                'Train_RMSE': np.nan,
                'Train_MAE': np.nan,
                'Train_R2': np.nan,
                'Validation_RMSE': np.nan,
                'Validation_MAE': np.nan,
                'Validation_R2': np.nan,
                'Status': f'Skipped: requires at least {min_samples} rows',
                'Hall_Of_Fame_File': detail_path,
            })
            continue

        X = domain_df[predictor_columns].to_numpy(dtype=float, copy=False)
        y = domain_df['Target'].to_numpy(dtype=float, copy=False)

        if test_fraction > 0.0 and domain_row_count >= max(min_samples, 5):
            X_train, X_test, y_train, y_test = train_test_split(
                X,
                y,
                test_size=test_fraction,
                random_state=42,
            )
        else:
            X_train, X_test, y_train, y_test = X, X, y, y

        model = PySRRegressor(
            model_selection='best',
            niterations=niterations,
            timeout_in_seconds=timeout_value,
            binary_operators=['+', '-', '*', '/'],
            unary_operators=['square'],
            maxsize=30,
            verbosity=0,
            progress=False,
            input_stream='devnull',
            temp_equation_file=True,
            warm_start=False,
        )

        try:
            model.fit(X_train, y_train)
        except Exception as exc:
            raise RuntimeError(
                f'Equation search failed for domain "{domain_name}". PySR requires a working Julia runtime/bootstrap on this machine. Original error: {exc}'
            ) from exc
        print(f'Equation finder: finished domain {domain_name}.')

        train_pred = np.asarray(model.predict(X_train), dtype=float)
        test_pred = np.asarray(model.predict(X_test), dtype=float)
        equations_df = getattr(model, 'equations_', None)
        best_equation = str(model.sympy()) if hasattr(model, 'sympy') else ''
        complexity_value = np.nan
        loss_value = np.nan

        if isinstance(equations_df, pd.DataFrame) and not equations_df.empty:
            equations_export = equations_df.copy()
            selected_mask = pd.Series(False, index=equations_export.index)
            if 'sympy_format' in equations_export.columns and best_equation:
                selected_mask = equations_export['sympy_format'].astype(str) == best_equation
            if not selected_mask.any() and 'pick' in equations_export.columns:
                pick_text = equations_export['pick'].astype(str)
                selected_mask = pick_text.str.contains('>', regex=False) | pick_text.isin(['True', 'true', '1'])
            equations_export.insert(0, 'Domain', domain_name)
            equations_export.insert(1, 'Selected_By_App', selected_mask.to_numpy(dtype=bool, copy=False))
            equations_export.to_csv(detail_path, index=False)

            if selected_mask.any():
                selected_row = equations_export.loc[selected_mask].iloc[0]
                complexity_value = selected_row.get('complexity', np.nan)
                loss_value = selected_row.get('loss', np.nan)
        else:
            pd.DataFrame([
                {
                    'Domain': domain_name,
                    'Selected_By_App': True,
                    'equation': best_equation,
                }
            ]).to_csv(detail_path, index=False)

        processed_domains += 1
        summary_rows.append({
            'Domain': domain_name,
            'Samples': domain_row_count,
            'Predictor_Count': len(predictor_columns),
            'Predictors': ', '.join(predictor_columns),
            'Skipped_Predictors': ', '.join(skipped_predictors),
            'Equation': best_equation,
            'Complexity': complexity_value,
            'Loss': loss_value,
            'Train_RMSE': float(np.sqrt(mean_squared_error(y_train, train_pred))),
            'Train_MAE': float(mean_absolute_error(y_train, train_pred)),
            'Train_R2': _safe_metric_r2(y_train, train_pred),
            'Validation_RMSE': float(np.sqrt(mean_squared_error(y_test, test_pred))),
            'Validation_MAE': float(mean_absolute_error(y_test, test_pred)),
            'Validation_R2': _safe_metric_r2(y_test, test_pred),
            'Status': 'OK',
            'Hall_Of_Fame_File': detail_path,
        })

    summary_df = pd.DataFrame(summary_rows)
    _emit_progress(progress_callback, 100, 100, 'Equation search complete.')
    return {
        'output_file': output_file,
        'details_directory': details_dir,
        'summary_dataframe': summary_df,
        'target_column': target_column,
        'domain_column': domain_column,
        'predictor_columns': list(predictor_columns),
        'skipped_predictor_columns': list(skipped_predictors),
        'input_rows': int(len(df_samples)),
        'valid_rows': filtered_row_count,
        'domain_count': len(domain_names),
        'processed_domains': processed_domains,
        'skipped_domains': skipped_domains,
    }


def load_full_samples_dataframe(samples_file, samples_delimiter=None, samples_header_line=1,
                                sample_filters=None, progress_label=None, progress_callback=None,
                                preserve_empty_columns=False):
    delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
    if samples_header_line and samples_header_line != 1:
        df, _ = read_csv_with_selected_header(
            samples_file,
            delimiter,
            samples_header_line,
            expected_min_cols=3,
            progress_label=progress_label,
            progress_callback=progress_callback,
            preserve_empty_columns=preserve_empty_columns,
        )
        if sample_filters:
            df, _ = apply_sample_filters(df, sample_filters=sample_filters)
        return df, delimiter

    df = read_autodetect_csv(
        samples_file,
        forced_delimiter=delimiter,
        progress_label=progress_label,
        progress_callback=progress_callback,
        preserve_empty_columns=preserve_empty_columns,
    )
    if sample_filters:
        df, _ = apply_sample_filters(df, sample_filters=sample_filters)
    return df, delimiter


def load_full_blocks_dataframe(blocks_file, blocks_delimiter=None, blocks_header_line=1,
                               block_filters=None, progress_label=None, progress_callback=None):
    if is_bmf_file(blocks_file):
        df = read_autodetect_csv(
            blocks_file,
            progress_label=progress_label,
            progress_callback=progress_callback,
        )
        if block_filters:
            print(
                f"Finished reading {len(df):,} rows from {os.path.basename(blocks_file)}. "
                "Preparing block filters..."
            )
            df, _ = apply_dataframe_filters(df, filters=block_filters, filter_subject='block', source_label='blocks file')
        return df, ','

    delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
    if blocks_header_line and blocks_header_line != 1:
        df, _ = read_csv_with_selected_header(
            blocks_file,
            delimiter,
            blocks_header_line,
            expected_min_cols=3,
            progress_label=progress_label,
            progress_callback=progress_callback,
        )
        if block_filters:
            print(
                f"Finished reading {len(df):,} rows from {os.path.basename(blocks_file)}. "
                "Preparing block filters..."
            )
            df, _ = apply_dataframe_filters(df, filters=block_filters, filter_subject='block', source_label='blocks file')
        return df, delimiter

    df = read_autodetect_csv(
        blocks_file,
        forced_delimiter=delimiter,
        progress_label=progress_label,
        progress_callback=progress_callback,
    )
    if block_filters:
        print(
            f"Finished reading {len(df):,} rows from {os.path.basename(blocks_file)}. "
            "Preparing block filters..."
        )
        df, _ = apply_dataframe_filters(df, filters=block_filters, filter_subject='block', source_label='blocks file')
    return df, delimiter

class FilterDataSource:
    def __init__(self, source, delimiter=None, header_line=1):
        self._delimiter = delimiter
        self._header_line = max(int(header_line or 1), 1)
        self._series_cache = {}
        self._categorical_value_cache = {}
        self._numeric_range_cache = {}

        if isinstance(source, pd.DataFrame):
            self._dataframe = source
            self._path = None
            display_columns = build_unique_column_names([str(column) for column in source.columns])
            self._columns = display_columns
            self._column_positions = {column: index for index, column in enumerate(display_columns)}
        elif isinstance(source, (str, os.PathLike)):
            self._dataframe = None
            self._path = os.fspath(source)
            if not os.path.isfile(self._path):
                raise ValueError(f'Filter data source file not found: {self._path}')
            self._delimiter = delimiter or detect_csv_delimiter(self._path)
            self._header_line = resolve_effective_csv_header_line(self._path, self._header_line)
            headers = parse_header_line(self._path, self._delimiter, self._header_line)
            self._columns = build_unique_column_names(headers)
            self._column_positions = {column: index for index, column in enumerate(self._columns)}
        else:
            raise TypeError('FilterDataSource source must be a pandas DataFrame or a CSV file path.')

    @property
    def columns(self):
        return list(self._columns)

    def has_field(self, field):
        return str(field or '').strip() in self._column_positions

    def get_series(self, field):
        field_name = str(field or '').strip()
        if not field_name:
            raise ValueError('Filter field name cannot be empty.')
        if field_name not in self._column_positions:
            raise ValueError(f'Filter field not found: {field_name}')
        if field_name not in self._series_cache:
            if self._dataframe is not None:
                column_index = self._column_positions[field_name]
                self._series_cache[field_name] = self._dataframe.iloc[:, column_index]
            else:
                df_column, _ = read_selected_columns_with_header(
                    self._path,
                    self._delimiter,
                    self._header_line,
                    [field_name],
                )
                if field_name not in df_column.columns:
                    raise ValueError(f'Filter field not found in file after loading: {field_name}')
                self._series_cache[field_name] = df_column[field_name]
        return self._series_cache[field_name]

    def get_categorical_values(self, field, max_values=1000):
        field_name = str(field or '').strip()
        cache_key = (field_name, int(max_values))
        if cache_key not in self._categorical_value_cache:
            series = self.get_series(field_name)
            unique_values = sorted({str(value) for value in series.dropna().astype(str)})
            total_unique_values = len(unique_values)
            truncated = total_unique_values > max_values
            if truncated:
                unique_values = unique_values[:max_values]
            self._categorical_value_cache[cache_key] = (unique_values, total_unique_values, truncated)
        return self._categorical_value_cache[cache_key]

    def get_numeric_range(self, field):
        field_name = str(field or '').strip()
        if field_name not in self._numeric_range_cache:
            numeric_values = pd.to_numeric(self.get_series(field_name), errors='coerce').dropna()
            if len(numeric_values) == 0:
                self._numeric_range_cache[field_name] = None
            else:
                self._numeric_range_cache[field_name] = (float(numeric_values.min()), float(numeric_values.max()))
        return self._numeric_range_cache[field_name]


def resolve_sample_coordinate_columns(header_names, sample_x_col=None, sample_y_col=None, sample_z_col=None):
    available_cols = [c for c in header_names if str(c).strip() != '']
    if sample_x_col and sample_y_col and sample_z_col:
        selected = [sample_x_col, sample_y_col, sample_z_col]
        missing = [c for c in selected if c not in header_names]
        if missing:
            raise ValueError(f"Selected sample coordinate columns not found in file: {missing}")
        if len(set(selected)) != 3:
            raise ValueError(
                "Sample coordinate columns must be three distinct columns. "
                "If you changed the samples file, reselect X, Y, and Z."
            )
        return sample_x_col, sample_y_col, sample_z_col

    if len(available_cols) < 3:
        raise ValueError(
            f"Samples file must provide at least 3 non-empty columns for coordinate mapping. Parsed headers: {header_names}"
        )

    return available_cols[0], available_cols[1], available_cols[2]


def resolve_block_coordinate_columns(header_names, block_x_col=None, block_y_col=None, block_z_col=None):
    available_cols = [c for c in header_names if str(c).strip() != '']
    if block_x_col and block_y_col and block_z_col:
        selected = [block_x_col, block_y_col, block_z_col]
        missing = [c for c in selected if c not in header_names]
        if missing:
            raise ValueError(f"Selected block coordinate columns not found in file: {missing}")
        if len(set(selected)) != 3:
            raise ValueError(
                "Block coordinate columns must be three distinct columns. "
                "If you changed the blocks file, reselect X, Y, and Z."
            )
        return block_x_col, block_y_col, block_z_col

    if len(available_cols) < 3:
        raise ValueError(
            f"Blocks file must provide at least 3 non-empty columns for coordinate mapping. Parsed headers: {header_names}"
        )

    return available_cols[0], available_cols[1], available_cols[2]


def summarize_sample_filter_spec(filter_spec):
    field = str(filter_spec.get('field', '')).strip()
    filter_type = str(filter_spec.get('type', '')).strip().lower()
    if filter_type == 'categorical':
        values = [str(value) for value in filter_spec.get('values', [])]
        preview = ', '.join(values[:5])
        if len(values) > 5:
            preview += ', ...'
        return f"{field} in [{preview}]"
    if filter_type == 'numeric':
        min_value = filter_spec.get('min', None)
        max_value = filter_spec.get('max', None)
        if min_value is None and max_value is None:
            return f"{field}: all numeric values"
        if min_value is None:
            return f"{field} <= {max_value}"
        if max_value is None:
            return f"{field} >= {min_value}"
        return f"{field} in [{min_value}, {max_value}]"
    return field or 'Invalid filter'


def collect_filter_fields(filters):
    fields = []
    for raw_filter in filters or []:
        field = str((raw_filter or {}).get('field', '')).strip()
        if field and field not in fields:
            fields.append(field)
    return fields


def get_configured_sample_filters(config):
    if not config:
        return []
    filters = config.get('sample_filters', None)
    if filters is None:
        filters = config.get('block_domain_sample_filters', [])
    return [dict(entry) for entry in (filters or [])]


def get_configured_block_filters(config):
    if not config:
        return []
    filters = config.get('block_filters', None)
    if filters is None:
        filters = config.get('block_volume_weighted_filters', [])
    return [dict(entry) for entry in (filters or [])]


def _normalize_categorical_filter_token(value):
    if pd.isna(value):
        return ''

    text = str(value).strip()
    if not text:
        return ''

    lower_text = text.lower()
    if '.' not in text and 'e' not in lower_text:
        return text

    try:
        numeric_value = Decimal(text)
    except (InvalidOperation, ValueError):
        return text

    if not numeric_value.is_finite():
        return text

    normalized = format(numeric_value.normalize(), 'f')
    if '.' in normalized:
        normalized = normalized.rstrip('0').rstrip('.')
    if normalized == '-0':
        normalized = '0'
    return normalized or '0'


def apply_dataframe_filters(df_source, filters=None, filter_subject='row', source_label='file',
                            progress_callback=None, progress_label=None, emit_logs=True):
    if not filters:
        return df_source.copy(), []

    filtered_df = df_source.copy()
    input_row_count = int(len(filtered_df))
    applied_filters = []
    total_filters = len(filters)
    progress_text = progress_label or f'Applying {filter_subject} filters...'

    if emit_logs:
        print(
            f"Applying {total_filters:,} {filter_subject} filter(s) to {len(filtered_df):,} rows from {source_label}..."
        )
    if progress_callback:
        progress_callback(0, total_filters, progress_text)

    for filter_number, raw_filter in enumerate(filters, start=1):
        filter_spec = dict(raw_filter or {})
        field = str(filter_spec.get('field', '')).strip()
        filter_type = str(filter_spec.get('type', '')).strip().lower()

        if not field:
            raise ValueError(f'Each {filter_subject} filter must define a field name.')
        if field not in filtered_df.columns:
            raise ValueError(f'{filter_subject.capitalize()} filter field not found in {source_label}: {field}')

        rows_before = int(len(filtered_df))
        series = filtered_df[field]
        if filter_type == 'categorical':
            raw_values = filter_spec.get('values', [])
            allowed_values = [str(value).strip() for value in raw_values]
            if not allowed_values:
                raise ValueError(f'Categorical filter for {field} must include at least one value.')
            series_text = series.fillna('').astype(str).str.strip()
            mask = series_text.isin(allowed_values)

            normalized_allowed_values = {
                _normalize_categorical_filter_token(value)
                for value in raw_values
            }
            normalized_allowed_values.discard('')
            if normalized_allowed_values:
                normalized_series = series.map(_normalize_categorical_filter_token)
                mask |= normalized_series.isin(normalized_allowed_values)
        elif filter_type == 'numeric':
            numeric_series = pd.to_numeric(series, errors='coerce')
            min_value = filter_spec.get('min', None)
            max_value = filter_spec.get('max', None)
            if min_value in ('', None):
                min_value = None
            else:
                min_value = float(min_value)
            if max_value in ('', None):
                max_value = None
            else:
                max_value = float(max_value)
            if min_value is not None and max_value is not None and min_value > max_value:
                raise ValueError(f'Numeric filter for {field} has min greater than max.')

            mask = numeric_series.notna()
            if min_value is not None:
                mask &= numeric_series >= min_value
            if max_value is not None:
                mask &= numeric_series <= max_value
            filter_spec['min'] = min_value
            filter_spec['max'] = max_value
        else:
            raise ValueError(f'Unsupported {filter_subject} filter type for {field}: {filter_type}')

        filtered_df = filtered_df.loc[mask].copy()
        summary = summarize_sample_filter_spec(filter_spec)
        applied_filters.append({
            'field': field,
            'type': filter_type,
            'summary': summary,
        })
        rows_after = int(len(filtered_df))
        if emit_logs:
            print(
                f"Applied {filter_subject} filter {filter_number:,}/{total_filters:,}: {summary}; "
                f"kept {rows_after:,} of {rows_before:,} rows"
            )
        if progress_callback:
            progress_callback(filter_number, total_filters, progress_text)

    if emit_logs:
        print(
            f"{filter_subject.capitalize()} filtering complete: {len(filtered_df):,} of {input_row_count:,} rows "
            f"satisfy all {filter_subject} filters."
        )

    return filtered_df, applied_filters


def apply_sample_filters(df_samples, sample_filters=None, progress_callback=None, progress_label=None):
    return apply_dataframe_filters(
        df_samples,
        filters=sample_filters,
        filter_subject='sample',
        source_label='samples file',
        progress_callback=progress_callback,
        progress_label=progress_label,
    )


def _compute_point_to_set_distance_stats(query_points, reference_points, query_chunk_size=1024,
                                         reference_chunk_size=4096, progress_callback=None,
                                         progress_label='Computing distance statistics',
                                         return_nearest_index=False,
                                         distance_band_edges=None):
    query_points = np.asarray(query_points, dtype=float)
    reference_points = np.asarray(reference_points, dtype=float)
    band_edges = None if distance_band_edges is None else np.asarray(distance_band_edges, dtype=float)
    return_band_counts = band_edges is not None and band_edges.size > 0
    if return_band_counts:
        if np.any(~np.isfinite(band_edges)):
            raise ValueError('Distance band edges must be finite numbers.')
        if np.any(band_edges <= 0):
            raise ValueError('Distance band edges must be greater than zero.')
        if np.any(np.diff(band_edges) <= 0):
            raise ValueError('Distance band edges must be strictly increasing.')

    if len(query_points) == 0:
        nearest_empty = np.empty(0, dtype=float)
        average_empty = np.empty(0, dtype=float)
        band_counts_empty = np.empty((0, len(band_edges) + 1), dtype=int) if return_band_counts else None
        if return_nearest_index:
            if return_band_counts:
                return nearest_empty, average_empty, np.empty(0, dtype=int), band_counts_empty
            return nearest_empty, average_empty, np.empty(0, dtype=int)
        if return_band_counts:
            return nearest_empty, average_empty, band_counts_empty
        return nearest_empty, average_empty
    if len(reference_points) == 0:
        nan_values = np.full(len(query_points), np.nan, dtype=float)
        band_counts_empty = np.zeros((len(query_points), len(band_edges) + 1), dtype=int) if return_band_counts else None
        if return_nearest_index:
            if return_band_counts:
                return nan_values.copy(), nan_values, np.full(len(query_points), -1, dtype=int), band_counts_empty
            return nan_values.copy(), nan_values, np.full(len(query_points), -1, dtype=int)
        if return_band_counts:
            return nan_values.copy(), nan_values, band_counts_empty
        return nan_values.copy(), nan_values

    nearest = np.empty(len(query_points), dtype=float)
    average = np.empty(len(query_points), dtype=float)
    nearest_indices = np.full(len(query_points), -1, dtype=int) if return_nearest_index else None
    band_counts = np.zeros((len(query_points), len(band_edges) + 1), dtype=int) if return_band_counts else None

    for query_start in range(0, len(query_points), query_chunk_size):
        query_end = min(query_start + query_chunk_size, len(query_points))
        query_chunk = query_points[query_start:query_end]
        chunk_nearest = np.full(len(query_chunk), np.inf, dtype=float)
        chunk_distance_sum = np.zeros(len(query_chunk), dtype=float)
        chunk_count = 0
        chunk_nearest_indices = np.full(len(query_chunk), -1, dtype=int) if return_nearest_index else None
        chunk_band_counts = np.zeros((len(query_chunk), len(band_edges) + 1), dtype=int) if return_band_counts else None

        for reference_start in range(0, len(reference_points), reference_chunk_size):
            reference_end = min(reference_start + reference_chunk_size, len(reference_points))
            reference_chunk = reference_points[reference_start:reference_end]
            deltas = query_chunk[:, None, :] - reference_chunk[None, :, :]
            distances = np.sqrt(np.sum(deltas * deltas, axis=2))
            chunk_min_indices = distances.argmin(axis=1)
            chunk_min_values = distances[np.arange(len(query_chunk)), chunk_min_indices]
            update_mask = chunk_min_values < chunk_nearest
            chunk_nearest = np.minimum(chunk_nearest, chunk_min_values)
            if return_nearest_index and np.any(update_mask):
                chunk_nearest_indices[update_mask] = reference_start + chunk_min_indices[update_mask]
            chunk_distance_sum += distances.sum(axis=1)
            chunk_count += distances.shape[1]
            if return_band_counts:
                bucket_indices = np.searchsorted(band_edges, distances, side='right')
                for bucket_index in range(len(band_edges) + 1):
                    chunk_band_counts[:, bucket_index] += np.count_nonzero(bucket_indices == bucket_index, axis=1)

        nearest[query_start:query_end] = chunk_nearest
        average[query_start:query_end] = chunk_distance_sum / max(chunk_count, 1)
        if return_nearest_index:
            nearest_indices[query_start:query_end] = chunk_nearest_indices
        if return_band_counts:
            band_counts[query_start:query_end] = chunk_band_counts
        if progress_callback:
            progress_callback(query_end, len(query_points), progress_label)

    if return_nearest_index:
        if return_band_counts:
            return nearest, average, nearest_indices, band_counts
        return nearest, average, nearest_indices
    if return_band_counts:
        return nearest, average, band_counts
    return nearest, average


def _compute_point_to_set_kdtree_metrics(query_points, reference_points, query_chunk_size=16384,
                                         progress_callback=None,
                                         progress_label='Computing nearest-neighbor statistics',
                                         return_nearest_index=False,
                                         knn_average_k=None):
    query_points = np.asarray(query_points, dtype=float)
    reference_points = np.asarray(reference_points, dtype=float)

    if len(query_points) == 0:
        nearest_empty = np.empty(0, dtype=float)
        average_empty = np.empty(0, dtype=float)
        if return_nearest_index:
            return nearest_empty, average_empty, np.empty(0, dtype=int)
        return nearest_empty, average_empty
    if len(reference_points) == 0:
        nan_values = np.full(len(query_points), np.nan, dtype=float)
        if return_nearest_index:
            return nan_values.copy(), nan_values.copy(), np.full(len(query_points), -1, dtype=int)
        return nan_values.copy(), nan_values.copy()

    try:
        from scipy.spatial import cKDTree
    except Exception as exc:
        raise RuntimeError('scipy.spatial.cKDTree is required for nearest-neighbor block metrics.') from exc

    max_query_k = 1 if knn_average_k is None else max(1, min(int(knn_average_k), len(reference_points)))
    tree = cKDTree(reference_points)
    nearest = np.empty(len(query_points), dtype=float)
    average_knn = np.empty(len(query_points), dtype=float)
    nearest_indices = np.full(len(query_points), -1, dtype=int) if return_nearest_index else None

    for query_start in range(0, len(query_points), query_chunk_size):
        query_end = min(query_start + query_chunk_size, len(query_points))
        query_chunk = query_points[query_start:query_end]
        try:
            chunk_distances, chunk_indices = tree.query(query_chunk, k=max_query_k, workers=-1)
        except TypeError:
            chunk_distances, chunk_indices = tree.query(query_chunk, k=max_query_k)

        chunk_distances = np.asarray(chunk_distances, dtype=float)
        chunk_indices = np.asarray(chunk_indices, dtype=int)
        if max_query_k == 1:
            chunk_distances = chunk_distances.reshape(-1, 1)
            chunk_indices = chunk_indices.reshape(-1, 1)

        nearest[query_start:query_end] = chunk_distances[:, 0]
        average_knn[query_start:query_end] = chunk_distances.mean(axis=1)
        if return_nearest_index:
            nearest_indices[query_start:query_end] = chunk_indices[:, 0]
        if progress_callback:
            progress_callback(query_end, len(query_points), progress_label)

    if return_nearest_index:
        return nearest, average_knn, nearest_indices
    return nearest, average_knn


def _format_metric_distance_token(distance_value):
    rounded_value = round(float(distance_value), 6)
    if math.isclose(rounded_value, round(rounded_value), abs_tol=1e-6):
        return str(int(round(rounded_value)))
    formatted = f'{rounded_value:.6f}'.rstrip('0').rstrip('.')
    return formatted.replace('-', 'neg').replace('.', 'p')


def _build_distance_band_column_names(domain_column_name, distance_step, distance_factor):
    if distance_step is None or distance_factor is None:
        return [], np.empty(0, dtype=float)

    step_value = float(distance_step)
    factor_value = int(distance_factor)
    if not np.isfinite(step_value) or step_value <= 0:
        raise ValueError('Distance count step must be greater than zero.')
    if factor_value <= 0:
        raise ValueError('Distance count max factor must be at least 1.')

    band_edges = np.arange(1, factor_value + 1, dtype=float) * step_value
    column_names = []
    previous_label = '0'
    for edge in band_edges:
        edge_label = _format_metric_distance_token(edge)
        column_names.append(f'{domain_column_name}_Sample_Count_{previous_label}_{edge_label}')
        previous_label = edge_label
    column_names.append(f'{domain_column_name}_Sample_Count_GE_{_format_metric_distance_token(band_edges[-1])}')
    return column_names, band_edges


def _format_distance_band_label(lower_bound, upper_bound):
    lower_label = _format_metric_distance_token(lower_bound)
    if np.isfinite(upper_bound):
        upper_label = _format_metric_distance_token(upper_bound)
        return f'{lower_label}-{upper_label}'
    return f'>={lower_label}'


def _build_block_domain_metrics_summary_dataframe(output_df, domain_column_name,
                                                  nearest_distance_column, distance_band_edges,
                                                  domain_sample_counts=None, block_row_volumes=None):
    if not nearest_distance_column or nearest_distance_column not in output_df.columns:
        return pd.DataFrame()

    if domain_sample_counts is None:
        domain_sample_counts = {}
    if block_row_volumes is None:
        block_row_volumes = np.full(len(output_df), np.nan, dtype=float)

    domain_series = output_df[domain_column_name].fillna('').astype(str).str.strip()
    distance_thresholds = np.asarray(distance_band_edges, dtype=float)
    if distance_thresholds.ndim != 1:
        raise ValueError('Distance thresholds must be a one-dimensional array.')
    summary_rows = []

    for domain in sorted(value for value in domain_series.unique() if value):
        domain_mask = domain_series == domain
        domain_row_volumes = np.asarray(block_row_volumes, dtype=float)[domain_mask.to_numpy()]
        has_domain_volume = bool(np.any(np.isfinite(domain_row_volumes)))
        domain_total_volume = float(np.nansum(domain_row_volumes)) if has_domain_volume else np.nan
        domain_sample_count = int(domain_sample_counts.get(domain, 0))
        domain_sample_density = (
            float(domain_sample_count) / float(domain_total_volume)
            if has_domain_volume and domain_total_volume > 0 else np.nan
        )
        domain_nearest_distances = pd.to_numeric(
            output_df.loc[domain_mask, nearest_distance_column],
            errors='coerce',
        ).to_numpy(dtype=float, copy=False)
        if domain_nearest_distances.size == 0:
            continue

        row = {
            'Domain': domain,
            'Domain_Sample_Count': domain_sample_count,
            'Domain_Block_Volume': domain_total_volume,
            'Domain_Sample_Density': domain_sample_density,
        }

        finite_distance_mask = np.isfinite(domain_nearest_distances)
        for threshold in distance_thresholds:
            threshold_label = _format_metric_distance_token(threshold)
            covered_mask = finite_distance_mask & (domain_nearest_distances <= threshold)
            covered_block_count = int(np.count_nonzero(covered_mask))
            covered_volume = (
                float(np.nansum(domain_row_volumes[covered_mask]))
                if has_domain_volume and covered_block_count > 0 else 0.0
            )
            coverage_fraction = (
                float(covered_volume) / float(domain_total_volume)
                if has_domain_volume and domain_total_volume > 0 else np.nan
            )
            coverage_density = (
                float(domain_sample_count) / float(covered_volume)
                if domain_sample_count > 0 and covered_volume > 0 else np.nan
            )

            row[f'Covered_Block_Count_LEQ_{threshold_label}'] = covered_block_count
            row[f'Covered_Block_Volume_LEQ_{threshold_label}'] = covered_volume
            row[f'Coverage_Fraction_LEQ_{threshold_label}'] = coverage_fraction
            row[f'Coverage_Density_LEQ_{threshold_label}'] = coverage_density

        summary_rows.append(row)

    return pd.DataFrame(summary_rows)


def _compute_average_pairwise_distance(points, query_chunk_size=512, reference_chunk_size=2048,
                                       progress_callback=None, progress_label='Computing pairwise distance statistics'):
    points = np.asarray(points, dtype=float)
    point_count = len(points)
    if point_count < 2:
        return np.nan

    total_distance = 0.0
    pair_count = 0

    for query_start in range(0, point_count, query_chunk_size):
        query_end = min(query_start + query_chunk_size, point_count)
        query_chunk = points[query_start:query_end]

        within_deltas = query_chunk[:, None, :] - query_chunk[None, :, :]
        within_distances = np.sqrt(np.sum(within_deltas * within_deltas, axis=2))
        upper_indices = np.triu_indices(len(query_chunk), k=1)
        if len(upper_indices[0]) > 0:
            total_distance += float(within_distances[upper_indices].sum())
            pair_count += int(len(upper_indices[0]))

        for reference_start in range(query_end, point_count, reference_chunk_size):
            reference_end = min(reference_start + reference_chunk_size, point_count)
            reference_chunk = points[reference_start:reference_end]
            deltas = query_chunk[:, None, :] - reference_chunk[None, :, :]
            distances = np.sqrt(np.sum(deltas * deltas, axis=2))
            total_distance += float(distances.sum())
            pair_count += int(distances.size)

        if progress_callback:
            progress_callback(query_end, point_count, progress_label)

    if pair_count == 0:
        return np.nan
    return total_distance / pair_count


def _compute_average_pairwise_axis_distances(points, query_chunk_size=512, reference_chunk_size=2048,
                                             progress_callback=None, progress_label='Computing pairwise axis distance statistics'):
    points = np.asarray(points, dtype=float)
    point_count = len(points)
    if point_count < 2:
        return np.full(3, np.nan, dtype=float)

    total_axis_distance = np.zeros(3, dtype=float)
    pair_count = 0

    for query_start in range(0, point_count, query_chunk_size):
        query_end = min(query_start + query_chunk_size, point_count)
        query_chunk = points[query_start:query_end]

        within_abs_deltas = np.abs(query_chunk[:, None, :] - query_chunk[None, :, :])
        upper_indices = np.triu_indices(len(query_chunk), k=1)
        if len(upper_indices[0]) > 0:
            total_axis_distance += within_abs_deltas[upper_indices].sum(axis=0)
            pair_count += int(len(upper_indices[0]))

        for reference_start in range(query_end, point_count, reference_chunk_size):
            reference_end = min(reference_start + reference_chunk_size, point_count)
            reference_chunk = points[reference_start:reference_end]
            abs_deltas = np.abs(query_chunk[:, None, :] - reference_chunk[None, :, :])
            total_axis_distance += abs_deltas.sum(axis=(0, 1))
            pair_count += int(abs_deltas.shape[0] * abs_deltas.shape[1])

        if progress_callback:
            progress_callback(query_end, point_count, progress_label)

    if pair_count == 0:
        return np.full(3, np.nan, dtype=float)
    return total_axis_distance / pair_count


def _compute_average_point_to_set_axis_distances(query_points, reference_points, query_chunk_size=1024,
                                                 reference_chunk_size=4096, progress_callback=None,
                                                 progress_label='Computing point-to-set axis distance statistics'):
    query_points = np.asarray(query_points, dtype=float)
    reference_points = np.asarray(reference_points, dtype=float)

    if len(query_points) == 0 or len(reference_points) == 0:
        return np.full(3, np.nan, dtype=float)

    total_axis_distance = np.zeros(3, dtype=float)
    pair_count = 0

    for query_start in range(0, len(query_points), query_chunk_size):
        query_end = min(query_start + query_chunk_size, len(query_points))
        query_chunk = query_points[query_start:query_end]

        for reference_start in range(0, len(reference_points), reference_chunk_size):
            reference_end = min(reference_start + reference_chunk_size, len(reference_points))
            reference_chunk = reference_points[reference_start:reference_end]
            abs_deltas = np.abs(query_chunk[:, None, :] - reference_chunk[None, :, :])
            total_axis_distance += abs_deltas.sum(axis=(0, 1))
            pair_count += int(abs_deltas.shape[0] * abs_deltas.shape[1])

        if progress_callback:
            progress_callback(query_end, len(query_points), progress_label)

    if pair_count == 0:
        return np.full(3, np.nan, dtype=float)
    return total_axis_distance / pair_count


def _cluster_axis_centers(values, tolerance):
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return np.empty(0, dtype=float)

    sorted_values = np.sort(values)
    clustered = [[float(sorted_values[0])]]
    for value in sorted_values[1:]:
        current = float(value)
        if abs(current - clustered[-1][-1]) <= tolerance:
            clustered[-1].append(current)
        else:
            clustered.append([current])

    return np.array([float(np.mean(group)) for group in clustered], dtype=float)


def _infer_subblock_axis_widths(local_axis_values, axis_extent, tolerance):
    axis_extent = float(axis_extent)
    if axis_extent <= 0:
        raise ValueError('Block size components must be greater than zero.')

    centers = _cluster_axis_centers(local_axis_values, tolerance)
    if len(centers) == 0:
        return np.empty(0, dtype=float)
    if len(centers) == 1:
        return np.full(len(local_axis_values), axis_extent, dtype=float)

    boundaries = np.empty(len(centers) + 1, dtype=float)
    boundaries[0] = 0.0
    boundaries[-1] = axis_extent
    boundaries[1:-1] = (centers[:-1] + centers[1:]) / 2.0

    widths = np.diff(boundaries)
    widths = np.clip(widths, 0.0, axis_extent)
    nearest_indices = np.abs(np.asarray(local_axis_values, dtype=float)[:, None] - centers[None, :]).argmin(axis=1)
    return widths[nearest_indices]


def infer_block_row_volumes(block_coords, block_size, group_indices=None,
                            progress_callback=None, progress_label='Inferring block volumes'):
    coords = np.asarray(block_coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError('Block coordinates must be an Nx3 array.')

    if isinstance(block_size, (list, tuple, np.ndarray)):
        unified_dims = np.asarray(block_size, dtype=float)
    else:
        unified_dims = np.array([block_size, block_size, block_size], dtype=float)

    if unified_dims.shape != (3,):
        raise ValueError('Block size must contain exactly three values.')
    if np.any(unified_dims <= 0):
        raise ValueError('Block size values must be greater than zero.')

    if len(coords) == 0:
        empty_indices = np.empty((0, 3), dtype=int)
        return np.empty(0, dtype=float), empty_indices

    all_min_bounds = np.floor(coords.min(axis=0) / unified_dims) * unified_dims
    if group_indices is None:
        base_indices = np.floor((coords - all_min_bounds) / unified_dims + 1e-6).astype(int)
    else:
        base_indices = np.asarray(group_indices, dtype=int)
        if base_indices.shape != coords.shape:
            raise ValueError('Group indices must match the shape of block coordinates.')

    local_coords = coords - (all_min_bounds + base_indices * unified_dims)
    local_coords = np.clip(local_coords, 0.0, unified_dims)
    volumes = np.empty(len(coords), dtype=float)

    unique_base_indices, inverse = np.unique(base_indices, axis=0, return_inverse=True)
    sorted_order = np.argsort(inverse, kind='stable')
    sorted_inverse = inverse[sorted_order]
    group_starts = np.flatnonzero(np.r_[True, sorted_inverse[1:] != sorted_inverse[:-1]])
    group_ends = np.r_[group_starts[1:], len(sorted_order)]
    total_groups = len(group_starts)
    started_at = time.perf_counter()
    next_terminal_report_at = started_at

    for group_number, (group_start, group_end) in enumerate(zip(group_starts, group_ends), start=1):
        group_indices_slice = sorted_order[group_start:group_end]
        group_local = local_coords[group_indices_slice]
        widths = []
        for axis in range(3):
            axis_tol = max(abs(float(unified_dims[axis])) * 1e-6, 1e-6)
            widths.append(
                _infer_subblock_axis_widths(
                    group_local[:, axis],
                    unified_dims[axis],
                    axis_tol,
                )
            )
        volumes[group_indices_slice] = widths[0] * widths[1] * widths[2]

        if progress_callback and (group_number == total_groups or group_number % 500 == 0):
            progress_callback(group_number, total_groups, progress_label)

        now = time.perf_counter()
        if now >= next_terminal_report_at:
            elapsed_seconds = max(int(now - started_at), 0)
            percent = int(round((group_number / max(total_groups, 1)) * 100)) if total_groups > 0 else 100
            current_base_idx = unique_base_indices[sorted_inverse[group_start]].tolist()
            print(
                f"Volume inference progress: {group_number:,}/{total_groups:,} base blocks ({percent}%) processed; "
                f"current base block={current_base_idx}; elapsed~{elapsed_seconds}s"
            )
            next_terminal_report_at = now + 5.0

    return volumes, base_indices


def build_concatenated_sample_ids(sample_df, id_columns):
    selected_columns = []
    for column_name in id_columns or []:
        normalized = str(column_name or '').strip()
        if not normalized or normalized == '(None)' or normalized in selected_columns:
            continue
        if normalized not in sample_df.columns:
            raise ValueError(f'Selected sample ID column not found in samples file: {normalized}')
        selected_columns.append(normalized)

    if not selected_columns:
        return None, []

    def stringify(value):
        if pd.isna(value):
            return ''
        return str(value).strip()

    concatenated = sample_df[selected_columns].apply(
        lambda row: ' | '.join(part for part in (stringify(value) for value in row) if part),
        axis=1,
    )
    return concatenated.to_numpy(dtype=object, copy=False), selected_columns


BLOCK_DOMAIN_METRIC_DEFINITIONS = (
    {
        'id': 'nearest_distance',
        'label': 'Nearest sample distance',
        'description': 'Exact distance from each block to its nearest same-domain sample with a valid value.',
        'cost_label': 'Low',
        'cost_note': 'Exact nearest valid-valued neighbor query using cKDTree; requires a sample value column.',
        'default_checked': True,
    },
    {
        'id': 'average_distance_knn',
        'label': 'Average distance to k nearest samples',
        'description': 'Mean distance from each block to its k nearest same-domain samples.',
        'cost_label': 'Medium',
        'cost_note': 'cKDTree query scales much better than averaging over every sample in the domain.',
        'default_checked': True,
    },
    {
        'id': 'average_distance_exact',
        'label': 'Average distance to all domain samples',
        'description': 'Exact mean distance from each block to every same-domain sample.',
        'cost_label': 'Very high',
        'cost_note': 'Requires all-pairs block-to-sample distances inside each domain.',
        'default_checked': False,
    },
    {
        'id': 'closest_sample_id',
        'label': 'Closest sample ID',
        'description': 'Export the configured sample ID columns from the nearest same-domain sample with a valid value.',
        'cost_label': 'Low',
        'cost_note': 'Reuses the nearest valid-valued sample query; requires a sample value column.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_value',
        'label': 'Nearest sample value',
        'description': 'Export the sample value attached to the nearest same-domain sample.',
        'cost_label': 'Low',
        'cost_note': 'Reuses the nearest-neighbor query; requires a sample value column.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_residual',
        'label': 'Nearest-sample residual',
        'description': 'Nearest same-domain sample value minus block value.',
        'cost_label': 'Low',
        'cost_note': 'Reuses the nearest-neighbor query; requires sample and block value columns.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_abs_residual',
        'label': 'Nearest-sample absolute residual',
        'description': 'Absolute value of the nearest-sample residual.',
        'cost_label': 'Low',
        'cost_note': 'Reuses the nearest-neighbor query; requires sample and block value columns.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_group_block_count',
        'label': 'Nearest-sample group block count',
        'description': 'Number of blocks that share the same nearest same-domain sample.',
        'cost_label': 'Low',
        'cost_note': 'Built from nearest-neighbor assignments; requires sample and block value columns.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_group_mean_residual',
        'label': 'Nearest-sample group mean residual',
        'description': 'Mean residual among blocks anchored to the same nearest sample.',
        'cost_label': 'Low',
        'cost_note': 'Grouped residual aggregation; requires sample and block value columns.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_group_rms_residual',
        'label': 'Nearest-sample group RMS residual',
        'description': 'Root-mean-square residual among blocks anchored to the same nearest sample.',
        'cost_label': 'Low',
        'cost_note': 'Grouped residual aggregation; requires sample and block value columns.',
        'default_checked': False,
    },
    {
        'id': 'nearest_sample_group_std_residual',
        'label': 'Nearest-sample group standardized residual',
        'description': 'Residual divided by the anchored group RMS residual.',
        'cost_label': 'Low',
        'cost_note': 'Grouped residual aggregation; requires sample and block value columns.',
        'default_checked': False,
    },
    {
        'id': 'distance_band_summary',
        'label': 'Distance-threshold summary CSV',
        'description': 'Write a per-domain summary with covered block volume and N/V(d) at each threshold.',
        'cost_label': 'Medium',
        'cost_note': 'Uses nearest distances plus block-volume inference for the summary file.',
        'default_checked': False,
    },
)
BLOCK_DOMAIN_METRIC_DEFINITIONS_BY_ID = {
    entry['id']: entry for entry in BLOCK_DOMAIN_METRIC_DEFINITIONS
}
BLOCK_DOMAIN_METRIC_DEFINITION_ORDER = [
    entry['id'] for entry in BLOCK_DOMAIN_METRIC_DEFINITIONS
]
BLOCK_DOMAIN_METRIC_DEFAULT_SELECTION = [
    entry['id'] for entry in BLOCK_DOMAIN_METRIC_DEFINITIONS if entry.get('default_checked')
]
BLOCK_DOMAIN_METRIC_LEGACY_SELECTION = [
    'nearest_distance',
    'average_distance_exact',
]
BLOCK_DOMAIN_SAMPLE_VALUE_METRICS = {
    'nearest_sample_value',
    'nearest_sample_residual',
    'nearest_sample_abs_residual',
    'nearest_sample_group_block_count',
    'nearest_sample_group_mean_residual',
    'nearest_sample_group_rms_residual',
    'nearest_sample_group_std_residual',
}
BLOCK_DOMAIN_BLOCK_VALUE_METRICS = {
    'nearest_sample_residual',
    'nearest_sample_abs_residual',
    'nearest_sample_group_block_count',
    'nearest_sample_group_mean_residual',
    'nearest_sample_group_rms_residual',
    'nearest_sample_group_std_residual',
}
BLOCK_DOMAIN_NEAREST_NEIGHBOR_METRICS = {
    'nearest_distance',
    'average_distance_knn',
    'closest_sample_id',
    'distance_band_summary',
} | BLOCK_DOMAIN_SAMPLE_VALUE_METRICS


def _normalize_block_domain_metric_selection(selected_metrics, include_legacy_defaults=False):
    if selected_metrics is None:
        metric_ids = (
            BLOCK_DOMAIN_METRIC_LEGACY_SELECTION
            if include_legacy_defaults else
            BLOCK_DOMAIN_METRIC_DEFAULT_SELECTION
        )
    elif isinstance(selected_metrics, str):
        metric_ids = [token.strip() for token in selected_metrics.split(',') if token.strip()]
    else:
        metric_ids = [str(token).strip() for token in (selected_metrics or []) if str(token).strip()]

    normalized = []
    seen = set()
    for metric_id in metric_ids:
        if metric_id in BLOCK_DOMAIN_METRIC_DEFINITIONS_BY_ID and metric_id not in seen:
            normalized.append(metric_id)
            seen.add(metric_id)
    return normalized


def _build_block_domain_metric_column_name(base_name, block_value_column_name='', use_block_value_prefix=True):
    base_label = str(base_name or '').strip()
    if not base_label:
        return ''

    if not use_block_value_prefix:
        return base_label

    prefix = str(block_value_column_name or '').strip()
    if not prefix or prefix == '(None)':
        return base_label

    return f'{prefix}_{base_label}'


def export_block_domain_sample_metrics(samples_file, blocks_file, output_file=None,
                                       samples_delimiter=None, blocks_delimiter=None,
                                       samples_header_line=1, blocks_header_line=1,
                                       sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                       sample_domain_col=None,
                                       sample_value_col=None,
                                       selected_metrics=None,
                                       average_distance_knn_k=8,
                                       closest_sample_id_cols=None,
                                       distance_count_step=None,
                                       distance_count_max_factor=None,
                                       use_block_value_prefix=True,
                                       block_x_col=None, block_y_col=None, block_z_col=None,
                                       block_domain_col=None, block_size=None,
                                       block_value_col=None,
                                       sample_filters=None, block_filters=None, progress_callback=None,
                                       blank_sample_domain_behavior='skip'):
    if not samples_file or not os.path.isfile(samples_file):
        raise ValueError('Please select a valid samples file.')
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')

    domain_column_name = str(block_domain_col or '').strip()
    if not domain_column_name or domain_column_name == '(None)':
        raise ValueError('Please select a domain column in "Blocks Columns" first.')
    sample_domain_column_name = str(sample_domain_col or '').strip()
    use_explicit_sample_domains = bool(sample_domain_column_name and sample_domain_column_name != '(None)')
    sample_value_column_name = str(sample_value_col or '').strip()
    if sample_value_column_name == '(None)':
        sample_value_column_name = ''
    block_value_column_name = str(block_value_col or '').strip()
    if block_value_column_name == '(None)':
        block_value_column_name = ''
    selected_metric_ids = _normalize_block_domain_metric_selection(
        selected_metrics,
        include_legacy_defaults=selected_metrics is None,
    )
    if selected_metrics is None:
        legacy_requested_metrics = []
        if any(str(column_name or '').strip() and str(column_name).strip() != '(None)' for column_name in (closest_sample_id_cols or [])):
            legacy_requested_metrics.append('closest_sample_id')
        if sample_value_column_name and block_value_column_name:
            legacy_requested_metrics.extend(
                metric_id
                for metric_id in BLOCK_DOMAIN_METRIC_DEFINITION_ORDER
                if metric_id in BLOCK_DOMAIN_SAMPLE_VALUE_METRICS
            )
        if distance_count_step is not None and distance_count_max_factor is not None:
            legacy_requested_metrics.append('distance_band_summary')
        for metric_id in legacy_requested_metrics:
            if metric_id not in selected_metric_ids:
                selected_metric_ids.append(metric_id)
    if not selected_metric_ids:
        raise ValueError('Please select at least one block-domain metric to export.')
    selected_metric_set = set(selected_metric_ids)

    wants_nearest_distance_metric = 'nearest_distance' in selected_metric_set
    wants_average_distance_exact = 'average_distance_exact' in selected_metric_set
    wants_average_distance_knn = 'average_distance_knn' in selected_metric_set
    wants_closest_sample_id = 'closest_sample_id' in selected_metric_set
    wants_distance_band_summary = 'distance_band_summary' in selected_metric_set
    wants_any_nearest_sample_metric = bool(selected_metric_set & BLOCK_DOMAIN_SAMPLE_VALUE_METRICS)
    wants_block_value_metrics = bool(selected_metric_set & BLOCK_DOMAIN_BLOCK_VALUE_METRICS)
    wants_value_aware_nearest_sample = wants_nearest_distance_metric or wants_closest_sample_id or wants_any_nearest_sample_metric

    if wants_value_aware_nearest_sample and not sample_value_column_name:
        if wants_closest_sample_id and not (wants_nearest_distance_metric or wants_any_nearest_sample_metric):
            raise ValueError('Select a sample value column or uncheck "Closest sample ID".')
        if wants_nearest_distance_metric and not (wants_closest_sample_id or wants_any_nearest_sample_metric):
            raise ValueError('Select a sample value column or uncheck "Nearest sample distance".')
        raise ValueError('Select a sample value column or uncheck nearest-sample metrics.')
    if wants_block_value_metrics and not block_value_column_name:
        raise ValueError('Select a block value column or uncheck residual-based nearest-sample metrics.')
    if wants_average_distance_knn:
        try:
            average_distance_knn_k = int(average_distance_knn_k)
        except Exception as exc:
            raise ValueError('Average Distance KNN requires an integer k value.') from exc
        if average_distance_knn_k <= 0:
            raise ValueError('Average Distance KNN requires k to be at least 1.')

    blocks_delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
    output_file = resolve_block_domain_metrics_export_path(
        output_file,
        blocks_file=blocks_file,
        domain_col=domain_column_name,
    )
    _emit_progress(progress_callback, 0, 100, 'Preparing block-domain sample metrics export...')
    block_metadata = None
    if not use_explicit_sample_domains:
        if block_size is None:
            raise ValueError('Block size must be specified for block domain metrics when no sample domain column is configured.')

        print(
            f"Using inferred sample-domain matching: samples will be mapped into block domains using block size {block_size}."
        )
        print(f"Loading block domain mapping from {blocks_file}...")
        block_metadata = load_large_blocks_metadata(
            blocks_file,
            blocks_delimiter,
            blocks_header_line or 1,
            block_size,
            None,
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
            block_domain_col=domain_column_name,
            block_filters=block_filters,
            config=None,
            progress_callback=_make_scaled_progress_callback(progress_callback, 0, 30, 'Loading block domain mapping...'),
        )
    else:
        print(
            f"Using explicit domain matching: sample column '{sample_domain_column_name}' -> block column '{domain_column_name}'. "
            f"Skipping block-size domain inference."
        )

    print(f"Loading samples from {samples_file}...")
    df_samples, _ = load_full_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        progress_label='Reading sample file',
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            10 if use_explicit_sample_domains else 30,
            40 if use_explicit_sample_domains else 55,
            'Reading sample file...',
        ),
    )

    sample_filter_progress_start = 41 if use_explicit_sample_domains else 56
    sample_filter_progress_end = 50 if use_explicit_sample_domains else 65
    _emit_progress(progress_callback, sample_filter_progress_start, 100, 'Applying sample filters...')
    filtered_samples_df, applied_filters = apply_sample_filters(
        df_samples,
        sample_filters=sample_filters,
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            sample_filter_progress_start,
            sample_filter_progress_end,
            'Applying sample filters...',
        ),
        progress_label='Applying sample filters...',
    )

    sample_x_col, sample_y_col, sample_z_col = resolve_sample_coordinate_columns(
        list(filtered_samples_df.columns),
        sample_x_col=sample_x_col,
        sample_y_col=sample_y_col,
        sample_z_col=sample_z_col,
    )

    sample_coord_frame = filtered_samples_df[[sample_x_col, sample_y_col, sample_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_sample_mask = sample_coord_frame.notna().all(axis=1)
    valid_sample_coords = sample_coord_frame.loc[valid_sample_mask].to_numpy(copy=False)
    candidate_sample_coords = valid_sample_coords
    candidate_sample_values = None
    candidate_sample_ids = None
    selected_id_columns = []
    if wants_closest_sample_id:
        candidate_sample_ids, selected_id_columns = build_concatenated_sample_ids(filtered_samples_df, closest_sample_id_cols)
        if not selected_id_columns:
            raise ValueError('Select one or more closest-sample ID columns or uncheck "Closest sample ID".')
        if candidate_sample_ids is not None:
            candidate_sample_ids = np.asarray(candidate_sample_ids, dtype=object)[valid_sample_mask.to_numpy()]
    if wants_value_aware_nearest_sample:
        if sample_value_column_name not in filtered_samples_df.columns:
            raise ValueError(f'Selected sample value column not found in samples file: {sample_value_column_name}')
        candidate_sample_values = pd.to_numeric(
            filtered_samples_df.loc[valid_sample_mask, sample_value_column_name],
            errors='coerce',
        ).to_numpy(copy=False)

    if use_explicit_sample_domains:
        if sample_domain_column_name not in filtered_samples_df.columns:
            raise ValueError(f'Selected sample domain column not found in samples file: {sample_domain_column_name}')

        filtered_samples_df, _ = apply_blank_sample_domain_behavior(
            filtered_samples_df,
            blank_domain_behavior=blank_sample_domain_behavior,
            domain_col=sample_domain_column_name,
            x_col=sample_x_col,
            y_col=sample_y_col,
            z_col=sample_z_col,
            blocks_file=blocks_file,
            blocks_delimiter=blocks_delimiter,
            blocks_header_line=blocks_header_line,
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
            block_domain_col=domain_column_name,
            block_size=block_size,
        )
        sample_coord_frame = filtered_samples_df[[sample_x_col, sample_y_col, sample_z_col]].apply(pd.to_numeric, errors='coerce')
        valid_sample_mask = sample_coord_frame.notna().all(axis=1)
        valid_sample_coords = sample_coord_frame.loc[valid_sample_mask].to_numpy(copy=False)
        candidate_sample_coords = valid_sample_coords
        if wants_closest_sample_id and candidate_sample_ids is not None:
            candidate_sample_ids = np.asarray(build_concatenated_sample_ids(filtered_samples_df, closest_sample_id_cols)[0], dtype=object)[valid_sample_mask.to_numpy()]
        if wants_value_aware_nearest_sample:
            candidate_sample_values = pd.to_numeric(
                filtered_samples_df.loc[valid_sample_mask, sample_value_column_name],
                errors='coerce',
            ).to_numpy(copy=False)

        _emit_progress(progress_callback, 44, 100, 'Grouping filtered samples by explicit domain...')
        explicit_domain_values = filtered_samples_df.loc[valid_sample_mask, sample_domain_column_name].fillna('').astype(str).str.strip()
        explicit_domain_values = explicit_domain_values.replace('nan', '')
        candidate_sample_domains = explicit_domain_values.to_numpy(dtype=object, copy=False)
    else:
        sample_coords_for_mapping = valid_sample_coords

        if block_metadata.get('is_rotated') and len(sample_coords_for_mapping) > 0:
            rotation_center = block_metadata['rotation_center']
            rotation_matrix = block_metadata['rotation_matrix']
            sample_coords_for_mapping = (sample_coords_for_mapping - rotation_center) @ rotation_matrix.T

        all_min_bounds = np.asarray(block_metadata.get('grid_index_origin', block_metadata['all_min_bounds']), dtype=float)
        unified_dims = np.asarray(block_metadata['unified_dims'], dtype=float)
        sample_block_indices = np.floor((sample_coords_for_mapping - all_min_bounds) / unified_dims + 1e-6).astype(int)
        domain_mapping = block_metadata['domain_mapping']

        _emit_progress(progress_callback, 58, 100, 'Mapping filtered samples to domains...')
        candidate_sample_domains = np.array(
            [domain_mapping.get((int(idx[0]), int(idx[1]), int(idx[2])), '') for idx in sample_block_indices],
            dtype=object,
        )

    print(f"Loading blocks from {blocks_file}...")
    df_blocks, output_delimiter = load_full_blocks_dataframe(
        blocks_file,
        blocks_delimiter=blocks_delimiter,
        blocks_header_line=blocks_header_line,
        block_filters=block_filters,
        progress_label='Reading blocks file',
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            45 if use_explicit_sample_domains else 60,
            75 if use_explicit_sample_domains else 85,
            'Reading blocks file...',
        ),
    )

    block_x_col, block_y_col, block_z_col = resolve_block_coordinate_columns(
        list(df_blocks.columns),
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
    )
    if domain_column_name not in df_blocks.columns:
        raise ValueError(f'Selected domain column not found in blocks file: {domain_column_name}')
    if wants_block_value_metrics and block_value_column_name not in df_blocks.columns:
        raise ValueError(f'Selected block metric column not found in blocks file: {block_value_column_name}')

    block_coord_frame = df_blocks[[block_x_col, block_y_col, block_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_block_mask = block_coord_frame.notna().all(axis=1)
    block_domain_series = df_blocks[domain_column_name].fillna('').astype(str).str.strip()
    block_domains = {value for value in block_domain_series.unique() if value}
    print(
        f"Finished reading blocks file. Computing domain distance statistics for {len(block_domains):,} domains and {int(valid_block_mask.sum()):,} valid blocks..."
    )

    candidate_domain_mask = np.array(
        [str(domain).strip() in block_domains for domain in candidate_sample_domains],
        dtype=bool,
    ) if len(candidate_sample_domains) else np.empty(0, dtype=bool)
    matched_sample_coords = candidate_sample_coords[candidate_domain_mask]
    matched_sample_domains = np.asarray(candidate_sample_domains, dtype=object)[candidate_domain_mask]
    matched_sample_ids = None if candidate_sample_ids is None else np.asarray(candidate_sample_ids, dtype=object)[candidate_domain_mask]
    matched_sample_values = None if candidate_sample_values is None else np.asarray(candidate_sample_values, dtype=float)[candidate_domain_mask]
    domain_sample_counts = {
        str(domain): int(np.count_nonzero(matched_sample_domains == domain))
        for domain in np.unique(matched_sample_domains)
    }

    distance_count_columns, distance_band_edges = _build_distance_band_column_names(
        domain_column_name,
        distance_count_step if wants_distance_band_summary else None,
        distance_count_max_factor if wants_distance_band_summary else None,
    )
    if distance_count_columns and block_size is None:
        raise ValueError('Block size must be specified when distance-band summary export is enabled.')

    samples_by_domain = {}
    for domain in np.unique(matched_sample_domains):
        domain_mask = matched_sample_domains == domain
        domain_values = None if matched_sample_values is None else matched_sample_values[domain_mask]
        valid_value_mask = None if domain_values is None else np.isfinite(domain_values)
        samples_by_domain[str(domain)] = {
            'coords': matched_sample_coords[domain_mask],
            'ids': None if matched_sample_ids is None else matched_sample_ids[domain_mask],
            'values': domain_values,
            'value_coords': (
                matched_sample_coords[domain_mask][valid_value_mask]
                if valid_value_mask is not None else None
            ),
            'value_values': (
                domain_values[valid_value_mask]
                if valid_value_mask is not None else None
            ),
            'value_ids': (
                matched_sample_ids[domain_mask][valid_value_mask]
                if valid_value_mask is not None and matched_sample_ids is not None else None
            ),
        }

    nn_column = _build_block_domain_metric_column_name(
        f'{domain_column_name}_NN_Distance',
        block_value_column_name,
        use_block_value_prefix,
    )
    avg_column = f'{domain_column_name}_Avg_Distance'
    avg_knn_column = f'{domain_column_name}_Avg_Distance_KNN'
    closest_id_column = (
        _build_block_domain_metric_column_name(
            f'{domain_column_name}_Closest_Sample_ID',
            block_value_column_name,
            use_block_value_prefix,
        )
        if wants_closest_sample_id and selected_id_columns else None
    )
    nearest_sample_value_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Value', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_value' in selected_metric_set else None
    )
    nearest_sample_residual_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Residual', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_residual' in selected_metric_set else None
    )
    nearest_sample_abs_residual_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Abs_Residual', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_abs_residual' in selected_metric_set else None
    )
    nearest_sample_group_block_count_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Group_Block_Count', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_group_block_count' in selected_metric_set else None
    )
    nearest_sample_group_mean_residual_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Group_Mean_Residual', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_group_mean_residual' in selected_metric_set else None
    )
    nearest_sample_group_rms_residual_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Group_RMS_Residual', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_group_rms_residual' in selected_metric_set else None
    )
    nearest_sample_group_std_residual_column = (
        _build_block_domain_metric_column_name('Nearest_Sample_Group_StdResidual', block_value_column_name, use_block_value_prefix)
        if 'nearest_sample_group_std_residual' in selected_metric_set else None
    )

    output_df = df_blocks.copy()
    if wants_nearest_distance_metric or wants_distance_band_summary:
        output_df[nn_column] = np.nan
    if wants_average_distance_exact:
        output_df[avg_column] = np.nan
    if wants_average_distance_knn:
        output_df[avg_knn_column] = np.nan
    if closest_id_column:
        output_df[closest_id_column] = ''
    if nearest_sample_value_column:
        output_df[nearest_sample_value_column] = np.nan
    if nearest_sample_residual_column:
        output_df[nearest_sample_residual_column] = np.nan
    if nearest_sample_abs_residual_column:
        output_df[nearest_sample_abs_residual_column] = np.nan
    if nearest_sample_group_block_count_column:
        output_df[nearest_sample_group_block_count_column] = np.nan
    if nearest_sample_group_mean_residual_column:
        output_df[nearest_sample_group_mean_residual_column] = np.nan
    if nearest_sample_group_rms_residual_column:
        output_df[nearest_sample_group_rms_residual_column] = np.nan
    if nearest_sample_group_std_residual_column:
        output_df[nearest_sample_group_std_residual_column] = np.nan

    processed_block_count = 0
    populated_block_count = 0
    domains_to_process = [value for value in sorted(block_domain_series.unique()) if value]
    total_domain_blocks = int(sum(int((valid_block_mask & (block_domain_series == domain)).sum()) for domain in domains_to_process))
    completed_domain_blocks = 0
    compute_started_at = time.perf_counter()
    next_terminal_report_at = compute_started_at
    _emit_progress(progress_callback, 80 if use_explicit_sample_domains else 85, 100, 'Computing domain distance statistics...')

    needs_exact_distance_query = wants_average_distance_exact
    needs_kdtree_query = wants_average_distance_knn or (
        not wants_average_distance_exact and bool(selected_metric_set & BLOCK_DOMAIN_NEAREST_NEIGHBOR_METRICS)
    )

    for domain in domains_to_process:
        domain_block_mask = valid_block_mask & (block_domain_series == domain)
        domain_block_indices = output_df.index[domain_block_mask]
        if len(domain_block_indices) == 0:
            continue

        processed_block_count += len(domain_block_indices)
        domain_samples = samples_by_domain.get(domain)
        if domain_samples is None or len(domain_samples['coords']) == 0:
            completed_domain_blocks += len(domain_block_indices)
            if total_domain_blocks > 0:
                progress_base = 80 if use_explicit_sample_domains else 85
                progress_span = 18 if use_explicit_sample_domains else 13
                progress_value = progress_base + int(round((completed_domain_blocks / total_domain_blocks) * progress_span))
                _emit_progress(progress_callback, progress_value, 100, f'Computing domain distance statistics... ({domain})')
            now = time.perf_counter()
            if now >= next_terminal_report_at:
                elapsed_seconds = max(int(now - compute_started_at), 0)
                percent = int(round((completed_domain_blocks / max(total_domain_blocks, 1)) * 100)) if total_domain_blocks > 0 else 100
                print(
                    f"Metrics compute progress: {completed_domain_blocks:,}/{total_domain_blocks:,} blocks ({percent}%) processed; "
                    f"current domain={domain}; no matching samples; elapsed~{elapsed_seconds}s"
                )
                next_terminal_report_at = now + 5.0
            continue

        query_points = block_coord_frame.loc[domain_block_indices].to_numpy(copy=False)
        nearest_distances = None
        average_distances = None
        average_knn_distances = None
        nearest_sample_indices = None
        nearest_value_sample_indices = None
        nearest_value_distances = None
        needs_nearest_index = bool(closest_id_column or wants_any_nearest_sample_metric)
        progress_fn = (
            None if total_domain_blocks <= 0 else
            lambda current, total, label, completed=completed_domain_blocks: _emit_progress(
                progress_callback,
                (80 if use_explicit_sample_domains else 85) + int(round(((completed + current) / total_domain_blocks) * (18 if use_explicit_sample_domains else 13))),
                100,
                f'{label}... ({domain})',
            )
        )

        if needs_exact_distance_query:
            distance_stats = _compute_point_to_set_distance_stats(
                query_points,
                domain_samples['coords'],
                progress_callback=progress_fn,
                progress_label='Computing distance statistics (exact average)',
                return_nearest_index=needs_nearest_index,
            )
            if needs_nearest_index:
                nearest_distances, average_distances, nearest_sample_indices = distance_stats
            else:
                nearest_distances, average_distances = distance_stats

        if needs_kdtree_query:
            kd_tree_stats = _compute_point_to_set_kdtree_metrics(
                query_points,
                domain_samples['coords'],
                progress_callback=progress_fn,
                progress_label='Computing distance statistics (nearest-neighbor)',
                return_nearest_index=needs_nearest_index and nearest_sample_indices is None,
                knn_average_k=average_distance_knn_k if wants_average_distance_knn else None,
            )
            if needs_nearest_index and nearest_sample_indices is None:
                kd_nearest_distances, average_knn_distances, nearest_sample_indices = kd_tree_stats
            else:
                kd_nearest_distances, average_knn_distances = kd_tree_stats
            if nearest_distances is None:
                nearest_distances = kd_nearest_distances

        if wants_value_aware_nearest_sample:
            domain_value_coords = domain_samples.get('value_coords')
            if domain_value_coords is not None and len(domain_value_coords) > 0:
                value_tree_stats = _compute_point_to_set_kdtree_metrics(
                    query_points,
                    domain_value_coords,
                    progress_callback=progress_fn,
                    progress_label='Computing nearest valued-sample statistics',
                    return_nearest_index=needs_nearest_index,
                )
                if needs_nearest_index:
                    nearest_value_distances, _, nearest_value_sample_indices = value_tree_stats
                else:
                    nearest_value_distances, _ = value_tree_stats
            else:
                nearest_value_distances = np.full(len(domain_block_indices), np.nan, dtype=float)
                nearest_value_sample_indices = np.full(len(domain_block_indices), -1, dtype=int)

        if wants_nearest_distance_metric or wants_distance_band_summary:
            if wants_nearest_distance_metric and wants_value_aware_nearest_sample:
                nearest_distances = nearest_value_distances
            output_df.loc[domain_block_indices, nn_column] = nearest_distances
        if wants_average_distance_exact:
            output_df.loc[domain_block_indices, avg_column] = average_distances
        if wants_average_distance_knn:
            output_df.loc[domain_block_indices, avg_knn_column] = average_knn_distances
        if closest_id_column:
            domain_ids = domain_samples.get('value_ids')
            closest_id_indices = nearest_value_sample_indices
            if domain_ids is None or closest_id_indices is None:
                domain_ids = domain_samples['ids']
                closest_id_indices = nearest_sample_indices
            if wants_value_aware_nearest_sample and nearest_value_sample_indices is not None:
                domain_ids = domain_samples.get('value_ids')
                closest_id_indices = nearest_value_sample_indices
            closest_ids = [domain_ids[index] if index >= 0 else '' for index in closest_id_indices]
            output_df.loc[domain_block_indices, closest_id_column] = closest_ids
        if wants_any_nearest_sample_metric:
            domain_sample_values = domain_samples.get('value_values')
            nearest_sample_values = np.full(len(domain_block_indices), np.nan, dtype=float)
            valid_nearest_mask = (
                (nearest_value_sample_indices >= 0)
                & (domain_sample_values is not None)
                & (nearest_value_sample_indices < len(domain_sample_values))
            )
            if np.any(valid_nearest_mask):
                nearest_sample_values[valid_nearest_mask] = domain_sample_values[nearest_value_sample_indices[valid_nearest_mask]]
            if nearest_sample_value_column:
                output_df.loc[domain_block_indices, nearest_sample_value_column] = nearest_sample_values

            if wants_block_value_metrics:
                block_values = pd.to_numeric(output_df.loc[domain_block_indices, block_value_column_name], errors='coerce').to_numpy(copy=False)
                residual_values = nearest_sample_values - block_values
                abs_residual_values = np.abs(residual_values)

                group_count_values = np.full(len(domain_block_indices), np.nan, dtype=float)
                group_mean_residual_values = np.full(len(domain_block_indices), np.nan, dtype=float)
                group_rms_residual_values = np.full(len(domain_block_indices), np.nan, dtype=float)
                group_std_residual_values = np.full(len(domain_block_indices), np.nan, dtype=float)

                valid_group_residual_mask = np.isfinite(residual_values) & (nearest_value_sample_indices >= 0)
                if np.any(valid_group_residual_mask):
                    group_counts = {}
                    group_residual_sums = {}
                    group_residual_square_sums = {}
                    for nearest_index, residual_value in zip(
                        nearest_value_sample_indices[valid_group_residual_mask],
                        residual_values[valid_group_residual_mask],
                    ):
                        nearest_index = int(nearest_index)
                        group_counts[nearest_index] = group_counts.get(nearest_index, 0) + 1
                        group_residual_sums[nearest_index] = group_residual_sums.get(nearest_index, 0.0) + float(residual_value)
                        group_residual_square_sums[nearest_index] = (
                            group_residual_square_sums.get(nearest_index, 0.0) + float(residual_value) * float(residual_value)
                        )

                    for nearest_index, group_count in group_counts.items():
                        group_mask = nearest_value_sample_indices == nearest_index
                        mean_residual = group_residual_sums[nearest_index] / float(group_count)
                        rms_residual = math.sqrt(group_residual_square_sums[nearest_index] / float(group_count))
                        group_count_values[group_mask] = float(group_count)
                        group_mean_residual_values[group_mask] = mean_residual
                        group_rms_residual_values[group_mask] = rms_residual
                        if rms_residual > 0.0:
                            group_std_residual_values[group_mask] = residual_values[group_mask] / rms_residual
                        else:
                            zero_mask = group_mask & np.isfinite(residual_values) & np.isclose(residual_values, 0.0)
                            group_std_residual_values[zero_mask] = 0.0

                if nearest_sample_residual_column:
                    output_df.loc[domain_block_indices, nearest_sample_residual_column] = residual_values
                if nearest_sample_abs_residual_column:
                    output_df.loc[domain_block_indices, nearest_sample_abs_residual_column] = abs_residual_values
                if nearest_sample_group_block_count_column:
                    output_df.loc[domain_block_indices, nearest_sample_group_block_count_column] = group_count_values
                if nearest_sample_group_mean_residual_column:
                    output_df.loc[domain_block_indices, nearest_sample_group_mean_residual_column] = group_mean_residual_values
                if nearest_sample_group_rms_residual_column:
                    output_df.loc[domain_block_indices, nearest_sample_group_rms_residual_column] = group_rms_residual_values
                if nearest_sample_group_std_residual_column:
                    output_df.loc[domain_block_indices, nearest_sample_group_std_residual_column] = group_std_residual_values
        populated_block_count += len(domain_block_indices)
        completed_domain_blocks += len(domain_block_indices)
        now = time.perf_counter()
        if now >= next_terminal_report_at:
            elapsed_seconds = max(int(now - compute_started_at), 0)
            percent = int(round((completed_domain_blocks / max(total_domain_blocks, 1)) * 100)) if total_domain_blocks > 0 else 100
            print(
                f"Metrics compute progress: {completed_domain_blocks:,}/{total_domain_blocks:,} blocks ({percent}%) processed; "
                f"current domain={domain}; elapsed~{elapsed_seconds}s"
            )
            next_terminal_report_at = now + 5.0

    elapsed_seconds = max(int(time.perf_counter() - compute_started_at), 0)
    print(
        f"Metrics compute complete: {completed_domain_blocks:,}/{total_domain_blocks:,} blocks processed; "
        f"blocks with matching domain samples={populated_block_count:,}; elapsed~{elapsed_seconds}s"
    )

    block_row_volumes = np.full(len(output_df), np.nan, dtype=float)
    if wants_distance_band_summary and distance_count_columns and block_size is not None:
        valid_block_coords = block_coord_frame.loc[valid_block_mask].to_numpy(copy=False)
        if len(valid_block_coords) > 0:
            _emit_progress(progress_callback, 98, 100, 'Inferring block volumes for summary...')
            volume_coords = valid_block_coords
            rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(volume_coords, block_size_hint=block_size)
            if is_rotated:
                volume_coords = (volume_coords - rotation_center) @ rotation_matrix.T
            row_volumes, _ = infer_block_row_volumes(volume_coords, block_size)
            block_row_volumes[valid_block_mask.to_numpy()] = row_volumes

    summary_output_file = None
    if wants_distance_band_summary:
        summary_df = _build_block_domain_metrics_summary_dataframe(
            output_df,
            domain_column_name,
            nn_column,
            distance_band_edges,
            domain_sample_counts=domain_sample_counts,
            block_row_volumes=block_row_volumes,
        )
    else:
        summary_df = pd.DataFrame()

    if nn_column in output_df.columns and not wants_nearest_distance_metric:
        output_df.drop(columns=[nn_column], inplace=True)

    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    _emit_progress(progress_callback, 99, 100, 'Writing metrics export...')
    output_df.to_csv(output_file, index=False, sep=output_delimiter)
    if not summary_df.empty:
        summary_output_file = resolve_block_domain_metrics_summary_export_path(
            metrics_file=output_file,
            blocks_file=blocks_file,
            domain_col=domain_column_name,
        )
        summary_df.to_csv(summary_output_file, index=False, sep=output_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Block-domain sample metrics export complete.')

    return {
        'output_file': output_file,
        'summary_output_file': summary_output_file,
        'domain_column': domain_column_name,
        'selected_metrics': list(selected_metric_ids),
        'nearest_distance_column': nn_column if wants_nearest_distance_metric else None,
        'average_distance_column': avg_column if wants_average_distance_exact else None,
        'average_distance_knn_column': avg_knn_column if wants_average_distance_knn else None,
        'average_distance_knn_k': int(average_distance_knn_k) if wants_average_distance_knn else None,
        'closest_sample_id_column': closest_id_column,
        'closest_sample_id_source_columns': list(selected_id_columns),
        'use_block_value_prefix': bool(use_block_value_prefix),
        'nearest_sample_value_column': nearest_sample_value_column,
        'nearest_sample_residual_column': nearest_sample_residual_column,
        'nearest_sample_abs_residual_column': nearest_sample_abs_residual_column,
        'nearest_sample_group_block_count_column': nearest_sample_group_block_count_column,
        'nearest_sample_group_mean_residual_column': nearest_sample_group_mean_residual_column,
        'nearest_sample_group_rms_residual_column': nearest_sample_group_rms_residual_column,
        'nearest_sample_group_std_residual_column': nearest_sample_group_std_residual_column,
        'distance_count_columns': [],
        'distance_summary_thresholds': [float(value) for value in np.asarray(distance_band_edges, dtype=float)] if wants_distance_band_summary else [],
        'summary_columns': list(summary_df.columns),
        'summary_row_count': int(len(summary_df)),
        'distance_count_step': None if not wants_distance_band_summary or distance_count_step is None else float(distance_count_step),
        'distance_count_max_factor': None if not wants_distance_band_summary or distance_count_max_factor is None else int(distance_count_max_factor),
        'input_samples': int(len(df_samples)),
        'filtered_samples': int(len(filtered_samples_df)),
        'filters_applied': applied_filters,
        'valid_coordinate_samples': int(valid_sample_mask.sum()),
        'matched_samples': int(len(matched_sample_domains)),
        'unmatched_samples': int(valid_sample_mask.sum() - len(matched_sample_domains)),
        'invalid_coordinate_samples': int((~valid_sample_mask).sum()),
        'processed_blocks': int(processed_block_count),
        'blocks_with_samples_in_domain': int(populated_block_count),
        'invalid_coordinate_blocks': int((~valid_block_mask).sum()),
    }


def export_domain_interpolation_confidence_metrics(samples_file, blocks_file, output_file=None,
                                                   samples_delimiter=None, blocks_delimiter=None,
                                                   samples_header_line=1, blocks_header_line=1,
                                                   sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                                   sample_domain_col=None,
                                                   block_x_col=None, block_y_col=None, block_z_col=None,
                                                   block_domain_col=None, block_size=None,
                                                   sample_filters=None, block_filters=None, progress_callback=None,
                                                   blank_sample_domain_behavior='skip'):
    if not samples_file or not os.path.isfile(samples_file):
        raise ValueError('Please select a valid samples file.')
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')

    domain_column_name = str(block_domain_col or '').strip()
    if not domain_column_name or domain_column_name == '(None)':
        raise ValueError('Please select a domain column in "Blocks Columns" first.')
    sample_domain_column_name = str(sample_domain_col or '').strip()
    use_explicit_sample_domains = bool(sample_domain_column_name and sample_domain_column_name != '(None)')

    blocks_delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
    output_file = resolve_domain_interpolation_confidence_export_path(
        output_file,
        blocks_file=blocks_file,
        domain_col=domain_column_name,
    )
    _emit_progress(progress_callback, 0, 100, 'Preparing domain interpolation confidence export...')

    block_metadata = None
    if not use_explicit_sample_domains:
        if block_size is None:
            raise ValueError('Block size must be specified for interpolation confidence metrics when no sample domain column is configured.')

        print(
            f"Using inferred sample-domain matching: samples will be mapped into block domains using block size {block_size}."
        )
        print(f"Loading block domain mapping from {blocks_file}...")
        block_metadata = load_large_blocks_metadata(
            blocks_file,
            blocks_delimiter,
            blocks_header_line or 1,
            block_size,
            None,
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
            block_domain_col=domain_column_name,
            block_filters=block_filters,
            config=None,
            progress_callback=_make_scaled_progress_callback(progress_callback, 0, 25, 'Loading block domain mapping...'),
        )
    else:
        print(
            f"Using explicit domain matching: sample column '{sample_domain_column_name}' -> block column '{domain_column_name}'. "
            f"Skipping block-size domain inference."
        )

    print(f"Loading samples from {samples_file}...")
    df_samples, _ = load_full_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        progress_label='Reading sample file',
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            10 if use_explicit_sample_domains else 25,
            35 if use_explicit_sample_domains else 50,
            'Reading sample file...',
        ),
    )

    _emit_progress(progress_callback, 36 if use_explicit_sample_domains else 51, 100, 'Applying sample filters...')
    filtered_samples_df, applied_filters = apply_sample_filters(
        df_samples,
        sample_filters=sample_filters,
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            36 if use_explicit_sample_domains else 51,
            45 if use_explicit_sample_domains else 60,
            'Applying sample filters...',
        ),
        progress_label='Applying sample filters...',
    )

    sample_x_col, sample_y_col, sample_z_col = resolve_sample_coordinate_columns(
        list(filtered_samples_df.columns),
        sample_x_col=sample_x_col,
        sample_y_col=sample_y_col,
        sample_z_col=sample_z_col,
    )

    if use_explicit_sample_domains:
        if sample_domain_column_name not in filtered_samples_df.columns:
            raise ValueError(f'Selected sample domain column not found in samples file: {sample_domain_column_name}')

        filtered_samples_df, _ = apply_blank_sample_domain_behavior(
            filtered_samples_df,
            blank_domain_behavior=blank_sample_domain_behavior,
            domain_col=sample_domain_column_name,
            x_col=sample_x_col,
            y_col=sample_y_col,
            z_col=sample_z_col,
            blocks_file=blocks_file,
            blocks_delimiter=blocks_delimiter,
            blocks_header_line=blocks_header_line,
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
            block_domain_col=domain_column_name,
            block_size=block_size,
        )

    sample_coord_frame = filtered_samples_df[[sample_x_col, sample_y_col, sample_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_sample_mask = sample_coord_frame.notna().all(axis=1)
    valid_sample_coords = sample_coord_frame.loc[valid_sample_mask].to_numpy(copy=False)

    if use_explicit_sample_domains:
        _emit_progress(progress_callback, 46, 100, 'Grouping filtered samples by explicit domain...')
        explicit_domain_values = filtered_samples_df.loc[valid_sample_mask, sample_domain_column_name].fillna('').astype(str).str.strip()
        explicit_domain_values = explicit_domain_values.replace('nan', '')
        candidate_sample_domains = explicit_domain_values.to_numpy(dtype=object, copy=False)
    else:
        sample_coords_for_mapping = valid_sample_coords
        if block_metadata.get('is_rotated') and len(sample_coords_for_mapping) > 0:
            rotation_center = block_metadata['rotation_center']
            rotation_matrix = block_metadata['rotation_matrix']
            sample_coords_for_mapping = (sample_coords_for_mapping - rotation_center) @ rotation_matrix.T

        all_min_bounds = np.asarray(block_metadata.get('grid_index_origin', block_metadata['all_min_bounds']), dtype=float)
        unified_dims = np.asarray(block_metadata['unified_dims'], dtype=float)
        sample_block_indices = np.floor((sample_coords_for_mapping - all_min_bounds) / unified_dims + 1e-6).astype(int)
        domain_mapping = block_metadata['domain_mapping']

        _emit_progress(progress_callback, 61, 100, 'Mapping filtered samples to domains...')
        candidate_sample_domains = np.array(
            [domain_mapping.get((int(idx[0]), int(idx[1]), int(idx[2])), '') for idx in sample_block_indices],
            dtype=object,
        )

    print(f"Loading blocks from {blocks_file}...")
    df_blocks, output_delimiter = load_full_blocks_dataframe(
        blocks_file,
        blocks_delimiter=blocks_delimiter,
        blocks_header_line=blocks_header_line,
        block_filters=block_filters,
        progress_label='Reading blocks file',
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            50 if use_explicit_sample_domains else 65,
            75 if use_explicit_sample_domains else 85,
            'Reading blocks file...',
        ),
    )

    block_x_col, block_y_col, block_z_col = resolve_block_coordinate_columns(
        list(df_blocks.columns),
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
    )
    if domain_column_name not in df_blocks.columns:
        raise ValueError(f'Selected domain column not found in blocks file: {domain_column_name}')

    block_coord_frame = df_blocks[[block_x_col, block_y_col, block_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_block_mask = block_coord_frame.notna().all(axis=1)
    block_domain_series = df_blocks[domain_column_name].fillna('').astype(str).str.strip()
    block_domains = [value for value in sorted(block_domain_series.unique()) if value]
    block_domain_set = set(block_domains)

    candidate_domain_mask = np.array(
        [str(domain).strip() in block_domain_set for domain in candidate_sample_domains],
        dtype=bool,
    ) if len(candidate_sample_domains) else np.empty(0, dtype=bool)
    matched_sample_coords = valid_sample_coords[candidate_domain_mask]
    matched_sample_domains = np.asarray(candidate_sample_domains, dtype=object)[candidate_domain_mask]

    valid_block_count = int(valid_block_mask.sum())
    print(
        f"Finished reading blocks file. Computing interpolation confidence metrics for {len(block_domains):,} domains, "
        f"{len(matched_sample_coords):,} matched samples, and {valid_block_count:,} valid blocks..."
    )

    domain_rows = []
    axis_labels = ('X', 'Y', 'Z')
    _emit_progress(progress_callback, 80 if use_explicit_sample_domains else 88, 100, 'Computing interpolation confidence metrics...')

    for domain_index, domain in enumerate(block_domains, start=1):
        progress_base = 80 if use_explicit_sample_domains else 88
        progress_span = 18 if use_explicit_sample_domains else 10
        progress_value = progress_base + int(round((domain_index / max(len(block_domains), 1)) * progress_span))
        _emit_progress(progress_callback, progress_value, 100, f'Computing interpolation confidence metrics... ({domain})')

        domain_sample_mask = matched_sample_domains == domain
        domain_sample_coords = matched_sample_coords[domain_sample_mask]
        domain_block_mask = valid_block_mask & (block_domain_series == domain)
        domain_block_coords = block_coord_frame.loc[domain_block_mask].to_numpy(copy=False)

        avg_source_sample_distance = _compute_average_pairwise_distance(domain_sample_coords)
        avg_source_sample_axis_distances = _compute_average_pairwise_axis_distances(domain_sample_coords)
        if len(domain_block_coords) > 0 and len(domain_sample_coords) > 0:
            _, block_average_distances = _compute_point_to_set_distance_stats(
                domain_block_coords,
                domain_sample_coords,
            )
            avg_block_to_source_sample_distance = float(np.nanmean(block_average_distances))
            avg_block_to_source_sample_axis_distances = _compute_average_point_to_set_axis_distances(
                domain_block_coords,
                domain_sample_coords,
            )
        else:
            avg_block_to_source_sample_distance = np.nan
            avg_block_to_source_sample_axis_distances = np.full(3, np.nan, dtype=float)

        if (
            np.isfinite(avg_source_sample_distance)
            and np.isfinite(avg_block_to_source_sample_distance)
            and avg_block_to_source_sample_distance != 0
        ):
            sample_to_block_ratio = avg_source_sample_distance / avg_block_to_source_sample_distance
        else:
            sample_to_block_ratio = np.nan

        row = {
            'Domain': domain,
            'Source_Sample_Count': int(len(domain_sample_coords)),
            'Domain_Block_Count': int(len(domain_block_coords)),
            'Avg_Source_Sample_Distance': avg_source_sample_distance,
            'Avg_Block_To_Source_Sample_Distance': avg_block_to_source_sample_distance,
            'Sample_To_Block_Distance_Ratio': sample_to_block_ratio,
        }
        for axis_index, axis_label in enumerate(axis_labels):
            row[f'Avg_Source_Sample_Distance_{axis_label}'] = avg_source_sample_axis_distances[axis_index]
            row[f'Avg_Block_To_Source_Sample_Distance_{axis_label}'] = avg_block_to_source_sample_axis_distances[axis_index]
        domain_rows.append(row)

    output_df = pd.DataFrame(domain_rows)
    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    _emit_progress(progress_callback, 99, 100, 'Writing interpolation confidence export...')
    output_df.to_csv(output_file, index=False, sep=output_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Domain interpolation confidence export complete.')

    return {
        'output_file': output_file,
        'domain_column': domain_column_name,
        'domain_count': int(len(domain_rows)),
        'input_samples': int(len(df_samples)),
        'filtered_samples': int(len(filtered_samples_df)),
        'filters_applied': applied_filters,
        'valid_coordinate_samples': int(valid_sample_mask.sum()),
        'matched_samples': int(len(matched_sample_domains)),
        'unmatched_samples': int(valid_sample_mask.sum() - len(matched_sample_domains)),
        'invalid_coordinate_samples': int((~valid_sample_mask).sum()),
        'processed_blocks': int(valid_block_count),
        'invalid_coordinate_blocks': int((~valid_block_mask).sum()),
        'columns': list(output_df.columns),
    }


def export_block_volume_weighted_average(blocks_file, value_col, output_file=None,
                                         blocks_delimiter=None, blocks_header_line=1,
                                         block_x_col=None, block_y_col=None, block_z_col=None,
                                         block_domain_col=None,
                                         weight_col=None, block_filters=None,
                                         block_size=None, progress_callback=None):
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')

    value_column_name = str(value_col or '').strip()
    if not value_column_name or value_column_name == '(None)':
        raise ValueError('Please select a value column in the blocks file.')

    selected_weight_column = str(weight_col or '').strip()
    use_volume_weights = selected_weight_column in ('', '(None)', '(Volume)')
    weight_column_name = None if use_volume_weights else selected_weight_column
    if use_volume_weights and block_size is None:
        raise ValueError('Block size must be specified to infer sub-block volumes.')

    output_file = resolve_block_volume_weighted_average_export_path(
        output_file,
        blocks_file=blocks_file,
        value_col=value_column_name,
    )
    _emit_progress(progress_callback, 0, 100, 'Preparing block volume export...')

    print(f"Loading blocks from {blocks_file}...")
    df_blocks, output_delimiter = load_full_blocks_dataframe(
        blocks_file,
        blocks_delimiter=blocks_delimiter,
        blocks_header_line=blocks_header_line,
        progress_label='Reading blocks file',
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 45, 'Reading blocks file...'),
    )

    input_block_count = int(len(df_blocks))
    block_filter_progress_start = 46
    block_filter_progress_end = 54
    _emit_progress(progress_callback, block_filter_progress_start, 100, 'Applying block filters...')
    df_blocks, applied_filters = apply_dataframe_filters(
        df_blocks,
        filters=block_filters,
        filter_subject='block',
        source_label='blocks file',
        progress_callback=_make_scaled_progress_callback(
            progress_callback,
            block_filter_progress_start,
            block_filter_progress_end,
            'Applying block filters...',
        ),
        progress_label='Applying block filters...',
    )
    filtered_block_count = int(len(df_blocks))
    print(
        f"Block filtering complete: {filtered_block_count:,} of {input_block_count:,} rows remain after filtering."
    )

    if value_column_name not in df_blocks.columns:
        raise ValueError(f'Selected value column not found in blocks file: {value_column_name}')
    if weight_column_name and weight_column_name not in df_blocks.columns:
        raise ValueError(f'Selected weight column not found in blocks file: {weight_column_name}')

    domain_column_name = str(block_domain_col or '').strip()
    use_domain_column = bool(domain_column_name and domain_column_name != '(None)')
    if use_domain_column and domain_column_name not in df_blocks.columns:
        raise ValueError(f'Selected domain column not found in blocks file: {domain_column_name}')

    _emit_progress(progress_callback, 55, 100, f"Preparing weighted column '{value_column_name}'...")
    print(f"Converting weighted column '{value_column_name}' to numeric...")
    value_series = pd.to_numeric(df_blocks[value_column_name], errors='coerce')
    if use_domain_column:
        domain_series = df_blocks[domain_column_name].fillna('').astype(str).str.strip()
    else:
        domain_series = pd.Series([''] * len(df_blocks), index=df_blocks.index, dtype=object)

    if use_volume_weights:
        _emit_progress(progress_callback, 58, 100, 'Preparing block coordinates...')
        block_x_col, block_y_col, block_z_col = resolve_block_coordinate_columns(
            list(df_blocks.columns),
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
        )
        print(
            f"Converting block coordinates to numeric using columns: {block_x_col}, {block_y_col}, {block_z_col}..."
        )
        block_coord_frame = df_blocks[[block_x_col, block_y_col, block_z_col]].apply(pd.to_numeric, errors='coerce')
        valid_coord_mask = block_coord_frame.notna().all(axis=1)
        valid_coords = block_coord_frame.loc[valid_coord_mask].to_numpy(copy=False)
        print(
            f"Coordinate preparation complete: {int(valid_coord_mask.sum()):,} valid rows, "
            f"{int((~valid_coord_mask).sum()):,} invalid rows."
        )

        if len(valid_coords) > 0:
            _emit_progress(progress_callback, 60, 100, f'Detecting grid rotation... ({len(valid_coords):,} valid rows)')
            print(f"Starting rotation detection for {len(valid_coords):,} valid block rows...")
            rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(valid_coords, block_size_hint=block_size)
            if is_rotated:
                valid_coords = (valid_coords - rotation_center) @ rotation_matrix.T

        _emit_progress(progress_callback, 62, 100, 'Inferring block volumes...')
        inferred_volumes = np.full(len(df_blocks), np.nan, dtype=float)

        if len(valid_coords) > 0:
            print(f"Inferring sub-block volumes for {len(valid_coords):,} valid block rows...")
            row_volumes, _ = infer_block_row_volumes(
                valid_coords,
                block_size,
                progress_callback=_make_scaled_progress_callback(progress_callback, 62, 92, 'Inferring block volumes...'),
                progress_label='Inferring block volumes...',
            )
            inferred_volumes[valid_coord_mask.to_numpy()] = row_volumes

        weight_values = inferred_volumes
        valid_weight_mask = (
            valid_coord_mask.to_numpy()
            & value_series.notna().to_numpy()
            & np.isfinite(inferred_volumes)
        )
        total_volume = float(np.sum(inferred_volumes[valid_weight_mask])) if np.any(valid_weight_mask) else 0.0
    else:
        _emit_progress(progress_callback, 58, 100, f"Preparing custom weight column '{weight_column_name}'...")
        print(f"Converting custom weight column '{weight_column_name}' to numeric...")
        weight_series = pd.to_numeric(df_blocks[weight_column_name], errors='coerce')
        weight_values = weight_series.to_numpy(dtype=float, copy=False)
        valid_coord_mask = pd.Series([True] * len(df_blocks), index=df_blocks.index)
        inferred_volumes = np.full(len(df_blocks), np.nan, dtype=float)
        valid_weight_mask = value_series.notna().to_numpy() & np.isfinite(weight_values)
        total_volume = float('nan')

    _emit_progress(progress_callback, 93, 100, 'Summarizing weighted averages...')

    weighted_components = weight_values * value_series.to_numpy(dtype=float, copy=False)
    total_weight = float(np.sum(weight_values[valid_weight_mask])) if np.any(valid_weight_mask) else 0.0
    weighted_sum = float(np.sum(weighted_components[valid_weight_mask])) if np.any(valid_weight_mask) else 0.0
    weighted_average = float(weighted_sum / total_weight) if total_weight > 0 else np.nan

    domain_summaries = {}

    if use_domain_column:
        domain_values = domain_series.to_numpy(dtype=object, copy=False)
        valid_domain_mask = valid_weight_mask & domain_series.ne('').to_numpy()
        for domain in sorted(domain_series.loc[valid_domain_mask].unique()):
            domain_mask = valid_domain_mask & (domain_values == domain)
            domain_weight = float(np.sum(weight_values[domain_mask])) if np.any(domain_mask) else 0.0
            domain_sum = float(np.sum(weighted_components[domain_mask])) if np.any(domain_mask) else 0.0
            domain_average = float(domain_sum / domain_weight) if domain_weight > 0 else np.nan
            domain_volume = float(np.sum(inferred_volumes[domain_mask])) if use_volume_weights and np.any(domain_mask) else float('nan')
            domain_summaries[str(domain)] = {
                'total_volume': domain_volume,
                'total_weight': domain_weight,
                'weighted_sum': domain_sum,
                'weighted_average': domain_average,
                'rows_with_numeric_value': int(np.count_nonzero(domain_mask)),
            }

    summary_rows = []
    if use_domain_column:
        for domain, summary in sorted(domain_summaries.items()):
            summary_rows.append(
                {
                    domain_column_name: domain,
                    'Weight_Column': weight_column_name if weight_column_name else 'Volume',
                    'Weighted_Average': summary['weighted_average'],
                    'Weighted_Sum': summary['weighted_sum'],
                    'Total_Weight': summary['total_weight'],
                    'Total_Volume': summary['total_volume'],
                    'Rows_With_Numeric_Value': summary['rows_with_numeric_value'],
                }
            )
    else:
        summary_rows.append(
            {
                'Weight_Column': weight_column_name if weight_column_name else 'Volume',
                'Weighted_Average': weighted_average,
                'Weighted_Sum': weighted_sum,
                'Total_Weight': total_weight,
                'Total_Volume': total_volume,
                'Rows_With_Numeric_Value': int(np.count_nonzero(valid_weight_mask)),
            }
        )

    output_df = pd.DataFrame(summary_rows)

    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    _emit_progress(progress_callback, 99, 100, 'Writing block volume export...')
    output_df.to_csv(output_file, index=False, sep=output_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Block volume export complete.')

    print(f"Exported {len(output_df):,} summary rows to {output_file}")
    print(
        f"Weighted average for '{value_column_name}' using '{weight_column_name or 'Volume'}': {weighted_average} "
        f"(weighted sum={weighted_sum}, total weight={total_weight})"
    )

    return {
        'output_file': output_file,
        'value_column': value_column_name,
        'weight_column': weight_column_name,
        'domain_column': domain_column_name if use_domain_column else None,
        'domain_summaries': domain_summaries,
        'filters_applied': applied_filters,
        'weighted_average': weighted_average,
        'weighted_sum': weighted_sum,
        'total_weight': total_weight,
        'total_volume': total_volume,
        'input_blocks': input_block_count,
        'filtered_blocks': filtered_block_count,
        'processed_rows': int(len(df_blocks)) if not use_volume_weights else int(valid_coord_mask.sum()),
        'rows_with_numeric_value': int(valid_weight_mask.sum()),
        'exported_rows': int(len(output_df)),
        'invalid_coordinate_rows': 0 if not use_volume_weights else int((~valid_coord_mask).sum()),
        'invalid_value_rows': int(value_series.isna().sum()),
        'invalid_weight_rows': int(np.count_nonzero(~np.isfinite(weight_values))) if len(df_blocks) > 0 else 0,
    }


def export_sample_blocks_from_samples_and_blocks(samples_file, blocks_file, output_file=None,
                                                 samples_delimiter=None, blocks_delimiter=None,
                                                 samples_header_line=1, blocks_header_line=1,
                                                 sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                                 sample_value_col=None, sample_domain_col=None,
                                                 sample_weight_col=None,
                                                 include_sample_ids=False,
                                                 sample_id_cols=None,
                                                 block_x_col=None, block_y_col=None, block_z_col=None,
                                                 block_domain_col=None, block_size=None,
                                                 sample_filters=None, block_filters=None,
                                                 blank_sample_domain_behavior='skip',
                                                 progress_callback=None):
    if not samples_file or not os.path.isfile(samples_file):
        raise ValueError('Please select a valid samples file.')
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')

    value_column_name = str(sample_value_col or '').strip()
    if not value_column_name or value_column_name == '(None)':
        raise ValueError('Please select a value column in "Samples Columns" first.')
    if block_size is None:
        raise ValueError('Block size must be specified for sample-block export.')

    output_file = resolve_sample_blocks_export_path(output_file, samples_file=samples_file)
    sample_delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
    selected_weight_column = str(sample_weight_col or '').strip()
    use_sample_weights = bool(selected_weight_column and selected_weight_column != '(None)')
    use_domain_matching = bool(
        str(sample_domain_col or '').strip()
        and str(sample_domain_col).strip() != '(None)'
        and str(block_domain_col or '').strip()
        and str(block_domain_col).strip() != '(None)'
    )

    _emit_progress(progress_callback, 0, 100, 'Preparing sample-block export...')
    print(f"Loading block assignment metadata from {blocks_file}...")
    block_metadata = _load_block_assignment_metadata(
        blocks_file,
        block_size,
        blocks_delimiter=blocks_delimiter,
        blocks_header_line=blocks_header_line,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=block_domain_col,
        block_filters=block_filters,
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 40, 'Loading block assignment metadata...'),
    )

    print(f"Loading samples from {samples_file}...")
    df_samples, parsed_cols, explicit_sample_map = load_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        sample_x_col=sample_x_col,
        sample_y_col=sample_y_col,
        sample_z_col=sample_z_col,
        sample_value_col=sample_value_col,
        sample_domain_col=sample_domain_col,
        sample_filters=sample_filters,
        progress_label='Reading sample file',
        extra_columns=[
            *([selected_weight_column] if use_sample_weights else []),
            *([column_name for column_name in (sample_id_cols or [])] if include_sample_ids else []),
        ],
    )
    if parsed_cols is not None:
        print(f"Samples file (custom header line {samples_header_line}) parsed columns: {parsed_cols}")
    elif hasattr(df_samples, '_detected_delimiter'):
        print(f"Samples file delimiter used: '{df_samples._detected_delimiter}'")

    df_samples = normalize_selected_sample_domain_column(df_samples, sample_domain_col=sample_domain_col)
    df_samples = normalize_selected_sample_weight_column(
        df_samples,
        sample_weight_col=sample_weight_col,
        sample_value_col=sample_value_col,
    )
    if explicit_sample_map:
        print(f"Applied user sample column mapping: {explicit_sample_map}")
    elif sample_x_col and sample_y_col and sample_z_col and sample_value_col:
        rename_map = {
            sample_x_col: 'x',
            sample_y_col: 'y',
            sample_z_col: 'z',
            sample_value_col: 'Value',
        }
        df_samples = df_samples.rename(columns=rename_map)
        print(f"Applied user sample column mapping: {rename_map}")
    else:
        expected = ['x', 'y', 'z', 'Value']
        if any(column_name not in df_samples.columns for column_name in expected):
            if 'Value' not in df_samples.columns:
                for column_name in df_samples.columns:
                    if str(column_name).strip().lower() == 'value':
                        df_samples = df_samples.rename(columns={column_name: 'Value'})
                        break
            if any(column_name not in df_samples.columns for column_name in expected):
                first_four = list(df_samples.columns[:4])
                if len(first_four) != 4:
                    raise ValueError('Samples file must have at least four columns for automatic mapping (x, y, z, Value).')
                rename_map = {original: mapped for original, mapped in zip(first_four, expected)}
                df_samples = df_samples.rename(columns=rename_map)
                print(f"Mapped first four sample columns to {expected}: {rename_map}")

    if use_domain_matching:
        df_samples, _ = apply_blank_sample_domain_behavior(
            df_samples,
            blank_domain_behavior=blank_sample_domain_behavior,
            domain_col='Domain',
            x_col='x',
            y_col='y',
            z_col='z',
            blocks_file=blocks_file,
            blocks_delimiter=blocks_delimiter,
            blocks_header_line=blocks_header_line,
            block_x_col=block_x_col,
            block_y_col=block_y_col,
            block_z_col=block_z_col,
            block_domain_col=block_domain_col,
            block_size=block_size,
            block_filters=block_filters,
        )

    coord_frame = df_samples[['x', 'y', 'z']].apply(pd.to_numeric, errors='coerce')
    value_series = pd.to_numeric(df_samples['Value'], errors='coerce')
    weight_series = None
    if use_sample_weights:
        weight_series = pd.to_numeric(df_samples['Weight'], errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1) & value_series.notna()
    if weight_series is not None:
        valid_mask &= weight_series.notna() & np.isfinite(weight_series) & (weight_series > 0.0)
    if not valid_mask.any():
        if use_sample_weights:
            raise ValueError('No samples remain with valid numeric coordinates, values, and positive weights after preprocessing.')
        raise ValueError('No samples remain with valid numeric coordinates and values after preprocessing.')

    valid_coords = coord_frame.loc[valid_mask].to_numpy(copy=False)
    valid_values = value_series.loc[valid_mask].to_numpy(dtype=float, copy=False)
    valid_weights = None if weight_series is None else weight_series.loc[valid_mask].to_numpy(dtype=float, copy=False)
    sample_ids = None
    selected_sample_id_columns = []
    if include_sample_ids:
        sample_ids, selected_sample_id_columns = build_concatenated_sample_ids(df_samples, sample_id_cols)
        if not selected_sample_id_columns:
            raise ValueError('Select one or more sample ID columns or disable sample ID export.')
        if sample_ids is not None:
            sample_ids = np.asarray(sample_ids, dtype=object)[valid_mask.to_numpy()]
    sample_domains = None
    if use_domain_matching:
        sample_domains = df_samples.loc[valid_mask, 'Domain'].fillna('').astype(str).str.strip().replace('nan', '')
        sample_domains = sample_domains.to_numpy(dtype=object, copy=False)

    coords_for_assignment = np.array(valid_coords, copy=True)
    if block_metadata.get('is_rotated') and len(coords_for_assignment) > 0:
        rotation_center = block_metadata['rotation_center']
        rotation_matrix = block_metadata['rotation_matrix']
        coords_for_assignment = (coords_for_assignment - rotation_center) @ rotation_matrix.T

    sample_block_data = aggregate_samples_into_blocks(
        coords_for_assignment,
        valid_values,
        block_metadata.get('grid_index_origin', block_metadata['all_min_bounds']),
        block_metadata['unified_dims'],
        allowed_grid=block_metadata.get('allowed_grid'),
        domain_mapping=block_metadata['domain_mapping'] if use_domain_matching else None,
        sample_domains=sample_domains if use_domain_matching else None,
        sample_ids=sample_ids,
        sample_weights=valid_weights,
        progress_label='Assigning samples to blocks',
    )

    sample_block_values = sample_block_data['sample_block_values']
    sample_block_counts = sample_block_data['sample_block_counts']
    if not sample_block_values:
        raise ValueError('No sample blocks were generated. Check block size, filters, and domain matching configuration.')

    row_domain_mapping = None
    selected_block_domain_col = str(block_domain_col or '').strip()
    if selected_block_domain_col and selected_block_domain_col != '(None)':
        row_domain_mapping = block_metadata['domain_mapping']
    elif sample_domains is not None:
        row_domain_mapping = _select_majority_sample_block_domains(sample_block_data['sample_block_domain_counts'])

    output_value_column_name = f'BLK_{value_column_name}' if value_column_name else 'BLK_Value'
    output_domain_column_name = 'Domain'
    selected_sample_domain_col = str(sample_domain_col or '').strip()
    if selected_sample_domain_col and selected_sample_domain_col != '(None)':
        output_domain_column_name = selected_sample_domain_col
    elif selected_block_domain_col and selected_block_domain_col != '(None)':
        output_domain_column_name = selected_block_domain_col

    sample_count_column_name = f'{value_column_name}_SampleCount' if value_column_name else 'Value_SampleCount'
    sample_ids_column_name = f'{value_column_name}_SampleIDs' if value_column_name else 'Value_SampleIDs'
    rows = _build_sample_block_rows(
        sample_block_values,
        sample_block_counts,
        block_metadata.get('grid_index_origin', block_metadata['all_min_bounds']),
        block_metadata['unified_dims'],
        rotation_matrix=block_metadata.get('rotation_matrix') if block_metadata.get('is_rotated') else None,
        rotation_center=block_metadata.get('rotation_center') if block_metadata.get('is_rotated') else None,
        domain_mapping=row_domain_mapping,
        value_column_name=output_value_column_name,
        domain_column_name=output_domain_column_name,
        sample_count_column_name=sample_count_column_name,
        sample_ids_by_block=sample_block_data.get('sample_block_ids') if include_sample_ids else None,
        sample_ids_column_name=sample_ids_column_name,
    )
    output_df = pd.DataFrame(rows)

    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    _emit_progress(progress_callback, 95, 100, 'Writing sample-block export...')
    output_df.to_csv(output_file, index=False, sep=sample_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Sample-block export complete.')

    assigned_sample_count = int(np.count_nonzero(sample_block_data['assigned_mask']))
    valid_sample_count = int(valid_mask.sum())
    invalid_coordinate_or_value_count = int((~valid_mask).sum())
    rejected_sample_count = valid_sample_count - assigned_sample_count
    print(f"Exported {len(output_df):,} sample blocks to {output_file}")
    print(
        f"Sample-block export summary: input rows={len(df_samples):,}; valid rows={valid_sample_count:,}; "
        f"assigned samples={assigned_sample_count:,}; sample blocks={len(output_df):,}."
    )
    if use_sample_weights:
        print(f"Sample-block values were computed as weighted averages using sample column '{selected_weight_column}'.")

    return {
        'output_file': output_file,
        'total_samples': int(len(df_samples)),
        'valid_samples': valid_sample_count,
        'assigned_samples': assigned_sample_count,
        'rejected_samples': int(rejected_sample_count),
        'invalid_coordinate_or_value_samples': invalid_coordinate_or_value_count,
        'domain_mismatch_samples': int(sample_block_data['domain_mismatch_count']),
        'sample_block_count': int(len(output_df)),
        'weight_column': selected_weight_column if use_sample_weights else None,
        'sample_ids_column': sample_ids_column_name if include_sample_ids else None,
        'sample_id_source_columns': list(selected_sample_id_columns),
        'columns': list(output_df.columns),
    }


def export_domained_samples_from_blocks(samples_file, blocks_file, output_file=None,
                                        samples_delimiter=None, blocks_delimiter=None,
                                        samples_header_line=1, blocks_header_line=1,
                                        sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                        sample_domain_col=None,
                                        block_x_col=None, block_y_col=None, block_z_col=None,
                                        block_domain_col=None, block_size=None,
                                        sample_filters=None, block_filters=None,
                                        progress_callback=None):
    if not samples_file or not os.path.isfile(samples_file):
        raise ValueError('Please select a valid samples file.')
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')

    domain_column_name = str(block_domain_col or '').strip()
    if not domain_column_name or domain_column_name == '(None)':
        raise ValueError('Please select a domain column in "Blocks Columns" first.')

    if block_size is None:
        raise ValueError('Block size must be specified for sample domaining.')

    blocks_delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
    output_file = resolve_domain_samples_export_path(
        output_file,
        samples_file=samples_file,
        domain_col=domain_column_name,
    )
    _emit_progress(progress_callback, 0, 100, 'Preparing sample domaining export...')

    print(f"Loading block domain mapping from {blocks_file}...")
    block_metadata = load_large_blocks_metadata(
        blocks_file,
        blocks_delimiter,
        blocks_header_line or 1,
        block_size,
        None,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=domain_column_name,
        block_filters=block_filters,
        config=None,
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 45, 'Loading block domain mapping...'),
    )

    print(f"Loading samples from {samples_file}...")
    df_samples, sample_delimiter = load_full_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        sample_filters=sample_filters,
        progress_label='Reading sample file',
        progress_callback=_make_scaled_progress_callback(progress_callback, 45, 80, 'Reading sample file...'),
        preserve_empty_columns=True,
    )

    sample_coordinate_columns = list(df_samples.columns)
    if not (sample_x_col and sample_y_col and sample_z_col):
        empty_sample_columns = set(_list_all_empty_columns(df_samples))
        sample_coordinate_columns = [
            column_name for column_name in sample_coordinate_columns
            if column_name not in empty_sample_columns
        ]

    sample_x_col, sample_y_col, sample_z_col = resolve_sample_coordinate_columns(
        sample_coordinate_columns,
        sample_x_col=sample_x_col,
        sample_y_col=sample_y_col,
        sample_z_col=sample_z_col,
    )

    coord_frame = df_samples[[sample_x_col, sample_y_col, sample_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1)
    valid_coords = coord_frame.loc[valid_mask].to_numpy(copy=False)

    if block_metadata.get('is_rotated') and len(valid_coords) > 0:
        rotation_center = block_metadata['rotation_center']
        rotation_matrix = block_metadata['rotation_matrix']
        valid_coords = (valid_coords - rotation_center) @ rotation_matrix.T

    all_min_bounds = np.asarray(block_metadata.get('grid_index_origin', block_metadata['all_min_bounds']), dtype=float)
    unified_dims = np.asarray(block_metadata['unified_dims'], dtype=float)
    block_indices = np.floor((valid_coords - all_min_bounds) / unified_dims + 1e-6).astype(int)
    domain_mapping = block_metadata['domain_mapping']

    assigned_domains = []
    matched_count = 0
    total_indices = len(block_indices)
    assign_started_at = time.perf_counter()
    next_terminal_report_at = assign_started_at
    _emit_progress(progress_callback, 80, 100, 'Assigning samples to block domains...')
    for index, idx in enumerate(block_indices, start=1):
        block_idx = (int(idx[0]), int(idx[1]), int(idx[2]))
        domain_value = domain_mapping.get(block_idx, '')
        if domain_value != '':
            matched_count += 1
        assigned_domains.append(domain_value)
        if progress_callback and (index == total_indices or index % 50_000 == 0):
            progress_value = 80 + int(round((index / max(total_indices, 1)) * 15))
            _emit_progress(progress_callback, progress_value, 100, 'Assigning samples to block domains...')
        now = time.perf_counter()
        if now >= next_terminal_report_at:
            elapsed_seconds = max(int(now - assign_started_at), 0)
            percent = int(round((index / max(total_indices, 1)) * 100)) if total_indices > 0 else 100
            print(
                f"Sample domaining progress: {index:,}/{total_indices:,} valid samples ({percent}%) assigned; "
                f"matched={matched_count:,}; elapsed~{elapsed_seconds}s"
            )
            next_terminal_report_at = now + 5.0

    elapsed_seconds = max(int(time.perf_counter() - assign_started_at), 0)
    print(
        f"Sample domaining complete: {total_indices:,} valid samples processed; "
        f"matched={matched_count:,}; elapsed~{elapsed_seconds}s"
    )

    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)

    output_df = df_samples.copy()
    domain_series = pd.Series([''] * len(output_df), index=output_df.index, dtype=object)
    domain_series.loc[valid_mask] = assigned_domains
    output_df[domain_column_name] = domain_series
    _emit_progress(progress_callback, 98, 100, 'Writing domained samples export...')
    output_df.to_csv(output_file, index=False, sep=sample_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Sample domaining export complete.')

    invalid_coordinate_count = int((~valid_mask).sum())
    unmatched_count = int(len(output_df) - matched_count)
    print(f"Exported {len(output_df):,} domained samples to {output_file}")
    print(
        f"Matched {matched_count:,} samples to block domains; "
        f"{unmatched_count:,} unmatched ({invalid_coordinate_count:,} invalid coordinates)."
    )

    return {
        'output_file': output_file,
        'total_samples': int(len(output_df)),
        'matched_samples': int(matched_count),
        'unmatched_samples': int(unmatched_count),
        'invalid_coordinate_samples': int(invalid_coordinate_count),
        'domain_column': domain_column_name,
    }


def export_samples_with_block_values_from_blocks(samples_file, blocks_file, output_file=None,
                                                 samples_delimiter=None, blocks_delimiter=None,
                                                 samples_header_line=1, blocks_header_line=1,
                                                 sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                                 block_x_col=None, block_y_col=None, block_z_col=None,
                                                 block_value_cols=None, block_size=None,
                                                 sample_filters=None, block_filters=None,
                                                 progress_callback=None):
    if not samples_file or not os.path.isfile(samples_file):
        raise ValueError('Please select a valid samples file.')
    if not blocks_file or not os.path.isfile(blocks_file):
        raise ValueError('Please select a valid blocks file.')
    if block_size is None:
        raise ValueError('Block size must be specified for block-to-sample value transfer.')

    selected_columns = _normalize_block_transfer_columns(
        block_value_cols,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
    )
    blocks_delimiter = blocks_delimiter or detect_csv_delimiter(blocks_file)
    output_file = resolve_block_value_transfer_export_path(
        output_file,
        samples_file=samples_file,
        block_value_cols=selected_columns,
    )
    _emit_progress(progress_callback, 0, 100, 'Preparing block value transfer export...')

    print(f"Loading block geometry from {blocks_file}...")
    block_metadata = load_large_blocks_metadata(
        blocks_file,
        blocks_delimiter,
        blocks_header_line or 1,
        block_size,
        None,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=None,
        block_filters=block_filters,
        config=None,
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 45, 'Loading block geometry...'),
    )

    print(f"Building block value mapping from {blocks_file}...")
    block_value_mapping_data = load_block_value_mappings(
        blocks_file,
        blocks_delimiter,
        blocks_header_line or 1,
        block_size,
        selected_columns,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_filters=block_filters,
        block_metadata=block_metadata,
        progress_callback=_make_scaled_progress_callback(progress_callback, 45, 75, 'Building block value mapping...'),
    )
    value_mappings = block_value_mapping_data['value_mappings']

    print(f"Loading samples from {samples_file}...")
    df_samples, sample_delimiter = load_full_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        sample_filters=sample_filters,
        progress_label='Reading sample file',
        progress_callback=_make_scaled_progress_callback(progress_callback, 75, 90, 'Reading sample file...'),
    )

    sample_x_col, sample_y_col, sample_z_col = resolve_sample_coordinate_columns(
        list(df_samples.columns),
        sample_x_col=sample_x_col,
        sample_y_col=sample_y_col,
        sample_z_col=sample_z_col,
    )

    coord_frame = df_samples[[sample_x_col, sample_y_col, sample_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1)
    valid_coords = coord_frame.loc[valid_mask].to_numpy(copy=False)
    block_indices = _compute_sample_block_indices_from_metadata(valid_coords, block_metadata)

    assigned_values = {column_name: [] for column_name in selected_columns}
    matched_count = 0
    total_indices = len(block_indices)
    assign_started_at = time.perf_counter()
    next_terminal_report_at = assign_started_at
    _emit_progress(progress_callback, 90, 100, 'Assigning block values to samples...')
    for index, idx in enumerate(block_indices, start=1):
        block_idx = (int(idx[0]), int(idx[1]), int(idx[2]))
        block_values = value_mappings.get(block_idx)
        if block_values:
            matched_count += 1
        for column_name in selected_columns:
            assigned_values[column_name].append(block_values.get(column_name, '') if block_values else '')
        if progress_callback and (index == total_indices or index % 50_000 == 0):
            progress_value = 90 + int(round((index / max(total_indices, 1)) * 8))
            _emit_progress(progress_callback, progress_value, 100, 'Assigning block values to samples...')
        now = time.perf_counter()
        if now >= next_terminal_report_at:
            elapsed_seconds = max(int(now - assign_started_at), 0)
            percent = int(round((index / max(total_indices, 1)) * 100)) if total_indices > 0 else 100
            print(
                f"Block transfer progress: {index:,}/{total_indices:,} valid samples ({percent}%) assigned; "
                f"matched={matched_count:,}; elapsed~{elapsed_seconds}s"
            )
            next_terminal_report_at = now + 5.0

    elapsed_seconds = max(int(time.perf_counter() - assign_started_at), 0)
    print(
        f"Block transfer complete: {total_indices:,} valid samples processed; "
        f"matched={matched_count:,}; elapsed~{elapsed_seconds}s"
    )

    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)

    output_df = df_samples.copy()
    for column_name in selected_columns:
        output_series = pd.Series([''] * len(output_df), index=output_df.index, dtype=object)
        output_series.loc[valid_mask] = assigned_values[column_name]
        output_df[column_name] = output_series
    _emit_progress(progress_callback, 99, 100, 'Writing block value transfer export...')
    output_df.to_csv(output_file, index=False, sep=sample_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Block value transfer export complete.')

    invalid_coordinate_count = int((~valid_mask).sum())
    unmatched_count = int(len(output_df) - matched_count)
    print(f"Exported {len(output_df):,} samples with transferred block values to {output_file}")
    print(
        f"Matched {matched_count:,} samples to block values; "
        f"{unmatched_count:,} unmatched ({invalid_coordinate_count:,} invalid coordinates)."
    )

    return {
        'output_file': output_file,
        'total_samples': int(len(output_df)),
        'matched_samples': int(matched_count),
        'unmatched_samples': int(unmatched_count),
        'invalid_coordinate_samples': int(invalid_coordinate_count),
        'transferred_columns': list(selected_columns),
        'column_modes': dict(block_value_mapping_data['column_modes']),
    }


def _normalize_optional_size_columns(size_columns):
    values = [str(value or '').strip() for value in (size_columns or ())]
    if len(values) != 3 or any(value in {'', '(None)', '(Infer)'} for value in values):
        return None
    return tuple(values)


def _infer_block_row_bounds(coords, base_block_size, progress_callback=None, progress_label='Resolving block geometry'):
    coords = np.asarray(coords, dtype=float)
    base_dims = np.asarray(base_block_size, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError('Block coordinates must be an Nx3 array.')
    if base_dims.shape != (3,) or np.any(~np.isfinite(base_dims)) or np.any(base_dims <= 0):
        raise ValueError('Base block size must contain three positive values.')
    if len(coords) == 0:
        empty = np.empty((0, 3), dtype=float)
        return empty, empty.copy(), empty.copy()

    grid_origin = np.floor(coords.min(axis=0) / base_dims) * base_dims
    parent_indices = np.floor((coords - grid_origin) / base_dims + 1e-6).astype(np.int64)
    parent_origins = grid_origin + parent_indices * base_dims
    local_coords = np.clip(coords - parent_origins, 0.0, base_dims)
    lower_bounds = np.empty_like(coords, dtype=float)
    upper_bounds = np.empty_like(coords, dtype=float)
    _, inverse = np.unique(parent_indices, axis=0, return_inverse=True)
    order = np.argsort(inverse, kind='stable')
    sorted_inverse = inverse[order]
    starts = np.flatnonzero(np.r_[True, sorted_inverse[1:] != sorted_inverse[:-1]])
    ends = np.r_[starts[1:], len(order)]
    total_groups = len(starts)

    for group_number, (start, end) in enumerate(zip(starts, ends), start=1):
        rows = order[start:end]
        group_local = local_coords[rows]
        group_origin = parent_origins[rows[0]]
        for axis in range(3):
            tolerance = max(base_dims[axis] * 1e-7, 1e-9)
            centers = _cluster_axis_centers(group_local[:, axis], tolerance)
            boundaries = np.empty(len(centers) + 1, dtype=float)
            boundaries[0], boundaries[-1] = 0.0, base_dims[axis]
            if len(centers) > 1:
                boundaries[1:-1] = (centers[:-1] + centers[1:]) / 2.0
            center_indices = np.abs(group_local[:, axis, None] - centers[None, :]).argmin(axis=1)
            lower_bounds[rows, axis] = group_origin[axis] + boundaries[center_indices]
            upper_bounds[rows, axis] = group_origin[axis] + boundaries[center_indices + 1]
        if progress_callback and (group_number == total_groups or group_number % 500 == 0):
            progress_callback(group_number, total_groups, progress_label)
    return lower_bounds, upper_bounds, upper_bounds - lower_bounds


def _resolve_block_row_geometry(df, coordinate_columns, base_block_size, size_columns=None,
                                progress_callback=None, progress_label='Resolving block geometry'):
    if len(coordinate_columns) != 3 or any(column not in df.columns for column in coordinate_columns):
        raise ValueError('Three valid block coordinate columns are required.')
    coord_frame = df[list(coordinate_columns)].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1).to_numpy()
    coords = coord_frame.loc[valid_mask].to_numpy(dtype=float, copy=False)
    explicit_size_columns = _normalize_optional_size_columns(size_columns)
    if explicit_size_columns:
        missing = [column for column in explicit_size_columns if column not in df.columns]
        if missing:
            raise ValueError(f"Block size column(s) not found: {', '.join(missing)}")
        sizes = df.loc[valid_mask, list(explicit_size_columns)].apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)
        valid_sizes = np.isfinite(sizes).all(axis=1) & (sizes > 0).all(axis=1)
        valid_positions = np.flatnonzero(valid_mask)
        valid_mask[valid_positions[~valid_sizes]] = False
        coords, sizes = coords[valid_sizes], sizes[valid_sizes]
        lower_bounds, upper_bounds = coords - sizes / 2.0, coords + sizes / 2.0
        mode = 'explicit-size-columns'
        if progress_callback:
            progress_callback(1, 1, progress_label)
    else:
        lower_bounds, upper_bounds, sizes = _infer_block_row_bounds(
            coords,
            base_block_size,
            progress_callback=progress_callback,
            progress_label=progress_label,
        )
        mode = 'inferred-from-base-grid'
    return {
        'row_indices': np.flatnonzero(valid_mask), 'centers': coords,
        'lower_bounds': lower_bounds, 'upper_bounds': upper_bounds,
        'sizes': sizes, 'mode': mode,
    }


def _prepare_exact_block_transfer_keys(df, coordinate_columns, size_columns=None):
    if len(coordinate_columns) != 3 or any(column not in df.columns for column in coordinate_columns):
        raise ValueError('Three valid block coordinate columns are required.')
    coord_frame = df[list(coordinate_columns)].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1).to_numpy()
    key_columns = ['__x__', '__y__', '__z__']
    key_frame = coord_frame.loc[valid_mask].copy()
    key_frame.columns = key_columns

    explicit_size_columns = _normalize_optional_size_columns(size_columns)
    if explicit_size_columns:
        missing = [column for column in explicit_size_columns if column not in df.columns]
        if missing:
            raise ValueError(f"Block size column(s) not found: {', '.join(missing)}")
        size_frame = df.loc[valid_mask, list(explicit_size_columns)].apply(pd.to_numeric, errors='coerce')
        valid_sizes = size_frame.notna().all(axis=1).to_numpy() & (size_frame.to_numpy(dtype=float, copy=False) > 0).all(axis=1)
        valid_positions = np.flatnonzero(valid_mask)
        valid_mask[valid_positions[~valid_sizes]] = False
        key_frame = key_frame.loc[valid_sizes].copy()
        size_frame = size_frame.loc[valid_sizes].copy()
        size_frame.columns = ['__dx__', '__dy__', '__dz__']
        key_frame = pd.concat([key_frame, size_frame], axis=1)
        key_columns.extend(['__dx__', '__dy__', '__dz__'])

    if len(key_frame):
        key_frame = key_frame.round(9)
    key_frame['__row_index__'] = np.flatnonzero(valid_mask)
    return {
        'row_indices': np.flatnonzero(valid_mask),
        'key_columns': key_columns,
        'key_frame': key_frame,
        'uses_explicit_sizes': bool(explicit_size_columns),
    }


def _try_exact_block_model_transfer(source_df, target_df, selected_columns,
                                    source_coordinate_columns, target_coordinate_columns,
                                    source_block_size, target_block_size,
                                    source_size_cols=None, target_size_cols=None,
                                    progress_callback=None):
    source_key_data = _prepare_exact_block_transfer_keys(source_df, source_coordinate_columns, source_size_cols)
    target_key_data = _prepare_exact_block_transfer_keys(target_df, target_coordinate_columns, target_size_cols)

    if source_key_data['uses_explicit_sizes'] != target_key_data['uses_explicit_sizes']:
        return None
    if not source_key_data['uses_explicit_sizes'] and not np.allclose(
        np.asarray(source_block_size, dtype=float),
        np.asarray(target_block_size, dtype=float),
        rtol=1e-9,
        atol=1e-9,
    ):
        return None

    source_keys = source_key_data['key_frame']
    target_keys = target_key_data['key_frame']
    key_columns = list(source_key_data['key_columns'])

    source_duplicates = source_keys.duplicated(subset=key_columns, keep=False)
    target_duplicates = target_keys.duplicated(subset=key_columns, keep=False)

    if progress_callback:
        _emit_progress(progress_callback, 46, 100, 'Matching exact blocks on common grid...')

    source_lookup = source_keys.loc[~source_duplicates, key_columns + ['__row_index__']].merge(
        source_df[selected_columns],
        left_on='__row_index__',
        right_index=True,
        how='left',
        sort=False,
    )
    source_lookup.rename(columns={'__row_index__': '__source_row_index__'}, inplace=True)
    merged = target_keys.loc[~target_duplicates, key_columns + ['__row_index__']].merge(
        source_lookup,
        on=key_columns,
        how='left',
        sort=False,
    )
    if progress_callback:
        _emit_progress(progress_callback, 56, 100, 'Matching exact blocks on common grid...')

    matched_mask = merged['__source_row_index__'].notna()
    source_values_df = source_df.iloc[source_key_data['row_indices']]
    column_modes = _detect_dataframe_transfer_column_modes(source_values_df, selected_columns)
    matched_target_row_indices = merged.loc[matched_mask, '__row_index__'].to_numpy(dtype=int, copy=False)
    matched_source_row_indices = merged.loc[matched_mask, '__source_row_index__'].to_numpy(dtype=np.int64, copy=False)
    valid_target_row_indices = target_key_data['row_indices']
    remaining_target_row_indices = valid_target_row_indices[
        ~np.isin(valid_target_row_indices, matched_target_row_indices, assume_unique=False)
    ]
    invalid_targets = len(target_df) - len(valid_target_row_indices)
    return {
        'column_modes': column_modes,
        'matched_target_row_indices': matched_target_row_indices,
        'matched_source_row_indices': matched_source_row_indices,
        'remaining_target_row_indices': remaining_target_row_indices,
        'total_target_blocks': int(len(target_df)),
        'overlap_matched_blocks': int(len(matched_target_row_indices)),
        'nearest_matched_blocks': 0,
        'unmatched_blocks': int(len(target_df) - len(matched_target_row_indices)),
        'invalid_target_blocks': int(invalid_targets),
        'source_geometry_mode': 'exact-grid',
        'target_geometry_mode': 'exact-grid',
        'all_valid_targets_matched': len(remaining_target_row_indices) == 0,
    }


def _detect_dataframe_transfer_column_modes(df, columns):
    modes = {}
    for column_name in columns:
        values = df[column_name]
        nonblank = values.loc[
            values.notna() & values.astype(str).str.strip().ne('') & values.astype(str).str.lower().ne('nan')
        ]
        modes[column_name] = (
            'numeric' if len(nonblank) and pd.to_numeric(nonblank, errors='coerce').notna().all() else 'categorical'
        )
    return modes


def _restore_list_widget_selection(list_widget, values):
    if list_widget is None:
        return
    desired = {str(value or '').strip() for value in (values or []) if str(value or '').strip()}
    blocker = QtCore.QSignalBlocker(list_widget)
    try:
        for index in range(list_widget.count()):
            item = list_widget.item(index)
            item.setSelected(item.text() in desired)
    finally:
        del blocker


def _normalize_table_attribute_columns(columns, selection_label):
    if isinstance(columns, str):
        raw_columns = [part.strip() for part in columns.split(',')]
    else:
        raw_columns = [str(part).strip() for part in (columns or [])]

    normalized_columns = []
    seen = set()
    for column_name in raw_columns:
        if not column_name or column_name == '(None)' or column_name in seen:
            continue
        normalized_columns.append(column_name)
        seen.add(column_name)

    if not normalized_columns:
        raise ValueError(f'Please select at least one {selection_label}.')

    return normalized_columns


def _normalize_table_attribute_key_series(series):
    normalized_values = []
    for value in series:
        if pd.isna(value):
            normalized_values.append(None)
            continue
        if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
            normalized_values.append(str(int(value)))
            continue
        if isinstance(value, (np.floating, float)):
            float_value = float(value)
            if not np.isfinite(float_value):
                normalized_values.append(None)
            elif float_value.is_integer():
                normalized_values.append(str(int(float_value)))
            else:
                normalized_values.append(format(float_value, '.15g'))
            continue

        text = str(value).strip()
        normalized_values.append(text if text and text.lower() != 'nan' else None)

    return pd.Series(normalized_values, index=series.index, dtype=object)


def export_block_model_with_table_attributes(block_model_file, table_file, output_file=None,
                                             block_model_delimiter=None, block_model_header_line=1,
                                             table_delimiter=None, table_header_line=1,
                                             key_columns=None, table_value_cols=None,
                                             progress_callback=None):
    if not block_model_file or not str(block_model_file).strip():
        raise ValueError('Block model file is required.')
    if not table_file or not str(table_file).strip():
        raise ValueError('Attribute table file is required.')

    block_model_file = str(block_model_file).strip()
    table_file = str(table_file).strip()
    if not os.path.isfile(block_model_file):
        raise ValueError(f'Block model file not found: {block_model_file}')
    if not os.path.isfile(table_file):
        raise ValueError(f'Attribute table file not found: {table_file}')
    if is_bmf_file(table_file):
        raise ValueError('Attribute table must be provided in CSV format.')

    key_columns = _normalize_table_attribute_columns(key_columns, 'match key column')
    table_value_cols = _normalize_table_attribute_columns(table_value_cols, 'attribute column to assign')
    key_column_set = set(key_columns)
    overlapping_columns = [column_name for column_name in table_value_cols if column_name in key_column_set]
    if overlapping_columns:
        raise ValueError(
            'Attribute columns cannot also be selected as match keys: '
            + ', '.join(overlapping_columns)
        )

    output_file = resolve_block_model_table_attribute_export_path(
        output_file,
        block_model_file=block_model_file,
        table_file=table_file,
    )

    _emit_progress(progress_callback, 0, 100, 'Preparing table attribute assignment...')
    block_df, resolved_block_delimiter = load_full_blocks_dataframe(
        block_model_file,
        blocks_delimiter=block_model_delimiter,
        blocks_header_line=block_model_header_line,
        progress_label='Reading block model',
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 45, 'Reading block model...'),
    )
    missing_block_keys = [column_name for column_name in key_columns if column_name not in block_df.columns]
    if missing_block_keys:
        raise ValueError(
            'Selected key columns were not found in the block model: '
            + ', '.join(missing_block_keys)
        )

    table_delimiter = table_delimiter or detect_csv_delimiter(table_file)
    selected_table_columns = list(dict.fromkeys(key_columns + table_value_cols))
    table_df, _ = read_selected_columns_with_header(
        table_file,
        table_delimiter,
        table_header_line,
        selected_table_columns,
        progress_label='Reading attribute table',
        progress_callback=_make_scaled_progress_callback(progress_callback, 45, 75, 'Reading attribute table...'),
    )
    missing_table_keys = [column_name for column_name in key_columns if column_name not in table_df.columns]
    if missing_table_keys:
        raise ValueError(
            'Selected key columns were not found in the attribute table: '
            + ', '.join(missing_table_keys)
        )

    block_work_df = block_df.copy()
    table_work_df = table_df.copy()
    helper_key_columns = []
    for key_index, column_name in enumerate(key_columns):
        helper_column = f'__table_attr_key_{key_index}'
        block_work_df[helper_column] = _normalize_table_attribute_key_series(block_work_df[column_name])
        table_work_df[helper_column] = _normalize_table_attribute_key_series(table_work_df[column_name])
        helper_key_columns.append(helper_column)

    duplicate_key_mask = table_work_df.duplicated(subset=helper_key_columns, keep=False)
    if duplicate_key_mask.any():
        duplicate_preview = table_df.loc[duplicate_key_mask, key_columns].drop_duplicates().head(5)
        preview_parts = []
        for _, row in duplicate_preview.iterrows():
            preview_parts.append(', '.join(f'{column_name}={row[column_name]}' for column_name in key_columns))
        preview_text = '; '.join(preview_parts)
        raise ValueError(
            'Attribute table contains duplicate key combinations for the selected keys'
            + (f': {preview_text}' if preview_text else '.')
        )

    _emit_progress(progress_callback, 76, 100, 'Matching table attributes...')
    helper_value_columns = {}
    table_lookup_df = table_work_df[helper_key_columns].copy()
    for value_index, column_name in enumerate(table_value_cols):
        helper_column = f'__table_attr_value_{value_index}'
        table_lookup_df[helper_column] = table_work_df[column_name]
        helper_value_columns[column_name] = helper_column

    merged_df = block_work_df.merge(
        table_lookup_df,
        how='left',
        on=helper_key_columns,
        sort=False,
        indicator='__table_attr_match',
    )
    matched_mask = merged_df['__table_attr_match'].eq('both')
    output_df = merged_df.drop(columns=helper_key_columns + ['__table_attr_match'])

    for column_name, helper_column in helper_value_columns.items():
        if column_name in block_df.columns:
            output_df.loc[matched_mask, column_name] = output_df.loc[matched_mask, helper_column]
            output_df = output_df.drop(columns=[helper_column])
        else:
            output_df[column_name] = output_df[helper_column]
            output_df = output_df.drop(columns=[helper_column])

    _emit_progress(progress_callback, 95, 100, 'Writing table attribute assignment...')
    output_df.to_csv(output_file, index=False)
    _emit_progress(progress_callback, 100, 100, 'Table attribute assignment complete.')

    return {
        'output_file': output_file,
        'block_model_file': block_model_file,
        'table_file': table_file,
        'block_model_delimiter': resolved_block_delimiter,
        'table_delimiter': table_delimiter,
        'key_columns': list(key_columns),
        'assigned_columns': list(table_value_cols),
        'matched_rows': int(matched_mask.sum()),
        'unmatched_rows': int((~matched_mask).sum()),
        'total_rows': int(len(output_df)),
    }


BLOCK_MODEL_TRANSFER_TARGET_CHUNK_SIZE = 100_000


def _prepare_source_block_transfer_dataframe(source_blocks_file, selected_columns,
                                             source_delimiter=None, source_header_line=1,
                                             source_x_col=None, source_y_col=None, source_z_col=None,
                                             source_size_cols=None, progress_callback=None):
    if is_bmf_file(source_blocks_file):
        source_df, _ = load_full_blocks_dataframe(
            source_blocks_file,
            source_delimiter,
            source_header_line,
            progress_label='Reading source block model',
            progress_callback=_make_scaled_progress_callback(progress_callback, 0, 25, 'Reading source block model...'),
        )
        source_x_col, source_y_col, source_z_col = resolve_block_coordinate_columns(
            list(source_df.columns),
            source_x_col,
            source_y_col,
            source_z_col,
        )
        return source_df, source_x_col, source_y_col, source_z_col

    delimiter = source_delimiter or detect_csv_delimiter(source_blocks_file)
    source_columns = parse_effective_header_line(source_blocks_file, delimiter, source_header_line)
    source_x_col, source_y_col, source_z_col = resolve_block_coordinate_columns(
        list(source_columns),
        source_x_col,
        source_y_col,
        source_z_col,
    )
    explicit_size_columns = _normalize_optional_size_columns(source_size_cols) or ()
    selected_source_columns = list(dict.fromkeys([
        source_x_col,
        source_y_col,
        source_z_col,
        *explicit_size_columns,
        *selected_columns,
    ]))
    source_df, _ = read_selected_columns_with_header(
        source_blocks_file,
        delimiter,
        source_header_line,
        selected_source_columns,
        progress_label='Reading source block model',
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 25, 'Reading source block model...'),
    )
    return source_df, source_x_col, source_y_col, source_z_col


def _prepare_exact_source_block_transfer_lookup(source_df, selected_columns,
                                                source_coordinate_columns,
                                                source_block_size, target_block_size,
                                                source_size_cols=None, target_size_cols=None):
    source_key_data = _prepare_exact_block_transfer_keys(source_df, source_coordinate_columns, source_size_cols)
    target_uses_explicit_sizes = bool(_normalize_optional_size_columns(target_size_cols))
    if source_key_data['uses_explicit_sizes'] != target_uses_explicit_sizes:
        return None
    if not source_key_data['uses_explicit_sizes'] and not np.allclose(
        np.asarray(source_block_size, dtype=float),
        np.asarray(target_block_size, dtype=float),
        rtol=1e-9,
        atol=1e-9,
    ):
        return None

    key_columns = list(source_key_data['key_columns'])
    source_keys = source_key_data['key_frame']
    source_duplicates = source_keys.duplicated(subset=key_columns, keep=False)
    source_lookup = source_keys.loc[~source_duplicates, key_columns + ['__row_index__']].merge(
        source_df[selected_columns],
        left_on='__row_index__',
        right_index=True,
        how='left',
        sort=False,
    )
    source_lookup.rename(columns={'__row_index__': '__source_row_index__'}, inplace=True)
    return {
        'key_columns': key_columns,
        'source_lookup': source_lookup,
        'row_indices': source_key_data['row_indices'],
    }


def _match_exact_block_transfer_chunk(target_chunk_df, exact_source_lookup,
                                      target_coordinate_columns, target_size_cols=None):
    if exact_source_lookup is None:
        return None
    target_key_data = _prepare_exact_block_transfer_keys(target_chunk_df, target_coordinate_columns, target_size_cols)
    key_columns = list(exact_source_lookup['key_columns'])
    target_keys = target_key_data['key_frame']
    target_duplicates = target_keys.duplicated(subset=key_columns, keep=False)
    matched_rows = target_keys.loc[~target_duplicates, key_columns + ['__row_index__']].merge(
        exact_source_lookup['source_lookup'],
        on=key_columns,
        how='left',
        sort=False,
    )
    matched_mask = matched_rows['__source_row_index__'].notna()
    matched_rows = matched_rows.loc[matched_mask].copy()
    matched_target_row_positions = matched_rows['__row_index__'].to_numpy(dtype=int, copy=False)
    valid_target_row_positions = target_key_data['row_indices']
    remaining_target_row_positions = valid_target_row_positions[
        ~np.isin(valid_target_row_positions, matched_target_row_positions, assume_unique=False)
    ]
    return {
        'matched_rows': matched_rows,
        'matched_target_row_positions': matched_target_row_positions,
        'remaining_target_row_positions': remaining_target_row_positions,
        'valid_target_row_positions': valid_target_row_positions,
    }


def _iter_block_model_transfer_target_chunks(target_blocks_file, target_delimiter=None, target_header_line=1,
                                             progress_callback=None, chunksize=BLOCK_MODEL_TRANSFER_TARGET_CHUNK_SIZE):
    delimiter = target_delimiter or detect_csv_delimiter(target_blocks_file)
    effective_header_line = resolve_effective_csv_header_line(target_blocks_file, target_header_line)
    target_columns = build_unique_column_names(parse_header_line(target_blocks_file, delimiter, effective_header_line))
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=target_columns,
        skiprows=effective_header_line,
        comment='#',
        chunksize=chunksize,
    )
    return delimiter, target_columns, iterate_csv_path_chunks_with_progress(
        target_blocks_file,
        'Reading target block model',
        progress_callback=progress_callback,
        header_line=effective_header_line,
        **read_kwargs,
    )


def _export_blocks_with_source_block_values_streaming_target(source_df, target_blocks_file, output_file,
                                                             selected_columns,
                                                             source_x_col=None, source_y_col=None, source_z_col=None,
                                                             target_delimiter=None, target_header_line=1,
                                                             target_x_col=None, target_y_col=None, target_z_col=None,
                                                             source_block_size=None, target_block_size=None,
                                                             source_size_cols=None, target_size_cols=None,
                                                             nearest_fallback=True, nearest_distance_limit=None,
                                                             progress_callback=None):
    target_delimiter, target_columns, target_chunks = _iter_block_model_transfer_target_chunks(
        target_blocks_file,
        target_delimiter=target_delimiter,
        target_header_line=target_header_line,
        progress_callback=_make_scaled_progress_callback(progress_callback, 25, 98, 'Reading target block model...'),
    )
    target_x_col, target_y_col, target_z_col = resolve_block_coordinate_columns(
        list(target_columns),
        target_x_col,
        target_y_col,
        target_z_col,
    )
    exact_source_lookup = _prepare_exact_source_block_transfer_lookup(
        source_df,
        selected_columns,
        (source_x_col, source_y_col, source_z_col),
        source_block_size,
        target_block_size,
        source_size_cols=source_size_cols,
        target_size_cols=target_size_cols,
    )

    if exact_source_lookup is not None:
        source_values_df = source_df.iloc[exact_source_lookup['row_indices']]
        column_modes = _detect_dataframe_transfer_column_modes(source_values_df, selected_columns)
    else:
        column_modes = None

    source_geometry = None
    values_by_column = None
    tree = None
    max_source_half_diagonal = None
    wrote_header = False
    total_target_blocks = 0
    overlap_matches = 0
    nearest_matches = 0
    invalid_target_blocks = 0
    target_geometry_mode = None

    def ensure_source_matching_context(use_progress=False):
        nonlocal source_geometry, values_by_column, tree, max_source_half_diagonal, column_modes
        if source_geometry is not None:
            return
        if use_progress:
            _emit_progress(progress_callback, 46, 100, 'Resolving source block geometry...')
        source_geometry_kwargs = {}
        if use_progress:
            source_geometry_kwargs = {
                'progress_callback': _make_scaled_progress_callback(progress_callback, 46, 56, 'Resolving source block geometry...'),
                'progress_label': 'Resolving source block geometry',
            }
        source_geometry = _resolve_block_row_geometry(
            source_df,
            (source_x_col, source_y_col, source_z_col),
            source_block_size,
            source_size_cols,
            **source_geometry_kwargs,
        )
        if not len(source_geometry['centers']):
            raise ValueError('The source block model has no rows with valid coordinates and dimensions.')

        from scipy.spatial import cKDTree

        source_values_df_local = source_df.iloc[source_geometry['row_indices']]
        if column_modes is None:
            column_modes = _detect_dataframe_transfer_column_modes(source_values_df_local, selected_columns)
        values_by_column = {
            column: (
                pd.to_numeric(source_values_df_local[column], errors='coerce').to_numpy(dtype=float)
                if mode == 'numeric' else
                source_values_df_local[column].fillna('').astype(str).str.strip().to_numpy(dtype=object)
            )
            for column, mode in column_modes.items()
        }
        tree = cKDTree(source_geometry['centers'])
        max_source_half_diagonal = float(np.linalg.norm(source_geometry['sizes'], axis=1).max() / 2.0)

    if exact_source_lookup is None:
        ensure_source_matching_context(use_progress=True)

    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)

    for chunk_number, target_chunk in enumerate(target_chunks, start=1):
        total_target_blocks += len(target_chunk)
        print(
            f"Target block chunk {chunk_number:,}: read {len(target_chunk):,} rows "
            f"({total_target_blocks:,} cumulative) from {os.path.basename(target_blocks_file)}."
        )
        exact_chunk = _match_exact_block_transfer_chunk(
            target_chunk,
            exact_source_lookup,
            (target_x_col, target_y_col, target_z_col),
            target_size_cols=target_size_cols,
        ) if exact_source_lookup is not None else None

        if exact_chunk is not None:
            invalid_target_blocks += int(len(target_chunk) - len(exact_chunk['valid_target_row_positions']))
            remaining_target_rows = exact_chunk['remaining_target_row_positions']
        else:
            remaining_target_rows = None

        if column_modes is None:
            ensure_source_matching_context(use_progress=(chunk_number == 1))

        assigned = {
            column: (
                np.full(len(target_chunk), np.nan, dtype=float)
                if mode == 'numeric' else
                np.full(len(target_chunk), '', dtype=object)
            )
            for column, mode in column_modes.items()
        }

        if exact_chunk is not None and len(exact_chunk['matched_target_row_positions']):
            overlap_matches += int(len(exact_chunk['matched_target_row_positions']))
            for column, mode in column_modes.items():
                values = exact_chunk['matched_rows'][column]
                if mode == 'numeric':
                    assigned[column][exact_chunk['matched_target_row_positions']] = pd.to_numeric(
                        values,
                        errors='coerce',
                    ).to_numpy(dtype=float)
                else:
                    normalized = values.fillna('').astype(str).str.strip()
                    normalized = normalized.mask(normalized.str.lower().eq('nan'), '')
                    assigned[column][exact_chunk['matched_target_row_positions']] = normalized.to_numpy(dtype=object)

        if remaining_target_rows is None or len(remaining_target_rows):
            ensure_source_matching_context(use_progress=(chunk_number == 1 and exact_source_lookup is not None))
            target_geometry_df = target_chunk if remaining_target_rows is None else target_chunk.iloc[remaining_target_rows]
            target_geometry = _resolve_block_row_geometry(
                target_geometry_df,
                (target_x_col, target_y_col, target_z_col),
                target_block_size,
                target_size_cols,
            )
            if exact_chunk is None:
                invalid_target_blocks += int(len(target_chunk) - len(target_geometry['row_indices']))
                target_row_positions = target_geometry['row_indices']
            else:
                target_row_positions = np.asarray(remaining_target_rows, dtype=int)[target_geometry['row_indices']]

            if len(target_geometry['centers']):
                target_geometry_mode = (
                    target_geometry['mode'] if exact_source_lookup is None else f'exact-prematch + {target_geometry["mode"]}'
                )

            for local_index, target_row in enumerate(target_row_positions):
                target_center = target_geometry['centers'][local_index]
                radius = float(np.linalg.norm(target_geometry['sizes'][local_index]) / 2.0 + max_source_half_diagonal)
                candidates = np.asarray(tree.query_ball_point(target_center, radius), dtype=int)
                overlaps = np.empty(0, dtype=int)
                volumes = np.empty(0, dtype=float)
                if len(candidates):
                    overlap_lengths = np.maximum(
                        0.0,
                        np.minimum(source_geometry['upper_bounds'][candidates], target_geometry['upper_bounds'][local_index])
                        - np.maximum(source_geometry['lower_bounds'][candidates], target_geometry['lower_bounds'][local_index]),
                    )
                    candidate_volumes = np.prod(overlap_lengths, axis=1)
                    positive = candidate_volumes > max(float(np.prod(target_geometry['sizes'][local_index])) * 1e-12, 1e-12)
                    overlaps, volumes = candidates[positive], candidate_volumes[positive]
                if len(overlaps):
                    overlap_matches += 1
                    for column, mode in column_modes.items():
                        values = values_by_column[column][overlaps]
                        if mode == 'numeric':
                            valid = np.isfinite(values)
                            if valid.any():
                                assigned[column][target_row] = float(np.average(values[valid], weights=volumes[valid]))
                        else:
                            totals = {}
                            for value, volume in zip(values, volumes):
                                value = str(value).strip()
                                if value and value.lower() != 'nan':
                                    totals[value] = totals.get(value, 0.0) + float(volume)
                            if totals:
                                assigned[column][target_row] = sorted(
                                    totals.items(),
                                    key=lambda item: (-item[1], item[0]),
                                )[0][0]
                elif nearest_fallback:
                    nearest_distance, nearest = tree.query(target_center, k=1)
                    if nearest_distance_limit is None or float(nearest_distance) <= nearest_distance_limit:
                        nearest = int(nearest)
                        nearest_matches += 1
                        for column, mode in column_modes.items():
                            value = values_by_column[column][nearest]
                            if mode == 'numeric' and np.isfinite(value):
                                assigned[column][target_row] = float(value)
                            elif mode != 'numeric' and str(value).strip() and str(value).strip().lower() != 'nan':
                                assigned[column][target_row] = str(value).strip()
        elif target_geometry_mode is None:
            target_geometry_mode = 'exact-grid'

        for column, mode in column_modes.items():
            if mode == 'numeric':
                target_chunk[column] = pd.Series(assigned[column], index=target_chunk.index, dtype=float)
            else:
                target_chunk[column] = pd.Series(assigned[column], index=target_chunk.index, dtype=object)
        target_chunk.to_csv(
            output_file,
            index=False,
            sep=target_delimiter,
            mode='w' if not wrote_header else 'a',
            header=not wrote_header,
        )
        wrote_header = True
        print(
            f"Target block chunk {chunk_number:,}: wrote {len(target_chunk):,} rows; "
            f"matches so far overlap/exact={overlap_matches:,}, nearest={nearest_matches:,}, "
            f"invalid={invalid_target_blocks:,}."
        )

    if not wrote_header:
        empty_columns = list(dict.fromkeys(list(target_columns) + list(selected_columns)))
        pd.DataFrame(columns=empty_columns).to_csv(output_file, index=False, sep=target_delimiter)

    _emit_progress(progress_callback, 99, 100, 'Writing block-model transfer export...')
    _emit_progress(progress_callback, 100, 100, 'Block-model transfer complete.')
    if target_geometry_mode is None:
        target_geometry_mode = 'exact-grid' if exact_source_lookup is not None else 'inferred-from-base-grid'
    source_geometry_mode = 'exact-grid' if source_geometry is None else source_geometry['mode']
    unmatched = total_target_blocks - overlap_matches - nearest_matches
    return {
        'output_file': output_file,
        'total_target_blocks': int(total_target_blocks),
        'overlap_matched_blocks': int(overlap_matches),
        'nearest_matched_blocks': int(nearest_matches),
        'unmatched_blocks': int(unmatched),
        'invalid_target_blocks': int(invalid_target_blocks),
        'transferred_columns': list(selected_columns),
        'column_modes': column_modes or {},
        'source_geometry_mode': source_geometry_mode,
        'target_geometry_mode': target_geometry_mode,
        'max_nearest_distance': nearest_distance_limit,
    }


def export_blocks_with_source_block_values(source_blocks_file, target_blocks_file, output_file=None,
                                            source_delimiter=None, target_delimiter=None,
                                            source_header_line=1, target_header_line=1,
                                            source_x_col=None, source_y_col=None, source_z_col=None,
                                            target_x_col=None, target_y_col=None, target_z_col=None,
                                            source_value_cols=None,
                                            source_block_size=None, target_block_size=None,
                                            source_size_cols=None, target_size_cols=None,
                                            nearest_fallback=True, max_nearest_distance=None,
                                            progress_callback=None):
    """Enrich existing target blocks from overlapping source blocks without creating rows."""
    if not source_blocks_file or not os.path.isfile(source_blocks_file):
        raise ValueError('Please select a valid source blocks file.')
    if not target_blocks_file or not os.path.isfile(target_blocks_file):
        raise ValueError('Please select a valid target blocks file.')
    if source_block_size is None or target_block_size is None:
        raise ValueError('Source and target base block sizes are required.')
    nearest_distance_limit = None
    if max_nearest_distance not in (None, ''):
        nearest_distance_limit = float(max_nearest_distance)
        if not np.isfinite(nearest_distance_limit) or nearest_distance_limit < 0:
            raise ValueError('Maximum nearest fallback distance must be a finite value greater than or equal to zero.')
        if nearest_distance_limit == 0:
            nearest_distance_limit = None
    selected_columns = _normalize_block_transfer_columns(
        source_value_cols, block_x_col=source_x_col, block_y_col=source_y_col, block_z_col=source_z_col,
    )
    output_file = resolve_block_model_transfer_export_path(output_file, target_blocks_file)
    source_df, source_x_col, source_y_col, source_z_col = _prepare_source_block_transfer_dataframe(
        source_blocks_file,
        selected_columns,
        source_delimiter=source_delimiter,
        source_header_line=source_header_line,
        source_x_col=source_x_col,
        source_y_col=source_y_col,
        source_z_col=source_z_col,
        source_size_cols=source_size_cols,
        progress_callback=progress_callback,
    )
    missing = [column for column in selected_columns if column not in source_df.columns]
    if missing:
        raise ValueError(f"Source block column(s) not found: {', '.join(missing)}")
    if not is_bmf_file(target_blocks_file):
        return _export_blocks_with_source_block_values_streaming_target(
            source_df,
            target_blocks_file,
            output_file,
            selected_columns,
            source_x_col=source_x_col,
            source_y_col=source_y_col,
            source_z_col=source_z_col,
            target_delimiter=target_delimiter,
            target_header_line=target_header_line,
            target_x_col=target_x_col,
            target_y_col=target_y_col,
            target_z_col=target_z_col,
            source_block_size=source_block_size,
            target_block_size=target_block_size,
            source_size_cols=source_size_cols,
            target_size_cols=target_size_cols,
            nearest_fallback=nearest_fallback,
            nearest_distance_limit=nearest_distance_limit,
            progress_callback=progress_callback,
        )
    target_df, output_delimiter = load_full_blocks_dataframe(
        target_blocks_file, target_delimiter, target_header_line,
        progress_label='Reading target block model',
        progress_callback=_make_scaled_progress_callback(progress_callback, 25, 45, 'Reading target block model...'),
    )
    target_x_col, target_y_col, target_z_col = resolve_block_coordinate_columns(
        list(target_df.columns), target_x_col, target_y_col, target_z_col,
    )
    exact_transfer = _try_exact_block_model_transfer(
        source_df,
        target_df,
        selected_columns,
        (source_x_col, source_y_col, source_z_col),
        (target_x_col, target_y_col, target_z_col),
        source_block_size,
        target_block_size,
        source_size_cols=source_size_cols,
        target_size_cols=target_size_cols,
        progress_callback=progress_callback,
    )
    if exact_transfer is not None and exact_transfer['all_valid_targets_matched']:
        output_df = target_df.copy()
        for column, mode in exact_transfer['column_modes'].items():
            matched_values = source_df.iloc[exact_transfer['matched_source_row_indices']][column]
            if mode == 'numeric':
                output_df[column] = pd.Series(np.nan, index=output_df.index, dtype=float)
                output_df.loc[exact_transfer['matched_target_row_indices'], column] = pd.to_numeric(
                    matched_values,
                    errors='coerce',
                ).to_numpy(dtype=float)
            else:
                output_df[column] = ''
                values = matched_values.fillna('').astype(str).str.strip()
                values = values.mask(values.str.lower().eq('nan'), '')
                output_df.loc[exact_transfer['matched_target_row_indices'], column] = values.to_numpy(dtype=object)
        os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
        _emit_progress(progress_callback, 99, 100, 'Writing block-model transfer export...')
        output_df.to_csv(output_file, index=False, sep=output_delimiter if output_delimiter != 'bmf' else ',')
        _emit_progress(progress_callback, 100, 100, 'Block-model transfer complete.')
        return {
            'output_file': output_file,
            'total_target_blocks': exact_transfer['total_target_blocks'],
            'overlap_matched_blocks': exact_transfer['overlap_matched_blocks'],
            'nearest_matched_blocks': exact_transfer['nearest_matched_blocks'],
            'unmatched_blocks': exact_transfer['unmatched_blocks'],
            'invalid_target_blocks': exact_transfer['invalid_target_blocks'],
            'transferred_columns': list(selected_columns),
            'column_modes': exact_transfer['column_modes'],
            'source_geometry_mode': exact_transfer['source_geometry_mode'],
            'target_geometry_mode': exact_transfer['target_geometry_mode'],
            'max_nearest_distance': nearest_distance_limit,
        }
    if exact_transfer is not None:
        column_modes = exact_transfer['column_modes']
        assigned = {column: np.full(len(target_df), '', dtype=object) for column in selected_columns}
        overlap_matches = int(exact_transfer['overlap_matched_blocks'])
        exact_target_rows = exact_transfer['matched_target_row_indices']
        exact_source_rows = exact_transfer['matched_source_row_indices']
        for column, mode in column_modes.items():
            matched_values = source_df.iloc[exact_source_rows][column]
            if mode == 'numeric':
                numeric_values = pd.to_numeric(matched_values, errors='coerce').to_numpy(dtype=float)
                valid = np.isfinite(numeric_values)
                assigned[column][exact_target_rows[valid]] = numeric_values[valid]
            else:
                values = matched_values.fillna('').astype(str).str.strip()
                values = values.mask(values.str.lower().eq('nan'), '')
                assigned[column][exact_target_rows] = values.to_numpy(dtype=object)
        remaining_target_rows = exact_transfer['remaining_target_row_indices']
    else:
        column_modes = None
        assigned = None
        overlap_matches = 0
        remaining_target_rows = None
    _emit_progress(progress_callback, 46, 100, 'Resolving source block geometry...')
    source_geometry = _resolve_block_row_geometry(
        source_df,
        (source_x_col, source_y_col, source_z_col),
        source_block_size,
        source_size_cols,
        progress_callback=_make_scaled_progress_callback(progress_callback, 46, 56, 'Resolving source block geometry...'),
        progress_label='Resolving source block geometry',
    )
    _emit_progress(progress_callback, 57, 100, 'Resolving target block geometry...')
    target_geometry_df = target_df if remaining_target_rows is None else target_df.iloc[remaining_target_rows]
    target_geometry = _resolve_block_row_geometry(
        target_geometry_df,
        (target_x_col, target_y_col, target_z_col),
        target_block_size,
        target_size_cols,
        progress_callback=_make_scaled_progress_callback(progress_callback, 57, 70, 'Resolving target block geometry...'),
        progress_label='Resolving target block geometry',
    )
    if not len(source_geometry['centers']):
        raise ValueError('The source block model has no rows with valid coordinates and dimensions.')

    from scipy.spatial import cKDTree
    source_values_df = source_df.iloc[source_geometry['row_indices']]
    if column_modes is None:
        column_modes = _detect_dataframe_transfer_column_modes(source_values_df, selected_columns)
    values_by_column = {
        column: (pd.to_numeric(source_values_df[column], errors='coerce').to_numpy(dtype=float)
                 if mode == 'numeric' else source_values_df[column].fillna('').astype(str).str.strip().to_numpy(dtype=object))
        for column, mode in column_modes.items()
    }
    if assigned is None:
        assigned = {column: np.full(len(target_df), '', dtype=object) for column in selected_columns}
    tree = cKDTree(source_geometry['centers'])
    max_source_half_diagonal = float(np.linalg.norm(source_geometry['sizes'], axis=1).max() / 2.0)
    nearest_matches = 0
    target_count = len(target_geometry['centers'])
    target_row_indices = (
        target_geometry['row_indices']
        if remaining_target_rows is None else
        np.asarray(remaining_target_rows, dtype=int)[target_geometry['row_indices']]
    )
    _emit_progress(progress_callback, 70, 100, 'Matching target blocks to source blocks...')

    for local_index, target_row in enumerate(target_row_indices):
        target_center = target_geometry['centers'][local_index]
        radius = float(np.linalg.norm(target_geometry['sizes'][local_index]) / 2.0 + max_source_half_diagonal)
        candidates = np.asarray(tree.query_ball_point(target_center, radius), dtype=int)
        overlaps = np.empty(0, dtype=int)
        volumes = np.empty(0, dtype=float)
        if len(candidates):
            overlap_lengths = np.maximum(
                0.0,
                np.minimum(source_geometry['upper_bounds'][candidates], target_geometry['upper_bounds'][local_index])
                - np.maximum(source_geometry['lower_bounds'][candidates], target_geometry['lower_bounds'][local_index]),
            )
            candidate_volumes = np.prod(overlap_lengths, axis=1)
            positive = candidate_volumes > max(float(np.prod(target_geometry['sizes'][local_index])) * 1e-12, 1e-12)
            overlaps, volumes = candidates[positive], candidate_volumes[positive]
        if len(overlaps):
            overlap_matches += 1
            for column, mode in column_modes.items():
                values = values_by_column[column][overlaps]
                if mode == 'numeric':
                    valid = np.isfinite(values)
                    if valid.any():
                        assigned[column][target_row] = float(np.average(values[valid], weights=volumes[valid]))
                else:
                    totals = {}
                    for value, volume in zip(values, volumes):
                        value = str(value).strip()
                        if value and value.lower() != 'nan':
                            totals[value] = totals.get(value, 0.0) + float(volume)
                    if totals:
                        assigned[column][target_row] = sorted(totals.items(), key=lambda item: (-item[1], item[0]))[0][0]
        elif nearest_fallback:
            nearest_distance, nearest = tree.query(target_center, k=1)
            if nearest_distance_limit is None or float(nearest_distance) <= nearest_distance_limit:
                nearest = int(nearest)
                nearest_matches += 1
                for column, mode in column_modes.items():
                    value = values_by_column[column][nearest]
                    if mode == 'numeric' and np.isfinite(value):
                        assigned[column][target_row] = float(value)
                    elif mode != 'numeric' and str(value).strip() and str(value).strip().lower() != 'nan':
                        assigned[column][target_row] = str(value).strip()
        if progress_callback and (local_index + 1 == target_count or (local_index + 1) % 10_000 == 0):
            _emit_progress(
                progress_callback, 70 + int(round(((local_index + 1) / max(target_count, 1)) * 28)), 100,
                'Matching target blocks to source blocks...',
            )

    output_df = target_df.copy()
    for column in selected_columns:
        output_df[column] = pd.Series(assigned[column], index=output_df.index)
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    _emit_progress(progress_callback, 99, 100, 'Writing block-model transfer export...')
    output_df.to_csv(output_file, index=False, sep=output_delimiter if output_delimiter != 'bmf' else ',')
    _emit_progress(progress_callback, 100, 100, 'Block-model transfer complete.')
    invalid_targets = (
        len(target_df) - target_count
        if exact_transfer is None else
        int(exact_transfer['invalid_target_blocks'])
    )
    target_geometry_mode = (
        target_geometry['mode']
        if exact_transfer is None else
        f"exact-prematch + {target_geometry['mode']}"
    )
    unmatched = len(target_df) - overlap_matches - nearest_matches
    return {
        'output_file': output_file, 'total_target_blocks': int(len(target_df)),
        'overlap_matched_blocks': int(overlap_matches), 'nearest_matched_blocks': int(nearest_matches),
        'unmatched_blocks': int(unmatched), 'invalid_target_blocks': int(invalid_targets),
        'transferred_columns': list(selected_columns), 'column_modes': column_modes,
        'source_geometry_mode': source_geometry['mode'], 'target_geometry_mode': target_geometry_mode,
        'max_nearest_distance': nearest_distance_limit,
    }


def _collect_export_block_data(blocks):
    data = []
    min_bounds = blocks._block_info['min_bounds']
    block_size = blocks._block_info['block_size']
    rotation_matrix = blocks._block_info.get('rotation_matrix')
    rotation_center = blocks._block_info.get('rotation_center')
    
    # Check if we have multiple interpolators (sequential domain processing)
    if hasattr(blocks, '_interpolators'):
        print(f"  Processing {len(blocks._interpolators)} domain interpolators...")
        for domain, interpolator_list in blocks._interpolators.items():
            # interpolator_list is [Pass1, (optional) Pass2]
            
            # We export the LAST interpolator, which contains the cumulative state
            last_interp = interpolator_list[-1]
            first_interp = interpolator_list[0]
            
            # Identify original samples to distinguish them from Pass 1 outputs
            original_sample_positions = set()
            if len(interpolator_list) > 0:
                interp1 = interpolator_list[0]
                for pos, block in interp1.blocks.items():
                    is_sample = False
                    if hasattr(block, 'is_sample'): is_sample = block.is_sample
                    elif isinstance(block, dict): is_sample = block.get('is_sample', False)
                    
                    if is_sample:
                        original_sample_positions.add(pos)
            
            print(f"    Exporting domain: {domain} (Final Pass)")
            _add_interpolator_blocks_to_data(last_interp, min_bounds, block_size, data, rotation_matrix, rotation_center, 
                                           original_samples=original_sample_positions,
                                           pass_count=len(interpolator_list),
                                           forced_domain=domain,
                                           first_pass_algorithm_name=first_interp.get_algorithm_name(),
                                           final_algorithm_name=last_interp.get_algorithm_name())
    else:
        # Single interpolator
        interpolator = blocks._ant_colony
        # Pass domain_mapping from block_info to restore original domains
        domain_mapping = blocks._block_info.get('domain_mapping')
        algo_name = interpolator.get_algorithm_name()
        _add_interpolator_blocks_to_data(
            interpolator,
            min_bounds,
            block_size,
            data,
            rotation_matrix,
            rotation_center,
            domain_mapping=domain_mapping,
            first_pass_algorithm_name=algo_name,
            final_algorithm_name=algo_name,
        )

    return data


def _build_base_block_export_dataframe(block_rows):
    return pd.DataFrame([{k: v for k, v in row.items() if k != '_Grid_Index'} for row in block_rows])


def _build_export_block_row_lookup(block_rows):
    export_columns = list(_build_base_block_export_dataframe(block_rows).columns)
    export_lookup = {}
    for row in block_rows:
        grid_index = tuple(row['_Grid_Index'])
        export_lookup[grid_index] = {k: v for k, v in row.items() if k != '_Grid_Index'}
    return export_lookup, export_columns


def _expand_export_block_rows_to_source_rows(block_rows, source_blocks_df, block_x_col, block_y_col, block_z_col,
                                             min_bounds, block_size, rotation_matrix=None, rotation_center=None,
                                             grid_index_origin=None):
    export_df = _build_base_block_export_dataframe(block_rows)
    export_df['_Grid_Index'] = [tuple(row['_Grid_Index']) for row in block_rows]

    expanded_index_df = pd.DataFrame(index=source_blocks_df.index)
    expanded_index_df['_Grid_Index'] = None

    coord_frame = source_blocks_df[[block_x_col, block_y_col, block_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1)
    if valid_mask.any():
        valid_coords = coord_frame.loc[valid_mask].to_numpy(copy=False)
        if rotation_matrix is not None and rotation_center is not None:
            valid_coords = (valid_coords - rotation_center) @ rotation_matrix.T
        index_origin = np.asarray(grid_index_origin if grid_index_origin is not None else min_bounds, dtype=float)
        block_indices = np.floor((valid_coords - index_origin) / np.asarray(block_size, dtype=float) + 1e-6).astype(int)
        expanded_index_df.loc[valid_mask, '_Grid_Index'] = [tuple(int(v) for v in idx) for idx in block_indices]

    expanded_df = expanded_index_df.merge(export_df, on='_Grid_Index', how='left')

    # Preserve the original source row coordinates so expanded output matches the input block rows.
    expanded_df['x'] = source_blocks_df[block_x_col].to_numpy(copy=False)
    expanded_df['y'] = source_blocks_df[block_y_col].to_numpy(copy=False)
    expanded_df['z'] = source_blocks_df[block_z_col].to_numpy(copy=False)

    return expanded_df.drop(columns=['_Grid_Index'])


def _expand_export_chunk_rows_to_source_rows(source_blocks_df, export_lookup, export_columns,
                                             block_x_col, block_y_col, block_z_col,
                                             min_bounds, block_size,
                                             rotation_matrix=None, rotation_center=None,
                                             grid_index_origin=None):
    expanded_records = [{} for _ in range(len(source_blocks_df))]

    coord_frame = source_blocks_df[[block_x_col, block_y_col, block_z_col]].apply(pd.to_numeric, errors='coerce')
    valid_mask = coord_frame.notna().all(axis=1)
    valid_positions = np.flatnonzero(valid_mask.to_numpy())

    if len(valid_positions):
        valid_coords = coord_frame.loc[valid_mask].to_numpy(copy=False)
        if rotation_matrix is not None and rotation_center is not None:
            valid_coords = (valid_coords - rotation_center) @ rotation_matrix.T
        index_origin = np.asarray(grid_index_origin if grid_index_origin is not None else min_bounds, dtype=float)
        block_indices = np.floor((valid_coords - index_origin) / np.asarray(block_size, dtype=float) + 1e-6).astype(int)
        for row_position, block_index in zip(valid_positions, block_indices):
            expanded_records[row_position] = export_lookup.get(tuple(int(v) for v in block_index), {})

    expanded_df = pd.DataFrame.from_records(expanded_records, columns=export_columns)
    expanded_df['x'] = source_blocks_df[block_x_col].to_numpy(copy=False)
    expanded_df['y'] = source_blocks_df[block_y_col].to_numpy(copy=False)
    expanded_df['z'] = source_blocks_df[block_z_col].to_numpy(copy=False)
    return expanded_df


def _iter_expanded_export_block_chunks(block_rows, source_blocks_file, blocks_delimiter=None,
                                       blocks_header_line=1, block_filters=None,
                                       block_x_col=None, block_y_col=None, block_z_col=None,
                                       min_bounds=None, block_size=None,
                                       rotation_matrix=None, rotation_center=None,
                                       grid_index_origin=None,
                                       progress_label='Reading source blocks for export expansion',
                                       chunksize=250_000):
    delimiter = blocks_delimiter or detect_csv_delimiter(source_blocks_file)
    header_line = resolve_effective_csv_header_line(source_blocks_file, blocks_header_line)
    headers = parse_header_line(source_blocks_file, delimiter, header_line)
    final_names = build_unique_column_names(headers)
    block_x_col, block_y_col, block_z_col = resolve_block_coordinate_columns(
        final_names,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
    )
    selected_columns = list(dict.fromkeys([
        block_x_col,
        block_y_col,
        block_z_col,
        *collect_filter_fields(block_filters or []),
    ]))
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=header_line,
        comment='#',
        usecols=selected_columns,
        chunksize=chunksize,
    )
    export_lookup, export_columns = _build_export_block_row_lookup(block_rows)

    for chunk_number, chunk in enumerate(
        iterate_csv_path_chunks_with_progress(
            source_blocks_file,
            progress_label,
            header_line=header_line,
            **read_kwargs,
        ),
        start=1,
    ):
        chunk = strip_leading_non_data_rows(chunk)
        if block_filters:
            chunk, _ = apply_dataframe_filters(
                chunk,
                filters=block_filters,
                filter_subject='block',
                source_label=f'blocks file chunk {chunk_number:,}',
                emit_logs=False,
            )
        if chunk.empty:
            continue
        yield _expand_export_chunk_rows_to_source_rows(
            chunk,
            export_lookup,
            export_columns,
            block_x_col,
            block_y_col,
            block_z_col,
            min_bounds,
            block_size,
            rotation_matrix=rotation_matrix,
            rotation_center=rotation_center,
            grid_index_origin=grid_index_origin,
        )


def _write_expanded_export_blocks_to_csv(block_rows, block_info, filepath):
    source_blocks_file = str(block_info.get('source_blocks_file') or '').strip()
    total_rows_written = 0
    wrote_header = False

    print(f"Streaming expanded export rows from {os.path.basename(source_blocks_file)}...")

    for chunk_number, expanded_chunk in enumerate(
        _iter_expanded_export_block_chunks(
            block_rows,
            source_blocks_file,
            blocks_delimiter=block_info.get('source_blocks_delimiter'),
            blocks_header_line=block_info.get('source_blocks_header_line', 1),
            block_filters=block_info.get('source_block_filters'),
            block_x_col=block_info.get('source_block_x_col'),
            block_y_col=block_info.get('source_block_y_col'),
            block_z_col=block_info.get('source_block_z_col'),
            min_bounds=block_info['min_bounds'],
            block_size=block_info['block_size'],
            rotation_matrix=block_info.get('rotation_matrix'),
            rotation_center=block_info.get('rotation_center'),
            grid_index_origin=block_info.get('grid_index_origin'),
        ),
        start=1,
    ):
        expanded_chunk.to_csv(filepath, index=False, mode='w' if not wrote_header else 'a', header=not wrote_header)
        wrote_header = True
        total_rows_written += len(expanded_chunk)
        if chunk_number == 1 or chunk_number % 10 == 0:
            print(
                f"Expanded export chunk {chunk_number:,}: wrote {len(expanded_chunk):,} rows "
                f"({total_rows_written:,} total)."
            )

    if not wrote_header:
        empty_columns = list(_build_base_block_export_dataframe(block_rows).columns)
        pd.DataFrame(columns=empty_columns).to_csv(filepath, index=False)

    print(
        f"Expanded {len(block_rows):,} base-block export rows to {total_rows_written:,} source block rows "
        f"from {os.path.basename(source_blocks_file)}."
    )
    return total_rows_written


def _build_export_blocks_dataframe(blocks, block_rows):
    base_df = _build_base_block_export_dataframe(block_rows)
    block_info = getattr(blocks, '_block_info', {}) or {}

    if not block_rows or base_df.empty:
        return base_df

    if not block_info.get('expand_interpolation_exports_to_subblocks', False) or _blocks_use_adaptive_leaf_cover(blocks):
        return base_df

    source_blocks_file = str(block_info.get('source_blocks_file') or '').strip()
    if not source_blocks_file or not os.path.isfile(source_blocks_file):
        return base_df

    try:
        source_blocks_df, _ = load_full_blocks_dataframe(
            source_blocks_file,
            blocks_delimiter=block_info.get('source_blocks_delimiter'),
            blocks_header_line=block_info.get('source_blocks_header_line', 1),
            block_filters=block_info.get('source_block_filters'),
            progress_label='Reading source blocks for export expansion',
        )
        block_x_col, block_y_col, block_z_col = resolve_block_coordinate_columns(
            list(source_blocks_df.columns),
            block_x_col=block_info.get('source_block_x_col'),
            block_y_col=block_info.get('source_block_y_col'),
            block_z_col=block_info.get('source_block_z_col'),
        )
        expanded_df = _expand_export_block_rows_to_source_rows(
            block_rows,
            source_blocks_df,
            block_x_col,
            block_y_col,
            block_z_col,
            block_info['min_bounds'],
            block_info['block_size'],
            rotation_matrix=block_info.get('rotation_matrix'),
            rotation_center=block_info.get('rotation_center'),
            grid_index_origin=block_info.get('grid_index_origin'),
        )
        print(
            f"Expanded {len(block_rows):,} base-block export rows to {len(expanded_df):,} source block rows "
            f"from {os.path.basename(source_blocks_file)}."
        )
        return expanded_df
    except Exception as exc:
        print(f"Failed to expand export back to source block rows, falling back to base-block export: {exc}")
        return base_df


def export_blocks_to_csv(blocks, filepath):
    output_dir = os.path.dirname(filepath) or '.'
    os.makedirs(output_dir, exist_ok=True)

    print(f"Exporting blocks to {filepath}...")
    data = _collect_export_block_data(blocks)
    block_info = getattr(blocks, '_block_info', {}) or {}
    source_blocks_file = str(block_info.get('source_blocks_file') or '').strip()
    should_stream_expanded_export = (
        bool(data)
        and block_info.get('expand_interpolation_exports_to_subblocks', False)
        and not _blocks_use_adaptive_leaf_cover(blocks)
        and source_blocks_file
        and os.path.isfile(source_blocks_file)
        and os.path.getsize(source_blocks_file) >= LARGE_BLOCK_FILE_THRESHOLD
    )

    if should_stream_expanded_export:
        rows_written = _write_expanded_export_blocks_to_csv(data, block_info, filepath)
        print(f"Exported {rows_written:,} rows to {filepath}")
        return

    df = _build_export_blocks_dataframe(blocks, data)
    df.to_csv(filepath, index=False)
    print(f"Exported {len(df):,} rows to {filepath}")


def export_dataframe_to_bmf(df, filepath, backend='tbms-config-text'):
    output_dir = os.path.dirname(filepath) or '.'
    os.makedirs(output_dir, exist_ok=True)

    if not {'x', 'y', 'z'}.issubset(df.columns):
        missing = sorted({'x', 'y', 'z'} - set(df.columns))
        raise ValueError(f"BMF export requires columns x, y, z. Missing: {missing}")

    value_cols = [
        column_name
        for column_name in df.columns
        if column_name not in {'x', 'y', 'z', 'grid_i', 'grid_j', 'grid_k'}
    ]
    with tempfile.NamedTemporaryFile('w', suffix='.csv', prefix='anterpolator_bmf_', delete=False, encoding='utf-8', newline='') as handle:
        temp_csv_path = handle.name
    try:
        df.to_csv(temp_csv_path, index=False)
        return _get_bmf_tools_module().export_bmf(
            temp_csv_path,
            filepath,
            backend=backend,
            x_col='x',
            y_col='y',
            z_col='z',
            value_cols=value_cols,
            dry_run=False,
        )
    finally:
        try:
            os.remove(temp_csv_path)
        except OSError:
            pass


def export_blocks_to_file(blocks, filepath, bmf_backend='tbms-config-text'):
    output_path = resolve_interpolation_csv_export_path(filepath)
    if output_path != str(filepath or '').strip():
        print(f"Interpolation BMF output is disabled; writing CSV instead: {output_path}")
    export_blocks_to_csv(blocks, output_path)
    return output_path


def parse_bmf_export_column_types(text):
    overrides = {}
    raw_text = str(text or '').strip()
    if not raw_text:
        return overrides

    for token in raw_text.split(','):
        entry = token.strip()
        if not entry:
            continue
        if '=' not in entry:
            raise ValueError(
                'Column type overrides must use the format column=type, separated by commas. '
                'Example: grade=double, domain=string'
            )
        column_name, field_type = entry.split('=', 1)
        column_name = column_name.strip()
        field_type = field_type.strip()
        if not column_name or not field_type:
            raise ValueError(
                'Column type overrides must use the format column=type, separated by commas. '
                'Example: grade=double, domain=string'
            )
        overrides[column_name] = normalize_bmf_export_field_type_name(field_type)
    return overrides


def normalize_bmf_export_field_type_name(field_type):
    normalized = str(field_type or '').strip().lower()
    alias_map = {
        '': '',
        'auto': '',
        'boolean': 'boolean',
        'bool': 'boolean',
        'int': 'int',
        'integer': 'int',
        'double': 'double',
        'float': 'double',
        'real': 'double',
        'string': 'string',
        'text': 'string',
        'namedshort': 'string',
    }
    if normalized not in alias_map:
        raise ValueError(
            'Unsupported BMF export field type. Use one of: auto, boolean, int, double, string.'
        )
    return alias_map[normalized]


def infer_bmf_export_field_types_from_preview(df, candidate_columns, delimiter=None):
    return _get_bmf_tools_module().infer_bmf_export_field_types_from_preview(
        df,
        candidate_columns,
        delimiter=delimiter,
    )


def export_csv_grid_to_bmf(input_csv, output_bmf, x_col='x', y_col='y', z_col='z',
                           value_cols=None, backend='tbms-config-text', delimiter=None,
                           header_line=1, column_types=None, value_exceptions=None, cell_size=None,
                           origin=None, null_float=-99.0, index_tolerance=1e-3,
                           size_cols=None, extent_cols=None,
                           regularize_to_base_block=False, summary_json=None, progress_callback=None):
    if progress_callback is not None:
        progress_callback(0, 100, 'Preparing BMF export...')
    result = _get_bmf_tools_module().export_bmf(
        input_csv=input_csv,
        output_bmf=output_bmf,
        backend=backend,
        x_col=x_col,
        y_col=y_col,
        z_col=z_col,
        delimiter=delimiter,
        header_line=header_line,
        cell_size=cell_size,
        origin=origin,
        size_cols=size_cols,
        extent_cols=extent_cols,
        value_cols=value_cols,
        column_types=column_types,
        value_exceptions=value_exceptions,
        null_float=null_float,
        index_tolerance=index_tolerance,
        regularize_to_base_block=regularize_to_base_block,
        dry_run=False,
        summary_json=summary_json,
        progress_callback=progress_callback,
    )
    if progress_callback is not None:
        progress_callback(100, 100, 'BMF export complete.')
    return result


def export_block_evaluated_samples_to_csv(blocks, filepath):
    output_dir = os.path.dirname(filepath) or '.'
    os.makedirs(output_dir, exist_ok=True)

    sample_assignment_data = getattr(blocks, '_sample_assignment_data', None)
    if not sample_assignment_data:
        raise ValueError("No sample assignment data available for evaluated sample export.")

    print(f"Exporting block-evaluated samples to {filepath}...")

    block_rows = _collect_export_block_data(blocks)
    final_block_lookup = {tuple(row['_Grid_Index']): row for row in block_rows}

    sample_points = np.asarray(sample_assignment_data.get('points', []))
    sample_values = np.asarray(sample_assignment_data.get('values', []))
    sample_domains = sample_assignment_data.get('domains')
    block_indices = np.asarray(sample_assignment_data.get('block_indices', []))
    assigned_mask = np.asarray(sample_assignment_data.get('assigned_mask', []), dtype=bool)

    if len(sample_points) != len(block_indices) or len(sample_points) != len(assigned_mask):
        raise ValueError("Sample assignment metadata is inconsistent; evaluated sample export aborted.")

    data = []
    matched_count = 0
    for i in range(len(sample_points)):
        sample_point = sample_points[i]
        sample_value = sample_values[i] if i < len(sample_values) else None
        sample_domain = None
        if sample_domains is not None and i < len(sample_domains):
            sample_domain = sample_domains[i]

        block_idx = tuple(int(v) for v in block_indices[i])
        final_block = final_block_lookup.get(block_idx) if assigned_mask[i] else None

        status = 'outside_allowed_grid'
        if assigned_mask[i]:
            status = 'matched' if final_block is not None else 'missing_final_block'
        if status == 'matched':
            matched_count += 1

        block_value = final_block.get('Value') if final_block is not None else None
        value_difference = None
        if block_value is not None and pd.notna(block_value) and sample_value is not None and pd.notna(sample_value):
            value_difference = float(block_value) - float(sample_value)

        data.append({
            'Sample_x': sample_point[0],
            'Sample_y': sample_point[1],
            'Sample_z': sample_point[2],
            'Sample_Value': sample_value,
            'Sample_Domain': sample_domain,
            'Block_ix': block_idx[0],
            'Block_iy': block_idx[1],
            'Block_iz': block_idx[2],
            'Block_x': final_block.get('x') if final_block is not None else None,
            'Block_y': final_block.get('y') if final_block is not None else None,
            'Block_z': final_block.get('z') if final_block is not None else None,
            'Block_Value': block_value,
            'Block_Domain': final_block.get('Domain') if final_block is not None else None,
            'Block_Source': final_block.get('Source') if final_block is not None else None,
            'Block_Algorithm': final_block.get('Algorithm') if final_block is not None else None,
            'Block_ID': final_block.get('Block_ID') if final_block is not None else None,
            'Value_Difference': value_difference,
            'Assignment_Status': status,
        })

    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Exported {len(data)} evaluated samples to {filepath} ({matched_count} matched to final blocks)")

def _normalize_export_algorithm_name(algo_name):
    return normalize_provenance_algorithm_name(algo_name)


def _add_interpolator_blocks_to_data(interpolator, min_bounds, block_size, data, rotation_matrix=None, rotation_center=None, domain_mapping=None, original_samples=None, pass_count=1, forced_domain=None, first_pass_algorithm_name=None, final_algorithm_name=None):
    """Process blocks from an interpolator and add to data list"""
    final_algorithm_label, final_algo_type = _normalize_export_algorithm_name(
        final_algorithm_name or interpolator.get_algorithm_name()
    )
    first_pass_algorithm_label, first_pass_algo_type = _normalize_export_algorithm_name(
        first_pass_algorithm_name or final_algorithm_name or interpolator.get_algorithm_name()
    )
    
    for pos, block in tqdm(interpolator.blocks.items(), desc="Processing blocks"):
        if hasattr(block, 'value'):
            relative_size = np.ones(3, dtype=float)
        else:
            relative_size = np.asarray(block.get('relative_size', (1, 1, 1)), dtype=float)
            if relative_size.shape != (3,):
                relative_size = np.ones(3, dtype=float)
        base_block_size = np.asarray(block_size, dtype=float)
        actual_block_size = base_block_size * relative_size
        centroid = min_bounds + (np.asarray(pos, dtype=float) + 0.5 * (relative_size - 1.0)) * base_block_size
        
        # Apply inverse rotation if needed
        if rotation_matrix is not None and rotation_center is not None:
            # P_orig = P_aligned @ R + Center
            centroid = centroid @ rotation_matrix + rotation_center
        
        # Determine Source / Algorithm provenance.
        export_provenance = get_export_provenance(interpolator, pos)
        if export_provenance is not None:
            source, algorithm_label, algo_type = export_provenance
        else:
            source = "Unknown"
            is_sample = False
            if hasattr(block, 'is_sample'): is_sample = block.is_sample
            elif isinstance(block, dict): is_sample = block.get('is_sample', False)
            
            if is_sample:
                if original_samples is None or pos in original_samples:
                    source = "Original Sample"
                    algorithm_label = "Sample"
                    algo_type = "sample"
                else:
                    # It's a sample in this pass, but not in original -> Must be from previous pass
                    source = "First Pass"
                    algorithm_label = first_pass_algorithm_label
                    algo_type = first_pass_algo_type
            else:
                # It's an interpolated block in this pass
                if pass_count == 1:
                    source = "First Pass"
                    algorithm_label = final_algorithm_label
                    algo_type = final_algo_type
                else:
                    source = "Second Pass"
                    algorithm_label = final_algorithm_label
                    algo_type = final_algo_type

        # Initialize common fields with None/NaN for all possible columns
        row = {
            '_Grid_Index': tuple(int(v) for v in pos),
            'x': centroid[0],
            'y': centroid[1],
            'z': centroid[2],
            'Algorithm': algorithm_label,
            'Algo_Type': algo_type,
            'Source': source,
            'Value': None,
            'Age': None,
            'Is_Sample': None,
            'Domain': None,
            # Ant Colony specific
            'Distance_To_Sample': None,
            'Nearest_Sample_Value': None,
            'Mark_Class': None,
            'Pheromone': None,
            'Visits': None,
            'Ant_Count': None,
            'Block_ID': None,
            # Molecular Clock specific
            'Event_ID': None,
            'Branch_ID': None,
            'Distance_To_Feeder': None,
            'Is_Feeder': None,
            'Block_Size_X': float(actual_block_size[0]),
            'Block_Size_Y': float(actual_block_size[1]),
            'Block_Size_Z': float(actual_block_size[2]),
            'AdaptiveLevel': None,
            'AdaptiveLeafID': None,
            'AdaptiveLeafSizeX': None,
            'AdaptiveLeafSizeY': None,
            'AdaptiveLeafSizeZ': None,
            'AdaptiveInherited': None,
            'AdaptiveOutputMode': None,
        }
        
        # Get block data - handle both AntColony Block dataclass and dict
        if hasattr(block, 'value'):  # AntColony Block dataclass
            # Age: negative distance from sample (0 at sample, more negative further away)
            age = -block.distance_to_sample if not block.is_sample else 0
            
            # Determine domain
            if forced_domain:
                domain = forced_domain
            elif domain_mapping and pos in domain_mapping:
                domain = domain_mapping[pos]
            else:
                domain = block.domain
            
            row.update({
                'Value': block.value,
                'Age': age,
                'Is_Sample': block.is_sample,
                'Domain': domain,
                'Distance_To_Sample': block.distance_to_sample,
                'Nearest_Sample_Value': block.nearest_sample_value,
                'Mark_Class': block.mark_class,
                'Pheromone': block.pheromone,
                'Visits': block.visit_count,
                'Ant_Count': block.ant_count,
                'Block_ID': block.block_id
            })
        else:  # BiochemicalClock dict format
            # Age: negative distance to feeder (0 at samples, more negative toward feeder/LUCA)
            # Samples have distance_to_feeder = 0, blocks closer to feeder have larger values
            # So we negate: Age = -distance_to_feeder (more negative = older/closer to LUCA)
            dist_to_feeder = block.get('distance_to_feeder', 0)
            is_sample = block.get('is_sample', False)
            is_feeder = block.get('is_feeder', False)
            
            if is_feeder:
                age = -dist_to_feeder  # Feeder is the LUCA (most negative)
            elif is_sample:
                age = 0  # Samples are the "present"
            else:
                age = -dist_to_feeder  # Intermediate blocks get negative ages
            
            # Determine domain
            if forced_domain:
                domain = forced_domain
            elif domain_mapping and pos in domain_mapping:
                domain = domain_mapping[pos]
            else:
                domain = block.get('domain', 'Undomained')
            
            row.update({
                'Value': block.get('value', 0),
                'Age': age,
                'Is_Sample': is_sample,
                'Domain': domain,
                'Event_ID': block.get('event_id', -1),
                'Branch_ID': block.get('branch_id', -1),
                'Distance_To_Feeder': dist_to_feeder,
                'Is_Feeder': is_feeder,
                'AdaptiveLevel': block.get('adaptive_level', block.get('level')),
                'AdaptiveLeafID': block.get('adaptive_leaf_id', block.get('leaf_id')),
                'AdaptiveLeafSizeX': float(np.asarray(block.get('adaptive_relative_size', block.get('relative_size', (1, 1, 1))), dtype=float)[0] * base_block_size[0]),
                'AdaptiveLeafSizeY': float(np.asarray(block.get('adaptive_relative_size', block.get('relative_size', (1, 1, 1))), dtype=float)[1] * base_block_size[1]),
                'AdaptiveLeafSizeZ': float(np.asarray(block.get('adaptive_relative_size', block.get('relative_size', (1, 1, 1))), dtype=float)[2] * base_block_size[2]),
                'AdaptiveInherited': block.get('is_inherited'),
                'AdaptiveOutputMode': interpolator.get_metadata().get('output_mode') if hasattr(interpolator, 'get_metadata') else None,
            })
            
        data.append(row)


def _iter_final_interpolators(blocks):
    if hasattr(blocks, '_interpolators') and blocks._interpolators:
        for interpolator_list in blocks._interpolators.values():
            if isinstance(interpolator_list, list) and interpolator_list:
                yield interpolator_list[-1]
            elif interpolator_list is not None:
                yield interpolator_list
        return
    interpolator = getattr(blocks, '_ant_colony', None)
    if interpolator is not None:
        yield interpolator


def _blocks_use_adaptive_leaf_cover(blocks):
    for interpolator in _iter_final_interpolators(blocks):
        try:
            metadata = interpolator.get_metadata()
        except Exception:
            metadata = {}
        if str(metadata.get('output_mode', '')).strip().lower() == 'adaptive_leaf_cover':
            return True
    return False

def _get_interpolator_run_profile(interpolator, dims, iterations):
    algo_name = interpolator.get_algorithm_name()
    if algo_name == 'Gaussian Kernel':
        allowed_positions = getattr(interpolator, 'allowed_grid_override', None)
        if allowed_positions is None:
            position_count = int(dims[0]) * int(dims[1]) * int(dims[2])
            scope_label = 'grid positions'
        else:
            position_count = len(allowed_positions)
            scope_label = 'allowed positions'

        sample_count = len(getattr(interpolator, 'sample_blocks', {}) or {})
        bandwidth = float(getattr(interpolator, 'bandwidth', 0.0))
        cutoff_sigma = float(getattr(interpolator, 'cutoff_sigma', 0.0))
        search_radius = bandwidth * cutoff_sigma
        message = (
            f"Running {algo_name} single pass over {position_count:,} {scope_label} "
            f"(samples={sample_count:,}, radius={search_radius:.1f} blocks)..."
        )
        return {
            'single_pass': True,
            'message': message,
            'desc_suffix': 'single pass',
            'total': position_count,
            'unit': 'pos',
        }

    if algo_name == 'String Theory':
        sample_count = len(getattr(interpolator, 'sorted_samples', []) or [])
        if sample_count == 0:
            sample_count = len(getattr(interpolator, 'sample_locations', set()) or set())
        return {
            'single_pass': True,
            'message': f"Running {algo_name} single pass over {sample_count:,} samples...",
            'desc_suffix': 'sample scan',
            'total': max(int(sample_count), 1),
            'unit': 'sample',
        }

    return {
        'single_pass': False,
        'message': (
            f"Running {algo_name} for {iterations} iterations"
            f"{_format_interpolator_run_counts(interpolator)}..."
        ),
        'desc_suffix': None,
        'total': max(int(iterations), 1),
        'unit': 'iter',
    }


def _format_interpolator_run_counts(interpolator):
    counts = []

    if hasattr(interpolator, 'ants'):
        sample_count = sum(1 for block in getattr(interpolator, 'blocks', {}).values() if getattr(block, 'is_sample', False))
        ant_count = len(getattr(interpolator, 'ants', []) or [])
        if sample_count:
            counts.append(f"samples={sample_count:,}")
        counts.append(f"ants={ant_count:,}")

    if not counts:
        return ""

    return f" ({', '.join(counts)})"

def _run_interpolator_with_progress(interpolator, dims, iterations, desc, on_first_iteration=None):
    profile = _get_interpolator_run_profile(interpolator, dims, iterations)
    print(profile['message'])

    if profile['single_pass']:
        pbar = tqdm(total=profile['total'], desc=f"{desc} ({profile['desc_suffix']})", unit=profile['unit'])
        previous_callback = getattr(interpolator, 'progress_callback', None)

        def update_progress(processed_positions, total_positions, interpolated_blocks):
            pbar.total = max(int(total_positions), 1)
            pbar.n = min(int(processed_positions), pbar.total)
            pbar.set_postfix_str(f"interpolated={int(interpolated_blocks):,}")
            pbar.refresh()

        interpolator.progress_callback = update_progress
        try:
            should_continue = interpolator.run_iteration(dims)
        finally:
            interpolator.progress_callback = previous_callback
        if on_first_iteration is not None:
            on_first_iteration()
        metadata = interpolator.get_metadata()
        interpolated_blocks = metadata.get('interpolated_blocks')
        if interpolated_blocks is not None:
            pbar.set_postfix_str(f"interpolated={interpolated_blocks:,}")
        pbar.n = pbar.total
        pbar.set_description(f"{desc} (completed)")
        pbar.close()
        print("Completed single pass")
        return should_continue

    pbar = tqdm(range(profile['total']), desc=desc, unit=profile['unit'])
    should_continue = True
    for i in pbar:
        should_continue = interpolator.run_iteration(dims)
        if i == 0 and on_first_iteration is not None:
            on_first_iteration()
        if not should_continue or interpolator.is_converged():
            pbar.set_description(f"{desc} (converged)")
            print(f"Converged at iteration {i+1}")
            break
    pbar.close()
    return should_continue


def _run_interpolator_statistics_with_retry(interpolator, output_dir, domain_name, parent=None):
    if not hasattr(interpolator, 'generate_statistics'):
        return True

    while True:
        try:
            interpolator.generate_statistics(output_dir, domain_name=domain_name)
            return True
        except PermissionError as exc:
            locked_path = getattr(exc, 'filename', None) or str(exc)
            print(f"Permission denied while saving interpolation statistics: {locked_path}")

            message_box = QtWidgets.QMessageBox(parent)
            message_box.setIcon(QtWidgets.QMessageBox.Warning)
            message_box.setWindowTitle('Statistics Save Failed')
            message_box.setText('Could not save the interpolation statistics file.')
            message_box.setInformativeText(
                f"Permission denied for:\n{locked_path}\n\n"
                'The file may be open in Excel or locked by OneDrive sync. '
                'Close the file and click Retry to resume saving statistics, or click Cancel to skip the statistics export and continue.'
            )
            message_box.setDetailedText(str(exc))
            message_box.setStandardButtons(QtWidgets.QMessageBox.Retry | QtWidgets.QMessageBox.Cancel)
            message_box.setDefaultButton(QtWidgets.QMessageBox.Retry)
            if message_box.exec_() == QtWidgets.QMessageBox.Retry:
                continue

            print(f"Skipped interpolation statistics for {domain_name}.")
            return False

def _normalize_domain_post_process_mode(mode):
    normalized = str(mode or 'skip').strip().lower()
    if normalized in ('fill_with_average', 'fill average', 'fill_average', 'fill'):
        return 'fill_with_average'
    return 'skip'

def _get_domain_post_process_mode(config, domain):
    if not config or not domain or str(domain).strip().lower() == 'global':
        return 'skip'
    domain_overrides = config.get('domain_algorithm_overrides') or {}
    domain_cfg = domain_overrides.get(domain) or {}
    if domain_cfg.get('skip', False):
        return 'skip'
    return _normalize_domain_post_process_mode(domain_cfg.get('post_process', 'skip'))

def _run_domain_post_process(interpolator, dims, mode, domain_label=None):
    if _normalize_domain_post_process_mode(mode) != 'fill_with_average':
        return 0, 0
    if not hasattr(interpolator, 'fill_unvisited_blocks_domainwise'):
        return 0, 0
    try:
        created, assigned = interpolator.fill_unvisited_blocks_domainwise(dims)
    except Exception as exc:
        label = str(domain_label or getattr(interpolator, 'get_algorithm_name', lambda: 'interpolator')())
        print(f"Post-process error for {label}: {exc}")
        return 0, 0
    return int(created or 0), int(assigned or 0)

def silent_interpolation(plotter, iterations, interpolation_file):
    blocks = plotter._blocks_data
    dims = tuple(blocks._block_info['dims'])
    block_evaluated_samples_file = getattr(plotter, '_block_evaluated_samples_file', None)
    post_process_config = {'domain_algorithm_overrides': getattr(plotter, '_domain_post_process_overrides', {})}
    
    # Check if we have multiple interpolators (sequential domain processing)
    if hasattr(blocks, '_interpolators'):
        print(f"Running sequential domain interpolation for {len(blocks._interpolators)} domains...")
        
        for domain_idx, (domain, interpolator_list) in enumerate(blocks._interpolators.items(), 1):
            # interpolator_list is [Pass1, (optional) Pass2]
            
            # --- Pass 1 ---
            interp1 = interpolator_list[0]
            algo_name1 = interp1.get_algorithm_name()
            print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} - Pass 1 ({algo_name1}) ===")
            seed_original_sample_provenance(interp1)
            pass1_snapshot = snapshot_interpolator_state(interp1)
            
            # Force verbose for first iteration if it's an AntColony
            if hasattr(interp1, 'verbose'):
                original_verbose = interp1.verbose
                interp1.verbose = True

            _run_interpolator_with_progress(
                interp1,
                dims,
                iterations,
                f"Domain {domain} - Pass 1",
                on_first_iteration=(lambda interp=interp1, verbose=original_verbose: setattr(interp, 'verbose', verbose)) if hasattr(interp1, 'verbose') else None,
            )
            
            # Generate stats for Pass 1 if String Theory
            output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
            _run_interpolator_statistics_with_retry(interp1, output_dir, f"{domain}_Pass1")
            finalize_phase_provenance(interp1, 'First Pass', algo_name1, pass1_snapshot)

            last_interp = interp1

            # --- Pass 2 ---
            if len(interpolator_list) > 1:
                interp2 = interpolator_list[1]
                algo_name2 = interp2.get_algorithm_name()
                print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} - Pass 2 ({algo_name2}) ===")
                
                print("  Transferring data from Pass 1 to Pass 2...")
                # Get Pass 1 results (pos -> value)
                pass1_values = interp1.get_interpolated_values()
                print(f"  Pass 1 generated {len(pass1_values)} blocks (samples + interpolated).")
                
                # Re-initialize Pass 2 with merged samples
                min_bounds = blocks._block_info['min_bounds']
                block_size = blocks._block_info['block_size']
                pass1_domain_mapping = {pos: domain for pos in pass1_values}
                
                # Determine if we should enforce domain mapping/grid restriction
                use_mapping = False
                if hasattr(interp2, 'allowed_grid_override') and interp2.allowed_grid_override is not None:
                    use_mapping = True
                
                # Double check for undomained special case where allowed_grid_override might be incorrectly set
                # or if we just want to be absolutely sure we don't trap the ants.
                if domain == "Undomained" and "Undomained" in blocks._interpolators and len(blocks._interpolators) == 1:
                     # If we are in the undomained special case, and we suspect the user didn't provide a blocks file
                     # (or provided an empty one which resulted in auto-blocks), we should ensure use_mapping is False
                     # UNLESS we are sure allowed_grid_override is valid and intended.
                     
                     # In the "Auto-blocks" path (create_blocks else branch), we explicitly do NOT set allowed_grid_override.
                     # So use_mapping should be False.
                    
                     # However, if the user provided a blocks file that was just samples (sparse), allowed_grid_override IS set.
                     # And use_mapping IS True.
                    
                     # And the ants ARE trapped.
                     
                     # The user claims "the blocks file was indeed empty".
                     # If the blocks file was empty, create_blocks should have gone to the 'else' branch (if my detection logic is correct).
                     # OR, it went to the 'if' branch but the dataframe was empty?
                     
                     # Let's look at create_blocks logic for empty file.
                     # if blocks_file is not None:
                     #    df_blocks = read_autodetect_csv(...)
                     
                     # If the file is empty (0 bytes), read_autodetect_csv might fail or return empty DF.
                     # If it returns empty DF, create_blocks raises ValueError("All block rows have non-numeric coordinates...").
                     
                     # So the user must have NOT provided a blocks file path in the UI?
                     # If the UI field is empty "", then blocks_file is None (if passed correctly).
                     
                     # In run_interpolation_cli:
                     # blocks_file = self.blocks_edit.text() if self.blocks_edit.text() else None
                     
                     # So if the field is empty, blocks_file is None.
                     # So we go to the 'else' branch of create_blocks.
                     
                     # In the 'else' branch:
                     # interp2 = create_interpolator(config2)
                     # interp2.initialize_blocks(..., use_domain_mapping=False)
                     
                     # create_interpolator does NOT set allowed_grid_override unless config has it.
                     # config comes from UI.
                     
                     # So interp2.allowed_grid_override should be None.
                     # So use_domain_mapping=True is fine?
                     # Wait, if use_domain_mapping=True, it tries to look up domains in self.domain_mapping.
                     # If self.domain_mapping is None (which it is in Auto-blocks), it defaults to "default" or "Undomained".
                     
                     # So why did the user get "Not in allowed_positions"?
                     # This error comes from AntColony.run_iteration -> move_ant.
                     # It checks: if npos not in self.allowed_positions:
                     
                     # In initialize_blocks (AntColony):
                     # if use_domain_mapping:
                     #    if hasattr(self, 'allowed_grid_override'): ...
                     #    else: self.allowed_positions = set(sample_blocks.keys())
                     # else:
                     #    self.allowed_positions = {full grid}
                     
                     # So if use_mapping is False, allowed_positions is full grid.
                     # If allowed_positions is full grid, "Not in allowed_positions" should only happen at the edges of the bounding box.
                     # 900k rejections suggests it's happening everywhere.
                     
                     # So allowed_positions MUST be restricted.
                     # So use_mapping MUST be True.
                     # So allowed_grid_override MUST be set.
                     
                     # How did allowed_grid_override get set in the 'else' branch?
                     # I checked the code, it's not set there.
                     
                     # Is it possible that 'interp2' in the loop is NOT the one created in the 'else' branch?
                     # multiblock._interpolators["Undomained"] = [interp1, interp2]
                     # In the loop: interp2 = interpolator_list[1]
                     
                     # It should be the same object.
                     
                     # Determine if we should enforce domain mapping/grid restriction
                     use_mapping = False
                     if hasattr(interp2, 'allowed_grid_override') and interp2.allowed_grid_override is not None:
                         use_mapping = True
                     
                     interp2.initialize_blocks(
                         pass1_values,
                         dims,
                         min_bounds,
                         block_size,
                         use_domain_mapping=use_mapping,
                         sample_domain_mapping=pass1_domain_mapping,
                     )
                     copy_interpolator_provenance(interp1, interp2)
                     pass2_snapshot = snapshot_interpolator_state(interp2)
                     
                     if hasattr(interp2, 'create_ants'):
                         interp2.create_ants()
                    
                     # Run Pass 2
                     if hasattr(interp2, 'verbose'):
                         original_verbose = interp2.verbose
                         interp2.verbose = True

                     _run_interpolator_with_progress(
                         interp2,
                         dims,
                         iterations,
                         f"Domain {domain} - Pass 2",
                         on_first_iteration=(lambda interp=interp2, verbose=original_verbose: setattr(interp, 'verbose', verbose)) if hasattr(interp2, 'verbose') else None,
                     )
                     finalize_phase_provenance(interp2, 'Second Pass', algo_name2, pass2_snapshot)
                
                # Generate stats for Pass 2 if String Theory
                output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                _run_interpolator_statistics_with_retry(interp2, output_dir, f"{domain}_Pass2")
                last_interp = interp2

            post_process_mode = _get_domain_post_process_mode(post_process_config, domain)
            post_process_snapshot = snapshot_interpolator_state(last_interp)
            created, assigned = _run_domain_post_process(last_interp, dims, post_process_mode, domain_label=domain)
            if created or assigned:
                finalize_phase_provenance(last_interp, 'Post-process', 'Fill with Average', post_process_snapshot)
                print(f"Applied post-process for {domain}: created={created}, assigned={assigned}")
            
            # Print domain summary (of the last pass)
            metadata = last_interp.get_metadata()
            print(f"\n=== Domain {domain} Summary ===")
            for key, value in metadata.items():
                print(f"{key}: {value}")
    else:
        # Single interpolator
        interpolator = blocks._ant_colony
        algo_name = interpolator.get_algorithm_name()
        seed_original_sample_provenance(interpolator)
        first_pass_snapshot = snapshot_interpolator_state(interpolator)
        _run_interpolator_with_progress(
            interpolator,
            dims,
            iterations,
            f"Interpolation ({algo_name})",
        )
        finalize_phase_provenance(interpolator, 'First Pass', algo_name, first_pass_snapshot)
        
        metadata = interpolator.get_metadata()
        print(f"\n=== Summary ===")
        for key, value in metadata.items():
            print(f"{key}: {value}")
            
        # Generate stats if String Theory
        output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
        _run_interpolator_statistics_with_retry(interpolator, output_dir, "Global")
    
    # Export results (handles both single and multiple interpolators)
    interpolation_file = export_blocks_to_file(blocks, interpolation_file)
    if block_evaluated_samples_file:
        export_block_evaluated_samples_to_csv(blocks, block_evaluated_samples_file)

def load_lfc_colormap(lfc_file):
    """Load a Leapfrog .lfc file returning (ListedColormap, boundaries, labels).
    boundaries: list of numeric boundary values (length n+1 for n colors) if present.
    labels: list of class labels (optional, fallback to range strings)."""
    colormap = []
    boundaries = []  # list of (min,max) when explicit ranges provided
    labels = []
    thresholds = []  # list of threshold values when file uses cumulative <value> entries
    try:
        if not os.path.exists(lfc_file):
            print(f"LFC file does not exist: {lfc_file}")
            return ListedColormap([]), [], []
        tree = ET.parse(lfc_file)
        root = tree.getroot()
        ranges = root.find('ranges')
        if ranges is None:
            print("No <ranges> element found in LFC file.")
            return ListedColormap([]), [], []
        entries = ranges.findall('entry')
        # Determine style: threshold style (<value>, <equal>) OR explicit min/max
        for entry in entries:
            colour_elem = entry.find('colour')
            if colour_elem is None:
                continue
            colour_text = colour_elem.text.strip()
            colour = tuple(map(float, colour_text.split()))
            colormap.append(colour)
            # Check for threshold style
            value_elem = entry.find('value')
            min_elem = entry.find('min_value')
            max_elem = entry.find('max_value')
            label_elem = entry.find('label')
            lbl = label_elem.text.strip() if label_elem is not None and label_elem.text else None
            labels.append(lbl)
            if value_elem is not None:
                # Accumulate threshold boundary
                try:
                    thresholds.append(float(value_elem.text.strip()))
                except Exception:
                    pass
            else:
                # Range style
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
        # Normalize colors: Leapfrog may export 0-255 or 0-1 floats; detect range
        flat = [c for tpl in colormap for c in tpl]
        max_val = max(flat) if flat else 1.0
        scale = 255.0 if max_val > 1.0 else 1.0
        normalized = []
        for rgba in colormap:
            # Ensure length 4
            if len(rgba) == 3:
                rgba = (*rgba, 1.0)
            elif len(rgba) > 4:
                rgba = rgba[:4]
            r,g,b,a = rgba
            normalized.append((r/scale, g/scale, b/scale, a if scale==1.0 else a/scale if a>1 else a))
        print(f"Successfully loaded LFC colormap with {len(normalized)} colors. First 3: {normalized[:3]}")
        # Build numeric_edges
        if thresholds:  # threshold style: thresholds define upper bounds of successive classes
            numeric_edges = sorted(set(thresholds))
        else:
            numeric_edges = []
            for idx,(vmin,vmax) in enumerate(boundaries):
                if vmin is not None:
                    if not numeric_edges:
                        numeric_edges.append(vmin)
                if vmax is not None:
                    numeric_edges.append(vmax)
            numeric_edges = sorted(set(numeric_edges))
        # Sync labels: replace None with generated
        final_labels = []
        for i,lbl in enumerate(labels):
            if lbl:
                final_labels.append(lbl)
            else:
                if i < len(boundaries):
                    vmin,vmax = boundaries[i]
                    final_labels.append(f"{vmin} - {vmax}")
                else:
                    final_labels.append(f"Class {i}")
        print(f"Successfully loaded LFC colormap with {len(normalized)} colors. Threshold style: {bool(thresholds)}. First labels: {final_labels[:3]}")
        # If threshold style used, we return numeric_edges as thresholds; caller will treat as cumulative
        return ListedColormap(normalized), numeric_edges, final_labels
    except Exception as e:
        print(f"Error loading colormap from {lfc_file}: {e}")
    return ListedColormap([]), [], []


def build_taichi_viewer_state_from_config(config):
    interpolation_file = str(config.get('interpolation_file') or '').strip()
    if config.get('_prefer_interpolation_file_for_viewer') and interpolation_file and os.path.isfile(interpolation_file):
        return build_taichi_viewer_state_from_interpolation_file(config)

    taichi_runtime = _load_taichi_runtime_module()

    samples_file = config['samples_file']
    print(f"Loading sample file from {samples_file}...")

    wants_st_domain = bool(
        config.get('algorithm') in ('string_theory', 'net_connector')
        and str(config.get('string_theory_params', {}).get('interpolate_target', 'value')).strip().lower() == 'domain'
    )
    wants_ant_domain = bool(
        config.get('algorithm') == 'ant_colony'
        and str(config.get('ant_colony_interpolate_target', 'value')).strip().lower() == 'domain'
    )
    wants_domain_any = wants_st_domain or wants_ant_domain
    sample_filters = get_configured_sample_filters(config)
    block_filters = get_configured_block_filters(config)
    needs_sample_domains = should_resolve_sample_domains_for_interpolation(
        wants_domain_any,
        blocks_file=config.get('blocks_file'),
        block_domain_col=config.get('block_domain_col'),
    )

    df, parsed_cols, explicit_sample_map = load_samples_dataframe(
        samples_file,
        samples_delimiter=config.get('samples_delimiter'),
        samples_header_line=config.get('samples_header_line', 1),
        sample_x_col=config.get('sample_x_col'),
        sample_y_col=config.get('sample_y_col'),
        sample_z_col=config.get('sample_z_col'),
        sample_value_col=config.get('sample_value_col'),
        sample_domain_col=config.get('sample_domain_col'),
        sample_filters=sample_filters,
        progress_label='Reading sample file',
        extra_columns=[config.get('sample_weight_col')],
    )
    if parsed_cols is not None:
        print(f"Samples file (custom header line {config.get('samples_header_line', 1)}) parsed columns: {parsed_cols}")
    elif hasattr(df, '_detected_delimiter'):
        print(f"Samples file delimiter used: '{df._detected_delimiter}'")

    if needs_sample_domains and explicit_sample_map:
        if config.get('samples_header_line', 1) and config.get('samples_header_line', 1) != 1 and config.get('samples_delimiter'):
            df, parsed_cols = read_csv_with_selected_header(
                samples_file,
                config.get('samples_delimiter'),
                config.get('samples_header_line', 1),
                expected_min_cols=4,
                progress_label='Reading sample file',
            )
            print(f"Samples file (custom header line {config.get('samples_header_line', 1)}) parsed columns: {parsed_cols}")
        else:
            df = read_autodetect_csv(
                samples_file,
                forced_delimiter=config.get('samples_delimiter'),
                progress_label='Reading sample file',
            )
            print(f"Samples file delimiter used: '{df._detected_delimiter}'")
        if sample_filters:
            df, _ = apply_sample_filters(df, sample_filters=sample_filters)
        explicit_sample_map = None

    if needs_sample_domains:
        df = normalize_selected_sample_domain_column(df, sample_domain_col=config.get('sample_domain_col'))
    df = normalize_selected_sample_weight_column(
        df,
        sample_weight_col=config.get('sample_weight_col'),
        sample_value_col=config.get('sample_value_col'),
    )

    if explicit_sample_map:
        print(f"Applied user sample column mapping: {explicit_sample_map}")
    elif config.get('sample_x_col') and config.get('sample_y_col') and config.get('sample_z_col') and config.get('sample_value_col'):
        rename_map = {
            config.get('sample_x_col'): 'x',
            config.get('sample_y_col'): 'y',
            config.get('sample_z_col'): 'z',
            config.get('sample_value_col'): 'Value',
        }
        df = df.rename(columns=rename_map)

    df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
    nan_before = int(df['Value'].isna().sum())

    if needs_sample_domains:
        df, domain_resolution = ensure_sample_domains_for_domain_operations(
            df,
            sample_domain_col=config.get('sample_domain_col'),
            blank_domain_behavior=config.get('blank_sample_domain_behavior', 'skip'),
            x_col='x',
            y_col='y',
            z_col='z',
            blocks_file=config.get('blocks_file'),
            blocks_delimiter=config.get('blocks_delimiter'),
            blocks_header_line=config.get('blocks_header_line', 1),
            block_x_col=config.get('block_x_col'),
            block_y_col=config.get('block_y_col'),
            block_z_col=config.get('block_z_col'),
            block_domain_col=config.get('block_domain_col'),
            block_size=config.get('block_size'),
            block_filters=block_filters,
        )
        df['Domain'] = df['Domain'].astype(str).str.strip()
        blank_domain = df['Domain'].isna() | (df['Domain'].str.strip() == '') | (df['Domain'].str.lower() == 'nan')
        blank_count = int(blank_domain.sum())
        if blank_count:
            print(f"Detected {blank_count} sample rows with blank Domain; these will be excluded from domain interpolation.")
            df = df.loc[~blank_domain].copy()

        if nan_before:
            print(f"Detected {nan_before} sample rows with blank/non-numeric Value; filling with 0.0 for domain interpolation.")
            df['Value'] = df['Value'].fillna(0.0)
    else:
        df = df.dropna(subset=['Value'])

    sample_weights = None
    selected_weight_column = str(config.get('sample_weight_col') or '').strip()
    if selected_weight_column and selected_weight_column != '(None)':
        df['Weight'] = pd.to_numeric(df['Weight'], errors='coerce')
        invalid_weight_mask = (~np.isfinite(df['Weight'])) | (df['Weight'] <= 0.0)
        invalid_weight_count = int(invalid_weight_mask.sum())
        if invalid_weight_count:
            print(f"Detected {invalid_weight_count} sample rows with blank/non-numeric/non-positive weights; these will be excluded from samples.")
            df = df.loc[~invalid_weight_mask].copy()
        if len(df) == 0:
            raise ValueError("After removing invalid sample weights, no samples remain. Check the configured 'sample_weight_col'.")
        sample_weights = df['Weight'].to_numpy(dtype=float, copy=False)

    points = df[['x', 'y', 'z']].values
    values = df['Value'].values
    sample_domains = df['Domain'].values if 'Domain' in df.columns else None
    print(f"Loaded {len(values)} samples from {samples_file}.")

    blocks = create_blocks(
        points,
        values,
        block_size=config['block_size'],
        verbose=config['verbose'],
        range_size=config['range_size'],
        max_pheromone=config['max_pheromone'],
        ants_per_sample=config['ants_per_sample'],
        blocks_file=config['blocks_file'],
        background_value=config['background_value'],
        background_distance=config['background_distance'],
        average_with_blocks=config['average_with_blocks'],
        blocks_delimiter=config.get('blocks_delimiter'),
        avoid_visited_threshold_enabled=config.get('avoid_visited_threshold_enabled', False),
        avoid_visited_threshold=config.get('avoid_visited_threshold', 100),
        blocks_header_line=config.get('blocks_header_line', 1),
        block_x_col=config.get('block_x_col'),
        block_y_col=config.get('block_y_col'),
        block_z_col=config.get('block_z_col'),
        block_domain_col=config.get('block_domain_col'),
        config=config,
        sample_domains=sample_domains,
        sample_weights=sample_weights,
    )

    lfc_colormap, lfc_bins, lfc_labels = load_lfc_colormap(config.get('color_file'))
    lfc_colors = [tuple(map(float, color)) for color in getattr(lfc_colormap, 'colors', [])]
    tick_labels = lfc_labels or taichi_runtime.build_lfc_tick_labels(lfc_colors, lfc_bins)

    return {
        'sample_points': np.asarray(points, dtype=np.float32),
        'sample_values': np.asarray(values, dtype=np.float32),
        'sample_domains': sample_domains,
        'blocks_data': blocks,
        'block_size': config['block_size'],
        'value_filter': config.get('value_filter', 0.0),
        'lfc_colors': lfc_colors,
        'lfc_bins': lfc_bins,
        'lfc_tick_labels': tick_labels,
        'config': config,
        'window_title': 'Anterpolator Taichi Viewer',
    }


def build_taichi_viewer_state_from_interpolation_file(config):
    taichi_runtime = _load_taichi_runtime_module()

    samples_file = config['samples_file']
    interpolation_file = str(config.get('interpolation_file') or '').strip()
    print(f"Loading interpolation viewer state from {interpolation_file}...")

    sample_filters = get_configured_sample_filters(config)
    df_samples, _, explicit_sample_map = load_samples_dataframe(
        samples_file,
        samples_delimiter=config.get('samples_delimiter'),
        samples_header_line=config.get('samples_header_line', 1),
        sample_x_col=config.get('sample_x_col'),
        sample_y_col=config.get('sample_y_col'),
        sample_z_col=config.get('sample_z_col'),
        sample_value_col=config.get('sample_value_col'),
        sample_domain_col=config.get('sample_domain_col'),
        sample_filters=sample_filters,
        progress_label='Reading sample file',
    )

    if explicit_sample_map:
        pass
    elif config.get('sample_x_col') and config.get('sample_y_col') and config.get('sample_z_col') and config.get('sample_value_col'):
        df_samples = df_samples.rename(columns={
            config.get('sample_x_col'): 'x',
            config.get('sample_y_col'): 'y',
            config.get('sample_z_col'): 'z',
            config.get('sample_value_col'): 'Value',
        })

    sample_coord_frame = df_samples[['x', 'y', 'z']].apply(pd.to_numeric, errors='coerce')
    sample_values_series = pd.to_numeric(df_samples['Value'], errors='coerce')
    sample_valid_mask = sample_coord_frame.notna().all(axis=1) & sample_values_series.notna()
    sample_points = sample_coord_frame.loc[sample_valid_mask].to_numpy(dtype=np.float32, copy=False)
    sample_values = sample_values_series.loc[sample_valid_mask].to_numpy(dtype=np.float32, copy=False)
    sample_domains = None
    if 'Domain' in df_samples.columns:
        sample_domains = df_samples.loc[sample_valid_mask, 'Domain'].astype(str).to_numpy(dtype=object, copy=False)

    df_blocks = read_autodetect_csv(
        interpolation_file,
        progress_label='Reading interpolation file',
    )
    required_columns = {'x', 'y', 'z', 'Value'}
    missing_columns = [column for column in required_columns if column not in df_blocks.columns]
    if missing_columns:
        raise ValueError(
            f"Interpolation file is missing required columns: {', '.join(missing_columns)}"
        )

    block_coord_frame = df_blocks[['x', 'y', 'z']].apply(pd.to_numeric, errors='coerce')
    block_values_series = pd.to_numeric(df_blocks['Value'], errors='coerce')
    block_valid_mask = block_coord_frame.notna().all(axis=1) & block_values_series.notna()
    df_blocks = df_blocks.loc[block_valid_mask].copy()
    block_points = block_coord_frame.loc[block_valid_mask].to_numpy(dtype=np.float32, copy=False)
    block_values = block_values_series.loc[block_valid_mask].to_numpy(dtype=np.float32, copy=False)

    sample_row_mask = np.zeros(len(df_blocks), dtype=bool)
    if 'Is_Sample' in df_blocks.columns:
        sample_tokens = df_blocks['Is_Sample'].astype(str).str.strip().str.lower()
        sample_row_mask = sample_tokens.isin({'true', '1', 'yes'}).to_numpy(dtype=bool, copy=False)

    block_positions = {
        tuple(point.astype(np.float32)): float(value)
        for point, value in zip(block_points, block_values)
    }
    block_domains = {}
    if 'Domain' in df_blocks.columns:
        for point, domain in zip(block_points, df_blocks['Domain'].fillna('Undomained').astype(str)):
            block_domains[tuple(point.astype(np.float32))] = str(domain).strip() or 'Undomained'

    sample_block_values = {
        tuple(point.astype(np.float32)): float(value)
        for point, value, is_sample in zip(block_points, block_values, sample_row_mask)
        if is_sample
    }

    class SimpleBlocks:
        pass

    class StaticInterpolator:
        def __init__(self, values, domains):
            self._values = values
            self.domain_mapping = domains
            self.interpolation_target = 'value'

        def get_interpolated_values(self):
            return self._values

        def run_iteration(self, dims):
            return False

        def is_converged(self):
            return True

    blocks_data = SimpleBlocks()
    blocks_data._sample_blocks = sample_block_values
    blocks_data._block_info = {
        'min_bounds': np.min(block_points, axis=0) if len(block_points) else (np.min(sample_points, axis=0) if len(sample_points) else np.zeros(3, dtype=np.float32)),
        'block_size': list(config['block_size']),
        'domain_mapping': block_domains,
        'positions_are_world': True,
    }
    blocks_data._ant_colony = StaticInterpolator(block_positions, block_domains)

    lfc_colormap, lfc_bins, lfc_labels = load_lfc_colormap(config.get('color_file'))
    lfc_colors = [tuple(map(float, color)) for color in getattr(lfc_colormap, 'colors', [])]
    tick_labels = lfc_labels or taichi_runtime.build_lfc_tick_labels(lfc_colors, lfc_bins)

    return {
        'sample_points': np.asarray(sample_points, dtype=np.float32),
        'sample_values': np.asarray(sample_values, dtype=np.float32),
        'sample_domains': None if sample_domains is None else np.asarray(sample_domains, dtype=object),
        'blocks_data': blocks_data,
        'block_size': config['block_size'],
        'value_filter': config.get('value_filter', 0.0),
        'lfc_colors': lfc_colors,
        'lfc_bins': lfc_bins,
        'lfc_tick_labels': tick_labels,
        'config': config,
        'window_title': 'Anterpolator Taichi Viewer',
    }


def build_taichi_viewer_state_from_existing_state(config, existing_state):
    taichi_runtime = _load_taichi_runtime_module()
    lfc_colormap, lfc_bins, lfc_labels = load_lfc_colormap(config.get('color_file'))
    lfc_colors = [tuple(map(float, color)) for color in getattr(lfc_colormap, 'colors', [])]
    tick_labels = lfc_labels or taichi_runtime.build_lfc_tick_labels(lfc_colors, lfc_bins)

    merged_config = dict(existing_state.get('config', {}) or {})
    merged_config.update(config)

    sample_values = existing_state.get('sample_values')
    sample_domains = existing_state.get('sample_domains')

    return {
        'sample_points': np.asarray(existing_state['sample_points'], dtype=np.float32),
        'sample_values': None if sample_values is None else np.asarray(sample_values, dtype=np.float32),
        'sample_domains': None if sample_domains is None else np.asarray(sample_domains, dtype=object),
        'blocks_data': existing_state['blocks_data'],
        'block_size': np.asarray(existing_state['block_size'], dtype=np.float32),
        'value_filter': config.get('value_filter', existing_state.get('value_filter', 0.0)),
        'lfc_colors': lfc_colors,
        'lfc_bins': lfc_bins,
        'lfc_tick_labels': tick_labels,
        'config': merged_config,
        'window_title': existing_state.get('window_title', 'Anterpolator Taichi Viewer'),
    }


def launch_taichi_viewer_with_state(viewer_state, external_state_callback=None):
    taichi_runtime = _load_taichi_runtime_module()
    viewer = taichi_runtime.TaichiInterpolationViewer(
        external_state_callback=external_state_callback,
        **viewer_state,
    )
    viewer.run()


def launch_taichi_viewer_from_config(config, external_state_callback=None):
    viewer_state = build_taichi_viewer_state_from_config(config)
    launch_taichi_viewer_with_state(viewer_state, external_state_callback=external_state_callback)

def load_and_visualize_samples(samples_file, block_size=10, value_filter=60, verbose=False, iterations=100, range_size=10, max_pheromone=150, ants_per_sample=3, blocks_file=None, color_file=None, background_value=0.0, background_distance=None, average_with_blocks=False,
                               samples_delimiter=None, blocks_delimiter=None, fill_unvisited_domainwise=False,
                               avoid_visited_threshold_enabled=False,
                               avoid_visited_threshold=100,
                               samples_header_line=1,
                               sample_x_col=None, sample_y_col=None, sample_z_col=None, sample_value_col=None,
                               blocks_header_line=1,
                               block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                               config=None):
    pv = _require_pyvista()
    try:
        print(f"Loading sample file from {samples_file}...")

        wants_st_domain = bool(
            config
            and config.get('algorithm') in ('string_theory', 'net_connector')
            and str(config.get('string_theory_params', {}).get('interpolate_target', 'value')).strip().lower() == 'domain'
        )
        wants_ant_domain = bool(
            config
            and config.get('algorithm') == 'ant_colony'
            and str(config.get('ant_colony_interpolate_target', 'value')).strip().lower() == 'domain'
        )
        wants_domain_any = wants_st_domain or wants_ant_domain
        sample_filters = get_configured_sample_filters(config)
        block_filters = get_configured_block_filters(config)
        needs_sample_domains = should_resolve_sample_domains_for_interpolation(
            wants_domain_any,
            blocks_file=blocks_file,
            block_domain_col=block_domain_col,
        )

        df, parsed_cols, explicit_sample_map = load_samples_dataframe(
            samples_file,
            samples_delimiter=samples_delimiter,
            samples_header_line=samples_header_line,
            sample_x_col=sample_x_col,
            sample_y_col=sample_y_col,
            sample_z_col=sample_z_col,
            sample_value_col=sample_value_col,
            sample_domain_col=config.get('sample_domain_col') if config else None,
            sample_filters=sample_filters,
            progress_label='Reading sample file',
            extra_columns=[config.get('sample_weight_col') if config else None],
        )
        if parsed_cols is not None:
            print(f"Samples file (custom header line {samples_header_line}) parsed columns: {parsed_cols}")
        elif hasattr(df, '_detected_delimiter'):
            print(f"Samples file delimiter used: '{df._detected_delimiter}'")

        if wants_domain_any and explicit_sample_map:
            if samples_header_line and samples_header_line != 1 and samples_delimiter:
                df, parsed_cols = read_csv_with_selected_header(
                    samples_file,
                    samples_delimiter,
                    samples_header_line,
                    expected_min_cols=4,
                    progress_label='Reading sample file',
                )
                print(f"Samples file (custom header line {samples_header_line}) parsed columns: {parsed_cols}")
            else:
                df = read_autodetect_csv(
                    samples_file,
                    forced_delimiter=samples_delimiter,
                    progress_label='Reading sample file',
                )
                print(f"Samples file delimiter used: '{df._detected_delimiter}'")
            if sample_filters:
                df, _ = apply_sample_filters(df, sample_filters=sample_filters)
            explicit_sample_map = None

        if wants_domain_any:
            df = normalize_selected_sample_domain_column(df, sample_domain_col=config.get('sample_domain_col') if config else None)
        df = normalize_selected_sample_weight_column(
            df,
            sample_weight_col=config.get('sample_weight_col') if config else None,
            sample_value_col=sample_value_col,
        )

        if explicit_sample_map:
            print(f"Applied user sample column mapping: {explicit_sample_map}")
        elif sample_x_col and sample_y_col and sample_z_col and sample_value_col:
            for chosen in [sample_x_col, sample_y_col, sample_z_col, sample_value_col]:
                if chosen not in df.columns:
                    raise ValueError(f"Selected samples column '{chosen}' not present in file.")
            rename_map = {sample_x_col: 'x', sample_y_col: 'y', sample_z_col: 'z', sample_value_col: 'Value'}
            df = df.rename(columns=rename_map)
            print(f"Applied user sample column mapping: {rename_map}")
        else:
            expected = ['x', 'y', 'z', 'Value']
            missing_any = any(col not in df.columns for col in expected)
            if missing_any:
                if 'Value' not in df.columns:
                    for c in df.columns:
                        if c.lower() == 'value':
                            df = df.rename(columns={c: 'Value'})
                            break
                if any(col not in df.columns for col in expected):
                    first_four = list(df.columns[:4])
                    if len(first_four) == 4:
                        rename_map = {orig: new for orig, new in zip(first_four, expected)}
                        df = df.rename(columns=rename_map)
                        print(f"Mapped first four sample columns to {expected}: {rename_map}")
                    else:
                        raise ValueError("Samples file must have at least four columns for automatic mapping (x,y,z,Value).")
        for col in ['x','y','z','Value']:
            if col not in df.columns:
                raise ValueError(f"Required column '{col}' not found after mapping.")

        # Optional Domain column (string). Used by String Theory domain interpolation.
        if needs_sample_domains:
            df, domain_resolution = ensure_sample_domains_for_domain_operations(
                df,
                sample_domain_col=config.get('sample_domain_col') if config else None,
                blank_domain_behavior=config.get('blank_sample_domain_behavior', 'skip') if config else 'skip',
                x_col='x',
                y_col='y',
                z_col='z',
                blocks_file=blocks_file,
                blocks_delimiter=blocks_delimiter,
                blocks_header_line=blocks_header_line,
                block_x_col=block_x_col,
                block_y_col=block_y_col,
                block_z_col=block_z_col,
                block_domain_col=block_domain_col,
                block_size=block_size,
                block_filters=block_filters,
            )
            df['Domain'] = df['Domain'].astype(str).str.strip()
            blank_domain = df['Domain'].isna() | (df['Domain'].str.strip() == '') | (df['Domain'].str.lower() == 'nan')
            blank_count = int(blank_domain.sum())
            if blank_count:
                print(f"Detected {blank_count} sample rows with blank Domain; these will be excluded from domain interpolation.")
                df = df.loc[~blank_domain].copy()
        print(f"Samples final headers: {list(df.columns)}")
        # --- Ensure the selected Value column is handled ---
        # For domain interpolation, Value is not used; keep rows even if Value is non-numeric.
        try:
            df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
        except Exception as e:
            print(f"Warning: could not coerce 'Value' column to numeric directly: {e}")

        nan_before = int(df['Value'].isna().sum())
        if wants_domain_any:
            if nan_before:
                print(f"Detected {nan_before} sample rows with blank/non-numeric Value; filling with 0.0 for domain interpolation.")
                df['Value'] = df['Value'].fillna(0.0)
        else:
            if nan_before:
                print(f"Detected {nan_before} sample rows with blank/non-numeric values in the selected value column; these will be excluded from samples.")
                df = df.dropna(subset=['Value'])
            if len(df) == 0:
                raise ValueError("After removing blank/non-numeric sample values, no samples remain. Check the configured 'sample_value_col'.")

        sample_weights = None
        selected_weight_column = str(config.get('sample_weight_col') if config else '').strip()
        if selected_weight_column and selected_weight_column != '(None)':
            df['Weight'] = pd.to_numeric(df['Weight'], errors='coerce')
            invalid_weight_mask = (~np.isfinite(df['Weight'])) | (df['Weight'] <= 0.0)
            invalid_weight_count = int(invalid_weight_mask.sum())
            if invalid_weight_count:
                print(f"Detected {invalid_weight_count} sample rows with blank/non-numeric/non-positive weights; these will be excluded from samples.")
                df = df.loc[~invalid_weight_mask].copy()
            if len(df) == 0:
                raise ValueError("After removing invalid sample weights, no samples remain. Check the configured 'sample_weight_col'.")
            sample_weights = df['Weight'].to_numpy(dtype=float, copy=False)

        print(f"Cleaned sample count: {len(df)} (removed {nan_before} blank value rows)")
        points = df[['x','y','z']].values
        values = df['Value'].values
        sample_domains = df['Domain'].values if 'Domain' in df.columns else None
        # Store scalar range early
        data_min, data_max = float(np.nanmin(values)), float(np.nanmax(values))
        print(f"Sample value range: min={data_min} max={data_max}")

        # Show the number of samples loaded
        print(f"Loaded {len(values)} samples from {samples_file}.")

        # Display the name of the blocks_file being loaded, if any
        if blocks_file:
            print(f"Loading blocks file from {blocks_file}...")

        # Load colormap from LFC file if provided
        lfc_boundaries = []
        lfc_labels = []
        if color_file and os.path.exists(color_file):
            print(f"Loading colormap from LFC file: {color_file}", flush=True)
            discrete_cmap, lfc_boundaries, lfc_labels = load_lfc_colormap(color_file)
            if isinstance(discrete_cmap, ListedColormap) and lfc_boundaries:
                print(f"Applied LFC discrete boundaries: {lfc_boundaries[:10]}")
        else:
            print("LFC file not provided or not found; using default 'rainbow'", flush=True)
            discrete_cmap = 'rainbow'

        plotter = pv.Plotter()
        if hasattr(pv, 'set_new_attribute'):
            if not hasattr(plotter, 'pickpoint'):
                pv.set_new_attribute(plotter, 'pickpoint', None)
            if not hasattr(plotter, 'picked_point'):
                pv.set_new_attribute(plotter, 'picked_point', None)
        plotter.set_background('white')

        # Store colormap in the plotter for reuse
        plotter._colormap = discrete_cmap

        # Create blocks and store in plotter using new parameters
        blocks = create_blocks(
            points,
            values,
            block_size,
            verbose,
            range_size,
            max_pheromone,
            ants_per_sample,
            blocks_file,
            background_value=background_value,  # Pass background_value
            background_distance=background_distance,  # Pass background_distance
            average_with_blocks=average_with_blocks,  # Pass average_with_blocks
            blocks_delimiter=blocks_delimiter,
            avoid_visited_threshold_enabled=avoid_visited_threshold_enabled,
            avoid_visited_threshold=avoid_visited_threshold,
            blocks_header_line=blocks_header_line,
            block_x_col=block_x_col, block_y_col=block_y_col, block_z_col=block_z_col, block_domain_col=block_domain_col,
            config=config,
            sample_domains=sample_domains,
            sample_weights=sample_weights,
        )

        # Store settings in plotter
        plotter._blocks_data = blocks
        plotter._verbose = verbose
        plotter._value_filter = value_filter  # Store value_filter
        plotter._domain_post_process_overrides = dict((config or {}).get('domain_algorithm_overrides', {}))
        plotter._avoid_visited_threshold_enabled = avoid_visited_threshold_enabled
        plotter._avoid_visited_threshold = avoid_visited_threshold
        plotter._block_lookup = None
        plotter._visible_blocks = None
        plotter._visible_positions = None
        plotter._blocks_actor = None
        plotter._blocks_display_prepared = False
        plotter._block_evaluated_samples_file = resolve_block_evaluated_samples_export_path(
            config.get('export_block_evaluated_samples', False) if config else False,
            config.get('block_evaluated_samples_file') if config else None,
            interpolation_file=config.get('interpolation_file') if config else None,
            samples_file=samples_file,
        )

        # Add point cloud with scalar bar settings
        cloud = pv.PolyData(points)
        # If we have discrete LFC boundaries treat values as classes via digitize
        if isinstance(plotter._colormap, ListedColormap) and lfc_boundaries:
            # Two possible styles:
            # 1) Threshold style: lfc_boundaries are ascending thresholds; colors = len(thresholds)+1 (end colour appended)
            # 2) Range edge style: lfc_boundaries are full set of edges (len = n+1 for n colors)
            bins = np.array(sorted(lfc_boundaries))
            n_colors = len(plotter._colormap.colors)
            threshold_style = (n_colors == len(bins) + 1)
            if threshold_style:
                # Class index determined by number of thresholds value exceeds
                # Use searchsorted to find insertion index; no minus 1 offset
                value_indices = np.searchsorted(bins, values, side='right')
                value_indices = np.clip(value_indices, 0, n_colors-1)
            else:
                # Treat bins as edges; need edges length = n_colors+1 ideally
                if len(bins) == n_colors + 1:
                    # np.digitize returns indices 1..n_colors for values within edges
                    value_indices = np.digitize(values, bins, right=False) - 1
                else:
                    # Fallback: linear mapping
                    value_indices = np.digitize(values, bins, right=False) - 1
                value_indices = np.clip(value_indices, 0, n_colors-1)
            cloud['Values'] = value_indices
            plotter._value_is_indexed = True
            plotter._lfc_bins = bins.tolist()
            plotter._lfc_labels = lfc_labels
            # Build tick labels strictly from thresholds/edges (numeric ranges)
            tick_labels = []
            if threshold_style:
                for i in range(n_colors):
                    if i == 0:
                        rng = f"< {bins[0]}"
                    elif i < len(bins):
                        rng = f"{bins[i-1]} - {bins[i]}"
                    else:
                        rng = f">= {bins[-1]}"
                    tick_labels.append(rng)
            else:
                if len(bins) == n_colors + 1:
                    for i in range(n_colors):
                        rng = f"{bins[i]} - {bins[i+1]}"
                        tick_labels.append(rng)
                else:
                    for i in range(n_colors):
                        tick_labels.append(f"Class {i}")
            plotter._lfc_tick_labels = tick_labels
            # Print textual legend
            print("\nLFC Discrete Legend:")
            colors = plotter._colormap.colors
            def fmt_rgb(c):
                r,g,b,a = c
                return f"RGB({r:.3f},{g:.3f},{b:.3f})"
            for i in range(n_colors):
                rng = tick_labels[i] if i < len(tick_labels) else f"Class {i}"
                print(f"Class {i}: {rng} | {fmt_rgb(colors[i])}")
        else:
            cloud['Values'] = values
        sargs = dict(
            vertical=True,
            position_x=0.93,
            position_y=0.2,
            height=0.6,
            label_font_size=10,
            title='Value'
        )
        if getattr(plotter, '_value_is_indexed', False):
            # Discrete scalar bar: category labels at integer indices
            sargs_indexed = dict(sargs)
            sargs_indexed['title'] = 'Value'
            sargs_indexed['n_labels'] = 0
            sargs_indexed['fmt'] = ''
            annotations = {i: plotter._lfc_tick_labels[i] for i in range(len(plotter._lfc_tick_labels))}
            plotter.add_mesh(
                cloud,
                render_points_as_spheres=True,
                point_size=5,
                scalars='Values',
                cmap=plotter._colormap,
                categories=True,
                annotations=annotations,
                clim=[-0.5, len(plotter._colormap.colors)-0.5],
                scalar_bar_args=sargs_indexed,
            )
        else:
            plotter.add_mesh(
                cloud,
                render_points_as_spheres=True,
                point_size=5,
                scalars='Values',
                cmap=plotter._colormap,
                clim=[data_min, data_max],
                scalar_bar_args=sargs,
            )

        plotter._scalar_range = (data_min, data_max)
        if len(blocks) <= INITIAL_BLOCK_RENDER_THRESHOLD:
            _ensure_blocks_actor(plotter, visible=True)
        else:
            print(
                f"Deferring initial block rendering for {len(blocks):,} blocks. "
                "Press 'b' after the window opens to build/show block meshes."
            )

        # Add controls
        plotter.add_key_event('b', lambda: toggle_blocks(plotter))
        plotter.add_key_event('i', lambda: update_interpolation(plotter))
        # Attach interpolation_file path if defined globally
        if 'interpolation_file' in globals():
            plotter._interpolation_file = globals()['interpolation_file']
        plotter.add_key_event(
            'I',
            lambda: silent_interpolation(
                plotter, iterations, getattr(plotter, '_interpolation_file', 'interpolation.csv')
            ),  # Shift+I
        )

        # Store iterations setting
        plotter._iterations = iterations

        # Show visualization
        plotter.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == '--launch-viewer-config':
        config_path = sys.argv[2]
        with open(config_path, 'r', encoding='utf-8') as f:
            cfg = json.load(f)
        cfg['viewer_backend'] = _normalize_viewer_backend(cfg.get('viewer_backend'))
        cfg.pop('_viewer_reload_mode', None)
        last_mtime = [os.path.getmtime(config_path) if os.path.isfile(config_path) else 0.0]
        current_state = [build_taichi_viewer_state_from_config(cfg)]

        def reload_callback():
            try:
                current_mtime = os.path.getmtime(config_path)
            except OSError:
                return None
            if current_mtime <= last_mtime[0]:
                return None
            try:
                with open(config_path, 'r', encoding='utf-8') as f:
                    updated_cfg = json.load(f)
                updated_cfg['viewer_backend'] = _normalize_viewer_backend(updated_cfg.get('viewer_backend'))
                reload_mode = str(updated_cfg.pop('_viewer_reload_mode', 'refresh')).strip().lower()
                if reload_mode == 'reload':
                    state = build_taichi_viewer_state_from_config(updated_cfg)
                else:
                    state = build_taichi_viewer_state_from_existing_state(updated_cfg, current_state[0])
                last_mtime[0] = current_mtime
                current_state[0] = state
                return state
            except Exception:
                print('Failed to live-reload viewer config:')
                traceback.print_exc()
                return None

        launch_taichi_viewer_with_state(current_state[0], external_state_callback=reload_callback)
        sys.exit(0)

    class DomainAlgorithmDialog(QtWidgets.QDialog):
        """Dialog for configuring algorithm per domain"""
        def __init__(self, domains, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Domain Algorithm Mapping")
            self.resize(800, 500)
            
            self.domain_configs = {}  # domain -> {'algorithm': str, 'second_pass_algorithm': str, 'post_process': str, 'skip': bool}
            
            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)
            
            # Info label
            info = QtWidgets.QLabel("Configure which algorithm to use for each domain.\n"
                                  "You can configure a second pass to run after the first one completes.\n"
                                  "The second pass uses the output of the first pass as input.\n"
                                  "Post-process runs once after the last enabled pass for that domain.")
            info.setWordWrap(True)
            layout.addWidget(info)
            
            # Table for domain mappings
            self.table = QtWidgets.QTableWidget()
            self.table.setColumnCount(4)
            self.table.setHorizontalHeaderLabels(['Domain', 'First Pass Algorithm', 'Second Pass Algorithm', 'Post-process'])
            self.table.horizontalHeader().setStretchLastSection(True)
            self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(3, QtWidgets.QHeaderView.Stretch)
            layout.addWidget(self.table)
            
            self.populate_domains(domains)
            
            # Buttons
            btn_layout = QtWidgets.QHBoxLayout()
            self.apply_all_btn = QtWidgets.QPushButton('Apply First Pass to All')
            self.apply_all_btn.clicked.connect(self.apply_to_all)
            btn_layout.addWidget(self.apply_all_btn)
            self.apply_second_pass_all_btn = QtWidgets.QPushButton('Apply Second Pass to All')
            self.apply_second_pass_all_btn.clicked.connect(self.apply_second_pass_to_all)
            btn_layout.addWidget(self.apply_second_pass_all_btn)
            self.apply_post_process_all_btn = QtWidgets.QPushButton('Apply Post Process to All')
            self.apply_post_process_all_btn.clicked.connect(self.apply_post_process_to_all)
            btn_layout.addWidget(self.apply_post_process_all_btn)
            btn_layout.addStretch()
            
            self.ok_btn = QtWidgets.QPushButton('OK')
            self.cancel_btn = QtWidgets.QPushButton('Cancel')
            self.ok_btn.clicked.connect(self.accept)
            self.cancel_btn.clicked.connect(self.reject)
            btn_layout.addWidget(self.ok_btn)
            btn_layout.addWidget(self.cancel_btn)
            layout.addLayout(btn_layout)
        
        def populate_domains(self, domains):
            """Populate the table with a preloaded domain catalog."""
            self.table.setRowCount(len(domains))
            for i, domain in enumerate(domains):
                domain_item = QtWidgets.QTableWidgetItem(domain)
                domain_item.setFlags(domain_item.flags() & ~QtCore.Qt.ItemIsEditable)
                self.table.setItem(i, 0, domain_item)

                algo1_combo = QtWidgets.QComboBox()
                algo1_combo.addItems(['(use default)', 'ant_colony', 'molecular_clock', 'gaussian_kernel', 'adaptive_octree', 'string_theory', 'skip'])
                algo1_combo.setCurrentText('(use default)')
                self.table.setCellWidget(i, 1, algo1_combo)

                algo2_combo = QtWidgets.QComboBox()
                algo2_combo.addItems(['skip', 'ant_colony', 'molecular_clock', 'gaussian_kernel', 'adaptive_octree', 'string_theory'])
                algo2_combo.setCurrentText('skip')
                self.table.setCellWidget(i, 2, algo2_combo)

                post_process_combo = QtWidgets.QComboBox()
                post_process_combo.addItems(['skip', 'fill_with_average'])
                post_process_combo.setCurrentText('skip')
                self.table.setCellWidget(i, 3, post_process_combo)

                algo1_combo.currentTextChanged.connect(
                    lambda text, row=i: self.on_first_pass_changed(row, text)
                )
                algo2_combo.currentTextChanged.connect(
                    lambda text, row=i: self.on_second_pass_changed(row, text)
                )

        def on_first_pass_changed(self, row, text):
            """Handle changes to first pass algorithm"""
            algo2_combo = self.table.cellWidget(row, 2)
            post_process_combo = self.table.cellWidget(row, 3)
            if not algo2_combo or not post_process_combo:
                return
                
            if text == 'skip':
                algo2_combo.setCurrentText('skip')
                algo2_combo.setEnabled(False)
                post_process_combo.setCurrentText('skip')
                post_process_combo.setEnabled(False)
            else:
                algo2_combo.setEnabled(True)
                post_process_combo.setEnabled(True)
                current_algo2 = algo2_combo.currentText()
                if current_algo2 == text and text != 'skip':
                    algo2_combo.setCurrentText('skip')
        
        def on_second_pass_changed(self, row, text):
            """Handle changes to second pass algorithm"""
            algo1_combo = self.table.cellWidget(row, 1)
            if not algo1_combo:
                return
                
            algo1 = algo1_combo.currentText()
            
            # Don't allow same algorithm (unless skip)
            if text != 'skip' and text == algo1:
                QtWidgets.QMessageBox.warning(self, "Invalid Selection", 
                    "Second pass algorithm cannot be the same as the first pass.")
                algo2_combo = self.table.cellWidget(row, 2)
                algo2_combo.setCurrentText('skip')

        def apply_to_all(self):
            """Apply same first pass algorithm to all domains"""
            algorithms = ['(use default)', 'ant_colony', 'molecular_clock', 'gaussian_kernel', 'adaptive_octree', 'string_theory', 'skip']
            algo, ok = QtWidgets.QInputDialog.getItem(
                self, 'Apply to All', 'Select first pass algorithm for all domains:', 
                algorithms, 0, False
            )
            if ok:
                for i in range(self.table.rowCount()):
                    combo = self.table.cellWidget(i, 1)
                    if combo:
                        combo.setCurrentText(algo)

        def apply_second_pass_to_all(self):
            """Apply same second pass algorithm to all eligible domains"""
            algorithms = ['skip', 'ant_colony', 'molecular_clock', 'gaussian_kernel', 'adaptive_octree', 'string_theory']
            algo, ok = QtWidgets.QInputDialog.getItem(
                self, 'Apply to All', 'Select second pass algorithm for all domains:',
                algorithms, 0, False
            )
            if ok:
                for i in range(self.table.rowCount()):
                    algo1_combo = self.table.cellWidget(i, 1)
                    algo2_combo = self.table.cellWidget(i, 2)
                    if not algo1_combo or not algo2_combo:
                        continue
                    if algo1_combo.currentText() == 'skip':
                        algo2_combo.setCurrentText('skip')
                        continue
                    algo2_combo.setCurrentText(algo)

        def apply_post_process_to_all(self):
            """Apply same post-process option to all eligible domains"""
            options = ['skip', 'fill_with_average']
            mode, ok = QtWidgets.QInputDialog.getItem(
                self, 'Apply to All', 'Select post-process for all domains:',
                options, 0, False
            )
            if ok:
                for i in range(self.table.rowCount()):
                    algo1_combo = self.table.cellWidget(i, 1)
                    post_process_combo = self.table.cellWidget(i, 3)
                    if not algo1_combo or not post_process_combo:
                        continue
                    if algo1_combo.currentText() == 'skip':
                        post_process_combo.setCurrentText('skip')
                        continue
                    post_process_combo.setCurrentText(mode)
        
        def get_domain_configs(self):
            """Get domain algorithm configurations"""
            configs = {}
            for i in range(self.table.rowCount()):
                domain_item = self.table.item(i, 0)
                algo1_combo = self.table.cellWidget(i, 1)
                algo2_combo = self.table.cellWidget(i, 2)
                post_process_combo = self.table.cellWidget(i, 3)
                
                if domain_item and algo1_combo and algo2_combo and post_process_combo:
                    domain = domain_item.text()
                    algo1 = algo1_combo.currentText()
                    algo2 = algo2_combo.currentText()
                    post_process = _normalize_domain_post_process_mode(post_process_combo.currentText())
                    
                    config = {}
                    
                    # First Pass
                    if algo1 == 'skip':
                        config['skip'] = True
                    elif algo1 != '(use default)':
                        config['algorithm'] = algo1
                    
                    # Second Pass
                    if algo2 != 'skip':
                        config['second_pass_algorithm'] = algo2

                    if algo1 != 'skip' and post_process != 'skip':
                        config['post_process'] = post_process
                    
                    if config:
                        configs[domain] = config
            
            return configs
        
        def set_domain_configs(self, configs):
            """Set domain algorithm configurations from loaded config"""
            for i in range(self.table.rowCount()):
                domain_item = self.table.item(i, 0)
                algo1_combo = self.table.cellWidget(i, 1)
                algo2_combo = self.table.cellWidget(i, 2)
                post_process_combo = self.table.cellWidget(i, 3)
                
                if domain_item and algo1_combo and algo2_combo and post_process_combo:
                    domain = domain_item.text()
                    if domain in configs:
                        config = configs[domain]
                        
                        # Set First Pass
                        if config.get('skip', False):
                            algo1_combo.setCurrentText('skip')
                        elif 'algorithm' in config:
                            algo1_combo.setCurrentText(config['algorithm'])
                        else:
                            algo1_combo.setCurrentText('(use default)')
                        
                        # Set Second Pass
                        if 'second_pass_algorithm' in config:
                            algo2_combo.setCurrentText(config['second_pass_algorithm'])
                        else:
                            algo2_combo.setCurrentText('skip')
                        post_process_combo.setCurrentText(_normalize_domain_post_process_mode(config.get('post_process', 'skip')))
                    else:
                        algo1_combo.setCurrentText('(use default)')
                        algo2_combo.setCurrentText('skip')
                        post_process_combo.setCurrentText('skip')

                    self.on_first_pass_changed(i, algo1_combo.currentText())

        def accept(self):
            super().accept()
        
        def reject(self):
            super().reject()

    class SampleFilterEditDialog(QtWidgets.QDialog):
        def __init__(self, filter_source, filter_spec=None, parent=None, subject_label='Sample'):
            super().__init__(parent)
            self.filter_source = filter_source if isinstance(filter_source, FilterDataSource) else FilterDataSource(filter_source)
            self.subject_label = str(subject_label or 'Sample').strip() or 'Sample'
            self.setWindowTitle(f'{self.subject_label} Filter')
            self.resize(520, 420)

            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)

            form = QtWidgets.QFormLayout()
            layout.addLayout(form)

            self.field_combo = QtWidgets.QComboBox()
            self.field_combo.addItem('(Select Field)')
            self.field_combo.addItems(self.filter_source.columns)
            form.addRow('Field', self.field_combo)

            self.type_combo = QtWidgets.QComboBox()
            self.type_combo.addItems(['(Select Type)', 'categorical', 'numeric'])
            form.addRow('Type', self.type_combo)

            self.criteria_stack = QtWidgets.QStackedWidget()

            categorical_page = QtWidgets.QWidget()
            categorical_layout = QtWidgets.QVBoxLayout()
            categorical_page.setLayout(categorical_layout)
            categorical_layout.addWidget(QtWidgets.QLabel('Select one or more values:'))
            self.value_list = QtWidgets.QListWidget()
            self.value_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            categorical_layout.addWidget(self.value_list)
            self.value_hint = QtWidgets.QLabel('')
            self.value_hint.setWordWrap(True)
            categorical_layout.addWidget(self.value_hint)
            self.criteria_stack.addWidget(categorical_page)

            numeric_page = QtWidgets.QWidget()
            numeric_form = QtWidgets.QFormLayout()
            numeric_page.setLayout(numeric_form)
            self.min_value_edit = QtWidgets.QLineEdit('')
            self.max_value_edit = QtWidgets.QLineEdit('')
            self.numeric_hint = QtWidgets.QLabel('')
            self.numeric_hint.setWordWrap(True)
            numeric_form.addRow('Minimum', self.min_value_edit)
            numeric_form.addRow('Maximum', self.max_value_edit)
            numeric_form.addRow('', self.numeric_hint)
            self.criteria_stack.addWidget(numeric_page)

            form.addRow('Criteria', self.criteria_stack)
            self.criteria_stack.setEnabled(False)

            button_row = QtWidgets.QHBoxLayout()
            button_row.addStretch(1)
            ok_btn = QtWidgets.QPushButton('OK')
            cancel_btn = QtWidgets.QPushButton('Cancel')
            ok_btn.clicked.connect(self.accept)
            cancel_btn.clicked.connect(self.reject)
            button_row.addWidget(ok_btn)
            button_row.addWidget(cancel_btn)
            layout.addLayout(button_row)

            self.field_combo.currentTextChanged.connect(self._refresh_criteria_ui)
            self.type_combo.currentTextChanged.connect(self._refresh_criteria_ui)

            if filter_spec:
                field = str(filter_spec.get('field', '')).strip()
                if field:
                    self.field_combo.setCurrentText(field)
                filter_type = str(filter_spec.get('type', 'categorical')).strip().lower() or 'categorical'
                self.type_combo.setCurrentText(filter_type)
            else:
                self.value_hint.setText('Select a field and choose categorical to load available values.')
                self.numeric_hint.setText('Select a field and choose numeric to inspect the available range.')

            self._refresh_criteria_ui()

            if filter_spec:
                if self.type_combo.currentText() == 'categorical':
                    selected_values = {str(value) for value in filter_spec.get('values', [])}
                    for idx in range(self.value_list.count()):
                        item = self.value_list.item(idx)
                        item.setSelected(item.text() in selected_values)
                else:
                    min_value = filter_spec.get('min', '')
                    max_value = filter_spec.get('max', '')
                    self.min_value_edit.setText('' if min_value is None else str(min_value))
                    self.max_value_edit.setText('' if max_value is None else str(max_value))

        def _refresh_criteria_ui(self):
            field = self.field_combo.currentText()
            filter_type = self.type_combo.currentText()
            has_field = bool(field and field != '(Select Field)' and self.filter_source.has_field(field))
            has_type = filter_type in ('categorical', 'numeric')

            if not has_field or not has_type:
                self.criteria_stack.setEnabled(False)
                self.value_list.clear()
                self.value_hint.setText('Select a field and choose categorical to load available values.')
                self.numeric_hint.setText('Select a field and choose numeric to inspect the available range.')
                self.criteria_stack.setCurrentIndex(0)
                return

            self.criteria_stack.setEnabled(True)

            if filter_type == 'categorical':
                self.criteria_stack.setCurrentIndex(0)
                unique_values, total_unique_values, truncated = self.filter_source.get_categorical_values(field)
                self.value_list.clear()
                self.value_list.addItems(unique_values)
                self.value_hint.setText(
                    f'Loaded {total_unique_values:,} distinct values.' +
                    (' Showing the first 1,000 values.' if truncated else '')
                )
            else:
                self.criteria_stack.setCurrentIndex(1)
                cached_range = self.filter_source.get_numeric_range(field)
                if cached_range is None:
                    self.numeric_hint.setText('No numeric values detected in this field.')
                else:
                    min_value, max_value = cached_range
                    self.numeric_hint.setText(f'Available numeric range: {min_value:g} to {max_value:g}')

        def get_filter_spec(self):
            field = self.field_combo.currentText().strip()
            filter_type = self.type_combo.currentText().strip().lower()
            if not field or field == '(Select Field)':
                raise ValueError('Please select a field.')
            if filter_type == '(select type)' or filter_type == '':
                raise ValueError('Please select a filter type.')

            if filter_type == 'categorical':
                selected_values = [item.text() for item in self.value_list.selectedItems()]
                if not selected_values:
                    raise ValueError('Select at least one value for a categorical filter.')
                return {
                    'field': field,
                    'type': 'categorical',
                    'values': selected_values,
                }

            min_text = self.min_value_edit.text().strip()
            max_text = self.max_value_edit.text().strip()
            min_value = None if min_text == '' else float(min_text)
            max_value = None if max_text == '' else float(max_text)
            if min_value is not None and max_value is not None and min_value > max_value:
                raise ValueError('Minimum value cannot be greater than maximum value.')
            if min_value is None and max_value is None:
                raise ValueError('Enter at least one numeric bound.')
            return {
                'field': field,
                'type': 'numeric',
                'min': min_value,
                'max': max_value,
            }

        def accept(self):
            try:
                self.get_filter_spec()
            except Exception as exc:
                QtWidgets.QMessageBox.warning(self, 'Invalid Filter', str(exc))
                return
            super().accept()

    class SampleFiltersDialog(QtWidgets.QDialog):
        def __init__(self, filter_source, filters=None, parent=None, subject_label='Sample'):
            super().__init__(parent)
            self.filter_source = filter_source if isinstance(filter_source, FilterDataSource) else FilterDataSource(filter_source)
            self.filters = [dict(entry) for entry in (filters or [])]
            self.subject_label = str(subject_label or 'Sample').strip() or 'Sample'
            subject_label_lower = self.subject_label.lower()
            self.setWindowTitle(f'{self.subject_label} Filters')
            self.resize(760, 420)

            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)

            info = QtWidgets.QLabel(
                f'Add one or more {subject_label_lower} filters. Categorical filters keep selected values. '
                'Numeric filters keep values inside the requested range.'
            )
            info.setWordWrap(True)
            layout.addWidget(info)

            self.table = QtWidgets.QTableWidget()
            self.table.setColumnCount(3)
            self.table.setHorizontalHeaderLabels(['Field', 'Type', 'Criteria'])
            self.table.horizontalHeader().setStretchLastSection(True)
            self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
            self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
            layout.addWidget(self.table)

            button_row = QtWidgets.QHBoxLayout()
            self.add_btn = QtWidgets.QPushButton('Add Filter')
            self.edit_btn = QtWidgets.QPushButton('Edit Filter')
            self.remove_btn = QtWidgets.QPushButton('Remove Filter')
            self.add_btn.clicked.connect(self.add_filter)
            self.edit_btn.clicked.connect(self.edit_selected_filter)
            self.remove_btn.clicked.connect(self.remove_selected_filter)
            button_row.addWidget(self.add_btn)
            button_row.addWidget(self.edit_btn)
            button_row.addWidget(self.remove_btn)
            button_row.addStretch(1)
            layout.addLayout(button_row)

            dialog_buttons = QtWidgets.QHBoxLayout()
            dialog_buttons.addStretch(1)
            ok_btn = QtWidgets.QPushButton('OK')
            cancel_btn = QtWidgets.QPushButton('Cancel')
            ok_btn.clicked.connect(self.accept)
            cancel_btn.clicked.connect(self.reject)
            dialog_buttons.addWidget(ok_btn)
            dialog_buttons.addWidget(cancel_btn)
            layout.addLayout(dialog_buttons)

            self._refresh_table()

        def _refresh_table(self):
            self.table.setRowCount(len(self.filters))
            for row, filter_spec in enumerate(self.filters):
                self.table.setItem(row, 0, QtWidgets.QTableWidgetItem(str(filter_spec.get('field', ''))))
                self.table.setItem(row, 1, QtWidgets.QTableWidgetItem(str(filter_spec.get('type', ''))))
                self.table.setItem(row, 2, QtWidgets.QTableWidgetItem(summarize_sample_filter_spec(filter_spec)))
            self.table.resizeRowsToContents()

        def add_filter(self):
            dialog = SampleFilterEditDialog(self.filter_source, parent=self, subject_label=self.subject_label)
            if dialog.exec_() == QtWidgets.QDialog.Accepted:
                self.filters.append(dialog.get_filter_spec())
                self._refresh_table()

        def edit_selected_filter(self):
            row = self.table.currentRow()
            if row < 0 or row >= len(self.filters):
                QtWidgets.QMessageBox.information(self, 'Edit Filter', 'Select a filter to edit.')
                return
            dialog = SampleFilterEditDialog(
                self.filter_source,
                filter_spec=self.filters[row],
                parent=self,
                subject_label=self.subject_label,
            )
            if dialog.exec_() == QtWidgets.QDialog.Accepted:
                self.filters[row] = dialog.get_filter_spec()
                self._refresh_table()

        def remove_selected_filter(self):
            row = self.table.currentRow()
            if row < 0 or row >= len(self.filters):
                QtWidgets.QMessageBox.information(self, 'Remove Filter', 'Select a filter to remove.')
                return
            del self.filters[row]
            self._refresh_table()

        def get_filters(self):
            return [dict(entry) for entry in self.filters]
        
    class ConfigDialog(QtWidgets.QDialog):
        def __init__(self):
            super().__init__()
            self.should_visualize = True
            self.viewer_backend = 'taichi'
            self._prefer_interpolation_file_for_viewer = False
            self._domain_catalog_cache = None
            self._viewer_process = None
            self._viewer_config_path = None
            self._viewer_render_mode = None
            self._active_operation_thread = None
            self._active_operation_worker = None
            self._active_operation_progress = None
            self._suspend_auto_viewer_refresh = False
            self._viewer_process_timer = QtCore.QTimer(self)
            self._viewer_process_timer.setInterval(2000)
            self._viewer_process_timer.timeout.connect(self._cleanup_finished_viewer)
            self._viewer_process_timer.start()
            self.taichi_block_render_mode_default = 'mesh'
            self.taichi_transparent_blocks_default = False
            self._legacy_fill_unvisited_domainwise = False
            self.setWindowTitle("Anterpolator 3D Viewer Configuration")
            self.setWindowFlags(
                self.windowFlags() | QtCore.Qt.WindowMinimizeButtonHint | QtCore.Qt.WindowMaximizeButtonHint
            )
            self.resize(1120, 680)
            
            # Main layout with tabs
            main_layout = QtWidgets.QVBoxLayout()
            self.setLayout(main_layout)
            
            # Create tab widget
            tabs = QtWidgets.QTabWidget()
            main_layout.addWidget(tabs)
            
            # Tab 1: Files & Data
            files_tab = QtWidgets.QWidget()
            files_form = QtWidgets.QFormLayout()
            files_tab.setLayout(files_form)
            tabs.addTab(files_tab, "Files & Data")

            operations_tab = QtWidgets.QWidget()
            operations_tab_layout = QtWidgets.QVBoxLayout()
            operations_tab_layout.setContentsMargins(0, 0, 0, 0)
            operations_tab.setLayout(operations_tab_layout)
            operations_scroll = QtWidgets.QScrollArea()
            operations_scroll.setWidgetResizable(True)
            operations_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
            operations_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
            operations_tab_layout.addWidget(operations_scroll)
            operations_content = QtWidgets.QWidget()
            operations_form = QtWidgets.QFormLayout()
            operations_content.setLayout(operations_form)
            operations_scroll.setWidget(operations_content)
            
            # Tab 2: Ant Colony Parameters
            ant_tab = QtWidgets.QWidget()
            ant_form = QtWidgets.QFormLayout()
            ant_tab.setLayout(ant_form)
            
            # Tab 3: Molecular Clock Parameters
            mc_tab = QtWidgets.QWidget()
            mc_form = QtWidgets.QFormLayout()
            mc_tab.setLayout(mc_form)

            gk_tab = QtWidgets.QWidget()
            gk_form = QtWidgets.QFormLayout()
            gk_tab.setLayout(gk_form)

            octree_tab = QtWidgets.QWidget()
            octree_form = QtWidgets.QFormLayout()
            octree_tab.setLayout(octree_form)

            st_tab = QtWidgets.QWidget()
            st_form = QtWidgets.QFormLayout()
            st_tab.setLayout(st_form)

            display_tab = QtWidgets.QWidget()
            display_form = QtWidgets.QFormLayout()
            display_tab.setLayout(display_form)

            # Tab 6: Advanced Options
            advanced_tab = QtWidgets.QWidget()
            advanced_form = QtWidgets.QFormLayout()
            advanced_tab.setLayout(advanced_form)

            tabs.addTab(operations_tab, "Operations")
            tabs.addTab(st_tab, "String Theory")
            tabs.addTab(ant_tab, "Ant Colony")
            tabs.addTab(mc_tab, "Molecular Clock")
            tabs.addTab(gk_tab, "Gaussian Kernel")
            tabs.addTab(octree_tab, "Adaptive Octree")
            tabs.addTab(display_tab, "Display")
            tabs.addTab(advanced_tab, "Advanced")

            # === FILES & DATA TAB ===
            self.samples_edit = QtWidgets.QLineEdit('Data/ANT-Samples.csv')
            self.blocks_edit = QtWidgets.QLineEdit('Data/ANT-Domains.csv')
            self.color_edit = QtWidgets.QLineEdit('Data/Value.lfc')
            self.interp_edit = QtWidgets.QLineEdit('')
            self.configure_block_domain_metrics_filters_btn = QtWidgets.QPushButton('Configure Sample Filters...')
            self.configure_block_volume_weighted_filters_btn = QtWidgets.QPushButton('Configure Block Filters...')
            self.block_domain_sample_filters = []
            self.block_volume_weighted_filters = []
            self.block_domain_metrics_filters_summary = QtWidgets.QLabel('No sample filters configured. These filters apply app-wide.')
            self.block_domain_metrics_filters_summary.setWordWrap(True)
            self.block_volume_weighted_filters_summary = QtWidgets.QLabel('No block filters configured. These filters apply app-wide.')
            self.block_volume_weighted_filters_summary.setWordWrap(True)

            def add_file_row(label, line_edit, filter_str, form_layout, extra_button=None):
                h = QtWidgets.QHBoxLayout()
                h.addWidget(line_edit)
                if extra_button is not None:
                    h.addWidget(extra_button)
                btn = QtWidgets.QPushButton('Browse')
                def pick():
                    if 'Interpolation' in label:
                        base = os.path.splitext(os.path.basename(self.samples_edit.text()))[0]
                        suggested = f"{base}_anterpolation.csv" if base else "interpolation.csv"
                        start_dir = os.path.dirname(self.samples_edit.text()) if os.path.isfile(self.samples_edit.text()) else '.'
                        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, label, os.path.join(start_dir, suggested), filter_str)
                    else:
                        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, label, '.', filter_str)
                    if path:
                        line_edit.setText(path)
                        if line_edit is self.samples_edit and not self.interp_edit.text():
                            base = os.path.splitext(os.path.basename(path))[0]
                            self.interp_edit.setText(os.path.join(os.path.dirname(path), f"{base}_anterpolation.csv"))
                btn.clicked.connect(pick)
                h.addWidget(btn)
                form_layout.addRow(label, h)
            
            add_file_row('Samples File', self.samples_edit, 'CSV Files (*.csv)', files_form, extra_button=self.configure_block_domain_metrics_filters_btn)
            files_form.addRow('', self.block_domain_metrics_filters_summary)
            add_file_row('Blocks File', self.blocks_edit, 'CSV Files (*.csv);;BMF Files (*.bmf);;All Files (*.*)', files_form, extra_button=self.configure_block_volume_weighted_filters_btn)
            files_form.addRow('', self.block_volume_weighted_filters_summary)
            add_file_row('Color File', self.color_edit, 'LFC Files (*.lfc);;All Files (*.*)', files_form)
            add_file_row('Interpolation File', self.interp_edit, 'CSV Files (*.csv);;All Files (*.*)', files_form)

            self.block_evaluated_samples_enabled = QtWidgets.QCheckBox('Block Evaluated Samples File')
            self.block_evaluated_samples_edit = QtWidgets.QLineEdit('')
            self.block_evaluated_samples_browse = QtWidgets.QPushButton('Browse')
            self.block_evaluated_samples_edit.setEnabled(False)
            self.block_evaluated_samples_browse.setEnabled(False)

            def suggested_block_evaluated_samples_path():
                sample_path = self.samples_edit.text().strip()
                base = os.path.splitext(os.path.basename(sample_path))[0] if sample_path else ''
                start_dir = os.path.dirname(sample_path) if sample_path and os.path.isfile(sample_path) else '.'
                filename = f"{base}_block_evaluated_samples.csv" if base else 'block_evaluated_samples.csv'
                return os.path.join(start_dir, filename)

            def update_block_evaluated_samples_controls(checked):
                self.block_evaluated_samples_edit.setEnabled(checked)
                self.block_evaluated_samples_browse.setEnabled(checked)
                if checked and not self.block_evaluated_samples_edit.text().strip():
                    self.block_evaluated_samples_edit.setText(suggested_block_evaluated_samples_path())

            def browse_block_evaluated_samples_file():
                initial_path = self.block_evaluated_samples_edit.text().strip() or suggested_block_evaluated_samples_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Block Evaluated Samples File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.block_evaluated_samples_edit.setText(path)

            self.block_evaluated_samples_enabled.toggled.connect(update_block_evaluated_samples_controls)
            self.block_evaluated_samples_browse.clicked.connect(browse_block_evaluated_samples_file)

            block_eval_layout = QtWidgets.QHBoxLayout()
            block_eval_layout.addWidget(self.block_evaluated_samples_edit)
            block_eval_layout.addWidget(self.block_evaluated_samples_browse)
            files_form.addRow(self.block_evaluated_samples_enabled, block_eval_layout)

            # Delimiter selectors
            delim_opts = [',',';','\t','|']
            self.samples_delim = QtWidgets.QComboBox(); self.samples_delim.addItems(delim_opts)
            self.blocks_delim = QtWidgets.QComboBox(); self.blocks_delim.addItems(delim_opts)
            # Detect initial if file exists
            if os.path.isfile(self.samples_edit.text()):
                det = detect_csv_delimiter(self.samples_edit.text())
                if det in delim_opts:
                    self.samples_delim.setCurrentText(det)
            if os.path.isfile(self.blocks_edit.text()):
                detb = detect_csv_delimiter(self.blocks_edit.text())
                if detb in delim_opts:
                    self.blocks_delim.setCurrentText(detb)
            files_form.addRow('Samples Delimiter', self.samples_delim)
            files_form.addRow('Blocks Delimiter', self.blocks_delim)

            # Header line selectors (1-based)
            self.samples_header_line = QtWidgets.QSpinBox(); self.samples_header_line.setRange(1, 1_000_000); self.samples_header_line.setValue(1)
            self.blocks_header_line = QtWidgets.QSpinBox(); self.blocks_header_line.setRange(1, 1_000_000); self.blocks_header_line.setValue(1)
            files_form.addRow('Samples Header Line', self.samples_header_line)
            files_form.addRow('Blocks Header Line', self.blocks_header_line)

            def configure_column_combo(combo_box, minimum_width=170, popup_width=320):
                combo_box.setEditable(False)
                combo_box.setMinimumWidth(minimum_width)
                combo_box.setMinimumContentsLength(16)
                combo_box.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
                combo_box.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
                combo_box.view().setMinimumWidth(popup_width)

            # Column mapping combo boxes for samples
            self.sample_x_col = QtWidgets.QComboBox(); self.sample_y_col = QtWidgets.QComboBox(); self.sample_z_col = QtWidgets.QComboBox(); self.sample_value_col = QtWidgets.QComboBox(); self.sample_domain_col = QtWidgets.QComboBox(); self.sample_weight_col = QtWidgets.QComboBox()
            self.sample_blocks_include_ids = QtWidgets.QCheckBox('Include concatenated sample IDs')
            self.sample_blocks_include_ids.setToolTip(
                'When enabled, export a text metadata column listing the contributing sample IDs for each sample block.'
            )
            self.sample_blocks_id_cols = [QtWidgets.QComboBox(), QtWidgets.QComboBox(), QtWidgets.QComboBox()]
            self.block_domain_metrics_id_cols = [QtWidgets.QComboBox(), QtWidgets.QComboBox(), QtWidgets.QComboBox()]
            self.block_domain_metrics_prefix_with_block_value = QtWidgets.QCheckBox('Prefix closest-sample ID and nearest-sample value metrics with block value column')
            self.block_domain_metrics_prefix_with_block_value.setChecked(True)
            self.block_domain_metrics_prefix_with_block_value.setToolTip(
                'When enabled, closest-sample ID and nearest-sample value/residual metrics are exported as <Block Value Column>_<Metric Name>.'
            )
            for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col, self.sample_domain_col, self.sample_weight_col]:
                configure_column_combo(cb)
            for cb in self.sample_blocks_id_cols:
                configure_column_combo(cb)
                cb.addItem('(None)')
            for cb in self.block_domain_metrics_id_cols:
                configure_column_combo(cb)
                cb.addItem('(None)')
            self.sample_domain_col.addItem('(None)')
            self.sample_weight_col.addItem('(None)')
            sample_map_layout = QtWidgets.QHBoxLayout(); sample_map_layout.addWidget(self.sample_x_col); sample_map_layout.addWidget(self.sample_y_col); sample_map_layout.addWidget(self.sample_z_col); sample_map_layout.addWidget(self.sample_value_col); sample_map_layout.addWidget(self.sample_domain_col)
            for index in range(5):
                sample_map_layout.setStretch(index, 1)
            files_form.addRow('Samples Columns (X Y Z Value Domain)', sample_map_layout)
            self.sample_weight_col.setToolTip('Optional weight column for block-level sample averaging. When selected, Sample Blocks export and interpolation preprocessing use weighted averages instead of plain means.')
            files_form.addRow('Sample Weight Column', self.sample_weight_col)

            # Column mapping combo boxes for blocks
            self.block_x_col = QtWidgets.QComboBox(); self.block_y_col = QtWidgets.QComboBox(); self.block_z_col = QtWidgets.QComboBox(); self.block_domain_col = QtWidgets.QComboBox(); self.block_domain_metrics_value_col = QtWidgets.QComboBox(); self.block_volume_weighted_value_col = QtWidgets.QComboBox(); self.block_weight_metric_col = QtWidgets.QComboBox()
            for cb in [self.block_x_col, self.block_y_col, self.block_z_col, self.block_domain_col, self.block_domain_metrics_value_col, self.block_volume_weighted_value_col, self.block_weight_metric_col]:
                configure_column_combo(cb)
            self.block_domain_col.addItem('(None)')
            self.block_domain_metrics_value_col.addItem('(None)')
            self.block_volume_weighted_value_col.addItem('(None)')
            self.block_weight_metric_col.addItem('(Volume)')
            block_map_layout = QtWidgets.QHBoxLayout(); block_map_layout.addWidget(self.block_x_col); block_map_layout.addWidget(self.block_y_col); block_map_layout.addWidget(self.block_z_col); block_map_layout.addWidget(self.block_domain_col)
            for index in range(4):
                block_map_layout.setStretch(index, 1)
            files_form.addRow('Blocks Columns (X Y Z Domain)', block_map_layout)

            # Block Size
            self.block_x = QtWidgets.QSpinBox(); self.block_x.setRange(1, 10000); self.block_x.setValue(10)
            self.block_y = QtWidgets.QSpinBox(); self.block_y.setRange(1, 10000); self.block_y.setValue(10)
            self.block_z = QtWidgets.QSpinBox(); self.block_z.setRange(1, 10000); self.block_z.setValue(10)
            block_size_tooltip = (
                'Set base block dimensions (x, y, z).\n'
                'You can change these values to override the implicit\n'
                'resolution in the source blocks file, or to\n'
                'explicitly confirm the base-block resolution\n'
                'defined by that file.'
            )
            self.block_x.setToolTip(block_size_tooltip)
            self.block_y.setToolTip(block_size_tooltip)
            self.block_z.setToolTip(block_size_tooltip)
            bx_layout = QtWidgets.QHBoxLayout(); bx_layout.addWidget(self.block_x); bx_layout.addWidget(self.block_y); bx_layout.addWidget(self.block_z)
            block_size_label = QtWidgets.QLabel('Block Size (x,y,z)')
            block_size_label.setToolTip(block_size_tooltip)
            files_form.addRow(block_size_label, bx_layout)

            def refresh_sample_columns():
                path = self.samples_edit.text().strip()
                delim = self.samples_delim.currentText()
                header_line = self.samples_header_line.value()
                current_domain = self.sample_domain_col.currentText()
                current_weight = self.sample_weight_col.currentText()
                current_sample_blocks_id_cols = [cb.currentText() for cb in self.sample_blocks_id_cols]
                current_metrics_id_cols = [cb.currentText() for cb in self.block_domain_metrics_id_cols]
                for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col]:
                    cb.clear()
                self.sample_domain_col.clear(); self.sample_domain_col.addItem('(None)')
                self.sample_weight_col.clear(); self.sample_weight_col.addItem('(None)')
                for cb in self.sample_blocks_id_cols:
                    cb.clear(); cb.addItem('(None)')
                for cb in self.block_domain_metrics_id_cols:
                    cb.clear(); cb.addItem('(None)')
                if not os.path.isfile(path):
                    return
                try:
                    cols = parse_effective_header_line(path, delim, header_line)
                    for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col, self.sample_domain_col, self.sample_weight_col]:
                        cb.addItems(cols)
                    for cb in self.sample_blocks_id_cols:
                        cb.addItems(cols)
                    for cb in self.block_domain_metrics_id_cols:
                        cb.addItems(cols)
                    # Attempt auto-suggest
                    def suggest(cb, keywords):
                        for k in keywords:
                            for i in range(cb.count()):
                                if cb.itemText(i).lower() == k:
                                    cb.setCurrentIndex(i); return
                    suggest(self.sample_x_col, ['x','easting'])
                    suggest(self.sample_y_col, ['y','northing'])
                    suggest(self.sample_z_col, ['z','elevation','rl'])
                    suggest(self.sample_value_col, ['value','grade'])
                    suggest(self.sample_domain_col, ['domain','dom','dg'])
                    suggest(self.sample_weight_col, ['length','interval_length','sample_length','weight'])
                    if current_domain and current_domain != '(None)':
                        idx = self.sample_domain_col.findText(current_domain)
                        if idx >= 0:
                            self.sample_domain_col.setCurrentIndex(idx)
                    if current_weight and current_weight != '(None)':
                        idx = self.sample_weight_col.findText(current_weight)
                        if idx >= 0:
                            self.sample_weight_col.setCurrentIndex(idx)
                    for cb, current_value in zip(self.sample_blocks_id_cols, current_sample_blocks_id_cols):
                        if current_value and current_value != '(None)':
                            idx = cb.findText(current_value)
                            if idx >= 0:
                                cb.setCurrentIndex(idx)
                    for cb, current_value in zip(self.block_domain_metrics_id_cols, current_metrics_id_cols):
                        if current_value and current_value != '(None)':
                            idx = cb.findText(current_value)
                            if idx >= 0:
                                cb.setCurrentIndex(idx)
                except Exception:
                    pass

            def refresh_block_columns():
                path = self.blocks_edit.text().strip()
                delim = self.blocks_delim.currentText()
                header_line = self.blocks_header_line.value()
                for cb in [self.block_x_col, self.block_y_col, self.block_z_col]:
                    cb.clear()
                # Domain retains (None) entry
                current_domain = self.block_domain_col.currentText()
                current_domain_metrics_value_col = self.block_domain_metrics_value_col.currentText()
                current_volume_weighted_value_col = self.block_volume_weighted_value_col.currentText()
                current_weight_col = self.block_weight_metric_col.currentText()
                self.block_domain_col.clear(); self.block_domain_col.addItem('(None)')
                self.block_domain_metrics_value_col.clear(); self.block_domain_metrics_value_col.addItem('(None)')
                self.block_volume_weighted_value_col.clear(); self.block_volume_weighted_value_col.addItem('(None)')
                self.block_weight_metric_col.clear(); self.block_weight_metric_col.addItem('(Volume)')
                if not os.path.isfile(path):
                    return
                try:
                    header_line = sync_csv_header_line_widget(self.blocks_header_line, path, header_line)
                    cols = parse_effective_header_line(path, delim, header_line)
                    for cb in [self.block_x_col, self.block_y_col, self.block_z_col, self.block_domain_col, self.block_domain_metrics_value_col, self.block_volume_weighted_value_col, self.block_weight_metric_col]:
                        for c in cols:
                            cb.addItem(c)
                    # Auto-suggest
                    def suggest(cb, keywords):
                        for k in keywords:
                            for i in range(cb.count()):
                                if cb.itemText(i).lower() == k:
                                    cb.setCurrentIndex(i); return
                    suggest(self.block_x_col, ['x','easting'])
                    suggest(self.block_y_col, ['y','northing'])
                    suggest(self.block_z_col, ['z','elevation','rl'])
                    suggest(self.block_domain_col, ['domain','dom'])
                    suggest(self.block_domain_metrics_value_col, ['value','grade'])
                    suggest(self.block_volume_weighted_value_col, ['value','grade'])
                    # Restore domain selection if possible
                    if current_domain and current_domain != '(None)':
                        idx = self.block_domain_col.findText(current_domain)
                        if idx >= 0: self.block_domain_col.setCurrentIndex(idx)
                    if current_domain_metrics_value_col and current_domain_metrics_value_col != '(None)':
                        idx = self.block_domain_metrics_value_col.findText(current_domain_metrics_value_col)
                        if idx >= 0:
                            self.block_domain_metrics_value_col.setCurrentIndex(idx)
                    if current_volume_weighted_value_col and current_volume_weighted_value_col != '(None)':
                        idx = self.block_volume_weighted_value_col.findText(current_volume_weighted_value_col)
                        if idx >= 0:
                            self.block_volume_weighted_value_col.setCurrentIndex(idx)
                    if current_weight_col and current_weight_col != '(Volume)':
                        idx = self.block_weight_metric_col.findText(current_weight_col)
                        if idx >= 0:
                            self.block_weight_metric_col.setCurrentIndex(idx)
                except Exception:
                    pass

            # Connect signals to refresh
            self.samples_delim.currentIndexChanged.connect(refresh_sample_columns)
            self.samples_header_line.valueChanged.connect(lambda _: refresh_sample_columns())
            self.samples_edit.textChanged.connect(lambda _: refresh_sample_columns())
            self.blocks_delim.currentIndexChanged.connect(refresh_block_columns)
            self.blocks_header_line.valueChanged.connect(lambda _: refresh_block_columns())
            self.blocks_edit.textChanged.connect(lambda _: refresh_block_columns())
            self.blocks_delim.currentIndexChanged.connect(lambda _: self._invalidate_domain_catalog_cache())
            self.blocks_header_line.valueChanged.connect(lambda _: self._invalidate_domain_catalog_cache())
            self.blocks_edit.textChanged.connect(lambda _: self._invalidate_domain_catalog_cache())
            self.block_domain_col.currentTextChanged.connect(lambda _: self._invalidate_domain_catalog_cache())

            # Initial refresh (silent if files missing)
            refresh_sample_columns()
            refresh_block_columns()

            def suggested_domain_samples_path():
                sample_path = self.samples_edit.text().strip()
                base_name = os.path.splitext(os.path.basename(sample_path))[0] if sample_path else ''
                output_dir = os.path.dirname(sample_path) if sample_path else '.'
                domain_name = self.block_domain_col.currentText() if hasattr(self, 'block_domain_col') else 'Domain'
                domain_suffix = _sanitize_filename_fragment(domain_name, fallback='Domain')
                filename = f"{base_name}+{domain_suffix}.csv" if base_name else f"samples+{domain_suffix}.csv"
                return os.path.join(output_dir or '.', filename)

            def suggested_sample_blocks_path():
                return resolve_sample_blocks_export_path(
                    None,
                    samples_file=self.samples_edit.text().strip(),
                )

            def suggested_block_domain_metrics_path():
                return resolve_block_domain_metrics_export_path(
                    None,
                    blocks_file=self.blocks_edit.text().strip(),
                    domain_col=self.block_domain_col.currentText(),
                )

            def suggested_block_value_transfer_path():
                return resolve_block_value_transfer_export_path(
                    None,
                    samples_file=self.samples_edit.text().strip(),
                    block_value_cols=self._get_selected_block_value_transfer_columns(),
                )

            def suggested_block_model_transfer_path():
                return resolve_block_model_transfer_export_path(
                    None,
                    target_blocks_file=self.block_model_target_edit.text().strip(),
                )

            def suggested_table_attribute_output_path():
                return resolve_block_model_table_attribute_export_path(
                    None,
                    block_model_file=self.table_attribute_block_model_edit.text().strip(),
                    table_file=self.table_attribute_table_edit.text().strip(),
                )

            def suggested_domain_interpolation_confidence_path():
                return resolve_domain_interpolation_confidence_export_path(
                    None,
                    blocks_file=self.blocks_edit.text().strip(),
                    domain_col=self.block_domain_col.currentText(),
                )

            def suggested_block_volume_weighted_path():
                return resolve_block_volume_weighted_average_export_path(
                    None,
                    blocks_file=self.blocks_edit.text().strip(),
                    value_col=self.block_volume_weighted_value_col.currentText(),
                )

            def suggested_equation_finder_path():
                return resolve_equation_finder_export_path(
                    None,
                    samples_file=self.samples_edit.text().strip(),
                    value_col=self.sample_value_col.currentText(),
                    domain_col=self.sample_domain_col.currentText(),
                )

            def suggested_bmf_export_input_path():
                interpolation_path = self.interp_edit.text().strip()
                blocks_path = self.blocks_edit.text().strip()
                candidates = [interpolation_path, blocks_path]
                for candidate in candidates:
                    if candidate and not is_bmf_file(candidate):
                        return candidate
                return interpolation_path or blocks_path

            def suggested_bmf_export_output_path():
                input_path = self.bmf_export_input_edit.text().strip() if hasattr(self, 'bmf_export_input_edit') else ''
                if not input_path:
                    input_path = suggested_bmf_export_input_path()
                base_name = os.path.splitext(os.path.basename(input_path))[0] if input_path else 'grid'
                output_dir = os.path.dirname(input_path) if input_path else '.'
                return os.path.join(output_dir or '.', f"{base_name}.bmf")

            self.sample_blocks_output_edit = QtWidgets.QLineEdit('')
            self.sample_blocks_browse = QtWidgets.QPushButton('Browse')
            self.start_sample_blocks_btn = QtWidgets.QPushButton('Export Sample Blocks')
            self.domain_samples_output_edit = QtWidgets.QLineEdit('')
            self.domain_samples_browse = QtWidgets.QPushButton('Browse')
            self.start_domaining_btn = QtWidgets.QPushButton('Start Domaining')
            self.block_value_transfer_output_edit = QtWidgets.QLineEdit('')
            self.block_value_transfer_browse = QtWidgets.QPushButton('Browse')
            self.start_block_value_transfer_btn = QtWidgets.QPushButton('Transfer Columns')
            self.block_domain_metrics_output_edit = QtWidgets.QLineEdit('')
            self.block_domain_metrics_browse = QtWidgets.QPushButton('Browse')
            self.start_block_domain_metrics_btn = QtWidgets.QPushButton('Export Metrics')
            self.block_value_transfer_cols = QtWidgets.QListWidget()
            self.block_value_transfer_cols.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.block_value_transfer_cols.setMinimumHeight(120)
            self.block_value_transfer_select_all_btn = QtWidgets.QPushButton('Select All')
            self.block_value_transfer_clear_btn = QtWidgets.QPushButton('Clear')
            self.block_value_transfer_summary = QtWidgets.QLabel('No block columns selected for transfer.')
            self.block_model_target_edit = QtWidgets.QLineEdit('')
            self.block_model_target_browse = QtWidgets.QPushButton('Browse')
            self.block_model_target_delim = QtWidgets.QComboBox()
            self.block_model_target_delim.addItems(delim_opts)
            self.block_model_target_header_line = QtWidgets.QSpinBox()
            self.block_model_target_header_line.setRange(1, 1_000_000)
            self.block_model_target_header_line.setValue(1)
            self.block_model_target_x_col = QtWidgets.QComboBox()
            self.block_model_target_y_col = QtWidgets.QComboBox()
            self.block_model_target_z_col = QtWidgets.QComboBox()
            self.block_model_source_size_cols = [QtWidgets.QComboBox() for _ in range(3)]
            self.block_model_target_size_cols = [QtWidgets.QComboBox() for _ in range(3)]
            self.block_model_target_size_spins = [TrimmedDisplayDoubleSpinBox() for _ in range(3)]
            for spin, source_spin in zip(self.block_model_target_size_spins, [self.block_x, self.block_y, self.block_z]):
                spin.setRange(0.001, 1_000_000_000.0)
                spin.setDecimals(6)
                spin.setValue(source_spin.value())
            self.block_model_transfer_cols = QtWidgets.QListWidget()
            self.block_model_transfer_cols.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.block_model_transfer_cols.setMinimumHeight(120)
            self.block_model_transfer_select_all_btn = QtWidgets.QPushButton('Select All')
            self.block_model_transfer_clear_btn = QtWidgets.QPushButton('Clear')
            self.block_model_transfer_summary = QtWidgets.QLabel('No source columns selected.')
            self.block_model_transfer_summary.setWordWrap(True)
            self.block_model_transfer_output_edit = QtWidgets.QLineEdit('')
            self.block_model_transfer_output_browse = QtWidgets.QPushButton('Browse')
            self.block_model_nearest_fallback = QtWidgets.QCheckBox('Use nearest source block when there is no overlap')
            self.block_model_nearest_fallback.setChecked(True)
            self.block_model_nearest_max_distance = TrimmedDisplayDoubleSpinBox()
            self.block_model_nearest_max_distance.setRange(0.0, 1_000_000_000.0)
            self.block_model_nearest_max_distance.setDecimals(3)
            self.block_model_nearest_max_distance.setValue(0.0)
            self.block_model_nearest_max_distance.setSpecialValueText('Unlimited')
            self.block_model_nearest_max_distance.setSuffix(' m')
            self.start_block_model_transfer_btn = QtWidgets.QPushButton('Transfer To Target Blocks')
            self.table_attribute_block_model_edit = QtWidgets.QLineEdit('')
            self.table_attribute_block_model_browse = QtWidgets.QPushButton('Browse')
            self.table_attribute_block_model_delim = QtWidgets.QComboBox()
            self.table_attribute_block_model_delim.addItems(delim_opts)
            self.table_attribute_block_model_header_line = QtWidgets.QSpinBox()
            self.table_attribute_block_model_header_line.setRange(1, 1_000_000)
            self.table_attribute_block_model_header_line.setValue(1)
            self.table_attribute_table_edit = QtWidgets.QLineEdit('')
            self.table_attribute_table_browse = QtWidgets.QPushButton('Browse')
            self.table_attribute_table_delim = QtWidgets.QComboBox()
            self.table_attribute_table_delim.addItems(delim_opts)
            self.table_attribute_table_header_line = QtWidgets.QSpinBox()
            self.table_attribute_table_header_line.setRange(1, 1_000_000)
            self.table_attribute_table_header_line.setValue(1)
            self.table_attribute_key_cols = QtWidgets.QListWidget()
            self.table_attribute_key_cols.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.table_attribute_key_cols.setMinimumHeight(100)
            self.table_attribute_key_select_all_btn = QtWidgets.QPushButton('Select All')
            self.table_attribute_key_clear_btn = QtWidgets.QPushButton('Clear')
            self.table_attribute_value_cols = QtWidgets.QListWidget()
            self.table_attribute_value_cols.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.table_attribute_value_cols.setMinimumHeight(120)
            self.table_attribute_value_select_all_btn = QtWidgets.QPushButton('Select All')
            self.table_attribute_value_clear_btn = QtWidgets.QPushButton('Clear')
            self.table_attribute_summary = QtWidgets.QLabel('No shared key columns are available yet.')
            self.table_attribute_summary.setWordWrap(True)
            self.table_attribute_output_edit = QtWidgets.QLineEdit('')
            self.table_attribute_output_browse = QtWidgets.QPushButton('Browse')
            self.start_table_attribute_assign_btn = QtWidgets.QPushButton('Assign Attributes')
            self.domain_interpolation_confidence_output_edit = QtWidgets.QLineEdit('')
            self.domain_interpolation_confidence_browse = QtWidgets.QPushButton('Browse')
            self.start_domain_interpolation_confidence_btn = QtWidgets.QPushButton('Export Confidence Metrics')
            self.block_volume_weighted_output_edit = QtWidgets.QLineEdit('')
            self.block_volume_weighted_browse = QtWidgets.QPushButton('Browse')
            self.start_block_volume_weighted_btn = QtWidgets.QPushButton('Export Volume Weighted Average')
            self.bmf_export_input_edit = QtWidgets.QLineEdit('')
            self.bmf_export_input_browse = QtWidgets.QPushButton('Browse')
            self.bmf_export_output_edit = QtWidgets.QLineEdit('')
            self.bmf_export_output_browse = QtWidgets.QPushButton('Browse')
            self.bmf_export_delim = QtWidgets.QComboBox()
            self.bmf_export_delim.addItems(delim_opts)
            self.bmf_export_header_line = QtWidgets.QSpinBox()
            self.bmf_export_header_line.setRange(1, 1_000_000)
            self.bmf_export_header_line.setValue(1)
            self.bmf_export_backend_combo = QtWidgets.QComboBox()
            self.bmf_export_backend_combo.addItems(['tbms-config-text', 'tbms-experimental', 'vulcan'])
            self.bmf_export_cell_x = QtWidgets.QDoubleSpinBox()
            self.bmf_export_cell_y = QtWidgets.QDoubleSpinBox()
            self.bmf_export_cell_z = QtWidgets.QDoubleSpinBox()
            for spin in [self.bmf_export_cell_x, self.bmf_export_cell_y, self.bmf_export_cell_z]:
                spin.setRange(0.0, 1_000_000_000.0)
                spin.setDecimals(6)
                spin.setSingleStep(1.0)
                spin.setValue(0.0)
                spin.setSpecialValueText('Auto')
            self.bmf_export_regularize_base_blocks = QtWidgets.QCheckBox('Regularize to base block size')
            self.bmf_export_regularize_base_blocks.setChecked(False)
            self.bmf_export_regularize_warning = QtWidgets.QLabel(
                'Info: unchecked preserves the source CSV block rows; checked aggregates rows to base blocks and may create a dense regular grid.'
            )
            self.bmf_export_regularize_warning.setWordWrap(True)
            self.bmf_export_regularize_warning.setStyleSheet('color: #b36b00;')
            self.bmf_export_x_col = QtWidgets.QComboBox()
            self.bmf_export_y_col = QtWidgets.QComboBox()
            self.bmf_export_z_col = QtWidgets.QComboBox()
            for cb in [self.bmf_export_x_col, self.bmf_export_y_col, self.bmf_export_z_col]:
                configure_column_combo(cb)
            self.bmf_export_value_cols = QtWidgets.QTableWidget(0, 2)
            self.bmf_export_value_cols.setHorizontalHeaderLabels(['Field', 'Type'])
            self.bmf_export_value_cols.verticalHeader().setVisible(False)
            self.bmf_export_value_cols.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
            self.bmf_export_value_cols.setMinimumHeight(160)
            self.bmf_export_value_cols.horizontalHeader().setStretchLastSection(False)
            self.bmf_export_value_cols.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
            self.bmf_export_value_cols.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
            self.bmf_export_select_all_btn = QtWidgets.QPushButton('Select All')
            self.bmf_export_clear_btn = QtWidgets.QPushButton('Clear')
            self.bmf_export_value_summary = QtWidgets.QLabel('No export value columns loaded.')
            self.bmf_export_value_summary.setWordWrap(True)
            self.bmf_export_value_exceptions = {}
            self.bmf_export_exception_table = QtWidgets.QTableWidget(0, 4)
            self.bmf_export_exception_table.setHorizontalHeaderLabels(['Field', 'CSV Value', 'Replacement', 'Use In Regularization'])
            self.bmf_export_exception_table.verticalHeader().setVisible(False)
            self.bmf_export_exception_table.setMinimumHeight(120)
            self.bmf_export_exception_table.horizontalHeader().setStretchLastSection(False)
            self.bmf_export_exception_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
            self.bmf_export_exception_table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
            self.bmf_export_exception_table.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
            self.bmf_export_exception_table.horizontalHeader().setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeToContents)
            self.bmf_export_exception_add_btn = QtWidgets.QPushButton('Add Exception')
            self.bmf_export_exception_remove_btn = QtWidgets.QPushButton('Remove Selected')
            self.bmf_export_exception_summary = QtWidgets.QLabel('No BMF value exceptions configured.')
            self.bmf_export_exception_summary.setWordWrap(True)
            self.start_bmf_export_btn = QtWidgets.QPushButton('Export CSV To BMF')
            self.equation_finder_output_edit = QtWidgets.QLineEdit('')
            self.equation_finder_browse = QtWidgets.QPushButton('Browse')
            self.start_equation_finder_btn = QtWidgets.QPushButton('Find Equations')
            self._sample_blocks_auto_path = ''
            self._domain_samples_auto_path = ''
            self._block_value_transfer_auto_path = ''
            self._block_model_transfer_auto_path = ''
            self._table_attribute_auto_path = ''
            self._block_domain_metrics_auto_path = ''
            self._domain_interpolation_confidence_auto_path = ''
            self._block_volume_weighted_auto_path = ''
            self._bmf_export_input_auto_path = ''
            self._bmf_export_output_auto_path = ''
            self._equation_finder_auto_path = ''
            self._pending_block_value_transfer_cols = []
            self._pending_block_model_transfer_cols = []
            self._pending_table_attribute_key_cols = []
            self._pending_table_attribute_value_cols = []
            self._pending_bmf_export_value_cols = None
            self._pending_bmf_export_column_types = None
            self._bmf_export_value_cols_initialized = False
            self._bmf_export_columns_source_path = ''

            def refresh_sample_blocks_output_path(force=False):
                suggested = suggested_sample_blocks_path()
                current = self.sample_blocks_output_edit.text().strip()
                if force or not current or current == self._sample_blocks_auto_path:
                    self.sample_blocks_output_edit.setText(suggested)
                self._sample_blocks_auto_path = suggested

            def refresh_domain_samples_output_path(force=False):
                suggested = suggested_domain_samples_path()
                current = self.domain_samples_output_edit.text().strip()
                if force or not current or current == self._domain_samples_auto_path:
                    self.domain_samples_output_edit.setText(suggested)
                self._domain_samples_auto_path = suggested

            def refresh_block_value_transfer_output_path(force=False):
                suggested = suggested_block_value_transfer_path()
                current = self.block_value_transfer_output_edit.text().strip()
                if force or not current or current == self._block_value_transfer_auto_path:
                    self.block_value_transfer_output_edit.setText(suggested)
                self._block_value_transfer_auto_path = suggested

            self._refresh_block_value_transfer_output_path = refresh_block_value_transfer_output_path

            def refresh_block_model_transfer_output_path(force=False):
                suggested = suggested_block_model_transfer_path()
                current = self.block_model_transfer_output_edit.text().strip()
                if force or not current or current == self._block_model_transfer_auto_path:
                    self.block_model_transfer_output_edit.setText(suggested)
                self._block_model_transfer_auto_path = suggested

            self._refresh_block_model_transfer_output_path = refresh_block_model_transfer_output_path

            def refresh_table_attribute_output_path(force=False):
                suggested = suggested_table_attribute_output_path()
                current = self.table_attribute_output_edit.text().strip()
                if force or not current or current == self._table_attribute_auto_path:
                    self.table_attribute_output_edit.setText(suggested)
                self._table_attribute_auto_path = suggested

            self._refresh_table_attribute_output_path = refresh_table_attribute_output_path

            def refresh_block_domain_metrics_output_path(force=False):
                suggested = suggested_block_domain_metrics_path()
                current = self.block_domain_metrics_output_edit.text().strip()
                if force or not current or current == self._block_domain_metrics_auto_path:
                    self.block_domain_metrics_output_edit.setText(suggested)
                self._block_domain_metrics_auto_path = suggested

            def refresh_domain_interpolation_confidence_output_path(force=False):
                suggested = suggested_domain_interpolation_confidence_path()
                current = self.domain_interpolation_confidence_output_edit.text().strip()
                if force or not current or current == self._domain_interpolation_confidence_auto_path:
                    self.domain_interpolation_confidence_output_edit.setText(suggested)
                self._domain_interpolation_confidence_auto_path = suggested

            def refresh_block_volume_weighted_output_path(force=False):
                suggested = suggested_block_volume_weighted_path()
                current = self.block_volume_weighted_output_edit.text().strip()
                if force or not current or current == self._block_volume_weighted_auto_path:
                    self.block_volume_weighted_output_edit.setText(suggested)
                self._block_volume_weighted_auto_path = suggested

            def refresh_bmf_export_input_path(force=False):
                suggested = suggested_bmf_export_input_path()
                current = self.bmf_export_input_edit.text().strip()
                if force or not current or current == self._bmf_export_input_auto_path:
                    self.bmf_export_input_edit.setText(suggested)
                self._bmf_export_input_auto_path = suggested

            def refresh_bmf_export_output_path(force=False):
                suggested = suggested_bmf_export_output_path()
                current = self.bmf_export_output_edit.text().strip()
                if force or not current or current == self._bmf_export_output_auto_path:
                    self.bmf_export_output_edit.setText(suggested)
                self._bmf_export_output_auto_path = suggested

            def sync_bmf_export_input_settings():
                path = self.bmf_export_input_edit.text().strip()
                if not path or not os.path.isfile(path):
                    return
                detected = detect_csv_delimiter(path)
                if detected in delim_opts and self.bmf_export_delim.currentText() != detected:
                    self.bmf_export_delim.setCurrentText(detected)

            def refresh_equation_finder_output_path(force=False):
                suggested = suggested_equation_finder_path()
                current = self.equation_finder_output_edit.text().strip()
                if force or not current or current == self._equation_finder_auto_path:
                    self.equation_finder_output_edit.setText(suggested)
                self._equation_finder_auto_path = suggested

            def browse_sample_blocks_output():
                initial_path = self.sample_blocks_output_edit.text().strip() or suggested_sample_blocks_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Sample Blocks Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.sample_blocks_output_edit.setText(path)

            def browse_domain_samples_output():
                initial_path = self.domain_samples_output_edit.text().strip() or suggested_domain_samples_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Domain Samples Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.domain_samples_output_edit.setText(path)

            def browse_block_value_transfer_output():
                initial_path = self.block_value_transfer_output_edit.text().strip() or suggested_block_value_transfer_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Block Value Transfer Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.block_value_transfer_output_edit.setText(path)

            def browse_block_model_target():
                path, _ = QtWidgets.QFileDialog.getOpenFileName(
                    self,
                    'Target Block Model',
                    self.block_model_target_edit.text().strip() or '.',
                    'Block Models (*.csv *.bmf);;All Files (*.*)',
                )
                if path:
                    self.block_model_target_edit.setText(path)

            def browse_block_model_transfer_output():
                initial_path = self.block_model_transfer_output_edit.text().strip() or suggested_block_model_transfer_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Block Model Transfer Output File',
                    initial_path,
                    'CSV Files (*.csv)',
                )
                if path:
                    self.block_model_transfer_output_edit.setText(path)

            def browse_table_attribute_block_model():
                path, _ = QtWidgets.QFileDialog.getOpenFileName(
                    self,
                    'Block Model To Enrich',
                    self.table_attribute_block_model_edit.text().strip() or '.',
                    'Block Models (*.csv *.bmf);;All Files (*.*)',
                )
                if path:
                    self.table_attribute_block_model_edit.setText(path)

            def browse_table_attribute_table():
                path, _ = QtWidgets.QFileDialog.getOpenFileName(
                    self,
                    'Attribute Table CSV',
                    self.table_attribute_table_edit.text().strip() or '.',
                    'CSV Files (*.csv);;All Files (*.*)',
                )
                if path:
                    self.table_attribute_table_edit.setText(path)

            def browse_table_attribute_output():
                initial_path = self.table_attribute_output_edit.text().strip() or suggested_table_attribute_output_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Assign Attributes Output File',
                    initial_path,
                    'CSV Files (*.csv)',
                )
                if path:
                    self.table_attribute_output_edit.setText(path)

            def browse_block_domain_metrics_output():
                initial_path = self.block_domain_metrics_output_edit.text().strip() or suggested_block_domain_metrics_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Block Domain Metrics Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.block_domain_metrics_output_edit.setText(path)

            def browse_domain_interpolation_confidence_output():
                initial_path = self.domain_interpolation_confidence_output_edit.text().strip() or suggested_domain_interpolation_confidence_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Domain Interpolation Confidence Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.domain_interpolation_confidence_output_edit.setText(path)

            def browse_block_volume_weighted_output():
                initial_path = self.block_volume_weighted_output_edit.text().strip() or suggested_block_volume_weighted_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Block Volume Weighted Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.block_volume_weighted_output_edit.setText(path)

            def browse_bmf_export_input():
                initial_path = self.bmf_export_input_edit.text().strip() or suggested_bmf_export_input_path()
                path, _ = QtWidgets.QFileDialog.getOpenFileName(
                    self,
                    'CSV Grid Input File',
                    initial_path,
                    'CSV Files (*.csv);;All Files (*.*)'
                )
                if path:
                    self.bmf_export_input_edit.setText(path)

            def browse_bmf_export_output():
                initial_path = self.bmf_export_output_edit.text().strip() or suggested_bmf_export_output_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'BMF Output File',
                    initial_path,
                    'BMF Files (*.bmf);;All Files (*.*)'
                )
                if path:
                    self.bmf_export_output_edit.setText(path)

            def browse_equation_finder_output():
                initial_path = self.equation_finder_output_edit.text().strip() or suggested_equation_finder_path()
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self,
                    'Equation Finder Output File',
                    initial_path,
                    'CSV Files (*.csv)'
                )
                if path:
                    self.equation_finder_output_edit.setText(path)

            def create_collapsible_operation_section(title, tooltip=''):
                section_widget = QtWidgets.QWidget()
                section_layout = QtWidgets.QVBoxLayout(section_widget)
                section_layout.setContentsMargins(0, 0, 0, 0)
                section_layout.setSpacing(4)

                header_widget = QtWidgets.QWidget()
                header_layout = QtWidgets.QHBoxLayout(header_widget)
                header_layout.setContentsMargins(0, 0, 0, 0)
                header_layout.setSpacing(6)

                toggle_button = QtWidgets.QToolButton()
                toggle_button.setText('+')
                toggle_button.setCheckable(True)
                toggle_button.setChecked(False)
                toggle_button.setFixedSize(20, 20)
                toggle_button.setToolTip(tooltip)

                title_label = QtWidgets.QLabel(title)
                title_font = title_label.font()
                title_font.setBold(True)
                title_label.setFont(title_font)
                title_label.setToolTip(tooltip)

                header_layout.addWidget(toggle_button)
                header_layout.addWidget(title_label)
                header_layout.addStretch(1)

                content_widget = QtWidgets.QWidget()
                content_widget.setToolTip(tooltip)
                content_layout = QtWidgets.QFormLayout()
                content_layout.setContentsMargins(28, 0, 0, 0)
                content_widget.setLayout(content_layout)
                content_widget.setVisible(False)

                def set_expanded(expanded):
                    toggle_button.setText('-' if expanded else '+')
                    content_widget.setVisible(bool(expanded))

                toggle_button.toggled.connect(set_expanded)
                set_expanded(False)

                section_widget.setToolTip(tooltip)
                section_layout.addWidget(header_widget)
                section_layout.addWidget(content_widget)
                return section_widget, content_layout

            sample_blocks_group, sample_blocks_form = create_collapsible_operation_section(
                'Sample Blocks',
                'Aggregate samples into their containing blocks and export one row per sampled block with the assigned mean sample value, or a weighted average when a sample weight column is configured.',
            )
            sample_blocks_output_layout = QtWidgets.QHBoxLayout()
            sample_blocks_output_layout.addWidget(self.sample_blocks_output_edit)
            sample_blocks_output_layout.addWidget(self.sample_blocks_browse)
            sample_blocks_form.addRow('Output File', sample_blocks_output_layout)
            sample_blocks_form.addRow('Sample ID Metadata', self.sample_blocks_include_ids)
            sample_blocks_id_layout = QtWidgets.QHBoxLayout()
            for index, cb in enumerate(self.sample_blocks_id_cols):
                sample_blocks_id_layout.addWidget(cb)
                sample_blocks_id_layout.setStretch(index, 1)
            sample_blocks_form.addRow('Sample ID Columns', sample_blocks_id_layout)
            self.start_sample_blocks_btn.setToolTip(
                'Export one row per grid block that contains samples, using the same sample-block aggregation created at the start of interpolation, including weighted averages when a sample weight column is configured.'
            )
            sample_blocks_form.addRow('', self.start_sample_blocks_btn)
            operations_form.addRow(sample_blocks_group)

            domain_samples_group, domain_samples_form = create_collapsible_operation_section(
                'Domain Samples',
                'Assign a domain to each sample row from the blocks model and export the domained samples to a new CSV.',
            )
            domain_output_layout = QtWidgets.QHBoxLayout()
            domain_output_layout.addWidget(self.domain_samples_output_edit)
            domain_output_layout.addWidget(self.domain_samples_browse)
            domain_samples_form.addRow('Output File', domain_output_layout)
            self.start_domaining_btn.setToolTip(
                'Export a copy of the samples file with one domain assigned to each sample from the configured blocks and block size.'
            )
            domain_samples_form.addRow('', self.start_domaining_btn)
            operations_form.addRow(domain_samples_group)

            block_value_transfer_group, block_value_transfer_form = create_collapsible_operation_section(
                'Assign Block Columns To Samples',
                'Copy selected block columns onto sample rows by matching each sample to its containing block.',
            )
            block_value_transfer_output_layout = QtWidgets.QHBoxLayout()
            block_value_transfer_output_layout.addWidget(self.block_value_transfer_output_edit)
            block_value_transfer_output_layout.addWidget(self.block_value_transfer_browse)
            block_value_transfer_form.addRow('Output File', block_value_transfer_output_layout)
            block_value_transfer_form.addRow('Block Columns', self.block_value_transfer_cols)
            block_value_transfer_buttons = QtWidgets.QHBoxLayout()
            block_value_transfer_buttons.addWidget(self.block_value_transfer_select_all_btn)
            block_value_transfer_buttons.addWidget(self.block_value_transfer_clear_btn)
            block_value_transfer_form.addRow('', block_value_transfer_buttons)
            block_value_transfer_form.addRow('', self.block_value_transfer_summary)
            self.start_block_value_transfer_btn.setToolTip(
                'Export the samples file with the selected block columns transferred onto each matched sample row.'
            )
            block_value_transfer_form.addRow('', self.start_block_value_transfer_btn)
            operations_form.addRow(block_value_transfer_group)

            block_model_transfer_group, block_model_transfer_form = create_collapsible_operation_section(
                'Assign Block Columns To Block Model',
                'Enrich the existing rows of a target block model from overlapping source blocks. Numeric fields are overlap-volume weighted; categorical fields use the category with the greatest overlap volume.',
            )
            self.block_model_target_edit.setToolTip('Target block model to enrich. The export preserves these rows and appends or replaces only the selected source columns.')
            self.block_model_target_browse.setToolTip('Choose the target block model CSV or BMF file.')
            self.block_model_target_delim.setToolTip('CSV delimiter for the target model. When a CSV is opened, the UI auto-detects the delimiter.')
            self.block_model_target_header_line.setToolTip('Header line for the target CSV. When a CSV is opened, this control is synchronized to the effective detected header line.')
            self.block_model_target_x_col.setToolTip('Target X coordinate column.')
            self.block_model_target_y_col.setToolTip('Target Y coordinate column.')
            self.block_model_target_z_col.setToolTip('Target Z coordinate column.')
            for combo, tip in zip(
                self.block_model_source_size_cols,
                ['Optional source DX column. Prefer explicit size columns for source sub-block models.',
                 'Optional source DY column. Prefer explicit size columns for source sub-block models.',
                 'Optional source DZ column. Prefer explicit size columns for source sub-block models.'],
            ):
                combo.setToolTip(tip)
            for combo, tip in zip(
                self.block_model_target_size_cols,
                ['Optional target DX column. Leave on Infer to use the target base size below.',
                 'Optional target DY column. Leave on Infer to use the target base size below.',
                 'Optional target DZ column. Leave on Infer to use the target base size below.'],
            ):
                combo.setToolTip(tip)
            for spin, axis in zip(self.block_model_target_size_spins, ['X', 'Y', 'Z']):
                spin.setToolTip(f'Fallback target base block size on {axis} used when no explicit target size column is selected.')
            self.block_model_transfer_cols.setToolTip('Source columns to transfer onto the target block rows. Numeric columns use overlap-volume-weighted averages; text columns use the category with the greatest overlap volume.')
            self.block_model_transfer_select_all_btn.setToolTip('Select every transferable source column.')
            self.block_model_transfer_clear_btn.setToolTip('Clear the source-column selection.')
            self.block_model_transfer_summary.setToolTip('Shows the current source-column selection and no-overlap fallback settings.')
            self.block_model_nearest_fallback.setToolTip('If checked, targets with no overlapping source block can inherit values from the nearest source-block center, subject to the distance limit.')
            self.block_model_nearest_max_distance.setToolTip('Maximum center-to-center distance allowed for the nearest fallback. Set to Unlimited to allow any distance.')
            self.block_model_transfer_output_edit.setToolTip('CSV file to write. The output keeps the target geometry rows and adds the selected source columns.')
            self.block_model_transfer_output_browse.setToolTip('Choose the CSV file to write for the enriched target model.')
            target_file_layout = QtWidgets.QHBoxLayout()
            target_file_layout.addWidget(self.block_model_target_edit)
            target_file_layout.addWidget(self.block_model_target_browse)
            block_model_transfer_form.addRow('Target Block Model', target_file_layout)
            target_format_layout = QtWidgets.QHBoxLayout()
            target_format_layout.addWidget(self.block_model_target_delim)
            target_format_layout.addWidget(self.block_model_target_header_line)
            block_model_transfer_form.addRow('Target Delimiter / Header', target_format_layout)
            target_coords_layout = QtWidgets.QHBoxLayout()
            for combo in [self.block_model_target_x_col, self.block_model_target_y_col, self.block_model_target_z_col]:
                target_coords_layout.addWidget(combo)
            block_model_transfer_form.addRow('Target X / Y / Z', target_coords_layout)
            source_sizes_layout = QtWidgets.QHBoxLayout()
            target_sizes_layout = QtWidgets.QHBoxLayout()
            for combo in self.block_model_source_size_cols:
                source_sizes_layout.addWidget(combo)
            for combo in self.block_model_target_size_cols:
                target_sizes_layout.addWidget(combo)
            block_model_transfer_form.addRow('Source DX / DY / DZ', source_sizes_layout)
            block_model_transfer_form.addRow('Target DX / DY / DZ', target_sizes_layout)
            target_base_size_layout = QtWidgets.QHBoxLayout()
            for spin in self.block_model_target_size_spins:
                target_base_size_layout.addWidget(spin)
            block_model_transfer_form.addRow('Target Base Size X / Y / Z', target_base_size_layout)
            block_model_transfer_form.addRow('Source Columns', self.block_model_transfer_cols)
            block_model_transfer_buttons = QtWidgets.QHBoxLayout()
            block_model_transfer_buttons.addWidget(self.block_model_transfer_select_all_btn)
            block_model_transfer_buttons.addWidget(self.block_model_transfer_clear_btn)
            block_model_transfer_form.addRow('', block_model_transfer_buttons)
            block_model_transfer_form.addRow('', self.block_model_transfer_summary)
            no_overlap_layout = QtWidgets.QHBoxLayout()
            no_overlap_layout.addWidget(self.block_model_nearest_fallback)
            no_overlap_layout.addWidget(self.block_model_nearest_max_distance)
            block_model_transfer_form.addRow('No-overlap Policy', no_overlap_layout)
            block_model_transfer_output_layout = QtWidgets.QHBoxLayout()
            block_model_transfer_output_layout.addWidget(self.block_model_transfer_output_edit)
            block_model_transfer_output_layout.addWidget(self.block_model_transfer_output_browse)
            block_model_transfer_form.addRow('Output File', block_model_transfer_output_layout)
            self.start_block_model_transfer_btn.setToolTip(
                'Preserve every target row and add the selected source columns using exact 3-D overlap. Explicit DX/DY/DZ columns are preferred for sub-block models; otherwise dimensions are inferred inside each configured base block.'
            )
            block_model_transfer_form.addRow('', self.start_block_model_transfer_btn)
            operations_form.addRow(block_model_transfer_group)

            table_attribute_group, table_attribute_form = create_collapsible_operation_section(
                'Assign Attributes From Table',
                'Join one or more CSV table columns onto a block model by matching one or more shared key columns.',
            )
            self.table_attribute_block_model_edit.setToolTip('Block model CSV or BMF whose rows will be preserved and enriched.')
            self.table_attribute_block_model_browse.setToolTip('Choose the block model to enrich.')
            self.table_attribute_block_model_delim.setToolTip('CSV delimiter for the block model. When a CSV is opened, the UI auto-detects the delimiter.')
            self.table_attribute_block_model_header_line.setToolTip('Header line for the block model CSV. When a CSV is opened, this control is synchronized to the effective detected header line.')
            self.table_attribute_table_edit.setToolTip('CSV table that provides key columns and attribute values to assign.')
            self.table_attribute_table_browse.setToolTip('Choose the CSV table to join onto the block model.')
            self.table_attribute_table_delim.setToolTip('CSV delimiter for the attribute table. The UI auto-detects it when the file is opened.')
            self.table_attribute_table_header_line.setToolTip('Header line for the attribute table CSV. The UI syncs this to the effective detected header line.')
            self.table_attribute_key_cols.setToolTip('Select one or more shared columns used to match block-model rows to table rows.')
            self.table_attribute_value_cols.setToolTip('Select one or more table columns to copy onto the block model.')
            self.table_attribute_summary.setToolTip('Shows how many keys and table columns are selected for assignment.')
            self.table_attribute_output_edit.setToolTip('CSV file to write. The output preserves all block-model rows and appends or updates the selected table columns.')
            self.table_attribute_output_browse.setToolTip('Choose the CSV file to write for the enriched block model.')
            table_attribute_block_model_layout = QtWidgets.QHBoxLayout()
            table_attribute_block_model_layout.addWidget(self.table_attribute_block_model_edit)
            table_attribute_block_model_layout.addWidget(self.table_attribute_block_model_browse)
            table_attribute_form.addRow('Block Model', table_attribute_block_model_layout)
            table_attribute_block_model_format_layout = QtWidgets.QHBoxLayout()
            table_attribute_block_model_format_layout.addWidget(self.table_attribute_block_model_delim)
            table_attribute_block_model_format_layout.addWidget(self.table_attribute_block_model_header_line)
            table_attribute_form.addRow('Block Delimiter / Header', table_attribute_block_model_format_layout)
            table_attribute_table_layout = QtWidgets.QHBoxLayout()
            table_attribute_table_layout.addWidget(self.table_attribute_table_edit)
            table_attribute_table_layout.addWidget(self.table_attribute_table_browse)
            table_attribute_form.addRow('Attribute Table', table_attribute_table_layout)
            table_attribute_table_format_layout = QtWidgets.QHBoxLayout()
            table_attribute_table_format_layout.addWidget(self.table_attribute_table_delim)
            table_attribute_table_format_layout.addWidget(self.table_attribute_table_header_line)
            table_attribute_form.addRow('Table Delimiter / Header', table_attribute_table_format_layout)
            table_attribute_form.addRow('Match Keys', self.table_attribute_key_cols)
            table_attribute_key_buttons = QtWidgets.QHBoxLayout()
            table_attribute_key_buttons.addWidget(self.table_attribute_key_select_all_btn)
            table_attribute_key_buttons.addWidget(self.table_attribute_key_clear_btn)
            table_attribute_form.addRow('', table_attribute_key_buttons)
            table_attribute_form.addRow('Assign Columns', self.table_attribute_value_cols)
            table_attribute_value_buttons = QtWidgets.QHBoxLayout()
            table_attribute_value_buttons.addWidget(self.table_attribute_value_select_all_btn)
            table_attribute_value_buttons.addWidget(self.table_attribute_value_clear_btn)
            table_attribute_form.addRow('', table_attribute_value_buttons)
            table_attribute_form.addRow('', self.table_attribute_summary)
            table_attribute_output_layout = QtWidgets.QHBoxLayout()
            table_attribute_output_layout.addWidget(self.table_attribute_output_edit)
            table_attribute_output_layout.addWidget(self.table_attribute_output_browse)
            table_attribute_form.addRow('Output File', table_attribute_output_layout)
            table_attribute_form.addRow('', self.start_table_attribute_assign_btn)
            operations_form.addRow(table_attribute_group)

            block_metrics_group, block_metrics_form = create_collapsible_operation_section(
                'Block Domain Sample Metrics',
                'Export selectable per-block and per-domain support metrics, including exact nearest distance, scalable k-nearest averages, optional nearest-sample residuals, and an optional distance-threshold summary.',
            )
            block_metrics_output_layout = QtWidgets.QHBoxLayout()
            block_metrics_output_layout.addWidget(self.block_domain_metrics_output_edit)
            block_metrics_output_layout.addWidget(self.block_domain_metrics_browse)
            block_metrics_form.addRow('Output File', block_metrics_output_layout)
            block_metrics_form.addRow('Block Value Column', self.block_domain_metrics_value_col)
            block_metrics_form.addRow('Metric Name Prefix', self.block_domain_metrics_prefix_with_block_value)
            self.block_domain_metrics_metric_list = QtWidgets.QListWidget()
            self.block_domain_metrics_metric_list.setAlternatingRowColors(True)
            self.block_domain_metrics_metric_list.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)
            self.block_domain_metrics_metric_list.setWordWrap(True)
            self.block_domain_metrics_metric_list.setMinimumHeight(240)
            for metric_definition in BLOCK_DOMAIN_METRIC_DEFINITIONS:
                item = QtWidgets.QListWidgetItem(
                    f"{metric_definition['label']}: {metric_definition['description']} Cost: {metric_definition['cost_label']}."
                )
                item.setData(QtCore.Qt.UserRole, metric_definition['id'])
                item.setToolTip(f"{metric_definition['description']}\nCost: {metric_definition['cost_note']}")
                item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
                item.setCheckState(
                    QtCore.Qt.Checked if metric_definition.get('default_checked') else QtCore.Qt.Unchecked
                )
                self.block_domain_metrics_metric_list.addItem(item)
            self.block_domain_metrics_select_all_btn = QtWidgets.QPushButton('Select All')
            self.block_domain_metrics_clear_btn = QtWidgets.QPushButton('Clear')
            self.block_domain_metrics_metric_summary = QtWidgets.QLabel('')
            self.block_domain_metrics_metric_summary.setWordWrap(True)
            block_metrics_selection_buttons = QtWidgets.QHBoxLayout()
            block_metrics_selection_buttons.addWidget(self.block_domain_metrics_select_all_btn)
            block_metrics_selection_buttons.addWidget(self.block_domain_metrics_clear_btn)
            block_metrics_selection_buttons.addStretch(1)
            block_metrics_form.addRow('Metrics', self.block_domain_metrics_metric_list)
            block_metrics_form.addRow('', block_metrics_selection_buttons)
            block_metrics_form.addRow('', self.block_domain_metrics_metric_summary)
            self.block_domain_metrics_knn_k = QtWidgets.QSpinBox()
            self.block_domain_metrics_knn_k.setRange(1, 9999)
            self.block_domain_metrics_knn_k.setValue(8)
            self.block_domain_metrics_knn_k.setToolTip(
                'For Avg Distance KNN, average the distances to the k nearest same-domain samples instead of all samples in the domain.'
            )
            block_metrics_form.addRow('Avg Distance KNN k', self.block_domain_metrics_knn_k)
            block_metrics_id_layout = QtWidgets.QHBoxLayout()
            for index, cb in enumerate(self.block_domain_metrics_id_cols):
                block_metrics_id_layout.addWidget(cb)
                block_metrics_id_layout.setStretch(index, 1)
            block_metrics_form.addRow('Closest Sample ID Columns', block_metrics_id_layout)
            self.block_domain_metrics_distance_step = TrimmedDisplayDoubleSpinBox()
            self.block_domain_metrics_distance_step.setDecimals(3)
            self.block_domain_metrics_distance_step.setRange(0.001, 1_000_000.0)
            self.block_domain_metrics_distance_step.setValue(50.0)
            self.block_domain_metrics_distance_step.setSuffix(' m')
            self.block_domain_metrics_distance_step.setToolTip('Distance step used to build the distance-threshold summary when that metric is selected.')
            self.block_domain_metrics_max_factor = QtWidgets.QSpinBox()
            self.block_domain_metrics_max_factor.setRange(1, 999)
            self.block_domain_metrics_max_factor.setValue(5)
            self.block_domain_metrics_max_factor.setToolTip('Maximum multiple of the distance step before the >= overflow bucket.')
            block_metrics_distance_layout = QtWidgets.QHBoxLayout()
            block_metrics_distance_layout.addWidget(self.block_domain_metrics_distance_step)
            block_metrics_distance_layout.addWidget(self.block_domain_metrics_max_factor)
            block_metrics_form.addRow('Distance Summary Thresholds', block_metrics_distance_layout)
            self.start_block_domain_metrics_btn.setToolTip(
                'Export the checked block-domain metrics. Nearest-neighbor metrics use cKDTree; exact Avg Distance to all domain samples remains available but is expensive on large domains.'
            )
            block_metrics_form.addRow('', self.start_block_domain_metrics_btn)
            operations_form.addRow(block_metrics_group)

            domain_confidence_group, domain_confidence_form = create_collapsible_operation_section(
                'Domain Interpolation Confidence',
                'Summarize each domain with average sample spacing, average distance from blocks to samples, and their ratio.',
            )
            domain_confidence_output_layout = QtWidgets.QHBoxLayout()
            domain_confidence_output_layout.addWidget(self.domain_interpolation_confidence_output_edit)
            domain_confidence_output_layout.addWidget(self.domain_interpolation_confidence_browse)
            domain_confidence_form.addRow('Output File', domain_confidence_output_layout)
            self.start_domain_interpolation_confidence_btn.setToolTip(
                'Export one row per domain with average distance between samples, average distance from blocks to samples, and the ratio between those two metrics.'
            )
            domain_confidence_form.addRow('', self.start_domain_interpolation_confidence_btn)
            operations_form.addRow(domain_confidence_group)

            block_volume_group, block_volume_form = create_collapsible_operation_section(
                'Block Volume Weighted Average',
                'Infer per-row block volumes when needed and export a volume-weighted average for the selected blocks column.',
            )
            block_volume_output_layout = QtWidgets.QHBoxLayout()
            block_volume_output_layout.addWidget(self.block_volume_weighted_output_edit)
            block_volume_output_layout.addWidget(self.block_volume_weighted_browse)
            block_volume_form.addRow('Output File', block_volume_output_layout)
            block_volume_form.addRow('Block Value Column', self.block_volume_weighted_value_col)
            block_volume_form.addRow('Weight Column', self.block_weight_metric_col)
            self.start_block_volume_weighted_btn.setToolTip(
                'Export a summary CSV containing inferred total volume, weighted sum, and weighted average for the selected blocks field.'
            )
            block_volume_form.addRow('', self.start_block_volume_weighted_btn)
            operations_form.addRow(block_volume_group)

            bmf_export_group, bmf_export_form = create_collapsible_operation_section(
                'CSV Grid To BMF',
                'Convert a regular CSV grid into the experimental TBMS2.0 BMF writer backend used by the standalone exporter.',
            )
            bmf_export_input_layout = QtWidgets.QHBoxLayout()
            bmf_export_input_layout.addWidget(self.bmf_export_input_edit)
            bmf_export_input_layout.addWidget(self.bmf_export_input_browse)
            bmf_export_form.addRow('Input CSV', bmf_export_input_layout)
            bmf_export_output_layout = QtWidgets.QHBoxLayout()
            bmf_export_output_layout.addWidget(self.bmf_export_output_edit)
            bmf_export_output_layout.addWidget(self.bmf_export_output_browse)
            bmf_export_form.addRow('Output BMF', bmf_export_output_layout)
            bmf_export_form.addRow('Input Delimiter', self.bmf_export_delim)
            bmf_export_form.addRow('Input Header Line', self.bmf_export_header_line)
            bmf_export_form.addRow('Backend', self.bmf_export_backend_combo)
            bmf_export_cell_layout = QtWidgets.QHBoxLayout()
            bmf_export_cell_layout.addWidget(self.bmf_export_cell_x)
            bmf_export_cell_layout.addWidget(self.bmf_export_cell_y)
            bmf_export_cell_layout.addWidget(self.bmf_export_cell_z)
            bmf_export_form.addRow('Cell Size X / Y / Z', bmf_export_cell_layout)
            bmf_export_regularize_layout = QtWidgets.QVBoxLayout()
            bmf_export_regularize_layout.addWidget(self.bmf_export_regularize_base_blocks)
            bmf_export_regularize_layout.addWidget(self.bmf_export_regularize_warning)
            bmf_export_form.addRow('', bmf_export_regularize_layout)
            bmf_export_coords_layout = QtWidgets.QHBoxLayout()
            bmf_export_coords_layout.addWidget(self.bmf_export_x_col)
            bmf_export_coords_layout.addWidget(self.bmf_export_y_col)
            bmf_export_coords_layout.addWidget(self.bmf_export_z_col)
            bmf_export_form.addRow('X / Y / Z Columns', bmf_export_coords_layout)
            bmf_export_form.addRow('Value Columns', self.bmf_export_value_cols)
            bmf_export_value_buttons = QtWidgets.QHBoxLayout()
            bmf_export_value_buttons.addWidget(self.bmf_export_select_all_btn)
            bmf_export_value_buttons.addWidget(self.bmf_export_clear_btn)
            bmf_export_form.addRow('', bmf_export_value_buttons)
            bmf_export_form.addRow('', self.bmf_export_value_summary)
            bmf_export_form.addRow('Value Exceptions', self.bmf_export_exception_table)
            bmf_export_exception_buttons = QtWidgets.QHBoxLayout()
            bmf_export_exception_buttons.addWidget(self.bmf_export_exception_add_btn)
            bmf_export_exception_buttons.addWidget(self.bmf_export_exception_remove_btn)
            bmf_export_form.addRow('', bmf_export_exception_buttons)
            bmf_export_form.addRow('', self.bmf_export_exception_summary)
            self.start_bmf_export_btn.setToolTip(
                'Export the selected CSV grid to BMF using the standalone reverse-engineered backend.'
            )
            bmf_export_form.addRow('', self.start_bmf_export_btn)
            operations_form.addRow(bmf_export_group)

            equation_finder_group, equation_finder_form = create_collapsible_operation_section(
                'Equation Finder By Domain',
                'Run symbolic regression independently on each domain to search for equations that explain the selected target from the chosen predictors.',
            )
            equation_output_layout = QtWidgets.QHBoxLayout()
            equation_output_layout.addWidget(self.equation_finder_output_edit)
            equation_output_layout.addWidget(self.equation_finder_browse)
            equation_finder_form.addRow('Summary Output', equation_output_layout)
            self.equation_include_coords = QtWidgets.QCheckBox('Include coordinate columns (x, y, z)')
            self.equation_include_coords.setChecked(True)
            self.equation_include_coords.setToolTip('When enabled, the configured sample coordinate columns are included in the selectable numeric predictor list and selected by default.')
            equation_finder_form.addRow('Coordinate Predictors', self.equation_include_coords)
            self.equation_predictor_list = QtWidgets.QListWidget()
            self.equation_predictor_list.setSelectionMode(QtWidgets.QAbstractItemView.MultiSelection)
            self.equation_predictor_list.setMinimumHeight(140)
            equation_finder_form.addRow('Predictor Columns', self.equation_predictor_list)
            equation_predictor_buttons = QtWidgets.QHBoxLayout()
            self.equation_refresh_predictors_btn = QtWidgets.QPushButton('Refresh')
            self.equation_select_all_predictors_btn = QtWidgets.QPushButton('Select All')
            self.equation_clear_predictors_btn = QtWidgets.QPushButton('Clear')
            equation_predictor_buttons.addWidget(self.equation_refresh_predictors_btn)
            equation_predictor_buttons.addWidget(self.equation_select_all_predictors_btn)
            equation_predictor_buttons.addWidget(self.equation_clear_predictors_btn)
            equation_finder_form.addRow('', equation_predictor_buttons)
            self.equation_predictor_summary = QtWidgets.QLabel('No numeric predictor columns loaded.')
            self.equation_predictor_summary.setWordWrap(True)
            equation_finder_form.addRow('', self.equation_predictor_summary)
            self.equation_min_samples_per_domain = QtWidgets.QSpinBox()
            self.equation_min_samples_per_domain.setRange(2, 1_000_000)
            self.equation_min_samples_per_domain.setValue(25)
            self.equation_min_samples_per_domain.setToolTip('Domains with fewer rows than this are skipped.')
            equation_finder_form.addRow('Min Samples Per Domain', self.equation_min_samples_per_domain)
            self.equation_validation_fraction = QtWidgets.QDoubleSpinBox()
            self.equation_validation_fraction.setRange(0.0, 0.5)
            self.equation_validation_fraction.setSingleStep(0.05)
            self.equation_validation_fraction.setDecimals(2)
            self.equation_validation_fraction.setValue(0.20)
            self.equation_validation_fraction.setToolTip('Fraction of each domain reserved for validation metrics. Set to 0 to score on the full domain only.')
            equation_finder_form.addRow('Validation Fraction', self.equation_validation_fraction)
            self.equation_max_iterations = QtWidgets.QSpinBox()
            self.equation_max_iterations.setRange(1, 10_000_000)
            self.equation_max_iterations.setValue(100)
            self.equation_max_iterations.setToolTip('Maximum PySR search iterations per domain. The best equation found so far is returned when this limit is reached.')
            equation_finder_form.addRow('Max Iterations', self.equation_max_iterations)
            self.equation_timeout_seconds = QtWidgets.QSpinBox()
            self.equation_timeout_seconds.setRange(0, 86_400)
            self.equation_timeout_seconds.setValue(60)
            self.equation_timeout_seconds.setSuffix(' s')
            self.equation_timeout_seconds.setToolTip('Maximum wall-clock time per domain. Set to 0 to disable the timeout. The best equation found so far is returned when the timeout is reached.')
            equation_finder_form.addRow('Timeout', self.equation_timeout_seconds)
            self.start_equation_finder_btn.setToolTip(
                'Fit candidate equations per domain and export both the per-domain winner summary and detailed equation search results.'
            )
            equation_finder_form.addRow('', self.start_equation_finder_btn)
            operations_form.addRow(equation_finder_group)

            self.sample_blocks_browse.clicked.connect(browse_sample_blocks_output)
            self.start_sample_blocks_btn.clicked.connect(self.run_sample_blocks_only)
            self.sample_blocks_include_ids.toggled.connect(lambda _: self._update_sample_blocks_id_controls())
            self.domain_samples_browse.clicked.connect(browse_domain_samples_output)
            self.start_domaining_btn.clicked.connect(self.run_domain_samples_only)
            self.block_value_transfer_browse.clicked.connect(browse_block_value_transfer_output)
            self.start_block_value_transfer_btn.clicked.connect(self.run_block_value_transfer_only)
            self.block_value_transfer_select_all_btn.clicked.connect(self.block_value_transfer_cols.selectAll)
            self.block_value_transfer_select_all_btn.clicked.connect(self._update_block_value_transfer_summary)
            self.block_value_transfer_select_all_btn.clicked.connect(lambda: refresh_block_value_transfer_output_path())
            self.block_value_transfer_clear_btn.clicked.connect(self.block_value_transfer_cols.clearSelection)
            self.block_value_transfer_clear_btn.clicked.connect(self._update_block_value_transfer_summary)
            self.block_value_transfer_clear_btn.clicked.connect(lambda: refresh_block_value_transfer_output_path())
            self.block_value_transfer_cols.itemSelectionChanged.connect(self._update_block_value_transfer_summary)
            self.block_value_transfer_cols.itemSelectionChanged.connect(lambda: refresh_block_value_transfer_output_path())
            self.block_model_target_browse.clicked.connect(browse_block_model_target)
            self.block_model_transfer_output_browse.clicked.connect(browse_block_model_transfer_output)
            self.start_block_model_transfer_btn.clicked.connect(self.run_block_model_transfer_only)
            self.block_model_transfer_select_all_btn.clicked.connect(self.block_model_transfer_cols.selectAll)
            self.block_model_transfer_clear_btn.clicked.connect(self.block_model_transfer_cols.clearSelection)
            self.block_model_transfer_cols.itemSelectionChanged.connect(self._update_block_model_transfer_summary)
            self.block_model_nearest_fallback.toggled.connect(self._update_block_model_transfer_fallback_controls)
            self.block_model_nearest_max_distance.valueChanged.connect(lambda _: self._update_block_model_transfer_summary())
            self.block_model_target_edit.textChanged.connect(lambda _: self._refresh_block_model_target_columns())
            self.block_model_target_edit.textChanged.connect(lambda _: refresh_block_model_transfer_output_path())
            self.block_model_target_delim.currentIndexChanged.connect(self._refresh_block_model_target_columns)
            self.block_model_target_header_line.valueChanged.connect(lambda _: self._refresh_block_model_target_columns())
            self.blocks_edit.textChanged.connect(lambda _: self._refresh_block_model_transfer_source_columns())
            self.blocks_delim.currentIndexChanged.connect(self._refresh_block_model_transfer_source_columns)
            self.blocks_header_line.valueChanged.connect(lambda _: self._refresh_block_model_transfer_source_columns())
            self.block_x_col.currentTextChanged.connect(lambda _: self._refresh_block_model_transfer_source_columns())
            self.block_y_col.currentTextChanged.connect(lambda _: self._refresh_block_model_transfer_source_columns())
            self.block_z_col.currentTextChanged.connect(lambda _: self._refresh_block_model_transfer_source_columns())
            self.table_attribute_block_model_browse.clicked.connect(browse_table_attribute_block_model)
            self.table_attribute_table_browse.clicked.connect(browse_table_attribute_table)
            self.table_attribute_output_browse.clicked.connect(browse_table_attribute_output)
            self.start_table_attribute_assign_btn.clicked.connect(self.run_table_attribute_assignment_only)
            self.table_attribute_key_select_all_btn.clicked.connect(self.table_attribute_key_cols.selectAll)
            self.table_attribute_key_select_all_btn.clicked.connect(self._update_table_attribute_summary)
            self.table_attribute_key_select_all_btn.clicked.connect(lambda: self._refresh_table_attribute_value_columns())
            self.table_attribute_key_clear_btn.clicked.connect(self.table_attribute_key_cols.clearSelection)
            self.table_attribute_key_clear_btn.clicked.connect(self._update_table_attribute_summary)
            self.table_attribute_key_clear_btn.clicked.connect(lambda: self._refresh_table_attribute_value_columns())
            self.table_attribute_value_select_all_btn.clicked.connect(self.table_attribute_value_cols.selectAll)
            self.table_attribute_value_select_all_btn.clicked.connect(self._update_table_attribute_summary)
            self.table_attribute_value_clear_btn.clicked.connect(self.table_attribute_value_cols.clearSelection)
            self.table_attribute_value_clear_btn.clicked.connect(self._update_table_attribute_summary)
            self.table_attribute_key_cols.itemSelectionChanged.connect(self._refresh_table_attribute_value_columns)
            self.table_attribute_key_cols.itemSelectionChanged.connect(self._update_table_attribute_summary)
            self.table_attribute_value_cols.itemSelectionChanged.connect(self._update_table_attribute_summary)
            self.table_attribute_block_model_edit.textChanged.connect(lambda _: self._refresh_table_attribute_shared_key_columns())
            self.table_attribute_block_model_edit.textChanged.connect(lambda _: refresh_table_attribute_output_path())
            self.table_attribute_block_model_delim.currentIndexChanged.connect(self._refresh_table_attribute_shared_key_columns)
            self.table_attribute_block_model_header_line.valueChanged.connect(lambda _: self._refresh_table_attribute_shared_key_columns())
            self.table_attribute_table_edit.textChanged.connect(lambda _: self._refresh_table_attribute_shared_key_columns())
            self.table_attribute_table_edit.textChanged.connect(lambda _: refresh_table_attribute_output_path())
            self.table_attribute_table_delim.currentIndexChanged.connect(self._refresh_table_attribute_shared_key_columns)
            self.table_attribute_table_header_line.valueChanged.connect(lambda _: self._refresh_table_attribute_shared_key_columns())
            self.block_domain_metrics_browse.clicked.connect(browse_block_domain_metrics_output)
            self.configure_block_domain_metrics_filters_btn.clicked.connect(self.configure_block_domain_metrics_filters)
            self.start_block_domain_metrics_btn.clicked.connect(self.run_block_domain_sample_metrics_only)
            self.block_domain_metrics_select_all_btn.clicked.connect(lambda: self._set_all_block_domain_metrics_checked(True))
            self.block_domain_metrics_clear_btn.clicked.connect(lambda: self._set_all_block_domain_metrics_checked(False))
            self.block_domain_metrics_metric_list.itemChanged.connect(lambda _item: self._update_block_domain_metrics_metric_summary())
            self.samples_edit.textChanged.connect(lambda _: refresh_block_value_transfer_output_path())
            self.blocks_edit.textChanged.connect(lambda _: self._refresh_block_value_transfer_columns())
            self.blocks_delim.currentIndexChanged.connect(lambda _: self._refresh_block_value_transfer_columns())
            self.blocks_header_line.valueChanged.connect(lambda _: self._refresh_block_value_transfer_columns())
            self.block_x_col.currentTextChanged.connect(lambda _: self._refresh_block_value_transfer_columns())
            self.block_y_col.currentTextChanged.connect(lambda _: self._refresh_block_value_transfer_columns())
            self.block_z_col.currentTextChanged.connect(lambda _: self._refresh_block_value_transfer_columns())
            self.domain_interpolation_confidence_browse.clicked.connect(browse_domain_interpolation_confidence_output)
            self.start_domain_interpolation_confidence_btn.clicked.connect(self.run_domain_interpolation_confidence_only)
            self.block_volume_weighted_browse.clicked.connect(browse_block_volume_weighted_output)
            self.configure_block_volume_weighted_filters_btn.clicked.connect(self.configure_block_volume_weighted_filters)
            self.start_block_volume_weighted_btn.clicked.connect(self.run_block_volume_weighted_average_only)
            self.bmf_export_input_browse.clicked.connect(browse_bmf_export_input)
            self.bmf_export_output_browse.clicked.connect(browse_bmf_export_output)
            self.bmf_export_select_all_btn.clicked.connect(lambda: self._set_all_bmf_export_value_columns_checked(True))
            self.bmf_export_clear_btn.clicked.connect(lambda: self._set_all_bmf_export_value_columns_checked(False))
            self.bmf_export_value_cols.itemChanged.connect(self._update_bmf_export_value_summary)
            self.bmf_export_exception_add_btn.clicked.connect(self._add_bmf_export_exception_row)
            self.bmf_export_exception_remove_btn.clicked.connect(self._remove_selected_bmf_export_exception_rows)
            self.bmf_export_exception_table.itemChanged.connect(lambda _item: self._update_bmf_export_exception_summary())
            self.start_bmf_export_btn.clicked.connect(self.run_csv_to_bmf_export_only)
            self.bmf_export_regularize_base_blocks.toggled.connect(self._update_bmf_export_regularize_warning)
            self._update_bmf_export_regularize_warning()
            self.equation_finder_browse.clicked.connect(browse_equation_finder_output)
            self.start_equation_finder_btn.clicked.connect(self.run_equation_finder_only)
            self.equation_refresh_predictors_btn.clicked.connect(self._refresh_equation_finder_predictor_columns)
            self.equation_select_all_predictors_btn.clicked.connect(self.equation_predictor_list.selectAll)
            self.equation_select_all_predictors_btn.clicked.connect(self._update_equation_finder_predictor_summary)
            self.equation_clear_predictors_btn.clicked.connect(self.equation_predictor_list.clearSelection)
            self.equation_clear_predictors_btn.clicked.connect(self._update_equation_finder_predictor_summary)
            self.equation_predictor_list.itemSelectionChanged.connect(self._update_equation_finder_predictor_summary)
            self.equation_include_coords.toggled.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.samples_edit.textChanged.connect(lambda _: refresh_sample_blocks_output_path())
            self.samples_edit.textChanged.connect(lambda _: refresh_domain_samples_output_path())
            self.blocks_edit.textChanged.connect(lambda _: refresh_block_domain_metrics_output_path())
            self.blocks_edit.textChanged.connect(lambda _: refresh_domain_interpolation_confidence_output_path())
            self.blocks_edit.textChanged.connect(lambda _: refresh_block_volume_weighted_output_path())
            self.interp_edit.textChanged.connect(lambda _: refresh_bmf_export_input_path())
            self.blocks_edit.textChanged.connect(lambda _: refresh_bmf_export_input_path())
            self.bmf_export_input_edit.textChanged.connect(lambda _: sync_bmf_export_input_settings())
            self.bmf_export_input_edit.textChanged.connect(lambda _: self._refresh_bmf_export_columns())
            self.bmf_export_input_edit.textChanged.connect(lambda _: refresh_bmf_export_output_path())
            self.bmf_export_delim.currentIndexChanged.connect(self._refresh_bmf_export_columns)
            self.bmf_export_header_line.valueChanged.connect(lambda _: self._refresh_bmf_export_columns())
            self.bmf_export_x_col.currentTextChanged.connect(lambda _: self._refresh_bmf_export_columns())
            self.bmf_export_y_col.currentTextChanged.connect(lambda _: self._refresh_bmf_export_columns())
            self.bmf_export_z_col.currentTextChanged.connect(lambda _: self._refresh_bmf_export_columns())
            self.samples_edit.textChanged.connect(lambda _: refresh_equation_finder_output_path())
            self.block_domain_col.currentTextChanged.connect(lambda _: refresh_domain_samples_output_path())
            self.block_domain_col.currentTextChanged.connect(lambda _: refresh_block_domain_metrics_output_path())
            self.block_domain_col.currentTextChanged.connect(lambda _: refresh_domain_interpolation_confidence_output_path())
            self.block_volume_weighted_value_col.currentTextChanged.connect(lambda _: refresh_block_volume_weighted_output_path())
            self.sample_value_col.currentTextChanged.connect(lambda _: refresh_equation_finder_output_path())
            self.sample_domain_col.currentTextChanged.connect(lambda _: refresh_equation_finder_output_path())
            self.samples_edit.textChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.samples_delim.currentIndexChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.samples_header_line.valueChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.sample_value_col.currentTextChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.sample_domain_col.currentTextChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.sample_x_col.currentTextChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.sample_y_col.currentTextChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.sample_z_col.currentTextChanged.connect(lambda _: self._refresh_equation_finder_predictor_columns())
            self.block_domain_metrics_value_col.currentTextChanged.connect(lambda _: self._update_block_domain_metrics_metric_summary())
            self.sample_value_col.currentTextChanged.connect(lambda _: self._update_block_domain_metrics_metric_summary())
            refresh_sample_blocks_output_path(force=True)
            self._update_sample_blocks_id_controls()
            refresh_domain_samples_output_path(force=True)
            self._refresh_block_value_transfer_columns()
            refresh_block_value_transfer_output_path(force=True)
            self._refresh_block_model_transfer_source_columns()
            self._refresh_block_model_target_columns()
            refresh_block_model_transfer_output_path(force=True)
            self._update_block_model_transfer_fallback_controls()
            self._refresh_table_attribute_shared_key_columns()
            refresh_table_attribute_output_path(force=True)
            self._update_block_domain_metrics_metric_summary()
            refresh_block_domain_metrics_output_path(force=True)
            refresh_domain_interpolation_confidence_output_path(force=True)
            refresh_block_volume_weighted_output_path(force=True)
            refresh_bmf_export_input_path(force=True)
            sync_bmf_export_input_settings()
            self._refresh_bmf_export_columns()
            refresh_bmf_export_output_path(force=True)
            refresh_equation_finder_output_path(force=True)
            self._refresh_equation_finder_predictor_columns()
            self._update_block_domain_metrics_filters_summary()
            self._update_block_volume_weighted_filters_summary()

            # === ANT COLONY TAB ===
            def dbl_spin(default, minv, maxv, step=0.1):
                s = QtWidgets.QDoubleSpinBox(); s.setRange(minv, maxv); s.setSingleStep(step); s.setValue(default); return s
            def int_spin(default, minv, maxv, step=1):
                s = QtWidgets.QSpinBox(); s.setRange(minv, maxv); s.setSingleStep(step); s.setValue(default); return s

            self.range_size = dbl_spin(0.2, 0.0001, 1e6, 0.1)
            self.range_size.setToolTip('Value tolerance / Bin size for grade classification.\nUnit: Grade units (e.g., %, ppm).\nLarger: Ants are more tolerant of grade differences (broader categories).\nSmaller: Ants are stricter about grade similarity (finer categories).')

            self.max_pheromone = int_spin(1000, 1, 10_000_000, 10)
            self.max_pheromone.setToolTip('Maximum pheromone level allowed in a single block.\nDetermines how far ants can travel and how long blocks persist.\nEquivalent to the maximum lifespan of a block in iterations.')

            self.ants_per_sample = int_spin(16, 1, 10000, 1)
            self.ants_per_sample.setToolTip('Number of ants spawned from each sample point per iteration.\nMore ants = better coverage but slower performance.')

            self.ants_sampling_percentage = dbl_spin(100.0, 0.001, 100.0, 1.0)
            self.ants_sampling_percentage.setToolTip('Percentage of samples that will spawn ants (0-100%).\nUseful for second pass interpolation where input blocks are numerous.')
            
            self.pheromone_decay_rate = int_spin(1, 0, 10, 1)
            self.pheromone_decay_rate.setToolTip('Amount of pheromone decay per iteration.\n0 = No decay (blocks persist indefinitely).\nHigher = Faster decay (blocks disappear quickly).')

            self.iterations = int_spin(500, 1, 10_000_000, 50)
            self.iterations.setToolTip('Number of iterations to run the simulation (if running in silent mode).\nNot used in visual mode.')

            self.background_value = dbl_spin(0.0, -1e9, 1e9, 0.1)
            self.background_value.setToolTip('Default value assigned to blocks that are not visited by any ants.')

            self.background_distance = dbl_spin(32.0, 0.0, 1e9, 1.0)
            self.background_distance.setToolTip('Distance from samples beyond which the background value is applied.\nUnit: Grid blocks.')

            self.value_filter = dbl_spin(0.0, -1e9, 1e9, 1.0)
            self.value_filter.setToolTip('Minimum value threshold for samples to spawn ants.\nSamples with values below this will be ignored.')

            self.avoid_visited_enabled = QtWidgets.QCheckBox(); self.avoid_visited_enabled.setChecked(False)
            self.avoid_visited_enabled.setToolTip('If checked, ants will avoid moving into blocks that have already been visited more than the threshold.')

            self.avoid_visited_threshold = QtWidgets.QSpinBox(); self.avoid_visited_threshold.setRange(1, 10_000_000); self.avoid_visited_threshold.setValue(100)
            self.avoid_visited_threshold.setToolTip('Maximum number of visits allowed per block before ants start avoiding it.\nRequires "Avoid Visited" to be checked.')

            self.average_with_blocks = QtWidgets.QCheckBox(); self.average_with_blocks.setChecked(True)
            self.average_with_blocks.setToolTip('Ant Colony only. When an ant revisits a non-sample block that already has a value, the new neighbor-based estimate is averaged 50/50 with the block\'s current value instead of replacing it outright. This smooths repeated updates and reduces oscillation.')

            ant_form.addRow('Range Size', self.range_size)
            ant_form.addRow('Max Pheromone', self.max_pheromone)
            ant_form.addRow('Ants per Sample', self.ants_per_sample)
            ant_form.addRow('Samples with Ants (%)', self.ants_sampling_percentage)
            ant_form.addRow('Pheromone Decay Rate', self.pheromone_decay_rate)
            ant_form.addRow('Iterations (silent)', self.iterations)
            ant_form.addRow('Background Value', self.background_value)
            ant_form.addRow('Background Distance', self.background_distance)
            ant_form.addRow('Value Filter', self.value_filter)
            ant_form.addRow('Avoid Heavily-Visited', self.avoid_visited_enabled)
            ant_form.addRow('Visited Threshold', self.avoid_visited_threshold)
            ant_form.addRow('Average With Blocks', self.average_with_blocks)

            # Ant Colony domain interpolation toggle
            self.ant_interpolate_target = QtWidgets.QComboBox(); self.ant_interpolate_target.addItems(['Value', 'Domain'])
            self.ant_interpolate_target.setToolTip('Select what Ant Colony should interpolate:\nValue: numeric grade/value (current behavior).\nDomain: categorical domain strings from samples (requires a Domain column in samples).')
            ant_form.addRow('Interpolate', self.ant_interpolate_target)

            # === MOLECULAR CLOCK TAB ===
            self.mc_spatial_weight = dbl_spin(1.0, 0.0, 100.0, 0.1)
            self.mc_spatial_weight.setToolTip('Relative importance of physical distance in calculating evolutionary divergence.\nUnit: Dimensionless multiplier.\n> 1.0: Spatial proximity dominates (samples must be physically close).\n< 1.0: Spatial distance matters less (distant samples can be connected if grades are similar).')
            
            self.mc_attr_weight = dbl_spin(1.0, 0.0, 100.0, 0.1)
            self.mc_attr_weight.setToolTip('Relative importance of grade/value difference in calculating evolutionary divergence.\nUnit: Dimensionless multiplier.\n> 1.0: Grade similarity dominates (samples must have very similar values).\n< 1.0: Grade differences are tolerated (chemically distinct samples can be connected).')
            
            self.mc_ancestor_depth_offset = dbl_spin(1.0, 0.0, 10.0, 0.1)
            self.mc_ancestor_depth_offset.setToolTip('Minimum vertical drop enforced when placing an ancestor node below its descendants.\nUnit: Grid blocks (1.0 = 1 block height).\nLarger: Forces ancestors deeper, creating steeper, more vertical tree structures.\nSmaller: Allows flatter branches (ancestors just slightly below deepest sample).')
            
            self.mc_branch_threshold = dbl_spin(2.0, 0.1, 10.0, 0.1)
            self.mc_branch_threshold.setToolTip('Sensitivity for splitting data into separate geological events (intrusions).\nUnit: Multiplier of "natural" nearest-neighbor distance.\nLarger: Lumps data together (connects distant groups into single large tree).\nSmaller: Splits data apart (identifies many small, isolated events/clusters).')
            
            self.mc_min_samples = int_spin(3, 1, 100, 1)
            self.mc_min_samples.setToolTip('Minimum number of samples required to form a valid cluster/event.\nUnit: Count (samples).\nLarger: Ignores small isolated groups (noise filtering).\nSmaller: Keeps small clusters, potentially including noise.')
            
            self.mc_max_samples = int_spin(1000, 0, 100000, 100)
            self.mc_max_samples.setToolTip('Performance limit. Maximum samples processed in a single event tree.\nUnit: Count (samples).\nLarger: More accurate for huge datasets but slower.\nSmaller: Faster, but randomly subsamples large events (may lose detail).')
            
            self.mc_detect_multiple = QtWidgets.QCheckBox(); self.mc_detect_multiple.setChecked(True)
            self.mc_detect_multiple.setToolTip('Whether to look for distinct clusters or treat everything as one system.\nChecked: Runs clustering first (can find multiple separate feeders).\nUnchecked: Forces all samples into one single tree with one common feeder (LUCA).')
            
            self.mc_interp_method = QtWidgets.QComboBox(); self.mc_interp_method.addItems(['linear', 'inverse_distance'])
            self.mc_interp_method.setToolTip('How block values are calculated along branches connecting ancestors and samples.\nlinear: Smooth gradient between nodes.\ninverse_distance: Weighted average based on distance (standard IDW).')
            
            mc_form.addRow('Spatial Weight', self.mc_spatial_weight)
            mc_form.addRow('Attribute Weight', self.mc_attr_weight)
            mc_form.addRow('Ancestor Depth Offset', self.mc_ancestor_depth_offset)
            mc_form.addRow('Branch Threshold', self.mc_branch_threshold)
            mc_form.addRow('Min Samples per Event', self.mc_min_samples)
            mc_form.addRow('Max Samples per Event', self.mc_max_samples)
            mc_form.addRow('Detect Multiple Events', self.mc_detect_multiple)
            mc_form.addRow('Interpolation Method', self.mc_interp_method)

            # === GAUSSIAN KERNEL TAB ===
            self.gk_bandwidth = dbl_spin(3.0, 0.001, 1_000_000.0, 0.1)
            self.gk_bandwidth.setToolTip('Kernel bandwidth in grid blocks. Larger values smooth farther across the grid; smaller values stay local.')

            self.gk_cutoff_sigma = dbl_spin(3.0, 0.5, 100.0, 0.1)
            self.gk_cutoff_sigma.setToolTip('Search radius expressed as multiples of the bandwidth. Actual radius = bandwidth x cutoff sigma.')

            self.gk_use_nearest_fallback = QtWidgets.QCheckBox()
            self.gk_use_nearest_fallback.setChecked(True)
            self.gk_use_nearest_fallback.setToolTip('If checked, uses the nearest compatible sample when no samples fall inside the kernel search radius.')

            self.gk_fill_background = QtWidgets.QCheckBox()
            self.gk_fill_background.setChecked(False)
            self.gk_fill_background.setToolTip('If checked, cells with no usable kernel estimate are filled with the background value instead of being left unset.')

            self.gk_background_value = dbl_spin(0.0, -1e9, 1e9, 1.0)
            self.gk_background_value.setToolTip('Fallback value assigned only when Fill Background is enabled and no kernel estimate is available.')
            self.gk_background_value.setEnabled(False)

            def update_gk_background_ui(checked):
                self.gk_background_value.setEnabled(checked)

            self.gk_fill_background.toggled.connect(update_gk_background_ui)

            gk_form.addRow('Bandwidth', self.gk_bandwidth)
            gk_form.addRow('Cutoff Sigma', self.gk_cutoff_sigma)
            gk_form.addRow('Use Nearest Fallback', self.gk_use_nearest_fallback)
            gk_form.addRow('Fill Background', self.gk_fill_background)
            gk_form.addRow('Background Value', self.gk_background_value)

            # === ADAPTIVE OCTREE TAB ===
            self.octree_output_mode = QtWidgets.QComboBox()
            self.octree_output_mode.addItems(['Dense Blocks Cover', 'Adaptive Leaf Cover'])
            self.octree_output_mode.setCurrentText('Dense Blocks Cover')
            self.octree_output_mode.setToolTip('Dense Blocks Cover populates every final block using the adaptive field. Adaptive Leaf Cover exports only the non-overlapping adaptive leaves.')

            self.octree_max_levels = QtWidgets.QSpinBox()
            self.octree_max_levels.setRange(0, 64)
            self.octree_max_levels.setValue(0)
            self.octree_max_levels.setToolTip('Maximum number of aggregation levels. Set to 0 to aggregate all the way to the domain root implied by the current block grid.')

            self.octree_support_density_alpha = dbl_spin(0.0, 0.0, 1.0, 0.05)
            self.octree_support_density_alpha.setDecimals(2)
            self.octree_support_density_alpha.setToolTip('Controls how strongly larger represented volumes are penalized when a parent octree node averages child values. Effective support is computed as s / (V^alpha), where s is sample support and V is the represented finest-cell volume. Set 0 for pure sample support or values closer to 1 to favor denser, smaller supported regions.')

            self.octree_include_dense_provenance = QtWidgets.QCheckBox()
            self.octree_include_dense_provenance.setChecked(True)
            self.octree_include_dense_provenance.setToolTip('When Dense Blocks Cover is selected, include adaptive leaf lineage columns on each dense output block for filtering and debugging.')

            def update_octree_ui():
                is_dense_output = self.octree_output_mode.currentText().strip().lower() == 'dense blocks cover'
                self.octree_include_dense_provenance.setEnabled(is_dense_output)

            self.octree_output_mode.currentTextChanged.connect(update_octree_ui)
            update_octree_ui()

            octree_form.addRow('Output Mode', self.octree_output_mode)
            octree_form.addRow('Max Levels (0 = auto)', self.octree_max_levels)
            octree_form.addRow('Support Density Alpha', self.octree_support_density_alpha)
            octree_form.addRow('Dense Provenance Columns', self.octree_include_dense_provenance)

            # === STRING THEORY TAB ===
            self.st_distance_threshold = dbl_spin(10.0, 0.1, 1000.0, 1.0)
            self.st_distance_threshold.setToolTip('Maximum distance to search for a connection (in blocks).')

            self.st_interpolate_target = QtWidgets.QComboBox()
            self.st_interpolate_target.addItems(['Value', 'Domain'])
            self.st_interpolate_target.setToolTip('Select what String Theory should interpolate:\nValue: numeric grade/value (current behavior).\nDomain: categorical domain strings from samples (requires a Domain column in samples).')
            
            self.st_grade_difference = dbl_spin(1.0, 0.0, 1000.0, 0.1)
            self.st_grade_difference.setToolTip('Maximum grade difference allowed for a connection.')

            def update_st_target_ui():
                is_domain = (self.st_interpolate_target.currentText().strip().lower() == 'domain')
                self.st_grade_difference.setEnabled(not is_domain)
                if is_domain:
                    self.st_grade_difference.setToolTip('Not used in Domain mode.')
                else:
                    self.st_grade_difference.setToolTip('Maximum grade difference allowed for a connection.')

            self.st_interpolate_target.currentIndexChanged.connect(update_st_target_ui)
            self._update_st_target_ui = update_st_target_ui
            self._update_st_target_ui()
            
            self.st_connect_to_all = QtWidgets.QCheckBox(); self.st_connect_to_all.setChecked(True)
            self.st_connect_to_all.setToolTip('If checked, connects to ALL valid samples within threshold.\nIf unchecked, connects only to the single best match (closest grade, then closest distance).')

            self.st_max_connections = int_spin(1, 1, 1000, 1)
            self.st_max_connections.setToolTip('Number of valid samples to connect to when "Connect to All Valid" is unchecked.\nMinimum: 1.')
            self.st_max_connections.setEnabled(False)
            
            self.st_min_connections = int_spin(1, 1, 1000, 1)
            self.st_min_connections.setToolTip('Minimum number of connections required. If a sample has fewer valid candidates than this, no connections are made.')
            self.st_min_connections.setEnabled(False)

            def update_connections_enabled(checked):
                self.st_max_connections.setEnabled(not checked)
                self.st_min_connections.setEnabled(not checked)

            self.st_connect_to_all.toggled.connect(update_connections_enabled)

            # Validation: Max >= Min
            def validate_connections():
                min_val = self.st_min_connections.value()
                max_val = self.st_max_connections.value()
                if max_val < min_val:
                    self.st_max_connections.setValue(min_val)

            self.st_min_connections.valueChanged.connect(validate_connections)
            self.st_max_connections.valueChanged.connect(validate_connections)

            conn_layout = QtWidgets.QHBoxLayout()
            conn_layout.addWidget(QtWidgets.QLabel("Min:"))
            conn_layout.addWidget(self.st_min_connections)
            conn_layout.addWidget(QtWidgets.QLabel("Max:"))
            conn_layout.addWidget(self.st_max_connections)

            self.st_collision_policy = QtWidgets.QComboBox()
            self.st_collision_policy.addItems(['overwrite', 'average', 'min', 'max'])
            self.st_collision_policy.setCurrentText('average')
            self.st_collision_policy.setToolTip('How to handle blocks where multiple strings overlap:\noverwrite: Use the latest value (order dependent).\naverage: Calculate running average of all strings.\nmin: Keep the minimum value.\nmax: Keep the maximum value.')

            self.st_processing_order = QtWidgets.QComboBox()
            self.st_processing_order.addItems(['ascending', 'random'])
            self.st_processing_order.setCurrentText('ascending')
            self.st_processing_order.setToolTip('Order in which samples are processed:\nascending: Process lowest grade samples first (builds network bottom-up).\nrandom: Process samples in random order.')

            self.st_filter_by_frequency = QtWidgets.QCheckBox(); self.st_filter_by_frequency.setChecked(False)
            self.st_filter_by_frequency.setToolTip('If checked, filters paths based on Azimuth and Dip frequency.')
            
            self.st_min_azimuth_freq = dbl_spin(10.0, 0.0, 100.0, 1.0)
            self.st_min_azimuth_freq.setToolTip('Minimum frequency (% of max) for Azimuth bin to be kept.')
            self.st_min_azimuth_freq.setEnabled(False)
            
            self.st_min_dip_freq = dbl_spin(10.0, 0.0, 100.0, 1.0)
            self.st_min_dip_freq.setToolTip('Minimum frequency (% of max) for Dip bin to be kept.')
            self.st_min_dip_freq.setEnabled(False)
            
            self.st_filter_by_frequency.toggled.connect(lambda checked: self.st_min_azimuth_freq.setEnabled(checked))
            self.st_filter_by_frequency.toggled.connect(lambda checked: self.st_min_dip_freq.setEnabled(checked))

            st_form.addRow('Interpolate', self.st_interpolate_target)
            st_form.addRow('Distance Threshold', self.st_distance_threshold)
            st_form.addRow('Grade Difference', self.st_grade_difference)
            st_form.addRow('Connect to All Valid', self.st_connect_to_all)
            st_form.addRow('Connections (Min/Max)', conn_layout)
            st_form.addRow('Collision Policy', self.st_collision_policy)
            st_form.addRow('Processing Order', self.st_processing_order)
            st_form.addRow('Filter by Frequency', self.st_filter_by_frequency)
            st_form.addRow('Min Azimuth Freq (%)', self.st_min_azimuth_freq)
            st_form.addRow('Min Dip Freq (%)', self.st_min_dip_freq)

            # === ADVANCED TAB ===
            self.verbose = QtWidgets.QCheckBox(); self.verbose.setChecked(False)
            self.process_domains_sequentially = QtWidgets.QCheckBox(); self.process_domains_sequentially.setChecked(True)
            self.expand_interpolation_exports_to_subblocks = QtWidgets.QCheckBox(); self.expand_interpolation_exports_to_subblocks.setChecked(True)
            self.expand_interpolation_exports_to_subblocks.setToolTip('If checked, interpolation CSV exports reuse the original blocks file rows and assign each sub-block the interpolated value of its parent base block. Interpolation still runs on the base-block grid.')
            self.force_rebuild_sample_blocks = QtWidgets.QCheckBox(); self.force_rebuild_sample_blocks.setChecked(False)
            self.force_rebuild_sample_blocks.setToolTip('If checked, ignore any valid sample-block cache and regenerate sample blocks from the source samples for this run.')
            self.blank_sample_domain_behavior = QtWidgets.QComboBox(); self.blank_sample_domain_behavior.addItems(['Skip', 'Infer From Blocks'])
            self.blank_sample_domain_behavior.setCurrentText('Skip')
            self.blank_sample_domain_behavior.setToolTip('How to handle blank sample domains during domain-based interpolation or metrics.\nSkip: exclude blank-domain rows.\nInfer From Blocks: infer missing sample domains from the blocks model when possible.')

            self.algorithm_combo = QtWidgets.QComboBox()
            self.algorithm_combo.addItems(['ant_colony', 'molecular_clock', 'gaussian_kernel', 'adaptive_octree', 'string_theory'])

            self.second_pass_combo = QtWidgets.QComboBox()
            self.second_pass_combo.addItems(['skip', 'ant_colony', 'molecular_clock', 'gaussian_kernel', 'adaptive_octree', 'string_theory'])
            self.second_pass_combo.setCurrentText('skip')

            def validate_algorithms():
                algo1 = self.algorithm_combo.currentText()
                algo2 = self.second_pass_combo.currentText()
                if algo2 != 'skip' and algo1 == algo2:
                    QtWidgets.QMessageBox.warning(self, "Invalid Selection",
                        f"Second pass algorithm cannot be the same as the first pass ({algo1}).")
                    self.second_pass_combo.setCurrentText('skip')

            self.algorithm_combo.currentTextChanged.connect(validate_algorithms)
            self.second_pass_combo.currentTextChanged.connect(validate_algorithms)

            default_algo_layout = QtWidgets.QHBoxLayout()
            default_algo_layout.addWidget(self.algorithm_combo)
            default_algo_layout.addWidget(QtWidgets.QLabel("Second Pass:"))
            default_algo_layout.addWidget(self.second_pass_combo)

            # === DISPLAY TAB ===
            self.taichi_sample_diameter = dbl_spin(1.0, 0.001, 1e6, 0.1)
            self.taichi_sample_diameter.setToolTip('Sample diameter in model units for the Taichi mesh viewer. Default is 1 unit.')
            display_form.addRow('Sample Diameter', self.taichi_sample_diameter)

            advanced_form.addRow('Default Pass Algorithms', default_algo_layout)

            # Domain Algorithm Mapping
            self.domain_overrides = {}  # Store domain -> algorithm mapping
            self.domain_mapping_btn = QtWidgets.QPushButton('Configure Domain Algorithms...')
            self.domain_mapping_btn.clicked.connect(self.open_domain_mapping)
            advanced_form.addRow('Domain-Specific Algorithms', self.domain_mapping_btn)

            advanced_form.addRow('Process Domains Sequentially', self.process_domains_sequentially)
            advanced_form.addRow('Expand CSV Export to Sub-Blocks', self.expand_interpolation_exports_to_subblocks)
            advanced_form.addRow('Force Rebuild Sample Blocks', self.force_rebuild_sample_blocks)
            advanced_form.addRow('Blank Sample Domains', self.blank_sample_domain_behavior)
            advanced_form.addRow('Verbose', self.verbose)

            # Buttons at bottom of main layout
            btn_box = QtWidgets.QHBoxLayout()
            self.load_btn = QtWidgets.QPushButton('Load Config')
            self.save_btn = QtWidgets.QPushButton('Save Config')
            self.run_only_btn = QtWidgets.QPushButton('Run Interpolation Only')
            self.start_viewer_btn = QtWidgets.QPushButton('Open Viewer')
            self.refresh_viewer_btn = QtWidgets.QPushButton('Refresh Viewer')
            self.reload_data_btn = QtWidgets.QPushButton('Reload Data')
            self.cancel_btn = QtWidgets.QPushButton('Cancel')
            self.viewer_status_label = QtWidgets.QLabel('Viewer: stopped')
            btn_box.addWidget(self.load_btn); btn_box.addWidget(self.save_btn)
            btn_box.addStretch(1)
            btn_box.addWidget(self.viewer_status_label)
            btn_box.addWidget(self.run_only_btn); btn_box.addWidget(self.start_viewer_btn); btn_box.addWidget(self.refresh_viewer_btn); btn_box.addWidget(self.reload_data_btn); btn_box.addWidget(self.cancel_btn)
            main_layout.addLayout(btn_box)

            self.load_btn.clicked.connect(self.load_config)
            self.save_btn.clicked.connect(self.save_config)
            self.run_only_btn.clicked.connect(self.run_interpolation_only)
            self.start_viewer_btn.clicked.connect(self.start_viewer)
            self.refresh_viewer_btn.clicked.connect(self.refresh_viewer)
            self.reload_data_btn.clicked.connect(self.reload_viewer_data)
            self.cancel_btn.clicked.connect(self.reject)

            self._update_viewer_status()

        def start_viewer(self):
            self.should_visualize = True
            self.viewer_backend = 'taichi'
            self.launch_viewer_process()

        def _is_viewer_running(self):
            return self._viewer_process is not None and self._viewer_process.poll() is None

        def _update_viewer_status(self):
            if self._is_viewer_running():
                self.viewer_status_label.setText('Viewer: running')
                self.start_viewer_btn.setEnabled(False)
                self.refresh_viewer_btn.setEnabled(True)
                self.reload_data_btn.setEnabled(True)
            else:
                self.viewer_status_label.setText('Viewer: stopped')
                self.start_viewer_btn.setEnabled(True)
                self.refresh_viewer_btn.setEnabled(False)
                self.reload_data_btn.setEnabled(False)

        def _cleanup_finished_viewer(self):
            if self._is_viewer_running():
                self._update_viewer_status()
                return
            if self._viewer_process is not None:
                self._viewer_process = None
            self._viewer_render_mode = None
            if self._viewer_config_path and os.path.isfile(self._viewer_config_path):
                try:
                    os.remove(self._viewer_config_path)
                except OSError:
                    pass
            self._viewer_config_path = None
            self._update_viewer_status()

        def _terminate_viewer_process(self):
            if not self._is_viewer_running():
                self._cleanup_finished_viewer()
                return
            try:
                self._viewer_process.terminate()
                self._viewer_process.wait(timeout=5)
            except Exception:
                try:
                    self._viewer_process.kill()
                except Exception:
                    pass
            self._cleanup_finished_viewer()

        def _prepare_for_shutdown(self):
            if getattr(self, '_viewer_process_timer', None) is not None and self._viewer_process_timer.isActive():
                self._viewer_process_timer.stop()
            self._terminate_viewer_process()

        def refresh_viewer(self):
            if not self._is_viewer_running():
                self._cleanup_finished_viewer()
                return
            cfg = self.to_dict(include_runtime_state=True)
            cfg['viewer_backend'] = _normalize_viewer_backend(cfg.get('viewer_backend'))
            if self._viewer_mode_requires_restart(cfg):
                self.launch_viewer_process(config_override=cfg, force_restart=True)
                return
            self._write_current_viewer_config(reload_mode='refresh', config_override=cfg)

        def reload_viewer_data(self):
            if not self._is_viewer_running():
                self._cleanup_finished_viewer()
                return
            cfg = self.to_dict(include_runtime_state=True)
            cfg['viewer_backend'] = _normalize_viewer_backend(cfg.get('viewer_backend'))
            if self._viewer_mode_requires_restart(cfg):
                self.launch_viewer_process(config_override=cfg, force_restart=True)
                return
            self._write_current_viewer_config(reload_mode='reload', config_override=cfg)

        def _viewer_mode_requires_restart(self, config):
            return False

        def _write_current_viewer_config(self, reload_mode='refresh', config_override=None):
            if not self._viewer_config_path:
                return
            cfg = dict(config_override or self.to_dict(include_runtime_state=True))
            cfg['viewer_backend'] = _normalize_viewer_backend(cfg.get('viewer_backend'))
            cfg['_viewer_reload_mode'] = reload_mode
            self._viewer_render_mode = _normalize_taichi_block_render_mode(cfg.get('taichi_block_render_mode'))
            _write_json_atomic(self._viewer_config_path, cfg)

        def launch_viewer_process(self, config_override=None, force_restart=False):
            try:
                cfg = dict(config_override or self.to_dict(include_runtime_state=True))
                cfg['viewer_backend'] = _normalize_viewer_backend(cfg.get('viewer_backend'))
                cfg['_viewer_reload_mode'] = 'reload'
                cfg['taichi_block_render_mode'] = _normalize_taichi_block_render_mode(cfg.get('taichi_block_render_mode'))
                if force_restart and self._is_viewer_running():
                    self._terminate_viewer_process()
                if self._is_viewer_running():
                    self._update_viewer_status()
                    return
                handle = tempfile.NamedTemporaryFile(
                    mode='w',
                    suffix='.json',
                    prefix='anterpolator_viewer_',
                    delete=False,
                    encoding='utf-8',
                )
                with handle:
                    json.dump(cfg, handle, indent=4)
                    config_path = handle.name

                script_path = os.path.abspath(__file__)
                popen_kwargs = {
                    'cwd': os.path.dirname(script_path),
                }
                if os.name == 'nt':
                    popen_kwargs['creationflags'] = getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP', 0)
                else:
                    popen_kwargs['start_new_session'] = True
                process = subprocess.Popen(
                    [sys.executable, script_path, '--launch-viewer-config', config_path],
                    **popen_kwargs,
                )
                self._viewer_process = process
                self._viewer_config_path = config_path
                self._viewer_render_mode = cfg['taichi_block_render_mode']
                self._update_viewer_status()
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Failed to launch viewer: {e}')

        def closeEvent(self, event):
            self._prepare_for_shutdown()
            super().closeEvent(event)

        def _invalidate_domain_catalog_cache(self):
            self._domain_catalog_cache = None

        def to_dict(self, include_runtime_state=False):
            config = {
                'samples_file': self.samples_edit.text(),
                'blocks_file': self.blocks_edit.text(),
                'color_file': self.color_edit.text(),
                'interpolation_file': resolve_interpolation_csv_export_path(self.interp_edit.text()),
                'export_block_evaluated_samples': self.block_evaluated_samples_enabled.isChecked(),
                'block_evaluated_samples_file': self.block_evaluated_samples_edit.text(),
                'sample_blocks_file': self.sample_blocks_output_edit.text(),
                'sample_blocks_include_ids': self.sample_blocks_include_ids.isChecked(),
                'sample_blocks_id_cols': [cb.currentText() for cb in self.sample_blocks_id_cols],
                'domain_samples_file': self.domain_samples_output_edit.text(),
                'block_value_transfer_file': self.block_value_transfer_output_edit.text(),
                'block_value_transfer_cols': self._get_selected_block_value_transfer_columns(),
                'block_model_target_file': self.block_model_target_edit.text(),
                'block_model_target_delimiter': self.block_model_target_delim.currentText(),
                'block_model_target_header_line': self.block_model_target_header_line.value(),
                'block_model_target_x_col': self.block_model_target_x_col.currentText(),
                'block_model_target_y_col': self.block_model_target_y_col.currentText(),
                'block_model_target_z_col': self.block_model_target_z_col.currentText(),
                'block_model_source_size_cols': [combo.currentText() for combo in self.block_model_source_size_cols],
                'block_model_target_size_cols': [combo.currentText() for combo in self.block_model_target_size_cols],
                'block_model_target_size': tuple(spin.value() for spin in self.block_model_target_size_spins),
                'block_model_transfer_cols': self._get_selected_block_model_transfer_columns(),
                'block_model_transfer_nearest_fallback': self.block_model_nearest_fallback.isChecked(),
                'block_model_transfer_max_nearest_distance': self.block_model_nearest_max_distance.value(),
                'block_model_transfer_file': self.block_model_transfer_output_edit.text(),
                'table_attribute_block_model_file': self.table_attribute_block_model_edit.text(),
                'table_attribute_block_model_delimiter': self.table_attribute_block_model_delim.currentText(),
                'table_attribute_block_model_header_line': self.table_attribute_block_model_header_line.value(),
                'table_attribute_table_file': self.table_attribute_table_edit.text(),
                'table_attribute_table_delimiter': self.table_attribute_table_delim.currentText(),
                'table_attribute_table_header_line': self.table_attribute_table_header_line.value(),
                'table_attribute_key_cols': self._get_selected_table_attribute_key_columns(),
                'table_attribute_value_cols': self._get_selected_table_attribute_value_columns(),
                'table_attribute_output_file': self.table_attribute_output_edit.text(),
                'block_domain_metrics_file': self.block_domain_metrics_output_edit.text(),
                'block_domain_metrics_selected_metrics': self._get_selected_block_domain_metrics(),
                'block_domain_metrics_avg_distance_knn_k': self.block_domain_metrics_knn_k.value(),
                'domain_interpolation_confidence_file': self.domain_interpolation_confidence_output_edit.text(),
                'block_volume_weighted_file': self.block_volume_weighted_output_edit.text(),
                'bmf_export_input_file': self.bmf_export_input_edit.text(),
                'bmf_export_output_file': self.bmf_export_output_edit.text(),
                'bmf_export_delimiter': self.bmf_export_delim.currentText(),
                'bmf_export_header_line': self.bmf_export_header_line.value(),
                'bmf_export_backend': self.bmf_export_backend_combo.currentText(),
                'bmf_export_cell_size': self._get_bmf_export_cell_size(),
                'bmf_export_regularize_to_base_block': self.bmf_export_regularize_base_blocks.isChecked(),
                'bmf_export_x_col': self.bmf_export_x_col.currentText(),
                'bmf_export_y_col': self.bmf_export_y_col.currentText(),
                'bmf_export_z_col': self.bmf_export_z_col.currentText(),
                'bmf_export_value_cols': self._get_selected_bmf_export_value_columns(),
                'bmf_export_column_types': self._get_bmf_export_column_type_overrides(include_unselected=True),
                'bmf_export_value_exceptions': self._get_bmf_export_value_exceptions(),
                'equation_finder_file': self.equation_finder_output_edit.text(),
                'block_domain_metrics_closest_sample_id_cols': [cb.currentText() for cb in self.block_domain_metrics_id_cols],
                'block_domain_metrics_use_block_value_prefix': self.block_domain_metrics_prefix_with_block_value.isChecked(),
                'block_domain_metrics_distance_counts_enabled': 'distance_band_summary' in self._get_selected_block_domain_metrics(),
                'block_domain_metrics_distance_step': self.block_domain_metrics_distance_step.value(),
                'block_domain_metrics_max_factor': self.block_domain_metrics_max_factor.value(),
                'block_domain_sample_filters': [dict(entry) for entry in self.block_domain_sample_filters],
                'block_volume_weighted_filters': [dict(entry) for entry in self.block_volume_weighted_filters],
                'equation_finder_predictor_cols': self._get_selected_equation_predictor_columns(),
                'equation_finder_include_coordinates': self.equation_include_coords.isChecked(),
                'equation_finder_min_samples_per_domain': self.equation_min_samples_per_domain.value(),
                'equation_finder_validation_fraction': self.equation_validation_fraction.value(),
                'equation_finder_max_iterations': self.equation_max_iterations.value(),
                'equation_finder_timeout_seconds': self.equation_timeout_seconds.value(),
                'algorithm': self.algorithm_combo.currentText(),
                'second_pass_algorithm': self.second_pass_combo.currentText(),
                'samples_delimiter': self.samples_delim.currentText(),
                'blocks_delimiter': self.blocks_delim.currentText(),
                'samples_header_line': self.samples_header_line.value(),
                'blocks_header_line': self.blocks_header_line.value(),
                'sample_x_col': self.sample_x_col.currentText(),
                'sample_y_col': self.sample_y_col.currentText(),
                'sample_z_col': self.sample_z_col.currentText(),
                'sample_value_col': self.sample_value_col.currentText(),
                'sample_domain_col': self.sample_domain_col.currentText(),
                'sample_weight_col': self.sample_weight_col.currentText(),
                'block_x_col': self.block_x_col.currentText(),
                'block_y_col': self.block_y_col.currentText(),
                'block_z_col': self.block_z_col.currentText(),
                'block_domain_col': self.block_domain_col.currentText(),
                'block_domain_metrics_value_col': self.block_domain_metrics_value_col.currentText(),
                'block_volume_weighted_value_col': self.block_volume_weighted_value_col.currentText(),
                'block_value_metric_col': self.block_volume_weighted_value_col.currentText(),
                'block_weight_metric_col': self.block_weight_metric_col.currentText(),
                'block_size': (self.block_x.value(), self.block_y.value(), self.block_z.value()),
                'range_size': self.range_size.value(),
                'max_pheromone': self.max_pheromone.value(),
                'ants_per_sample': self.ants_per_sample.value(),
                'ants_sampling_percentage': self.ants_sampling_percentage.value(),
                'pheromone_decay_rate': self.pheromone_decay_rate.value(),
                'iterations': self.iterations.value(),
                'background_value': self.background_value.value(),
                'background_distance': self.background_distance.value(),
                'value_filter': self.value_filter.value(),
                'avoid_visited_threshold_enabled': self.avoid_visited_enabled.isChecked(),
                'avoid_visited_threshold': self.avoid_visited_threshold.value(),
                'ant_colony_interpolate_target': self.ant_interpolate_target.currentText().strip().lower(),
                'molecular_clock_params': {
                    'spatial_weight': self.mc_spatial_weight.value(),
                    'attr_weight': self.mc_attr_weight.value(),
                    'ancestor_depth_offset': self.mc_ancestor_depth_offset.value(),
                    'branch_threshold': self.mc_branch_threshold.value(),
                    'min_samples': self.mc_min_samples.value(),
                    'max_samples': self.mc_max_samples.value(),
                    'detect_multiple': self.mc_detect_multiple.isChecked(),
                    'interp_method': self.mc_interp_method.currentText()
                },
                'gaussian_kernel_params': {
                    'bandwidth': self.gk_bandwidth.value(),
                    'cutoff_sigma': self.gk_cutoff_sigma.value(),
                    'use_nearest_fallback': self.gk_use_nearest_fallback.isChecked(),
                    'fill_background': self.gk_fill_background.isChecked(),
                    'background_value': self.gk_background_value.value(),
                },
                'adaptive_octree_params': {
                    'output_mode': 'dense_blocks_cover' if self.octree_output_mode.currentText().strip().lower() == 'dense blocks cover' else 'adaptive_leaf_cover',
                    'max_levels': self.octree_max_levels.value(),
                    'support_density_alpha': self.octree_support_density_alpha.value(),
                    'include_dense_provenance': self.octree_include_dense_provenance.isChecked(),
                },
                'string_theory_params': {
                    'interpolate_target': self.st_interpolate_target.currentText().strip().lower(),
                    'distance_threshold': self.st_distance_threshold.value(),
                    'grade_difference': self.st_grade_difference.value(),
                    'connect_to_all': self.st_connect_to_all.isChecked(),
                    'max_connections': self.st_max_connections.value(),
                    'min_connections': self.st_min_connections.value(),
                    'collision_policy': self.st_collision_policy.currentText(),
                    'processing_order': self.st_processing_order.currentText(),
                    'filter_by_frequency': self.st_filter_by_frequency.isChecked(),
                    'min_azimuth_freq': self.st_min_azimuth_freq.value(),
                    'min_dip_freq': self.st_min_dip_freq.value()
                },
                'average_with_blocks': self.average_with_blocks.isChecked(),
                'process_domains_sequentially': self.process_domains_sequentially.isChecked(),
                'expand_interpolation_exports_to_subblocks': self.expand_interpolation_exports_to_subblocks.isChecked(),
                'force_rebuild_sample_blocks': self.force_rebuild_sample_blocks.isChecked(),
                'blank_sample_domain_behavior': 'infer_from_blocks' if self.blank_sample_domain_behavior.currentText() == 'Infer From Blocks' else 'skip',
                'verbose': self.verbose.isChecked(),
                'viewer_backend': self.viewer_backend,
                'taichi_block_render_mode': 'mesh',
                'taichi_transparent_blocks': False,
                'taichi_sample_diameter': self.taichi_sample_diameter.value(),
                'domain_algorithm_overrides': self.domain_overrides
            }
            if include_runtime_state and self._domain_catalog_cache:
                config['_domain_catalog_cache'] = {
                    'signature': dict(self._domain_catalog_cache.get('signature', {})),
                    'domains': list(self._domain_catalog_cache.get('domains', [])),
                }
            if include_runtime_state:
                config['_prefer_interpolation_file_for_viewer'] = bool(self._prefer_interpolation_file_for_viewer)
            return config

        def from_dict(self, config):
            was_running = self._is_viewer_running()
            self._suspend_auto_viewer_refresh = True
            try:
                self.block_domain_sample_filters = []
                self.block_volume_weighted_filters = []
                self._update_block_domain_metrics_filters_summary()
                self._update_block_volume_weighted_filters_summary()
                if 'samples_file' in config: self.samples_edit.setText(config['samples_file'])
                if 'blocks_file' in config: self.blocks_edit.setText(config['blocks_file'])
                if 'color_file' in config: self.color_edit.setText(config['color_file'])
                if 'interpolation_file' in config: self.interp_edit.setText(resolve_interpolation_csv_export_path(config['interpolation_file']))
                if 'export_block_evaluated_samples' in config: self.block_evaluated_samples_enabled.setChecked(bool(config['export_block_evaluated_samples']))
                if 'block_evaluated_samples_file' in config: self.block_evaluated_samples_edit.setText(config['block_evaluated_samples_file'])
                if 'sample_blocks_file' in config: self.sample_blocks_output_edit.setText(config['sample_blocks_file'])
                if 'sample_blocks_include_ids' in config: self.sample_blocks_include_ids.setChecked(bool(config['sample_blocks_include_ids']))
                self.block_domain_metrics_output_edit.setText('')
                if 'algorithm' in config: self.algorithm_combo.setCurrentText(config['algorithm'])
                if 'second_pass_algorithm' in config: self.second_pass_combo.setCurrentText(config['second_pass_algorithm'])
                if 'samples_delimiter' in config: self.samples_delim.setCurrentText(config['samples_delimiter'])
                if 'blocks_delimiter' in config: self.blocks_delim.setCurrentText(config['blocks_delimiter'])
                if 'samples_header_line' in config: self.samples_header_line.setValue(config['samples_header_line'])
                if 'blocks_header_line' in config: self.blocks_header_line.setValue(config['blocks_header_line'])
            # Columns might need refresh first, but we can try setting text
                if 'sample_x_col' in config: self.sample_x_col.setCurrentText(config['sample_x_col'])
                if 'sample_y_col' in config: self.sample_y_col.setCurrentText(config['sample_y_col'])
                if 'sample_z_col' in config: self.sample_z_col.setCurrentText(config['sample_z_col'])
                if 'sample_value_col' in config: self.sample_value_col.setCurrentText(config['sample_value_col'])
                if 'sample_domain_col' in config: self.sample_domain_col.setCurrentText(config['sample_domain_col'])
                if 'sample_weight_col' in config: self.sample_weight_col.setCurrentText(config['sample_weight_col'])
                if 'block_x_col' in config: self.block_x_col.setCurrentText(config['block_x_col'])
                if 'block_y_col' in config: self.block_y_col.setCurrentText(config['block_y_col'])
                if 'block_z_col' in config: self.block_z_col.setCurrentText(config['block_z_col'])
                if 'block_domain_col' in config: self.block_domain_col.setCurrentText(config['block_domain_col'])
                if 'domain_samples_file' in config: self.domain_samples_output_edit.setText(config['domain_samples_file'])
                if 'block_value_transfer_file' in config: self.block_value_transfer_output_edit.setText(config['block_value_transfer_file'])
                if 'block_model_target_file' in config: self.block_model_target_edit.setText(config['block_model_target_file'])
                if 'block_model_target_delimiter' in config: self.block_model_target_delim.setCurrentText(config['block_model_target_delimiter'])
                if 'block_model_target_header_line' in config: self.block_model_target_header_line.setValue(int(config['block_model_target_header_line']))
                if 'block_model_transfer_file' in config: self.block_model_transfer_output_edit.setText(config['block_model_transfer_file'])
                if 'table_attribute_block_model_file' in config: self.table_attribute_block_model_edit.setText(config['table_attribute_block_model_file'])
                if 'table_attribute_block_model_delimiter' in config: self.table_attribute_block_model_delim.setCurrentText(config['table_attribute_block_model_delimiter'])
                if 'table_attribute_block_model_header_line' in config: self.table_attribute_block_model_header_line.setValue(int(config['table_attribute_block_model_header_line']))
                if 'table_attribute_table_file' in config: self.table_attribute_table_edit.setText(config['table_attribute_table_file'])
                if 'table_attribute_table_delimiter' in config: self.table_attribute_table_delim.setCurrentText(config['table_attribute_table_delimiter'])
                if 'table_attribute_table_header_line' in config: self.table_attribute_table_header_line.setValue(int(config['table_attribute_table_header_line']))
                if 'table_attribute_output_file' in config: self.table_attribute_output_edit.setText(config['table_attribute_output_file'])
                if 'block_domain_metrics_file' in config: self.block_domain_metrics_output_edit.setText(config['block_domain_metrics_file'])
                if 'domain_interpolation_confidence_file' in config: self.domain_interpolation_confidence_output_edit.setText(config['domain_interpolation_confidence_file'])
                if 'block_volume_weighted_file' in config: self.block_volume_weighted_output_edit.setText(config['block_volume_weighted_file'])
                if 'bmf_export_input_file' in config: self.bmf_export_input_edit.setText(config['bmf_export_input_file'])
                if 'bmf_export_output_file' in config: self.bmf_export_output_edit.setText(config['bmf_export_output_file'])
                if 'bmf_export_delimiter' in config: self.bmf_export_delim.setCurrentText(config['bmf_export_delimiter'])
                if 'bmf_export_header_line' in config: self.bmf_export_header_line.setValue(int(config['bmf_export_header_line']))
                if 'bmf_export_backend' in config: self.bmf_export_backend_combo.setCurrentText(config['bmf_export_backend'])
                if 'bmf_export_cell_size' in config: self._set_bmf_export_cell_size(config.get('bmf_export_cell_size'))
                if 'bmf_export_regularize_to_base_block' in config:
                    self.bmf_export_regularize_base_blocks.setChecked(bool(config.get('bmf_export_regularize_to_base_block')))
                self._update_bmf_export_regularize_warning()
                if 'bmf_export_value_exceptions' in config:
                    self._set_bmf_export_value_exceptions(config.get('bmf_export_value_exceptions'))
                else:
                    self._set_bmf_export_value_exceptions({})
                pending_bmf_x_col = str(config.get('bmf_export_x_col') or '').strip()
                pending_bmf_y_col = str(config.get('bmf_export_y_col') or '').strip()
                pending_bmf_z_col = str(config.get('bmf_export_z_col') or '').strip()
                pending_bmf_value_cols = None
                pending_bmf_column_types = {}
                if 'bmf_export_value_cols' in config:
                    raw_bmf_value_cols = config['bmf_export_value_cols']
                    if isinstance(raw_bmf_value_cols, str):
                        pending_bmf_value_cols = [
                            token.strip()
                            for token in raw_bmf_value_cols.split(',')
                            if token.strip()
                        ]
                    else:
                        pending_bmf_value_cols = [
                            str(token).strip()
                            for token in (raw_bmf_value_cols or [])
                            if str(token).strip()
                        ]
                raw_bmf_column_types = config.get('bmf_export_column_types')
                if isinstance(raw_bmf_column_types, dict):
                    pending_bmf_column_types = {
                        str(column_name).strip(): normalize_bmf_export_field_type_name(field_type)
                        for column_name, field_type in raw_bmf_column_types.items()
                        if str(column_name).strip()
                    }
                elif raw_bmf_column_types not in (None, ''):
                    pending_bmf_column_types = parse_bmf_export_column_types(raw_bmf_column_types)
                if hasattr(self, '_refresh_bmf_export_columns'):
                    self._pending_bmf_export_value_cols = pending_bmf_value_cols
                    self._pending_bmf_export_column_types = pending_bmf_column_types
                    self._bmf_export_value_cols_initialized = False
                    self._refresh_bmf_export_columns()
                    if pending_bmf_x_col:
                        self.bmf_export_x_col.setCurrentText(pending_bmf_x_col)
                    if pending_bmf_y_col:
                        self.bmf_export_y_col.setCurrentText(pending_bmf_y_col)
                    if pending_bmf_z_col:
                        self.bmf_export_z_col.setCurrentText(pending_bmf_z_col)
                    self._pending_bmf_export_value_cols = pending_bmf_value_cols
                    self._pending_bmf_export_column_types = pending_bmf_column_types
                    self._refresh_bmf_export_columns()
                if 'equation_finder_file' in config: self.equation_finder_output_edit.setText(config['equation_finder_file'])
                self._pending_block_value_transfer_cols = list(config.get('block_value_transfer_cols', []))
                if hasattr(self, '_refresh_block_value_transfer_columns'):
                    self._refresh_block_value_transfer_columns()
                self._pending_block_model_transfer_cols = list(config.get('block_model_transfer_cols', []))
                if hasattr(self, '_refresh_block_model_transfer_source_columns'):
                    self._refresh_block_model_transfer_source_columns()
                if hasattr(self, '_refresh_block_model_target_columns'):
                    self._refresh_block_model_target_columns()
                self._pending_table_attribute_key_cols = list(config.get('table_attribute_key_cols', []))
                self._pending_table_attribute_value_cols = list(config.get('table_attribute_value_cols', []))
                if hasattr(self, '_refresh_table_attribute_shared_key_columns'):
                    self._refresh_table_attribute_shared_key_columns()
                for combo, column_name in zip(
                    [self.block_model_target_x_col, self.block_model_target_y_col, self.block_model_target_z_col],
                    [config.get('block_model_target_x_col'), config.get('block_model_target_y_col'), config.get('block_model_target_z_col')],
                ):
                    if column_name: combo.setCurrentText(str(column_name))
                for combo, column_name in zip(self.block_model_source_size_cols, config.get('block_model_source_size_cols', [])):
                    combo.setCurrentText(str(column_name))
                for combo, column_name in zip(self.block_model_target_size_cols, config.get('block_model_target_size_cols', [])):
                    combo.setCurrentText(str(column_name))
                for spin, value in zip(self.block_model_target_size_spins, config.get('block_model_target_size', [])):
                    spin.setValue(float(value))
                if 'block_model_transfer_nearest_fallback' in config:
                    self.block_model_nearest_fallback.setChecked(bool(config['block_model_transfer_nearest_fallback']))
                if 'block_model_transfer_max_nearest_distance' in config:
                    self.block_model_nearest_max_distance.setValue(float(config['block_model_transfer_max_nearest_distance']))
                if 'sample_blocks_id_cols' in config:
                    for cb, column_name in zip(self.sample_blocks_id_cols, config['sample_blocks_id_cols']):
                        cb.setCurrentText(column_name)
                if 'block_domain_metrics_closest_sample_id_cols' in config:
                    for cb, column_name in zip(self.block_domain_metrics_id_cols, config['block_domain_metrics_closest_sample_id_cols']):
                        cb.setCurrentText(column_name)
                if 'block_domain_metrics_avg_distance_knn_k' in config:
                    self.block_domain_metrics_knn_k.setValue(int(config['block_domain_metrics_avg_distance_knn_k']))
                if 'block_domain_metrics_use_block_value_prefix' in config:
                    self.block_domain_metrics_prefix_with_block_value.setChecked(bool(config['block_domain_metrics_use_block_value_prefix']))
                if 'block_domain_metrics_distance_step' in config:
                    self.block_domain_metrics_distance_step.setValue(float(config['block_domain_metrics_distance_step']))
                if 'block_domain_metrics_max_factor' in config:
                    self.block_domain_metrics_max_factor.setValue(int(config['block_domain_metrics_max_factor']))
                legacy_block_value_metric_col = config.get('block_value_metric_col')
                if 'block_domain_metrics_value_col' in config:
                    self.block_domain_metrics_value_col.setCurrentText(config['block_domain_metrics_value_col'])
                elif legacy_block_value_metric_col is not None:
                    self.block_domain_metrics_value_col.setCurrentText(legacy_block_value_metric_col)
                if 'block_domain_metrics_selected_metrics' in config:
                    self._set_selected_block_domain_metrics(config.get('block_domain_metrics_selected_metrics'))
                else:
                    legacy_metric_selection = list(BLOCK_DOMAIN_METRIC_LEGACY_SELECTION)
                    if any(str(name or '').strip() and str(name).strip() != '(None)' for name in config.get('block_domain_metrics_closest_sample_id_cols', []) or []):
                        legacy_metric_selection.append('closest_sample_id')
                    effective_block_value_col = config.get('block_domain_metrics_value_col', legacy_block_value_metric_col)
                    sample_value_name = str(config.get('sample_value_col') or '').strip()
                    block_value_name = str(effective_block_value_col or '').strip()
                    if sample_value_name and sample_value_name != '(None)' and block_value_name and block_value_name != '(None)':
                        for metric_id in BLOCK_DOMAIN_METRIC_DEFINITION_ORDER:
                            if metric_id in BLOCK_DOMAIN_SAMPLE_VALUE_METRICS and metric_id not in legacy_metric_selection:
                                legacy_metric_selection.append(metric_id)
                    if bool(config.get('block_domain_metrics_distance_counts_enabled')):
                        legacy_metric_selection.append('distance_band_summary')
                    self._set_selected_block_domain_metrics(legacy_metric_selection)
                if 'block_volume_weighted_value_col' in config:
                    self.block_volume_weighted_value_col.setCurrentText(config['block_volume_weighted_value_col'])
                elif legacy_block_value_metric_col is not None:
                    self.block_volume_weighted_value_col.setCurrentText(legacy_block_value_metric_col)
                if 'block_weight_metric_col' in config: self.block_weight_metric_col.setCurrentText(config['block_weight_metric_col'])
                if 'equation_finder_include_coordinates' in config:
                    self.equation_include_coords.setChecked(bool(config['equation_finder_include_coordinates']))
                if 'equation_finder_min_samples_per_domain' in config:
                    self.equation_min_samples_per_domain.setValue(int(config['equation_finder_min_samples_per_domain']))
                if 'equation_finder_validation_fraction' in config:
                    self.equation_validation_fraction.setValue(float(config['equation_finder_validation_fraction']))
                if 'equation_finder_max_iterations' in config:
                    self.equation_max_iterations.setValue(int(config['equation_finder_max_iterations']))
                if 'equation_finder_timeout_seconds' in config:
                    self.equation_timeout_seconds.setValue(int(config['equation_finder_timeout_seconds']))
                self._pending_equation_predictor_selection = list(config.get('equation_finder_predictor_cols', []))
                if hasattr(self, '_refresh_equation_finder_predictor_columns'):
                    self._refresh_equation_finder_predictor_columns()
                if 'block_domain_sample_filters' in config:
                    self.block_domain_sample_filters = [dict(entry) for entry in config['block_domain_sample_filters']]
                    self._update_block_domain_metrics_filters_summary()
                if 'block_volume_weighted_filters' in config:
                    self.block_volume_weighted_filters = [dict(entry) for entry in config['block_volume_weighted_filters']]
                    self._update_block_volume_weighted_filters_summary()
                if 'block_size' in config:
                    bs = config['block_size']
                    if isinstance(bs, (list, tuple)) and len(bs) == 3:
                        self.block_x.setValue(bs[0])
                        self.block_y.setValue(bs[1])
                        self.block_z.setValue(bs[2])
                    elif isinstance(bs, (int, float)):
                        self.block_x.setValue(int(bs))
                        self.block_y.setValue(int(bs))
                        self.block_z.setValue(int(bs))
                if 'range_size' in config: self.range_size.setValue(config['range_size'])
                if 'max_pheromone' in config: self.max_pheromone.setValue(config['max_pheromone'])
                if 'ants_per_sample' in config: self.ants_per_sample.setValue(config['ants_per_sample'])
                if 'ants_sampling_percentage' in config: self.ants_sampling_percentage.setValue(config['ants_sampling_percentage'])
                if 'pheromone_decay_rate' in config: self.pheromone_decay_rate.setValue(config['pheromone_decay_rate'])
                if 'iterations' in config: self.iterations.setValue(config['iterations'])
                if 'background_value' in config: self.background_value.setValue(config['background_value'])
                if 'background_distance' in config: self.background_distance.setValue(config['background_distance'])
                if 'value_filter' in config: self.value_filter.setValue(config['value_filter'])
                if 'avoid_visited_threshold_enabled' in config: self.avoid_visited_enabled.setChecked(config['avoid_visited_threshold_enabled'])
                if 'avoid_visited_threshold' in config: self.avoid_visited_threshold.setValue(config['avoid_visited_threshold'])
                if 'ant_colony_interpolate_target' in config:
                    tgt = str(config['ant_colony_interpolate_target']).strip().lower()
                    self.ant_interpolate_target.setCurrentText('Domain' if tgt == 'domain' else 'Value')
                if 'molecular_clock_params' in config:
                    mc = config['molecular_clock_params']
                    if 'spatial_weight' in mc: self.mc_spatial_weight.setValue(mc['spatial_weight'])
                    if 'attr_weight' in mc: self.mc_attr_weight.setValue(mc['attr_weight'])
                    if 'ancestor_depth_offset' in mc: self.mc_ancestor_depth_offset.setValue(mc['ancestor_depth_offset'])
                    if 'branch_threshold' in mc: self.mc_branch_threshold.setValue(mc['branch_threshold'])
                    if 'min_samples' in mc: self.mc_min_samples.setValue(mc['min_samples'])
                    if 'max_samples' in mc: self.mc_max_samples.setValue(mc['max_samples'])
                    if 'detect_multiple' in mc: self.mc_detect_multiple.setChecked(mc['detect_multiple'])
                    if 'interp_method' in mc: self.mc_interp_method.setCurrentText(mc['interp_method'])
                if 'gaussian_kernel_params' in config:
                    gk = config['gaussian_kernel_params']
                    if 'bandwidth' in gk: self.gk_bandwidth.setValue(gk['bandwidth'])
                    if 'cutoff_sigma' in gk: self.gk_cutoff_sigma.setValue(gk['cutoff_sigma'])
                    if 'use_nearest_fallback' in gk: self.gk_use_nearest_fallback.setChecked(bool(gk['use_nearest_fallback']))
                    if 'fill_background' in gk:
                        self.gk_fill_background.setChecked(bool(gk['fill_background']))
                        self.gk_background_value.setEnabled(bool(gk['fill_background']))
                    if 'background_value' in gk: self.gk_background_value.setValue(gk['background_value'])
                if 'adaptive_octree_params' in config:
                    octree = config['adaptive_octree_params']
                    if 'output_mode' in octree:
                        mode = str(octree['output_mode']).strip().lower()
                        self.octree_output_mode.setCurrentText('Adaptive Leaf Cover' if mode == 'adaptive_leaf_cover' else 'Dense Blocks Cover')
                    if 'max_levels' in octree: self.octree_max_levels.setValue(int(octree['max_levels']))
                    if 'support_density_alpha' in octree: self.octree_support_density_alpha.setValue(float(octree['support_density_alpha']))
                    if 'include_dense_provenance' in octree: self.octree_include_dense_provenance.setChecked(bool(octree['include_dense_provenance']))
                if 'string_theory_params' in config:
                    st = config['string_theory_params']
                    if 'interpolate_target' in st:
                        tgt = str(st['interpolate_target']).strip().lower()
                        self.st_interpolate_target.setCurrentText('Domain' if tgt == 'domain' else 'Value')
                        if hasattr(self, '_update_st_target_ui'):
                            self._update_st_target_ui()
                    if 'distance_threshold' in st: self.st_distance_threshold.setValue(st['distance_threshold'])
                    if 'grade_difference' in st: self.st_grade_difference.setValue(st['grade_difference'])
                    if 'connect_to_all' in st: 
                        self.st_connect_to_all.setChecked(st['connect_to_all'])
                        self.st_max_connections.setEnabled(not st['connect_to_all'])
                        self.st_min_connections.setEnabled(not st['connect_to_all'])
                    if 'max_connections' in st: self.st_max_connections.setValue(st['max_connections'])
                    if 'min_connections' in st: self.st_min_connections.setValue(st['min_connections'])
                    if 'collision_policy' in st: self.st_collision_policy.setCurrentText(st['collision_policy'])
                    if 'processing_order' in st: self.st_processing_order.setCurrentText(st['processing_order'])
                    if 'filter_by_frequency' in st: 
                        self.st_filter_by_frequency.setChecked(st['filter_by_frequency'])
                        self.st_min_azimuth_freq.setEnabled(st['filter_by_frequency'])
                        self.st_min_dip_freq.setEnabled(st['filter_by_frequency'])
                    if 'min_azimuth_freq' in st: self.st_min_azimuth_freq.setValue(st['min_azimuth_freq'])
                    if 'min_dip_freq' in st: self.st_min_dip_freq.setValue(st['min_dip_freq'])
                # Backward compatibility for tolerance params (ignore or convert?)
                # If old params exist but new ones don't, we could try to map them, but they are different concepts.
                # We'll just ignore them for now.
                if 'average_with_blocks' in config: self.average_with_blocks.setChecked(config['average_with_blocks'])
                if 'process_domains_sequentially' in config: self.process_domains_sequentially.setChecked(config['process_domains_sequentially'])
                self.expand_interpolation_exports_to_subblocks.setChecked(bool(config.get('expand_interpolation_exports_to_subblocks', True)))
                self.force_rebuild_sample_blocks.setChecked(bool(config.get('force_rebuild_sample_blocks', False)))
                if 'blank_sample_domain_behavior' in config:
                    behavior = str(config['blank_sample_domain_behavior']).strip().lower()
                    self.blank_sample_domain_behavior.setCurrentText('Infer From Blocks' if behavior == 'infer_from_blocks' else 'Skip')
                if 'verbose' in config: self.verbose.setChecked(config['verbose'])
                if 'viewer_backend' in config: self.viewer_backend = _normalize_viewer_backend(config['viewer_backend'])
                if 'taichi_sample_diameter' in config: self.taichi_sample_diameter.setValue(config['taichi_sample_diameter'])
                self._legacy_fill_unvisited_domainwise = bool(config.get('fill_unvisited_domainwise', False))
                if 'domain_algorithm_overrides' in config:
                    self.domain_overrides = {
                        str(domain): dict(settings)
                        for domain, settings in (config['domain_algorithm_overrides'] or {}).items()
                    }
                    if self._legacy_fill_unvisited_domainwise:
                        for settings in self.domain_overrides.values():
                            if not settings.get('skip', False):
                                settings.setdefault('post_process', 'fill_with_average')
                count = len(self.domain_overrides)
                if count > 0:
                    self.domain_mapping_btn.setText(f'Configure Domain Algorithms... ({count} configured)')
                else:
                    self.domain_mapping_btn.setText('Configure Domain Algorithms...')
            finally:
                self._suspend_auto_viewer_refresh = False
        def save_config(self):
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Configuration", ".", "JSON Files (*.json)")
            if path:
                try:
                    config = self.to_dict()
                    with open(path, 'w') as f:
                        json.dump(config, f, indent=4)
                    QtWidgets.QMessageBox.information(self, "Success", "Configuration saved successfully.")
                except Exception as e:
                    QtWidgets.QMessageBox.critical(self, "Error", f"Failed to save configuration: {e}")

        def load_config(self):
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load Configuration", ".", "JSON Files (*.json)")
            if path:
                try:
                    with open(path, 'r') as f:
                        config = json.load(f)
                    self.from_dict(config)
                    QtWidgets.QMessageBox.information(self, "Success", "Configuration loaded successfully.")
                except Exception as e:
                    QtWidgets.QMessageBox.critical(self, "Error", f"Failed to load configuration: {e}")

        def open_domain_mapping(self):
            """Open dialog to configure domain-specific algorithms"""
            blocks_file = self.blocks_edit.text().strip()
            blocks_delimiter = self.blocks_delim.currentText() or detect_csv_delimiter(blocks_file)
            blocks_header_line = self.blocks_header_line.value()
            block_domain_col = self.block_domain_col.currentText()
            
            if not blocks_file or not os.path.isfile(blocks_file):
                QtWidgets.QMessageBox.warning(self, 'Warning', 
                    'Please select a valid blocks file first.')
                return
            
            if not block_domain_col or block_domain_col == '(None)':
                QtWidgets.QMessageBox.warning(self, 'Warning',
                    'Please select a domain column in "Blocks Columns" first.')
                return

            signature = build_domain_catalog_cache_signature(
                blocks_file,
                blocks_delimiter,
                blocks_header_line,
                block_domain_col,
                block_filters=self.block_volume_weighted_filters,
            )
            cached_catalog = self._domain_catalog_cache or {}
            domains = None
            if cached_catalog.get('signature') == signature:
                domains = list(cached_catalog.get('domains', []))

            if domains is None:
                progress = QtWidgets.QProgressDialog('Preparing domain catalog...', None, 0, 100, self)
                progress.setWindowTitle('Loading Domains')
                progress.setWindowModality(QtCore.Qt.WindowModal)
                progress.setMinimumDuration(0)
                progress.setAutoClose(False)
                progress.setAutoReset(False)
                progress.setValue(0)
                progress.show()
                QtWidgets.QApplication.processEvents()

                def update_progress(bytes_read, total_bytes, label):
                    total_bytes = max(int(total_bytes), 1)
                    bytes_read = max(0, min(int(bytes_read), total_bytes))
                    percent = int((bytes_read / total_bytes) * 100)
                    progress.setValue(percent)
                    progress.setLabelText(
                        f'{label}...\n{bytes_read / 1024**2:.1f} / {total_bytes / 1024**2:.1f} MiB'
                    )
                    QtWidgets.QApplication.processEvents()

                try:
                    domains = load_block_domain_catalog(
                        blocks_file,
                        blocks_delimiter,
                        blocks_header_line,
                        block_domain_col,
                        block_filters=self.block_volume_weighted_filters,
                        progress_callback=update_progress,
                    )
                    self._domain_catalog_cache = {
                        'signature': signature,
                        'domains': list(domains),
                    }
                except Exception as e:
                    progress.close()
                    QtWidgets.QMessageBox.critical(self, 'Error', f'Failed to load domains: {e}')
                    return
                finally:
                    if progress.value() < 100:
                        progress.setValue(100)
                    progress.close()

            if not domains:
                QtWidgets.QMessageBox.warning(self, 'Warning', 'No domains found in blocks file.')
                return
            
            dialog = DomainAlgorithmDialog(domains, self)
            
            # Load existing configuration
            if self.domain_overrides:
                dialog_configs = {
                    str(domain): dict(settings)
                    for domain, settings in self.domain_overrides.items()
                }
                if self._legacy_fill_unvisited_domainwise:
                    for domain in domains:
                        settings = dialog_configs.setdefault(domain, {})
                        if not settings.get('skip', False):
                            settings.setdefault('post_process', 'fill_with_average')
                dialog.set_domain_configs(dialog_configs)
            elif self._legacy_fill_unvisited_domainwise:
                dialog.set_domain_configs({
                    domain: {'post_process': 'fill_with_average'}
                    for domain in domains
                })
            
            if dialog.exec_() == QtWidgets.QDialog.Accepted:
                self.domain_overrides = dialog.get_domain_configs()
                self._legacy_fill_unvisited_domainwise = False
                count = len(self.domain_overrides)
                if count > 0:
                    self.domain_mapping_btn.setText(f'Configure Domain Algorithms... ({count} configured)')
                else:
                    self.domain_mapping_btn.setText('Configure Domain Algorithms...')

        def _build_filter_data_source(self, csv_file, delimiter=None, header_line=1):
            if not csv_file or not os.path.isfile(csv_file):
                raise ValueError('Please select a valid file first.')
            return FilterDataSource(
                csv_file,
                delimiter=delimiter,
                header_line=header_line,
            )

        def _get_selected_block_value_transfer_columns(self):
            if not hasattr(self, 'block_value_transfer_cols'):
                return []
            return [item.text() for item in self.block_value_transfer_cols.selectedItems()]

        def _set_selected_block_value_transfer_columns(self, column_names):
            if not hasattr(self, 'block_value_transfer_cols'):
                return
            _restore_list_widget_selection(self.block_value_transfer_cols, column_names)
            self._update_block_value_transfer_summary()

        def _update_block_value_transfer_summary(self):
            if not hasattr(self, 'block_value_transfer_summary'):
                return
            total = self.block_value_transfer_cols.count() if hasattr(self, 'block_value_transfer_cols') else 0
            selected = len(self._get_selected_block_value_transfer_columns())
            if total == 0:
                self.block_value_transfer_summary.setText('No transferable block columns are available.')
                return
            if selected == 0:
                self.block_value_transfer_summary.setText(
                    f'0 of {total} available block columns selected. Choose one or more columns to copy onto the samples file.'
                )
                return
            self.block_value_transfer_summary.setText(
                f'{selected} of {total} available block columns selected for transfer.'
            )

        def _refresh_block_value_transfer_columns(self):
            if not hasattr(self, 'block_value_transfer_cols'):
                return

            current_selection = set(self._get_selected_block_value_transfer_columns())
            pending_selection = set(getattr(self, '_pending_block_value_transfer_cols', []) or [])
            desired_selection = pending_selection or current_selection
            self.block_value_transfer_cols.clear()

            path = self.blocks_edit.text().strip()
            if not os.path.isfile(path):
                self._pending_block_value_transfer_cols = []
                self._update_block_value_transfer_summary()
                refresher = getattr(self, '_refresh_block_value_transfer_output_path', None)
                if callable(refresher):
                    refresher()
                return

            try:
                cols = parse_effective_header_line(path, self.blocks_delim.currentText(), self.blocks_header_line.value())
                excluded = {
                    str(self.block_x_col.currentText() or '').strip(),
                    str(self.block_y_col.currentText() or '').strip(),
                    str(self.block_z_col.currentText() or '').strip(),
                }
                excluded = {value for value in excluded if value and value != '(None)'}
                for column_name in cols:
                    if column_name in excluded:
                        continue
                    item = QtWidgets.QListWidgetItem(column_name)
                    item.setSelected(column_name in desired_selection)
                    self.block_value_transfer_cols.addItem(item)
            except Exception:
                pass

            self._pending_block_value_transfer_cols = []
            self._update_block_value_transfer_summary()
            refresher = getattr(self, '_refresh_block_value_transfer_output_path', None)
            if callable(refresher):
                refresher()

        def _get_selected_block_model_transfer_columns(self):
            if not hasattr(self, 'block_model_transfer_cols'):
                return []
            return [item.text() for item in self.block_model_transfer_cols.selectedItems()]

        def _set_selected_block_model_transfer_columns(self, column_names):
            if not hasattr(self, 'block_model_transfer_cols'):
                return
            _restore_list_widget_selection(self.block_model_transfer_cols, column_names)
            self._update_block_model_transfer_summary()

        def _update_block_model_transfer_summary(self):
            if not hasattr(self, 'block_model_transfer_summary'):
                return
            total = self.block_model_transfer_cols.count()
            selected = len(self._get_selected_block_model_transfer_columns())
            if not total:
                self.block_model_transfer_summary.setText('No transferable source columns are available.')
                return
            fallback_state = 'disabled'
            if self.block_model_nearest_fallback.isChecked():
                distance_limit = float(self.block_model_nearest_max_distance.value())
                fallback_state = (
                    'enabled with unlimited distance'
                    if distance_limit <= 0 else
                    f'enabled up to {distance_limit:g} m'
                )
            self.block_model_transfer_summary.setText(
                f'{selected} of {total} source columns selected. No-overlap fallback is {fallback_state}.'
            )

        def _update_block_model_transfer_fallback_controls(self):
            enabled = bool(self.block_model_nearest_fallback.isChecked())
            self.block_model_nearest_max_distance.setEnabled(enabled)
            self._update_block_model_transfer_summary()

        def _populate_block_model_column_combo(self, combo, columns, first_label=None, suggestions=()):
            previous = combo.currentText()
            blocker = QtCore.QSignalBlocker(combo)
            try:
                combo.clear()
                if first_label:
                    combo.addItem(first_label)
                combo.addItems([str(column) for column in columns])
                if previous and combo.findText(previous) >= 0:
                    combo.setCurrentText(previous)
                else:
                    lowered = {str(value).lower() for value in suggestions}
                    for index in range(combo.count()):
                        if combo.itemText(index).lower() in lowered:
                            combo.setCurrentIndex(index)
                            break
            finally:
                del blocker

        def _refresh_block_model_transfer_source_columns(self):
            if not hasattr(self, 'block_model_transfer_cols'):
                return
            path = self.blocks_edit.text().strip()
            columns = []
            if os.path.isfile(path):
                try:
                    columns = parse_effective_header_line(path, self.blocks_delim.currentText(), self.blocks_header_line.value())
                except Exception:
                    columns = []
            for combo, suggestions in zip(
                self.block_model_source_size_cols,
                [('dx', 'dim_x', 'size_x', 'x_size'), ('dy', 'dim_y', 'size_y', 'y_size'), ('dz', 'dim_z', 'size_z', 'z_size')],
            ):
                self._populate_block_model_column_combo(combo, columns, '(Infer)', suggestions)

            selected = set(self._get_selected_block_model_transfer_columns())
            selected.update(getattr(self, '_pending_block_model_transfer_cols', []) or [])
            excluded = {
                self.block_x_col.currentText(), self.block_y_col.currentText(), self.block_z_col.currentText(),
            }
            self.block_model_transfer_cols.clear()
            for column in columns:
                if column in excluded:
                    continue
                item = QtWidgets.QListWidgetItem(column)
                self.block_model_transfer_cols.addItem(item)
            self._set_selected_block_model_transfer_columns(selected)
            self._pending_block_model_transfer_cols = []
            self._update_block_model_transfer_summary()

        def _refresh_block_model_target_columns(self):
            if not hasattr(self, 'block_model_target_edit'):
                return
            path = self.block_model_target_edit.text().strip()
            columns = []
            if os.path.isfile(path):
                try:
                    if not is_bmf_file(path):
                        detected = detect_csv_delimiter(path)
                        if self.block_model_target_delim.findText(detected) >= 0:
                            self.block_model_target_delim.setCurrentText(detected)
                        sync_csv_header_line_widget(
                            self.block_model_target_header_line,
                            path,
                            self.block_model_target_header_line.value(),
                        )
                        columns = parse_effective_header_line(
                            path, self.block_model_target_delim.currentText(), self.block_model_target_header_line.value(),
                        )
                    else:
                        blocker = QtCore.QSignalBlocker(self.block_model_target_header_line)
                        self.block_model_target_header_line.setValue(1)
                        del blocker
                        preview, _ = _load_bmf_dataframe(path, row_limit=1)
                        columns = list(preview.columns)
                except Exception:
                    columns = []
            for combo, suggestions in zip(
                [self.block_model_target_x_col, self.block_model_target_y_col, self.block_model_target_z_col],
                [('x', 'easting'), ('y', 'northing'), ('z', 'elevation', 'rl')],
            ):
                self._populate_block_model_column_combo(combo, columns, None, suggestions)
            for combo, suggestions in zip(
                self.block_model_target_size_cols,
                [('dx', 'dim_x', 'size_x', 'x_size'), ('dy', 'dim_y', 'size_y', 'y_size'), ('dz', 'dim_z', 'size_z', 'z_size')],
            ):
                self._populate_block_model_column_combo(combo, columns, '(Infer)', suggestions)

        def _get_selected_table_attribute_key_columns(self):
            if not hasattr(self, 'table_attribute_key_cols'):
                return []
            return [item.text() for item in self.table_attribute_key_cols.selectedItems()]

        def _set_selected_table_attribute_key_columns(self, column_names):
            if not hasattr(self, 'table_attribute_key_cols'):
                return
            _restore_list_widget_selection(self.table_attribute_key_cols, column_names)

        def _get_selected_table_attribute_value_columns(self):
            if not hasattr(self, 'table_attribute_value_cols'):
                return []
            return [item.text() for item in self.table_attribute_value_cols.selectedItems()]

        def _set_selected_table_attribute_value_columns(self, column_names):
            if not hasattr(self, 'table_attribute_value_cols'):
                return
            _restore_list_widget_selection(self.table_attribute_value_cols, column_names)

        def _update_table_attribute_summary(self):
            if not hasattr(self, 'table_attribute_summary'):
                return
            total_keys = self.table_attribute_key_cols.count()
            selected_keys = len(self._get_selected_table_attribute_key_columns())
            total_values = self.table_attribute_value_cols.count()
            selected_values = len(self._get_selected_table_attribute_value_columns())
            if total_keys == 0:
                self.table_attribute_summary.setText('No shared key columns are available between the block model and the table.')
                return
            self.table_attribute_summary.setText(
                f'{selected_keys} of {total_keys} match keys selected. '
                f'{selected_values} of {total_values} table columns selected for assignment.'
            )

        def _read_table_attribute_source_columns(self, path, delimiter_combo, header_line_spin):
            columns = []
            if not path or not os.path.isfile(path):
                return columns
            try:
                if not is_bmf_file(path):
                    detected = detect_csv_delimiter(path)
                    if delimiter_combo.findText(detected) >= 0:
                        delimiter_combo.setCurrentText(detected)
                    sync_csv_header_line_widget(header_line_spin, path, header_line_spin.value())
                    columns = parse_effective_header_line(path, delimiter_combo.currentText(), header_line_spin.value())
                else:
                    blocker = QtCore.QSignalBlocker(header_line_spin)
                    header_line_spin.setValue(1)
                    del blocker
                    preview, _ = _load_bmf_dataframe(path, row_limit=1)
                    columns = list(preview.columns)
            except Exception:
                columns = []
            return columns

        def _refresh_table_attribute_shared_key_columns(self):
            if not hasattr(self, 'table_attribute_key_cols'):
                return
            block_columns = self._read_table_attribute_source_columns(
                self.table_attribute_block_model_edit.text().strip(),
                self.table_attribute_block_model_delim,
                self.table_attribute_block_model_header_line,
            )
            table_columns = self._read_table_attribute_source_columns(
                self.table_attribute_table_edit.text().strip(),
                self.table_attribute_table_delim,
                self.table_attribute_table_header_line,
            )
            table_column_set = set(table_columns)
            shared_columns = [column_name for column_name in block_columns if column_name in table_column_set]
            selected_keys = set(self._get_selected_table_attribute_key_columns())
            selected_keys.update(getattr(self, '_pending_table_attribute_key_cols', []) or [])
            self.table_attribute_key_cols.clear()
            for column_name in shared_columns:
                self.table_attribute_key_cols.addItem(QtWidgets.QListWidgetItem(column_name))
            self._set_selected_table_attribute_key_columns(selected_keys)
            self._pending_table_attribute_key_cols = []
            self._refresh_table_attribute_value_columns()

        def _refresh_table_attribute_value_columns(self):
            if not hasattr(self, 'table_attribute_value_cols'):
                return
            table_columns = self._read_table_attribute_source_columns(
                self.table_attribute_table_edit.text().strip(),
                self.table_attribute_table_delim,
                self.table_attribute_table_header_line,
            )
            selected_values = set(self._get_selected_table_attribute_value_columns())
            selected_values.update(getattr(self, '_pending_table_attribute_value_cols', []) or [])
            selected_keys = set(self._get_selected_table_attribute_key_columns())
            self.table_attribute_value_cols.clear()
            for column_name in table_columns:
                if column_name in selected_keys:
                    continue
                self.table_attribute_value_cols.addItem(QtWidgets.QListWidgetItem(column_name))
            self._set_selected_table_attribute_value_columns(selected_values)
            self._pending_table_attribute_value_cols = []
            self._update_table_attribute_summary()

        def _get_selected_block_domain_metrics(self):
            if not hasattr(self, 'block_domain_metrics_metric_list'):
                return []
            selected_metrics = []
            for index in range(self.block_domain_metrics_metric_list.count()):
                item = self.block_domain_metrics_metric_list.item(index)
                if item is not None and item.checkState() == QtCore.Qt.Checked:
                    metric_id = str(item.data(QtCore.Qt.UserRole) or '').strip()
                    if metric_id:
                        selected_metrics.append(metric_id)
            return selected_metrics

        def _set_selected_block_domain_metrics(self, metric_ids):
            if not hasattr(self, 'block_domain_metrics_metric_list'):
                return
            desired = set(_normalize_block_domain_metric_selection(metric_ids))
            blocker = QtCore.QSignalBlocker(self.block_domain_metrics_metric_list)
            try:
                for index in range(self.block_domain_metrics_metric_list.count()):
                    item = self.block_domain_metrics_metric_list.item(index)
                    metric_id = str(item.data(QtCore.Qt.UserRole) or '').strip()
                    item.setCheckState(QtCore.Qt.Checked if metric_id in desired else QtCore.Qt.Unchecked)
            finally:
                del blocker
            self._update_block_domain_metrics_metric_summary()

        def _set_all_block_domain_metrics_checked(self, checked):
            self._set_selected_block_domain_metrics(
                BLOCK_DOMAIN_METRIC_DEFINITION_ORDER if checked else []
            )

        def _update_sample_blocks_id_controls(self):
            enabled = bool(getattr(self, 'sample_blocks_include_ids', None) and self.sample_blocks_include_ids.isChecked())
            for cb in getattr(self, 'sample_blocks_id_cols', []):
                cb.setEnabled(enabled)

        def _update_block_domain_metrics_metric_summary(self):
            if not hasattr(self, 'block_domain_metrics_metric_summary'):
                return
            selected_metrics = self._get_selected_block_domain_metrics()
            selected_set = set(selected_metrics)
            if not selected_metrics:
                self.block_domain_metrics_metric_summary.setText(
                    'No metrics selected. Choose at least one metric before exporting.'
                )
            else:
                selected_labels = [
                    BLOCK_DOMAIN_METRIC_DEFINITIONS_BY_ID[metric_id]['label']
                    for metric_id in selected_metrics
                    if metric_id in BLOCK_DOMAIN_METRIC_DEFINITIONS_BY_ID
                ]
                self.block_domain_metrics_metric_summary.setText(
                    f"{len(selected_metrics)} metrics enabled: {', '.join(selected_labels)}"
                )
            if hasattr(self, 'block_domain_metrics_knn_k'):
                self.block_domain_metrics_knn_k.setEnabled('average_distance_knn' in selected_set)
            distance_summary_enabled = 'distance_band_summary' in selected_set
            if hasattr(self, 'block_domain_metrics_distance_step'):
                self.block_domain_metrics_distance_step.setEnabled(distance_summary_enabled)
            if hasattr(self, 'block_domain_metrics_max_factor'):
                self.block_domain_metrics_max_factor.setEnabled(distance_summary_enabled)
            closest_id_enabled = 'closest_sample_id' in selected_set
            for cb in getattr(self, 'block_domain_metrics_id_cols', []):
                cb.setEnabled(closest_id_enabled)
            if hasattr(self, 'block_domain_metrics_value_col'):
                self.block_domain_metrics_value_col.setEnabled(bool(selected_set & BLOCK_DOMAIN_BLOCK_VALUE_METRICS))
            if hasattr(self, 'block_domain_metrics_prefix_with_block_value'):
                self.block_domain_metrics_prefix_with_block_value.setEnabled(bool(selected_set & BLOCK_DOMAIN_SAMPLE_VALUE_METRICS))

        def _get_bmf_export_cell_size(self):
            if not all(hasattr(self, name) for name in ['bmf_export_cell_x', 'bmf_export_cell_y', 'bmf_export_cell_z']):
                return None
            values = [
                float(self.bmf_export_cell_x.value()),
                float(self.bmf_export_cell_y.value()),
                float(self.bmf_export_cell_z.value()),
            ]
            positive = [value for value in values if value > 0]
            if not positive:
                return None
            if len(positive) != 3:
                raise ValueError('BMF cell size must be set for X, Y, and Z, or left as Auto for all three axes.')
            return values

        def _set_bmf_export_cell_size(self, cell_size):
            if not all(hasattr(self, name) for name in ['bmf_export_cell_x', 'bmf_export_cell_y', 'bmf_export_cell_z']):
                return
            spins = [self.bmf_export_cell_x, self.bmf_export_cell_y, self.bmf_export_cell_z]
            if not cell_size:
                for spin in spins:
                    spin.setValue(0.0)
                return
            try:
                values = [float(value) for value in cell_size]
            except Exception:
                values = []
            if len(values) != 3 or any(value <= 0 for value in values):
                for spin in spins:
                    spin.setValue(0.0)
                return
            for spin, value in zip(spins, values):
                spin.setValue(value)

        def _update_bmf_export_regularize_warning(self):
            if not hasattr(self, 'bmf_export_regularize_warning') or not hasattr(self, 'bmf_export_regularize_base_blocks'):
                return
            if self.bmf_export_regularize_base_blocks.isChecked():
                self.bmf_export_regularize_warning.setText(
                    'Warning: regularization aggregates source CSV rows to base blocks and may create a dense regular grid.'
                )
            else:
                self.bmf_export_regularize_warning.setText(
                    'Info: tbms-config-text preserves source CSV rows as BMF rows instead of regularizing to a dense grid.'
                )
            self.bmf_export_regularize_warning.setVisible(True)

        def _get_selected_bmf_export_value_columns(self):
            if not hasattr(self, 'bmf_export_value_cols'):
                return []
            selected = []
            for row_index in range(self.bmf_export_value_cols.rowCount()):
                item = self.bmf_export_value_cols.item(row_index, 0)
                if item is not None and item.checkState() == QtCore.Qt.Checked:
                    selected.append(item.text())
            return selected

        def _normalize_bmf_export_value_exceptions(self, raw_exceptions):
            normalized = {}
            if isinstance(raw_exceptions, list):
                for entry in raw_exceptions:
                    if not isinstance(entry, dict):
                        continue
                    column_name = str(entry.get('column') or '').strip()
                    bad_value = str(entry.get('value') or '')
                    replacement = '' if entry.get('replacement') is None else str(entry.get('replacement'))
                    include_in_regularization = self._coerce_bmf_exception_include_flag(
                        entry.get('include_in_regularization', entry.get('include_regularization', False))
                    )
                    if column_name and bad_value:
                        normalized.setdefault(column_name, {})[bad_value] = {
                            'replacement': replacement,
                            'include_in_regularization': include_in_regularization,
                        }
                return normalized
            for raw_column, raw_rules in dict(raw_exceptions or {}).items():
                column_name = str(raw_column or '').strip()
                if not column_name or not isinstance(raw_rules, dict):
                    continue
                for raw_value, raw_rule in raw_rules.items():
                    bad_value = str(raw_value)
                    if not bad_value:
                        continue
                    if isinstance(raw_rule, dict):
                        replacement = raw_rule.get('replacement', '')
                        include_in_regularization = self._coerce_bmf_exception_include_flag(
                            raw_rule.get('include_in_regularization', raw_rule.get('include_regularization', False))
                        )
                    else:
                        replacement = raw_rule
                        include_in_regularization = False
                    normalized.setdefault(column_name, {})[bad_value] = {
                        'replacement': '' if replacement is None else str(replacement),
                        'include_in_regularization': include_in_regularization,
                    }
            return normalized

        def _coerce_bmf_exception_include_flag(self, value):
            if isinstance(value, bool):
                return value
            if value is None:
                return False
            text = str(value).strip().lower()
            return text in {'1', 'true', 'yes', 'y', 'on', 'checked'}

        def _make_bmf_exception_include_item(self, include_in_regularization=False):
            item = QtWidgets.QTableWidgetItem('')
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Checked if include_in_regularization else QtCore.Qt.Unchecked)
            item.setToolTip('When checked, the replacement is applied before base-block regularization and can contribute to numeric averages.')
            return item

        def _extract_bmf_exception_rule_fields(self, raw_rule):
            if isinstance(raw_rule, dict):
                replacement = raw_rule.get('replacement', '')
                include_in_regularization = self._coerce_bmf_exception_include_flag(
                    raw_rule.get('include_in_regularization', raw_rule.get('include_regularization', False))
                )
            else:
                replacement = raw_rule
                include_in_regularization = False
            return '' if replacement is None else str(replacement), include_in_regularization

        def _get_bmf_export_value_exceptions(self):
            if not hasattr(self, 'bmf_export_exception_table'):
                return {}
            rules = {}
            for row_index in range(self.bmf_export_exception_table.rowCount()):
                column_item = self.bmf_export_exception_table.item(row_index, 0)
                value_item = self.bmf_export_exception_table.item(row_index, 1)
                replacement_item = self.bmf_export_exception_table.item(row_index, 2)
                include_item = self.bmf_export_exception_table.item(row_index, 3)
                column_name = str(column_item.text() if column_item is not None else '').strip()
                bad_value = str(value_item.text() if value_item is not None else '')
                replacement = str(replacement_item.text() if replacement_item is not None else '')
                include_in_regularization = bool(
                    include_item is not None and include_item.checkState() == QtCore.Qt.Checked
                )
                if column_name and bad_value:
                    if include_in_regularization:
                        rules.setdefault(column_name, {})[bad_value] = {
                            'replacement': replacement,
                            'include_in_regularization': True,
                        }
                    else:
                        rules.setdefault(column_name, {})[bad_value] = replacement
            return rules

            def _normalize_optional_size_columns(size_columns):
                values = list(size_columns or ())
                if len(values) != 3:
                    return None
                normalized = [str(value or '').strip() for value in values]
                if any(not value or value in {'(None)', '(Infer)'} for value in normalized):
                    return None
                return tuple(normalized)

            def _infer_block_row_bounds(coords, base_block_size):
                coords = np.asarray(coords, dtype=float)
                base_dims = np.asarray(base_block_size, dtype=float)
                if coords.ndim != 2 or coords.shape[1] != 3:
                    raise ValueError('Block coordinates must be an Nx3 array.')
                if base_dims.shape != (3,) or np.any(~np.isfinite(base_dims)) or np.any(base_dims <= 0):
                    raise ValueError('Base block size must contain three positive values.')
                if len(coords) == 0:
                    return np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=float)

                grid_origin = np.floor(coords.min(axis=0) / base_dims) * base_dims
                parent_indices = np.floor((coords - grid_origin) / base_dims + 1e-6).astype(np.int64)
                parent_origins = grid_origin + parent_indices * base_dims
                local_coords = np.clip(coords - parent_origins, 0.0, base_dims)
                lower_bounds = np.empty_like(coords, dtype=float)
                upper_bounds = np.empty_like(coords, dtype=float)

                _, inverse = np.unique(parent_indices, axis=0, return_inverse=True)
                order = np.argsort(inverse, kind='stable')
                sorted_inverse = inverse[order]
                starts = np.flatnonzero(np.r_[True, sorted_inverse[1:] != sorted_inverse[:-1]])
                ends = np.r_[starts[1:], len(order)]
                for start, end in zip(starts, ends):
                    row_indices = order[start:end]
                    group_local = local_coords[row_indices]
                    group_origin = parent_origins[row_indices[0]]
                    for axis in range(3):
                        tolerance = max(base_dims[axis] * 1e-7, 1e-9)
                        centers = _cluster_axis_centers(group_local[:, axis], tolerance)
                        boundaries = np.empty(len(centers) + 1, dtype=float)
                        boundaries[0] = 0.0
                        boundaries[-1] = base_dims[axis]
                        if len(centers) > 1:
                            boundaries[1:-1] = (centers[:-1] + centers[1:]) / 2.0
                        center_indices = np.abs(group_local[:, axis, None] - centers[None, :]).argmin(axis=1)
                        lower_bounds[row_indices, axis] = group_origin[axis] + boundaries[center_indices]
                        upper_bounds[row_indices, axis] = group_origin[axis] + boundaries[center_indices + 1]

                return lower_bounds, upper_bounds, upper_bounds - lower_bounds

            def _resolve_block_row_geometry(df, coordinate_columns, base_block_size, size_columns=None):
                coordinate_columns = tuple(coordinate_columns or ())
                if len(coordinate_columns) != 3 or any(column not in df.columns for column in coordinate_columns):
                    raise ValueError('Three valid block coordinate columns are required.')
                coord_frame = df[list(coordinate_columns)].apply(pd.to_numeric, errors='coerce')
                valid_mask = coord_frame.notna().all(axis=1).to_numpy()
                coords = coord_frame.loc[valid_mask].to_numpy(dtype=float, copy=False)

                explicit_size_columns = _normalize_optional_size_columns(size_columns)
                if explicit_size_columns is not None:
                    missing = [column for column in explicit_size_columns if column not in df.columns]
                    if missing:
                        raise ValueError(f"Block size column(s) not found: {', '.join(missing)}")
                    size_frame = df.loc[valid_mask, list(explicit_size_columns)].apply(pd.to_numeric, errors='coerce')
                    sizes = size_frame.to_numpy(dtype=float, copy=False)
                    valid_sizes = np.isfinite(sizes).all(axis=1) & (sizes > 0).all(axis=1)
                    valid_positions = np.flatnonzero(valid_mask)
                    valid_mask[valid_positions[~valid_sizes]] = False
                    coords = coords[valid_sizes]
                    sizes = sizes[valid_sizes]
                    lower_bounds = coords - sizes / 2.0
                    upper_bounds = coords + sizes / 2.0
                    geometry_mode = 'explicit-size-columns'
                else:
                    lower_bounds, upper_bounds, sizes = _infer_block_row_bounds(coords, base_block_size)
                    geometry_mode = 'inferred-from-base-grid'

                return {
                    'valid_mask': valid_mask,
                    'row_indices': np.flatnonzero(valid_mask),
                    'centers': coords,
                    'lower_bounds': lower_bounds,
                    'upper_bounds': upper_bounds,
                    'sizes': sizes,
                    'mode': geometry_mode,
                }

            def _detect_dataframe_transfer_column_modes(df, columns):
                modes = {}
                for column_name in columns:
                    values = df[column_name]
                    nonblank_mask = values.notna() & values.astype(str).str.strip().ne('') & values.astype(str).str.lower().ne('nan')
                    nonblank = values.loc[nonblank_mask]
                    modes[column_name] = (
                        'numeric'
                        if len(nonblank) > 0 and pd.to_numeric(nonblank, errors='coerce').notna().all()
                        else 'categorical'
                    )
                return modes

            def export_blocks_with_source_block_values(source_blocks_file, target_blocks_file, output_file=None,
                                                        source_delimiter=None, target_delimiter=None,
                                                        source_header_line=1, target_header_line=1,
                                                        source_x_col=None, source_y_col=None, source_z_col=None,
                                                        target_x_col=None, target_y_col=None, target_z_col=None,
                                                        source_value_cols=None,
                                                        source_block_size=None, target_block_size=None,
                                                        source_size_cols=None, target_size_cols=None,
                                                        nearest_fallback=True, progress_callback=None):
                """Transfer source block attributes to existing target rows using 3-D overlap volume.

                Numeric fields use an overlap-volume-weighted mean. Categorical fields use the
                category with the greatest total overlap volume. Targets without overlap can use
                the nearest source block center as a deterministic fallback.
                """
                if not source_blocks_file or not os.path.isfile(source_blocks_file):
                    raise ValueError('Please select a valid source blocks file.')
                if not target_blocks_file or not os.path.isfile(target_blocks_file):
                    raise ValueError('Please select a valid target blocks file.')
                if source_block_size is None or target_block_size is None:
                    raise ValueError('Source and target base block sizes are required.')

                selected_columns = _normalize_block_transfer_columns(
                    source_value_cols,
                    block_x_col=source_x_col,
                    block_y_col=source_y_col,
                    block_z_col=source_z_col,
                )
                output_file = resolve_block_model_transfer_export_path(output_file, target_blocks_file)
                _emit_progress(progress_callback, 0, 100, 'Reading source block model...')
                source_df, _ = load_full_blocks_dataframe(
                    source_blocks_file,
                    blocks_delimiter=source_delimiter,
                    blocks_header_line=source_header_line,
                    progress_label='Reading source block model',
                    progress_callback=_make_scaled_progress_callback(progress_callback, 0, 25, 'Reading source block model...'),
                )
                source_x_col, source_y_col, source_z_col = resolve_block_coordinate_columns(
                    list(source_df.columns), source_x_col, source_y_col, source_z_col,
                )
                missing_columns = [column for column in selected_columns if column not in source_df.columns]
                if missing_columns:
                    raise ValueError(f"Source block column(s) not found: {', '.join(missing_columns)}")

                _emit_progress(progress_callback, 25, 100, 'Reading target block model...')
                target_df, target_output_delimiter = load_full_blocks_dataframe(
                    target_blocks_file,
                    blocks_delimiter=target_delimiter,
                    blocks_header_line=target_header_line,
                    progress_label='Reading target block model',
                    progress_callback=_make_scaled_progress_callback(progress_callback, 25, 45, 'Reading target block model...'),
                )
                target_x_col, target_y_col, target_z_col = resolve_block_coordinate_columns(
                    list(target_df.columns), target_x_col, target_y_col, target_z_col,
                )

                _emit_progress(progress_callback, 45, 100, 'Resolving source and target block geometry...')
                source_geometry = _resolve_block_row_geometry(
                    source_df,
                    (source_x_col, source_y_col, source_z_col),
                    source_block_size,
                    source_size_cols,
                )
                target_geometry = _resolve_block_row_geometry(
                    target_df,
                    (target_x_col, target_y_col, target_z_col),
                    target_block_size,
                    target_size_cols,
                )
                if len(source_geometry['centers']) == 0:
                    raise ValueError('The source block model has no rows with valid coordinates and dimensions.')

                from scipy.spatial import cKDTree

                source_valid_df = source_df.iloc[source_geometry['row_indices']]
                column_modes = _detect_dataframe_transfer_column_modes(source_valid_df, selected_columns)
                source_column_values = {
                    column: (
                        pd.to_numeric(source_valid_df[column], errors='coerce').to_numpy(dtype=float, copy=False)
                        if mode == 'numeric'
                        else source_valid_df[column].fillna('').astype(str).str.strip().to_numpy(dtype=object, copy=False)
                    )
                    for column, mode in column_modes.items()
                }
                assigned_values = {column: [''] * len(target_df) for column in selected_columns}
                tree = cKDTree(source_geometry['centers'])
                source_half_diagonal_max = float(np.linalg.norm(source_geometry['sizes'], axis=1).max() / 2.0)
                overlap_matches = 0
                nearest_matches = 0
                total_valid_targets = len(target_geometry['centers'])
                _emit_progress(progress_callback, 50, 100, 'Matching target blocks to source blocks...')

                for local_index, target_row_index in enumerate(target_geometry['row_indices']):
                    target_center = target_geometry['centers'][local_index]
                    target_lower = target_geometry['lower_bounds'][local_index]
                    target_upper = target_geometry['upper_bounds'][local_index]
                    query_radius = float(np.linalg.norm(target_geometry['sizes'][local_index]) / 2.0 + source_half_diagonal_max)
                    candidate_indices = tree.query_ball_point(target_center, query_radius)
                    overlap_indices = np.empty(0, dtype=int)
                    overlap_volumes = np.empty(0, dtype=float)
                    if candidate_indices:
                        candidate_indices = np.asarray(candidate_indices, dtype=int)
                        overlap_lengths = np.maximum(
                            0.0,
                            np.minimum(source_geometry['upper_bounds'][candidate_indices], target_upper)
                            - np.maximum(source_geometry['lower_bounds'][candidate_indices], target_lower),
                        )
                        volumes = np.prod(overlap_lengths, axis=1)
                        positive_mask = volumes > max(float(np.prod(target_geometry['sizes'][local_index])) * 1e-12, 1e-12)
                        overlap_indices = candidate_indices[positive_mask]
                        overlap_volumes = volumes[positive_mask]

                    if len(overlap_indices) > 0:
                        overlap_matches += 1
                        for column_name, mode in column_modes.items():
                            values = source_column_values[column_name][overlap_indices]
                            if mode == 'numeric':
                                valid_values = np.isfinite(values)
                                if valid_values.any():
                                    assigned_values[column_name][target_row_index] = float(
                                        np.average(values[valid_values], weights=overlap_volumes[valid_values])
                                    )
                            else:
                                category_volumes = {}
                                for value, volume in zip(values, overlap_volumes):
                                    normalized = str(value).strip()
                                    if not normalized or normalized.lower() == 'nan':
                                        continue
                                    category_volumes[normalized] = category_volumes.get(normalized, 0.0) + float(volume)
                                if category_volumes:
                                    assigned_values[column_name][target_row_index] = sorted(
                                        category_volumes.items(), key=lambda item: (-item[1], item[0])
                                    )[0][0]
                    elif nearest_fallback:
                        _, nearest_index = tree.query(target_center, k=1)
                        nearest_index = int(nearest_index)
                        nearest_matches += 1
                        for column_name, mode in column_modes.items():
                            value = source_column_values[column_name][nearest_index]
                            if mode == 'numeric':
                                if np.isfinite(value):
                                    assigned_values[column_name][target_row_index] = float(value)
                            else:
                                normalized = str(value).strip()
                                if normalized and normalized.lower() != 'nan':
                                    assigned_values[column_name][target_row_index] = normalized

                    if progress_callback and (local_index + 1 == total_valid_targets or (local_index + 1) % 10_000 == 0):
                        percent = 50 + int(round(((local_index + 1) / max(total_valid_targets, 1)) * 47))
                        _emit_progress(progress_callback, percent, 100, 'Matching target blocks to source blocks...')

                output_df = target_df.copy()
                for column_name in selected_columns:
                    output_df[column_name] = pd.Series(assigned_values[column_name], index=output_df.index)
                output_dir = os.path.dirname(output_file) or '.'
                os.makedirs(output_dir, exist_ok=True)
                output_separator = target_output_delimiter if target_output_delimiter and target_output_delimiter != 'bmf' else ','
                _emit_progress(progress_callback, 98, 100, 'Writing block-model transfer export...')
                output_df.to_csv(output_file, index=False, sep=output_separator)
                _emit_progress(progress_callback, 100, 100, 'Block-model transfer complete.')

                invalid_target_count = int(len(target_df) - total_valid_targets)
                unmatched_count = int(total_valid_targets - overlap_matches - nearest_matches + invalid_target_count)
                print(
                    f"Block-model transfer complete: targets={len(target_df):,}; overlap matches={overlap_matches:,}; "
                    f"nearest fallback matches={nearest_matches:,}; unmatched={unmatched_count:,}."
                )
                return {
                    'output_file': output_file,
                    'total_target_blocks': int(len(target_df)),
                    'overlap_matched_blocks': int(overlap_matches),
                    'nearest_matched_blocks': int(nearest_matches),
                    'unmatched_blocks': unmatched_count,
                    'invalid_target_blocks': invalid_target_count,
                    'transferred_columns': list(selected_columns),
                    'column_modes': dict(column_modes),
                    'source_geometry_mode': source_geometry['mode'],
                    'target_geometry_mode': target_geometry['mode'],
                }
                value_item = self.bmf_export_exception_table.item(row_index, 1)
                replacement_item = self.bmf_export_exception_table.item(row_index, 2)
                include_item = self.bmf_export_exception_table.item(row_index, 3)
                column_name = str(column_item.text() if column_item is not None else '').strip()
                bad_value = str(value_item.text() if value_item is not None else '')
                replacement = str(replacement_item.text() if replacement_item is not None else '')
                include_in_regularization = bool(
                    include_item is not None and include_item.checkState() == QtCore.Qt.Checked
                )
                if column_name and bad_value:
                    if include_in_regularization:
                        rules.setdefault(column_name, {})[bad_value] = {
                            'replacement': replacement,
                            'include_in_regularization': True,
                        }
                    else:
                        rules.setdefault(column_name, {})[bad_value] = replacement
            return rules

        def _set_bmf_export_value_exceptions(self, raw_exceptions):
            if not hasattr(self, 'bmf_export_exception_table'):
                return
            rules = self._normalize_bmf_export_value_exceptions(raw_exceptions)
            blocker = QtCore.QSignalBlocker(self.bmf_export_exception_table)
            try:
                rows = []
                for column_name, column_rules in rules.items():
                    for bad_value, raw_rule in column_rules.items():
                        replacement, include_in_regularization = self._extract_bmf_exception_rule_fields(raw_rule)
                        rows.append((column_name, bad_value, replacement, include_in_regularization))
                self.bmf_export_exception_table.setRowCount(len(rows))
                for row_index, (column_name, bad_value, replacement, include_in_regularization) in enumerate(rows):
                    self.bmf_export_exception_table.setItem(row_index, 0, QtWidgets.QTableWidgetItem(column_name))
                    self.bmf_export_exception_table.setItem(row_index, 1, QtWidgets.QTableWidgetItem(bad_value))
                    self.bmf_export_exception_table.setItem(row_index, 2, QtWidgets.QTableWidgetItem(replacement))
                    self.bmf_export_exception_table.setItem(
                        row_index,
                        3,
                        self._make_bmf_exception_include_item(include_in_regularization),
                    )
            finally:
                del blocker
            self.bmf_export_value_exceptions = rules
            self._update_bmf_export_exception_summary()

        def _add_bmf_export_exception_row(self, column_name='', bad_value='', replacement='', include_in_regularization=False):
            if not hasattr(self, 'bmf_export_exception_table'):
                return
            if not column_name:
                selected_columns = self._get_selected_bmf_export_value_columns()
                column_name = selected_columns[0] if selected_columns else ''
            row_index = self.bmf_export_exception_table.rowCount()
            self.bmf_export_exception_table.insertRow(row_index)
            self.bmf_export_exception_table.setItem(row_index, 0, QtWidgets.QTableWidgetItem(str(column_name or '')))
            self.bmf_export_exception_table.setItem(row_index, 1, QtWidgets.QTableWidgetItem(str(bad_value or '')))
            self.bmf_export_exception_table.setItem(row_index, 2, QtWidgets.QTableWidgetItem('' if replacement is None else str(replacement)))
            self.bmf_export_exception_table.setItem(row_index, 3, self._make_bmf_exception_include_item(include_in_regularization))
            self.bmf_export_exception_table.setCurrentCell(row_index, 1)
            self._update_bmf_export_exception_summary()

        def _remove_selected_bmf_export_exception_rows(self):
            if not hasattr(self, 'bmf_export_exception_table'):
                return
            rows = sorted({index.row() for index in self.bmf_export_exception_table.selectedIndexes()}, reverse=True)
            if not rows and self.bmf_export_exception_table.currentRow() >= 0:
                rows = [self.bmf_export_exception_table.currentRow()]
            for row_index in rows:
                self.bmf_export_exception_table.removeRow(row_index)
            self._update_bmf_export_exception_summary()

        def _update_bmf_export_exception_summary(self):
            if not hasattr(self, 'bmf_export_exception_summary'):
                return
            rules = self._get_bmf_export_value_exceptions()
            self.bmf_export_value_exceptions = rules
            count = sum(len(column_rules) for column_rules in rules.values())
            if count == 0:
                self.bmf_export_exception_summary.setText('No BMF value exceptions configured.')
                return
            self.bmf_export_exception_summary.setText(
                f'{count} BMF value exception rule(s) configured. Checked rules are applied before regularization; blank replacements are exported as null/default values.'
            )

        def _parse_bmf_numeric_type_error(self, error_text):
            text = str(error_text or '')
            match = re.search(
                r"Column (?P<column>.+?) cannot be exported as (?P<field_type>int|double).*?Invalid values include: (?P<values>\[[^\]]*\])",
                text,
                flags=re.IGNORECASE | re.DOTALL,
            )
            if not match:
                return None
            try:
                values = ast.literal_eval(match.group('values'))
            except Exception:
                values = []
            if not values:
                return None
            column_name = str(match.group('column')).strip().strip('"\'')
            return {
                'column': column_name,
                'field_type': match.group('field_type').lower(),
                'value': str(values[0]),
                'values': [str(value) for value in values],
            }

        def _prompt_bmf_value_exception(self, column_name, bad_value, field_type='double'):
            dialog = QtWidgets.QDialog(self)
            dialog.setWindowTitle('BMF Value Exception')
            layout = QtWidgets.QVBoxLayout(dialog)
            label = QtWidgets.QLabel(
                f"Column '{column_name}' is being exported as {field_type}, but value '{bad_value}' is not numeric."
            )
            label.setWordWrap(True)
            layout.addWidget(label)
            replacement_edit = QtWidgets.QLineEdit()
            replacement_edit.setPlaceholderText('Numeric replacement')
            blank_check = QtWidgets.QCheckBox('Replace with blank/null')
            blank_check.setChecked(True)
            replacement_edit.setEnabled(False)

            def sync_replacement_enabled(checked):
                replacement_edit.setEnabled(not checked)

            blank_check.toggled.connect(sync_replacement_enabled)
            layout.addWidget(blank_check)
            layout.addWidget(replacement_edit)
            buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
            buttons.accepted.connect(dialog.accept)
            buttons.rejected.connect(dialog.reject)
            layout.addWidget(buttons)
            if dialog.exec_() != QtWidgets.QDialog.Accepted:
                return None
            if blank_check.isChecked():
                return ''
            replacement = replacement_edit.text().strip()
            if not replacement:
                return ''
            try:
                float(replacement)
            except ValueError:
                QtWidgets.QMessageBox.warning(self, 'Invalid Replacement', 'Replacement must be numeric or blank.')
                return self._prompt_bmf_value_exception(column_name, bad_value, field_type=field_type)
            return replacement

        def _set_selected_bmf_export_value_columns(self, column_names):
            if not hasattr(self, 'bmf_export_value_cols'):
                return
            desired = {str(name or '').strip() for name in (column_names or []) if str(name or '').strip()}
            blocker = QtCore.QSignalBlocker(self.bmf_export_value_cols)
            try:
                for row_index in range(self.bmf_export_value_cols.rowCount()):
                    item = self.bmf_export_value_cols.item(row_index, 0)
                    if item is not None:
                        item.setCheckState(QtCore.Qt.Checked if item.text() in desired else QtCore.Qt.Unchecked)
            finally:
                del blocker
            self._update_bmf_export_value_summary()

        def _set_all_bmf_export_value_columns_checked(self, checked):
            if not hasattr(self, 'bmf_export_value_cols'):
                return
            blocker = QtCore.QSignalBlocker(self.bmf_export_value_cols)
            try:
                for row_index in range(self.bmf_export_value_cols.rowCount()):
                    item = self.bmf_export_value_cols.item(row_index, 0)
                    if item is not None:
                        item.setCheckState(QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
            finally:
                del blocker
            self._update_bmf_export_value_summary()

        def _get_bmf_export_column_type_overrides(self, include_unselected=False):
            if not hasattr(self, 'bmf_export_value_cols'):
                return {}
            overrides = {}
            for row_index in range(self.bmf_export_value_cols.rowCount()):
                item = self.bmf_export_value_cols.item(row_index, 0)
                combo = self.bmf_export_value_cols.cellWidget(row_index, 1)
                if item is None or combo is None:
                    continue
                if not include_unselected and item.checkState() != QtCore.Qt.Checked:
                    continue
                field_type = normalize_bmf_export_field_type_name(combo.currentData() or combo.currentText())
                inferred_type = normalize_bmf_export_field_type_name(combo.property('bmf_inferred_type') or '')
                is_explicit = bool(combo.property('bmf_type_is_explicit'))
                if field_type and inferred_type == field_type and not is_explicit:
                    continue
                if field_type:
                    overrides[item.text()] = field_type
            return overrides

        def _set_bmf_export_column_type_overrides(self, column_types, include_unselected=True):
            if not hasattr(self, 'bmf_export_value_cols'):
                return
            normalized_types = {
                str(column_name).strip(): normalize_bmf_export_field_type_name(field_type)
                for column_name, field_type in dict(column_types or {}).items()
                if str(column_name).strip()
            }
            for row_index in range(self.bmf_export_value_cols.rowCount()):
                item = self.bmf_export_value_cols.item(row_index, 0)
                combo = self.bmf_export_value_cols.cellWidget(row_index, 1)
                if item is None or combo is None:
                    continue
                if not include_unselected and item.checkState() != QtCore.Qt.Checked:
                    continue
                combo.setProperty('bmf_type_is_explicit', True)
                combo.setCurrentIndex(max(combo.findData(normalized_types.get(item.text(), '')), 0))

        def _update_bmf_export_value_summary(self):
            if not hasattr(self, 'bmf_export_value_summary'):
                return
            total = self.bmf_export_value_cols.rowCount() if hasattr(self, 'bmf_export_value_cols') else 0
            selected = len(self._get_selected_bmf_export_value_columns())
            if total == 0:
                self.bmf_export_value_summary.setText('No export value columns loaded.')
                return
            if selected == 0:
                self.bmf_export_value_summary.setText(
                    f'0 of {total} value columns selected. Select one or more columns to export.'
                )
                return
            self.bmf_export_value_summary.setText(
                f'{selected} of {total} eligible non-coordinate columns selected for export.'
            )

        def _refresh_bmf_export_columns(self):
            if not all(hasattr(self, name) for name in [
                'bmf_export_input_edit',
                'bmf_export_x_col',
                'bmf_export_y_col',
                'bmf_export_z_col',
                'bmf_export_value_cols',
            ]):
                return

            path = self.bmf_export_input_edit.text().strip()
            delimiter = self.bmf_export_delim.currentText() if hasattr(self, 'bmf_export_delim') else ','
            header_line = self.bmf_export_header_line.value() if hasattr(self, 'bmf_export_header_line') else 1
            pending_selection = getattr(self, '_pending_bmf_export_value_cols', None)
            pending_column_types = getattr(self, '_pending_bmf_export_column_types', None)
            current_selection = set(self._get_selected_bmf_export_value_columns())
            current_column_types = self._get_bmf_export_column_type_overrides(include_unselected=True)
            source_path = getattr(self, '_bmf_export_columns_source_path', '')
            path_changed = path != source_path
            if path_changed and pending_selection is None:
                current_selection = set()
                current_column_types = {}
                self._bmf_export_value_cols_initialized = False

            current_x = self.bmf_export_x_col.currentText()
            current_y = self.bmf_export_y_col.currentText()
            current_z = self.bmf_export_z_col.currentText()

            blockers = [
                QtCore.QSignalBlocker(self.bmf_export_x_col),
                QtCore.QSignalBlocker(self.bmf_export_y_col),
                QtCore.QSignalBlocker(self.bmf_export_z_col),
                QtCore.QSignalBlocker(self.bmf_export_value_cols),
            ]
            try:
                for cb in [self.bmf_export_x_col, self.bmf_export_y_col, self.bmf_export_z_col]:
                    cb.clear()
                self.bmf_export_value_cols.setRowCount(0)

                if not path or not os.path.isfile(path):
                    self._bmf_export_columns_source_path = path
                    self._update_bmf_export_value_summary()
                    return

                try:
                    cols = parse_effective_header_line(path, delimiter, header_line)
                except Exception:
                    self._bmf_export_columns_source_path = path
                    self._update_bmf_export_value_summary()
                    return

                for cb in [self.bmf_export_x_col, self.bmf_export_y_col, self.bmf_export_z_col]:
                    cb.addItems(cols)

                def restore_or_suggest(cb, current_value, keywords):
                    if current_value:
                        idx = cb.findText(current_value)
                        if idx >= 0:
                            cb.setCurrentIndex(idx)
                            return
                    for keyword in keywords:
                        for idx in range(cb.count()):
                            if cb.itemText(idx).lower() == keyword:
                                cb.setCurrentIndex(idx)
                                return
                    if cb.count() > 0 and cb.currentIndex() < 0:
                        cb.setCurrentIndex(0)

                restore_or_suggest(self.bmf_export_x_col, current_x, ['x', 'easting'])
                restore_or_suggest(self.bmf_export_y_col, current_y, ['y', 'northing'])
                restore_or_suggest(self.bmf_export_z_col, current_z, ['z', 'elevation', 'rl'])

                excluded = {
                    str(self.bmf_export_x_col.currentText() or '').strip(),
                    str(self.bmf_export_y_col.currentText() or '').strip(),
                    str(self.bmf_export_z_col.currentText() or '').strip(),
                }
                excluded = {value for value in excluded if value}
                eligible_columns = [column_name for column_name in cols if column_name not in excluded]

                desired_selection = set(pending_selection) if pending_selection is not None else current_selection
                desired_column_types = dict(pending_column_types) if pending_column_types is not None else current_column_types
                inferred_column_types = {}
                try:
                    preview_df, _ = load_samples_preview_dataframe(
                        path,
                        samples_delimiter=delimiter,
                        samples_header_line=header_line,
                        max_rows=1000,
                    )
                    inferred_column_types = infer_bmf_export_field_types_from_preview(
                        preview_df,
                        eligible_columns,
                        delimiter=delimiter,
                    )
                except Exception:
                    inferred_column_types = {}
                desired_selection = {name for name in desired_selection if name in eligible_columns}
                if not desired_selection and not getattr(self, '_bmf_export_value_cols_initialized', False):
                    desired_selection = set(eligible_columns)

                type_options = [
                    ('Auto', ''),
                    ('Boolean', 'boolean'),
                    ('Integer', 'int'),
                    ('Double', 'double'),
                    ('String', 'string'),
                ]
                self.bmf_export_value_cols.setRowCount(len(eligible_columns))
                for row_index, column_name in enumerate(eligible_columns):
                    item = QtWidgets.QTableWidgetItem(column_name)
                    item.setFlags(
                        QtCore.Qt.ItemIsEnabled
                        | QtCore.Qt.ItemIsSelectable
                        | QtCore.Qt.ItemIsUserCheckable
                    )
                    item.setCheckState(QtCore.Qt.Checked if column_name in desired_selection else QtCore.Qt.Unchecked)
                    self.bmf_export_value_cols.setItem(row_index, 0, item)

                    combo = QtWidgets.QComboBox()
                    for label, value in type_options:
                        combo.addItem(label, value)
                    has_explicit_type = column_name in desired_column_types
                    inferred_type = inferred_column_types.get(column_name, '')
                    default_type = desired_column_types.get(column_name, inferred_type)
                    combo.setCurrentIndex(max(combo.findData(default_type), 0))
                    combo.setProperty('bmf_inferred_type', inferred_type if not has_explicit_type else '')
                    combo.setProperty('bmf_type_is_explicit', has_explicit_type)
                    combo.activated.connect(lambda _index, cb=combo: cb.setProperty('bmf_type_is_explicit', True))
                    self.bmf_export_value_cols.setCellWidget(row_index, 1, combo)

                self._pending_bmf_export_value_cols = None
                self._pending_bmf_export_column_types = None
                self._bmf_export_value_cols_initialized = True
                self._bmf_export_columns_source_path = path
            finally:
                del blockers

            self._update_bmf_export_value_summary()

        def _get_selected_equation_predictor_columns(self):
            if not hasattr(self, 'equation_predictor_list'):
                return []
            return [item.text() for item in self.equation_predictor_list.selectedItems()]

        def _set_selected_equation_predictor_columns(self, column_names):
            desired = {str(name or '').strip() for name in (column_names or []) if str(name or '').strip()}
            for index in range(self.equation_predictor_list.count()):
                item = self.equation_predictor_list.item(index)
                item.setSelected(item.text() in desired)
            self._update_equation_finder_predictor_summary()

        def _update_equation_finder_predictor_summary(self):
            if not hasattr(self, 'equation_predictor_list'):
                return
            total = self.equation_predictor_list.count()
            selected = len(self.equation_predictor_list.selectedItems())
            if total == 0:
                self.equation_predictor_summary.setText('No numeric predictor columns loaded.')
                return
            self.equation_predictor_summary.setText(
                f'{selected} of {total} available predictor columns selected. '
                'Columns that do not look numeric in the preview are still shown, but the equation search will reject non-numeric selections at runtime.'
            )

        def _refresh_equation_finder_predictor_columns(self):
            if not hasattr(self, 'equation_predictor_list'):
                return

            current_selection = set(self._get_selected_equation_predictor_columns())
            pending_selection = set(getattr(self, '_pending_equation_predictor_selection', []) or [])
            self.equation_predictor_list.clear()

            samples_file = self.samples_edit.text().strip()
            if not samples_file or not os.path.isfile(samples_file):
                self._update_equation_finder_predictor_summary()
                return

            try:
                available_columns = parse_effective_header_line(
                    samples_file,
                    self.samples_delim.currentText(),
                    self.samples_header_line.value(),
                )
            except Exception:
                self._update_equation_finder_predictor_summary()
                return

            preview_df = None
            try:
                preview_df, _ = load_samples_preview_dataframe(
                    samples_file,
                    samples_delimiter=self.samples_delim.currentText(),
                    samples_header_line=self.samples_header_line.value(),
                )
            except Exception:
                preview_df = None

            target_column = self.sample_value_col.currentText().strip()
            domain_column = self.sample_domain_col.currentText().strip()
            coordinate_columns = [
                self.sample_x_col.currentText().strip(),
                self.sample_y_col.currentText().strip(),
                self.sample_z_col.currentText().strip(),
            ]
            available_columns = [str(column_name) for column_name in available_columns]
            available_columns = [col for col in available_columns if col != target_column and col != domain_column]
            if not self.equation_include_coords.isChecked():
                available_columns = [col for col in available_columns if col not in coordinate_columns]

            numeric_columns = infer_numeric_sample_columns(
                preview_df,
                delimiter=self.samples_delim.currentText(),
            ) if preview_df is not None else []
            numeric_column_set = {
                col for col in numeric_columns
                if col != target_column and col != domain_column and (self.equation_include_coords.isChecked() or col not in coordinate_columns)
            }

            for column_name in available_columns:
                item = QtWidgets.QListWidgetItem(column_name)
                if column_name not in numeric_column_set:
                    item.setToolTip('This column did not look numeric in the preview. The equation search may reject it if selected.')
                self.equation_predictor_list.addItem(item)

            desired_selection = pending_selection or current_selection
            if not desired_selection:
                desired_selection = set(numeric_column_set)
                if self.equation_include_coords.isChecked():
                    desired_selection.update({col for col in coordinate_columns if col in numeric_column_set})

            self._set_selected_equation_predictor_columns(desired_selection)
            self._pending_equation_predictor_selection = []

        def _set_operation_buttons_enabled(self, enabled):
            self.start_domaining_btn.setEnabled(enabled)
            self.start_block_value_transfer_btn.setEnabled(enabled)
            self.start_table_attribute_assign_btn.setEnabled(enabled)
            self.start_block_domain_metrics_btn.setEnabled(enabled)
            self.start_domain_interpolation_confidence_btn.setEnabled(enabled)
            self.start_block_volume_weighted_btn.setEnabled(enabled)
            self.start_bmf_export_btn.setEnabled(enabled)
            self.start_equation_finder_btn.setEnabled(enabled)

        def _run_operation_with_progress(self, title, initial_message, operation, kwargs,
                                         success_handler, error_context, failure_handler=None):
            active_thread = getattr(self, '_active_operation_thread', None)
            if active_thread is not None and active_thread.isRunning():
                QtWidgets.QMessageBox.warning(self, 'Operation Running', 'Another long-running operation is already in progress.')
                return

            progress = QtWidgets.QProgressDialog(initial_message, None, 0, 100, self)
            progress.setWindowTitle(title)
            progress.setWindowModality(QtCore.Qt.WindowModal)
            progress.setMinimumDuration(0)
            progress.setAutoClose(False)
            progress.setAutoReset(False)
            progress.setCancelButton(None)
            progress.setValue(0)
            progress.show()

            thread = QtCore.QThread(self)
            worker = BackgroundOperationWorker(operation, kwargs)
            worker.moveToThread(thread)

            self._active_operation_thread = thread
            self._active_operation_worker = worker
            self._active_operation_progress = progress
            self._set_operation_buttons_enabled(False)

            def handle_progress(value, maximum, message):
                progress.setMaximum(max(int(maximum or 0), 1))
                progress.setValue(max(0, min(int(value or 0), progress.maximum())))
                if message:
                    progress.setLabelText(_format_progress_message(message, value, progress.maximum()))

            def cleanup():
                if self._active_operation_progress is progress:
                    progress.close()
                    self._active_operation_progress = None
                self._set_operation_buttons_enabled(True)
                if self._active_operation_worker is worker:
                    self._active_operation_worker = None
                if self._active_operation_thread is thread:
                    self._active_operation_thread = None
                thread.quit()
                thread.wait(2000)
                worker.deleteLater()
                thread.deleteLater()

            def handle_finished(result):
                cleanup()
                success_handler(result)

            def handle_failed(error_text, stack_text):
                print(f"{error_context}: {error_text}")
                print(stack_text)
                cleanup()
                if callable(failure_handler) and failure_handler(error_text, stack_text):
                    return
                QtWidgets.QMessageBox.critical(self, 'Error', f'{error_context}:\n{error_text}')

            worker.progress.connect(handle_progress, QtCore.Qt.QueuedConnection)
            worker.finished.connect(handle_finished, QtCore.Qt.QueuedConnection)
            worker.failed.connect(handle_failed, QtCore.Qt.QueuedConnection)
            thread.started.connect(worker.run)
            thread.start()

        def _update_block_domain_metrics_filters_summary(self):
            count = len(self.block_domain_sample_filters)
            if count == 0:
                self.block_domain_metrics_filters_summary.setText('No sample filters configured. These filters apply app-wide.')
                self.configure_block_domain_metrics_filters_btn.setText('Configure Sample Filters...')
                return

            summaries = [summarize_sample_filter_spec(spec) for spec in self.block_domain_sample_filters]
            self.block_domain_metrics_filters_summary.setText('\n'.join(summaries[:4]))
            self.configure_block_domain_metrics_filters_btn.setText(f'Configure Sample Filters... ({count} active)')

        def _update_block_volume_weighted_filters_summary(self):
            count = len(self.block_volume_weighted_filters)
            if count == 0:
                self.block_volume_weighted_filters_summary.setText('No block filters configured. These filters apply app-wide.')
                self.configure_block_volume_weighted_filters_btn.setText('Configure Block Filters...')
                return

            summaries = [summarize_sample_filter_spec(spec) for spec in self.block_volume_weighted_filters]
            self.block_volume_weighted_filters_summary.setText('\n'.join(summaries[:4]))
            self.configure_block_volume_weighted_filters_btn.setText(f'Configure Block Filters... ({count} active)')

        def configure_block_domain_metrics_filters(self):
            cursor_set = False
            try:
                cfg = self.to_dict()
                sample_filter_source = self._build_filter_data_source(
                    cfg.get('samples_file'),
                    delimiter=cfg.get('samples_delimiter'),
                    header_line=cfg.get('samples_header_line', 1),
                )

                dialog = SampleFiltersDialog(sample_filter_source, filters=self.block_domain_sample_filters, parent=self, subject_label='Sample')
                if dialog.exec_() == QtWidgets.QDialog.Accepted:
                    self.block_domain_sample_filters = dialog.get_filters()
                    self._update_block_domain_metrics_filters_summary()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Could not configure sample filters:\n{exc}')
            finally:
                if cursor_set:
                    QtWidgets.QApplication.restoreOverrideCursor()

        def configure_block_volume_weighted_filters(self):
            cursor_set = False
            try:
                cfg = self.to_dict()
                block_filter_source = self._build_filter_data_source(
                    cfg.get('blocks_file'),
                    delimiter=cfg.get('blocks_delimiter'),
                    header_line=cfg.get('blocks_header_line', 1),
                )

                dialog = SampleFiltersDialog(block_filter_source, filters=self.block_volume_weighted_filters, parent=self, subject_label='Block')
                if dialog.exec_() == QtWidgets.QDialog.Accepted:
                    self.block_volume_weighted_filters = dialog.get_filters()
                    self._update_block_volume_weighted_filters_summary()
                    self._invalidate_domain_catalog_cache()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Could not configure block filters:\n{exc}')
            finally:
                if cursor_set:
                    QtWidgets.QApplication.restoreOverrideCursor()

        def _prompt_to_replace_samples_file(self, samples_file_path, title='Use Generated Samples File'):
            target_path = str(samples_file_path or '').strip()
            if not target_path:
                return False

            response = QtWidgets.QMessageBox.question(
                self,
                title,
                (
                    f"Do you want to set this file as the current Samples File?\n\n"
                    f"{target_path}"
                ),
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.Yes,
            )
            if response == QtWidgets.QMessageBox.Yes:
                self.samples_edit.setText(target_path)
                return True
            return False

        def run_sample_blocks_only(self):
            """Aggregate samples into their containing blocks and export the sampled block subset."""
            try:
                cfg = self.to_dict()
                output_file = resolve_sample_blocks_export_path(
                    cfg.get('sample_blocks_file'),
                    samples_file=cfg.get('samples_file'),
                )
                self.sample_blocks_output_edit.setText(output_file)

                print("=" * 60)
                print("Exporting sample blocks from samples + blocks...")
                print("=" * 60)

                def handle_success(result):
                    print("=" * 60)
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Sample-block export complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Average mode: {'Weighted by ' + result['weight_column'] if result.get('weight_column') else 'Unweighted mean'}\n"
                            f"Sample blocks: {result['sample_block_count']:,}\n"
                            f"Assigned samples: {result['assigned_samples']:,}\n"
                            f"Rejected samples: {result['rejected_samples']:,}\n"
                            f"Invalid coordinate/value rows: {result['invalid_coordinate_or_value_samples']:,}"
                        ),
                    )
                    self._prompt_to_replace_samples_file(result['output_file'], title='Use Sample Blocks As Samples File')

                self._run_operation_with_progress(
                    'Sample Blocks Export',
                    'Preparing sample-block export...',
                    export_sample_blocks_from_samples_and_blocks,
                    {
                        'samples_file': cfg.get('samples_file'),
                        'blocks_file': cfg.get('blocks_file'),
                        'output_file': output_file,
                        'samples_delimiter': cfg.get('samples_delimiter'),
                        'blocks_delimiter': cfg.get('blocks_delimiter'),
                        'samples_header_line': cfg.get('samples_header_line', 1),
                        'blocks_header_line': cfg.get('blocks_header_line', 1),
                        'sample_x_col': cfg.get('sample_x_col'),
                        'sample_y_col': cfg.get('sample_y_col'),
                        'sample_z_col': cfg.get('sample_z_col'),
                        'sample_value_col': cfg.get('sample_value_col'),
                        'sample_domain_col': cfg.get('sample_domain_col'),
                        'sample_weight_col': cfg.get('sample_weight_col'),
                        'include_sample_ids': cfg.get('sample_blocks_include_ids', False),
                        'sample_id_cols': cfg.get('sample_blocks_id_cols'),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_domain_col': cfg.get('block_domain_col'),
                        'block_size': cfg.get('block_size'),
                        'sample_filters': cfg.get('block_domain_sample_filters'),
                        'block_filters': cfg.get('block_volume_weighted_filters'),
                        'blank_sample_domain_behavior': cfg.get('blank_sample_domain_behavior', 'skip'),
                    },
                    handle_success,
                    'An error occurred during sample-block export',
                )
            except Exception as e:
                print(f"Error during sample-block export: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during sample-block export:\n{str(e)}')

        def run_domain_samples_only(self):
            """Assign block domains directly to the samples file and export the result."""
            try:
                cfg = self.to_dict()
                output_file = resolve_domain_samples_export_path(
                    cfg.get('domain_samples_file'),
                    samples_file=cfg.get('samples_file'),
                    domain_col=cfg.get('block_domain_col'),
                )
                self.domain_samples_output_edit.setText(output_file)

                print("=" * 60)
                print("Assigning sample domains from blocks...")
                print("=" * 60)

                def handle_success(result):
                    print("=" * 60)
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Domaining complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Domain column: {result['domain_column']}\n"
                            f"Matched samples: {result['matched_samples']:,}\n"
                            f"Unmatched samples: {result['unmatched_samples']:,}\n"
                            f"Invalid coordinates: {result['invalid_coordinate_samples']:,}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Sample Domaining',
                    'Preparing sample domaining export...',
                    export_domained_samples_from_blocks,
                    {
                        'samples_file': cfg.get('samples_file'),
                        'blocks_file': cfg.get('blocks_file'),
                        'output_file': output_file,
                        'samples_delimiter': cfg.get('samples_delimiter'),
                        'blocks_delimiter': cfg.get('blocks_delimiter'),
                        'samples_header_line': cfg.get('samples_header_line', 1),
                        'blocks_header_line': cfg.get('blocks_header_line', 1),
                        'sample_x_col': cfg.get('sample_x_col'),
                        'sample_y_col': cfg.get('sample_y_col'),
                        'sample_z_col': cfg.get('sample_z_col'),
                        'sample_domain_col': cfg.get('sample_domain_col'),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_domain_col': cfg.get('block_domain_col'),
                        'block_size': cfg.get('block_size'),
                        'sample_filters': cfg.get('block_domain_sample_filters'),
                        'block_filters': cfg.get('block_volume_weighted_filters'),
                    },
                    handle_success,
                    'An error occurred during sample domaining',
                )
            except Exception as e:
                print(f"Error during sample domaining: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during sample domaining:\n{str(e)}')

        def run_block_value_transfer_only(self):
            """Assign selected block columns directly to the samples file and export the result."""
            try:
                cfg = self.to_dict()
                selected_columns = cfg.get('block_value_transfer_cols') or []
                output_file = resolve_block_value_transfer_export_path(
                    cfg.get('block_value_transfer_file'),
                    samples_file=cfg.get('samples_file'),
                    block_value_cols=selected_columns,
                )
                self.block_value_transfer_output_edit.setText(output_file)

                print("=" * 60)
                print("Assigning selected block columns to samples...")
                print("=" * 60)

                def handle_success(result):
                    print("=" * 60)
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Block value transfer complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Transferred columns: {', '.join(result['transferred_columns'])}\n"
                            f"Matched samples: {result['matched_samples']:,}\n"
                            f"Unmatched samples: {result['unmatched_samples']:,}\n"
                            f"Invalid coordinates: {result['invalid_coordinate_samples']:,}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Block Value Transfer',
                    'Preparing block value transfer export...',
                    export_samples_with_block_values_from_blocks,
                    {
                        'samples_file': cfg.get('samples_file'),
                        'blocks_file': cfg.get('blocks_file'),
                        'output_file': output_file,
                        'samples_delimiter': cfg.get('samples_delimiter'),
                        'blocks_delimiter': cfg.get('blocks_delimiter'),
                        'samples_header_line': cfg.get('samples_header_line', 1),
                        'blocks_header_line': cfg.get('blocks_header_line', 1),
                        'sample_x_col': cfg.get('sample_x_col'),
                        'sample_y_col': cfg.get('sample_y_col'),
                        'sample_z_col': cfg.get('sample_z_col'),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_value_cols': selected_columns,
                        'block_size': cfg.get('block_size'),
                    },
                    handle_success,
                    'An error occurred during block value transfer',
                )
            except Exception as e:
                print(f"Error during block value transfer: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during block value transfer:\n{str(e)}')

        def run_block_model_transfer_only(self):
            """Transfer selected source block fields onto the existing target block rows."""
            try:
                cfg = self.to_dict()
                selected_columns = cfg.get('block_model_transfer_cols') or []
                output_file = resolve_block_model_transfer_export_path(
                    cfg.get('block_model_transfer_file'),
                    cfg.get('block_model_target_file'),
                )
                self.block_model_transfer_output_edit.setText(output_file)

                def handle_success(result):
                    max_nearest_distance = result.get('max_nearest_distance')
                    if max_nearest_distance in (None, ''):
                        nearest_distance_display = 'Unlimited'
                    else:
                        nearest_distance_display = f"{max_nearest_distance:g} m"
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Block-model transfer complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Transferred columns: {', '.join(result['transferred_columns'])}\n"
                            f"Target blocks preserved: {result['total_target_blocks']:,}\n"
                            f"Overlap matches: {result['overlap_matched_blocks']:,}\n"
                            f"Nearest fallback matches: {result['nearest_matched_blocks']:,}\n"
                            f"Unmatched blocks: {result['unmatched_blocks']:,}\n"
                            f"Nearest max distance: {nearest_distance_display}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Block Model Transfer',
                    'Preparing block-model transfer...',
                    export_blocks_with_source_block_values,
                    {
                        'source_blocks_file': cfg.get('blocks_file'),
                        'target_blocks_file': cfg.get('block_model_target_file'),
                        'output_file': output_file,
                        'source_delimiter': cfg.get('blocks_delimiter'),
                        'target_delimiter': cfg.get('block_model_target_delimiter'),
                        'source_header_line': cfg.get('blocks_header_line', 1),
                        'target_header_line': cfg.get('block_model_target_header_line', 1),
                        'source_x_col': cfg.get('block_x_col'),
                        'source_y_col': cfg.get('block_y_col'),
                        'source_z_col': cfg.get('block_z_col'),
                        'target_x_col': cfg.get('block_model_target_x_col'),
                        'target_y_col': cfg.get('block_model_target_y_col'),
                        'target_z_col': cfg.get('block_model_target_z_col'),
                        'source_value_cols': selected_columns,
                        'source_block_size': cfg.get('block_size'),
                        'target_block_size': cfg.get('block_model_target_size'),
                        'source_size_cols': cfg.get('block_model_source_size_cols'),
                        'target_size_cols': cfg.get('block_model_target_size_cols'),
                        'nearest_fallback': cfg.get('block_model_transfer_nearest_fallback', True),
                        'max_nearest_distance': cfg.get('block_model_transfer_max_nearest_distance'),
                    },
                    handle_success,
                    'An error occurred during block-model transfer',
                )
            except Exception as e:
                print(f"Error during block-model transfer: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during block-model transfer:\n{str(e)}')

        def run_table_attribute_assignment_only(self):
            """Assign selected table attributes onto a block model by matching one or more key columns."""
            try:
                cfg = self.to_dict()
                output_file = resolve_block_model_table_attribute_export_path(
                    cfg.get('table_attribute_output_file'),
                    block_model_file=cfg.get('table_attribute_block_model_file'),
                    table_file=cfg.get('table_attribute_table_file'),
                )
                self.table_attribute_output_edit.setText(output_file)

                def handle_success(result):
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Table attribute assignment complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Match keys: {', '.join(result['key_columns'])}\n"
                            f"Assigned columns: {', '.join(result['assigned_columns'])}\n"
                            f"Matched rows: {result['matched_rows']:,}\n"
                            f"Unmatched rows: {result['unmatched_rows']:,}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Assign Attributes From Table',
                    'Preparing table attribute assignment...',
                    export_block_model_with_table_attributes,
                    {
                        'block_model_file': cfg.get('table_attribute_block_model_file'),
                        'table_file': cfg.get('table_attribute_table_file'),
                        'output_file': output_file,
                        'block_model_delimiter': cfg.get('table_attribute_block_model_delimiter'),
                        'block_model_header_line': cfg.get('table_attribute_block_model_header_line', 1),
                        'table_delimiter': cfg.get('table_attribute_table_delimiter'),
                        'table_header_line': cfg.get('table_attribute_table_header_line', 1),
                        'key_columns': cfg.get('table_attribute_key_cols'),
                        'table_value_cols': cfg.get('table_attribute_value_cols'),
                    },
                    handle_success,
                    'An error occurred during table attribute assignment',
                )
            except Exception as e:
                print(f"Error during table attribute assignment: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during table attribute assignment:\n{str(e)}')

        def run_block_domain_sample_metrics_only(self):
            """Export block rows with distance statistics to filtered samples inside the same domain."""
            try:
                cfg = self.to_dict()
                output_file = resolve_block_domain_metrics_export_path(
                    cfg.get('block_domain_metrics_file'),
                    blocks_file=cfg.get('blocks_file'),
                    domain_col=cfg.get('block_domain_col'),
                )
                self.block_domain_metrics_output_edit.setText(output_file)

                print("=" * 60)
                print("Calculating block-domain sample metrics...")
                print("=" * 60)

                def handle_success(result):
                    filters_text = 'None'
                    if result['filters_applied']:
                        filters_text = '\n'.join(entry['summary'] for entry in result['filters_applied'])
                    selected_metric_text = ', '.join(
                        BLOCK_DOMAIN_METRIC_DEFINITIONS_BY_ID[metric_id]['label']
                        for metric_id in result.get('selected_metrics', [])
                        if metric_id in BLOCK_DOMAIN_METRIC_DEFINITIONS_BY_ID
                    ) or 'None'
                    closest_id_text = 'None'
                    if result.get('closest_sample_id_source_columns'):
                        closest_id_text = ', '.join(result['closest_sample_id_source_columns'])
                    nearest_sample_value_metrics_text = 'Disabled'
                    nearest_sample_value_columns = [
                        result.get('nearest_sample_value_column'),
                        result.get('nearest_sample_residual_column'),
                        result.get('nearest_sample_abs_residual_column'),
                        result.get('nearest_sample_group_block_count_column'),
                        result.get('nearest_sample_group_mean_residual_column'),
                        result.get('nearest_sample_group_rms_residual_column'),
                        result.get('nearest_sample_group_std_residual_column'),
                    ]
                    if any(nearest_sample_value_columns):
                        nearest_sample_value_metrics_text = ', '.join(
                            column_name for column_name in nearest_sample_value_columns if column_name
                        )
                    prefix_text = 'Enabled' if result.get('use_block_value_prefix', False) else 'Disabled'
                    distance_bands_text = 'Disabled'
                    if result.get('distance_summary_thresholds'):
                        threshold_tokens = ', '.join(
                            _format_metric_distance_token(value)
                            for value in result['distance_summary_thresholds']
                        )
                        distance_bands_text = (
                            f"step={result['distance_count_step']:,g}, max factor={result['distance_count_max_factor']:,}; "
                            f"summary thresholds={threshold_tokens}"
                        )
                    avg_distance_text = 'Not exported'
                    if result.get('average_distance_column'):
                        avg_distance_text = f"Exact all-sample average -> {result['average_distance_column']}"
                    if result.get('average_distance_knn_column'):
                        knn_text = (
                            f"KNN average (k={result.get('average_distance_knn_k')}) -> {result['average_distance_knn_column']}"
                        )
                        avg_distance_text = f"{avg_distance_text}; {knn_text}" if avg_distance_text != 'Not exported' else knn_text
                    summary_text = 'Not generated'
                    if result.get('summary_output_file'):
                        summary_text = (
                            f"{result['summary_output_file']} "
                            f"({result.get('summary_row_count', 0):,} rows)"
                        )

                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Block metrics export complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Summary saved to:\n{summary_text}\n\n"
                            f"Domain column: {result['domain_column']}\n"
                            f"Filtered samples: {result['filtered_samples']:,} of {result['input_samples']:,}\n"
                            f"Matched samples: {result['matched_samples']:,}\n"
                            f"Processed blocks: {result['processed_blocks']:,}\n"
                            f"Blocks with samples in domain: {result['blocks_with_samples_in_domain']:,}\n"
                            f"Invalid sample coordinates: {result['invalid_coordinate_samples']:,}\n"
                            f"Selected metrics: {selected_metric_text}\n"
                            f"Average-distance outputs: {avg_distance_text}\n"
                            f"Closest sample ID columns: {closest_id_text}\n"
                            f"Block-value column prefix: {prefix_text}\n"
                            f"Nearest-sample value metrics: {nearest_sample_value_metrics_text}\n"
                            f"Distance bands: {distance_bands_text}\n\n"
                            f"Filters:\n{filters_text}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Block Metrics Export',
                    'Preparing block-domain sample metrics export...',
                    export_block_domain_sample_metrics,
                    {
                        'samples_file': cfg.get('samples_file'),
                        'blocks_file': cfg.get('blocks_file'),
                        'output_file': output_file,
                        'samples_delimiter': cfg.get('samples_delimiter'),
                        'blocks_delimiter': cfg.get('blocks_delimiter'),
                        'samples_header_line': cfg.get('samples_header_line', 1),
                        'blocks_header_line': cfg.get('blocks_header_line', 1),
                        'sample_x_col': cfg.get('sample_x_col'),
                        'sample_y_col': cfg.get('sample_y_col'),
                        'sample_z_col': cfg.get('sample_z_col'),
                        'sample_domain_col': cfg.get('sample_domain_col'),
                        'sample_value_col': cfg.get('sample_value_col'),
                        'selected_metrics': cfg.get('block_domain_metrics_selected_metrics'),
                        'average_distance_knn_k': cfg.get('block_domain_metrics_avg_distance_knn_k', 8),
                        'closest_sample_id_cols': cfg.get('block_domain_metrics_closest_sample_id_cols'),
                        'distance_count_step': cfg.get('block_domain_metrics_distance_step') if cfg.get('block_domain_metrics_distance_counts_enabled') else None,
                        'distance_count_max_factor': cfg.get('block_domain_metrics_max_factor') if cfg.get('block_domain_metrics_distance_counts_enabled') else None,
                        'use_block_value_prefix': cfg.get('block_domain_metrics_use_block_value_prefix', True),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_domain_col': cfg.get('block_domain_col'),
                        'block_size': cfg.get('block_size'),
                        'block_value_col': cfg.get('block_domain_metrics_value_col'),
                        'sample_filters': cfg.get('block_domain_sample_filters'),
                        'block_filters': cfg.get('block_volume_weighted_filters'),
                        'blank_sample_domain_behavior': cfg.get('blank_sample_domain_behavior', 'skip'),
                    },
                    handle_success,
                    'An error occurred during block domain metrics export',
                )
            except Exception as e:
                print(f"Error during block domain metrics export: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during block domain metrics export:\n{str(e)}')
        def run_domain_interpolation_confidence_only(self):
            """Export one row per domain with source-sample, block-to-sample, and block-to-block spacing metrics."""
            try:
                cfg = self.to_dict()
                output_file = resolve_domain_interpolation_confidence_export_path(
                    cfg.get('domain_interpolation_confidence_file'),
                    blocks_file=cfg.get('blocks_file'),
                    domain_col=cfg.get('block_domain_col'),
                )
                self.domain_interpolation_confidence_output_edit.setText(output_file)

                print("=" * 60)
                print("Calculating domain interpolation confidence metrics...")
                print("=" * 60)

                def handle_success(result):
                    filters_text = 'None'
                    if result['filters_applied']:
                        filters_text = '\n'.join(entry['summary'] for entry in result['filters_applied'])

                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Domain confidence export complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Domain column: {result['domain_column']}\n"
                            f"Domains summarized: {result['domain_count']:,}\n"
                            f"Filtered samples: {result['filtered_samples']:,} of {result['input_samples']:,}\n"
                            f"Matched samples: {result['matched_samples']:,}\n"
                            f"Processed blocks: {result['processed_blocks']:,}\n"
                            f"Invalid sample coordinates: {result['invalid_coordinate_samples']:,}\n\n"
                            f"Ratio: Avg_Source_Sample_Distance / Avg_Block_To_Source_Sample_Distance\n"
                            f"Includes axis columns for sample and block-to-sample distances: *_X, *_Y, *_Z\n\n"
                            f"Filters:\n{filters_text}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Domain Interpolation Confidence Export',
                    'Preparing domain interpolation confidence export...',
                    export_domain_interpolation_confidence_metrics,
                    {
                        'samples_file': cfg.get('samples_file'),
                        'blocks_file': cfg.get('blocks_file'),
                        'output_file': output_file,
                        'samples_delimiter': cfg.get('samples_delimiter'),
                        'blocks_delimiter': cfg.get('blocks_delimiter'),
                        'samples_header_line': cfg.get('samples_header_line', 1),
                        'blocks_header_line': cfg.get('blocks_header_line', 1),
                        'sample_x_col': cfg.get('sample_x_col'),
                        'sample_y_col': cfg.get('sample_y_col'),
                        'sample_z_col': cfg.get('sample_z_col'),
                        'sample_domain_col': cfg.get('sample_domain_col'),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_domain_col': cfg.get('block_domain_col'),
                        'block_size': cfg.get('block_size'),
                        'sample_filters': cfg.get('block_domain_sample_filters'),
                        'block_filters': cfg.get('block_volume_weighted_filters'),
                        'blank_sample_domain_behavior': cfg.get('blank_sample_domain_behavior', 'skip'),
                    },
                    handle_success,
                    'An error occurred during domain interpolation confidence export',
                )
            except Exception as e:
                print(f"Error during domain interpolation confidence export: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during domain interpolation confidence export:\n{str(e)}')

        def run_block_volume_weighted_average_only(self):
            """Export inferred per-row block volumes and compute a volume-weighted average for a selected blocks column."""
            try:
                cfg = self.to_dict()
                output_file = resolve_block_volume_weighted_average_export_path(
                    cfg.get('block_volume_weighted_file'),
                    blocks_file=cfg.get('blocks_file'),
                    value_col=cfg.get('block_volume_weighted_value_col'),
                )
                self.block_volume_weighted_output_edit.setText(output_file)

                print("=" * 60)
                print("Calculating block volume weighted average...")
                print("=" * 60)

                def handle_success(result):
                    weighted_average = result['weighted_average']
                    weight_label = result.get('weight_column') or 'Volume'
                    weighted_average_text = 'NaN' if pd.isna(weighted_average) else f"{weighted_average:.12g}"
                    filters_text = 'None'
                    if result['filters_applied']:
                        filters_text = '\n'.join(entry['summary'] for entry in result['filters_applied'])
                    domain_lines = ''
                    if result.get('domain_column') and result.get('domain_summaries'):
                        preview = []
                        for domain, summary in sorted(result['domain_summaries'].items()):
                            avg_text = 'NaN' if pd.isna(summary['weighted_average']) else f"{summary['weighted_average']:.12g}"
                            detail = f"{domain}: avg={avg_text}, weight={summary['total_weight']:.12g}, rows={summary['rows_with_numeric_value']:,}"
                            if not pd.isna(summary['total_volume']):
                                detail += f", vol={summary['total_volume']:.12g}"
                            preview.append(detail)
                        domain_lines = (
                            f"\nDomain column: {result['domain_column']}\n"
                            f"Domain summaries:\n" + '\n'.join(preview[:8])
                        )
                        if len(preview) > 8:
                            domain_lines += f"\n... and {len(preview) - 8} more domains"
                    total_volume_line = ''
                    if not pd.isna(result['total_volume']):
                        total_volume_line = f"Total volume: {result['total_volume']:.12g}\n"
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Block volume export complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Weighted column: {result['value_column']}\n"
                            f"Weight column: {weight_label}\n"
                            f"Weighted average: {weighted_average_text}\n"
                            f"Total weight: {result['total_weight']:.12g}\n"
                            f"{total_volume_line}"
                            f"Weighted sum: {result['weighted_sum']:.12g}\n"
                            f"Filtered blocks: {result['filtered_blocks']:,} of {result['input_blocks']:,}\n"
                            f"Processed rows: {result['processed_rows']:,}\n"
                            f"Rows with numeric value: {result['rows_with_numeric_value']:,}\n"
                            f"Invalid coordinates: {result['invalid_coordinate_rows']:,}\n"
                            f"Invalid values: {result['invalid_value_rows']:,}\n"
                            f"Invalid weights: {result['invalid_weight_rows']:,}\n\n"
                            f"Filters:\n{filters_text}"
                            f"{domain_lines}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Block Volume Weighted Average',
                    'Preparing block volume export...',
                    export_block_volume_weighted_average,
                    {
                        'blocks_file': cfg.get('blocks_file'),
                        'value_col': cfg.get('block_volume_weighted_value_col'),
                        'output_file': output_file,
                        'blocks_delimiter': cfg.get('blocks_delimiter'),
                        'blocks_header_line': cfg.get('blocks_header_line', 1),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_domain_col': cfg.get('block_domain_col'),
                        'weight_col': cfg.get('block_weight_metric_col'),
                        'block_filters': cfg.get('block_volume_weighted_filters'),
                        'block_size': cfg.get('block_size'),
                    },
                    handle_success,
                    'An error occurred during block volume export',
                )
            except Exception as e:
                print(f"Error during block volume export: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during block volume export:\n{str(e)}')

        def run_csv_to_bmf_export_only(self):
            """Export a CSV grid to BMF using the standalone reverse-engineered backend."""
            try:
                cfg = self.to_dict()
                input_csv = str(cfg.get('bmf_export_input_file') or '').strip()
                output_bmf = str(cfg.get('bmf_export_output_file') or '').strip()
                if not input_csv or not os.path.isfile(input_csv):
                    raise ValueError('Please select a valid input CSV file.')
                if not output_bmf:
                    raise ValueError('Please select an output BMF file.')
                if not output_bmf.lower().endswith('.bmf'):
                    output_bmf = f"{output_bmf}.bmf"
                    self.bmf_export_output_edit.setText(output_bmf)

                value_cols = [
                    str(token).strip()
                    for token in (cfg.get('bmf_export_value_cols') or [])
                    if str(token).strip()
                ]
                if not value_cols:
                    raise ValueError('Please select at least one value column to export.')
                column_types = self._get_bmf_export_column_type_overrides(include_unselected=False)
                unknown_typed_columns = [column_name for column_name in column_types if column_name not in value_cols]
                if unknown_typed_columns:
                    raise ValueError(
                        'Column type overrides must refer only to selected value columns. '
                        f'Unknown or unselected columns: {unknown_typed_columns}'
                    )
                regularize_to_base_block = bool(cfg.get('bmf_export_regularize_to_base_block', False))
                bmf_cell_size = cfg.get('bmf_export_cell_size')
                if regularize_to_base_block and not bmf_cell_size:
                    block_size = cfg.get('block_size')
                    try:
                        block_size_values = [float(value) for value in block_size]
                    except Exception:
                        block_size_values = []
                    if len(block_size_values) == 3 and all(value > 0 for value in block_size_values):
                        bmf_cell_size = block_size_values
                    else:
                        raise ValueError(
                            'Regularize to base block size requires Cell Size X/Y/Z values or valid configured block dimensions.'
                        )
                print("=" * 60)
                print("Exporting CSV grid to BMF...")
                print("=" * 60)

                def handle_success(result):
                    summary = result.get('summary', {}) or {}
                    grid = summary.get('grid', {}) or {}
                    value_text = ', '.join(summary.get('value_columns') or []) or 'Auto-detected non-coordinate columns'
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"BMF export complete!\nResults saved to:\n{summary.get('output_bmf', output_bmf)}\n\n"
                            f"Backend: {summary.get('backend', cfg.get('bmf_export_backend'))}\n"
                            f"Rows exported: {summary.get('rows', 0):,}\n"
                            f"Grid dimensions: {tuple(grid.get('dimensions', []))}\n"
                            f"Cell size: {tuple(grid.get('cell_size', []))}\n"
                            f"Value columns: {value_text}"
                        ),
                    )

                def make_export_kwargs():
                    return {
                        'input_csv': input_csv,
                        'output_bmf': output_bmf,
                        'x_col': str(cfg.get('bmf_export_x_col') or 'x').strip() or 'x',
                        'y_col': str(cfg.get('bmf_export_y_col') or 'y').strip() or 'y',
                        'z_col': str(cfg.get('bmf_export_z_col') or 'z').strip() or 'z',
                        'value_cols': value_cols,
                        'column_types': column_types,
                        'value_exceptions': self._get_bmf_export_value_exceptions(),
                        'backend': cfg.get('bmf_export_backend') or 'tbms-config-text',
                        'delimiter': cfg.get('bmf_export_delimiter'),
                        'header_line': cfg.get('bmf_export_header_line', 1),
                        'cell_size': bmf_cell_size,
                        'regularize_to_base_block': regularize_to_base_block,
                    }

                def handle_failure(error_text, _stack_text):
                    parsed = self._parse_bmf_numeric_type_error(error_text)
                    if not parsed:
                        return False
                    column_name = parsed['column']
                    bad_value = parsed['value']
                    replacement = self._prompt_bmf_value_exception(
                        column_name,
                        bad_value,
                        field_type=parsed.get('field_type', 'double'),
                    )
                    if replacement is None:
                        return False
                    current_rules = self._get_bmf_export_value_exceptions()
                    current_rules.setdefault(column_name, {})[bad_value] = replacement
                    self._set_bmf_export_value_exceptions(current_rules)
                    print(
                        f"Added BMF value exception: column={column_name!r}, value={bad_value!r}, "
                        f"replacement={replacement!r}. Retrying export..."
                    )
                    QtCore.QTimer.singleShot(0, start_export)
                    return True

                def start_export():
                    self._run_operation_with_progress(
                        'CSV To BMF Export',
                        'Preparing BMF export...',
                        export_csv_grid_to_bmf,
                        make_export_kwargs(),
                        handle_success,
                        'An error occurred during CSV to BMF export',
                        failure_handler=handle_failure,
                    )

                start_export()
            except Exception as e:
                print(f"Error during CSV to BMF export: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during CSV to BMF export:\n{str(e)}')

        def run_equation_finder_only(self):
            """Run domain-wise symbolic regression on the configured samples file."""
            try:
                cfg = self.to_dict()
                output_file = resolve_equation_finder_export_path(
                    cfg.get('equation_finder_file'),
                    samples_file=cfg.get('samples_file'),
                    value_col=cfg.get('sample_value_col'),
                    domain_col=cfg.get('sample_domain_col'),
                )
                self.equation_finder_output_edit.setText(output_file)

                predictor_columns = self._get_selected_equation_predictor_columns()
                if not predictor_columns:
                    raise ValueError('Please select at least one numeric predictor column.')

                print("=" * 60)
                print("Finding equations by domain...")
                print("=" * 60)

                def handle_success(result):
                    summary_df = result.get('summary_dataframe')
                    while summary_df is not None:
                        try:
                            summary_df.to_csv(result['output_file'], index=False)
                            break
                        except PermissionError as exc:
                            message_box = QtWidgets.QMessageBox(self)
                            message_box.setIcon(QtWidgets.QMessageBox.Warning)
                            message_box.setWindowTitle('Equation Summary Save Failed')
                            message_box.setText('Could not save the equation summary CSV.')
                            message_box.setInformativeText(
                                f"Permission denied for:\n{result['output_file']}\n\n"
                                'The file may be open in Excel or locked by OneDrive sync. '
                                'Close the file and click Retry to save again, or click Cancel to leave only the per-domain detail files.'
                            )
                            message_box.setDetailedText(str(exc))
                            message_box.setStandardButtons(QtWidgets.QMessageBox.Retry | QtWidgets.QMessageBox.Cancel)
                            message_box.setDefaultButton(QtWidgets.QMessageBox.Retry)
                            if message_box.exec_() == QtWidgets.QMessageBox.Retry:
                                continue
                            QtWidgets.QMessageBox.warning(
                                self,
                                'Summary Not Saved',
                                (
                                    'Equation search completed, but the summary CSV was not saved.\n\n'
                                    f"Per-domain detail files are available in:\n{result['details_directory']}"
                                ),
                            )
                            return

                    predictor_text = ', '.join(result['predictor_columns']) if result['predictor_columns'] else 'None'
                    skipped_predictor_text = ', '.join(result.get('skipped_predictor_columns') or []) or 'None'
                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Equation finder complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Details directory:\n{result['details_directory']}\n\n"
                            f"Target column: {result['target_column']}\n"
                            f"Domain column: {result['domain_column']}\n"
                            f"Predictors: {predictor_text}\n"
                            f"Skipped predictors: {skipped_predictor_text}\n"
                            f"Input rows: {result['input_rows']:,}\n"
                            f"Valid rows: {result['valid_rows']:,}\n"
                            f"Domains found: {result['domain_count']:,}\n"
                            f"Processed domains: {result['processed_domains']:,}\n"
                            f"Skipped domains: {result['skipped_domains']:,}"
                        ),
                    )

                self._run_operation_with_progress(
                    'Equation Finder',
                    'Preparing equation search...',
                    export_domain_symbolic_regression_equations,
                    {
                        'samples_file': cfg.get('samples_file'),
                        'output_file': output_file,
                        'samples_delimiter': cfg.get('samples_delimiter'),
                        'samples_header_line': cfg.get('samples_header_line', 1),
                        'sample_value_col': cfg.get('sample_value_col'),
                        'sample_domain_col': cfg.get('sample_domain_col'),
                        'predictor_cols': predictor_columns,
                        'sample_filters': get_configured_sample_filters(cfg),
                        'min_samples_per_domain': cfg.get('equation_finder_min_samples_per_domain', 25),
                        'validation_fraction': cfg.get('equation_finder_validation_fraction', 0.2),
                        'max_iterations': cfg.get('equation_finder_max_iterations', 100),
                        'timeout_seconds': cfg.get('equation_finder_timeout_seconds', 60),
                    },
                    handle_success,
                    'An error occurred during equation search',
                )
            except Exception as e:
                print(f"Error during equation search: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during equation search:\n{str(e)}')

        def run_interpolation_only(self):
            """Run interpolation without visualization"""
            try:
                cfg = self.to_dict(include_runtime_state=True)
                interpolation_file = resolve_interpolation_csv_export_path(cfg['interpolation_file'])
                cfg['interpolation_file'] = interpolation_file
                block_evaluated_samples_file = resolve_block_evaluated_samples_export_path(
                    cfg.get('export_block_evaluated_samples', False),
                    cfg.get('block_evaluated_samples_file'),
                    interpolation_file=interpolation_file,
                    samples_file=cfg.get('samples_file'),
                )
                
                # Run interpolation directly without visualization
                print("=" * 60)
                print("Running interpolation without visualization...")
                print("=" * 60)
                
                # Load samples first
                samples_file = cfg['samples_file']
                samples_delimiter = cfg.get('samples_delimiter')
                samples_header_line = cfg.get('samples_header_line', 1)
                
                print(f"Loading sample file from {samples_file}...")
                sample_x_col = cfg.get('sample_x_col')
                sample_y_col = cfg.get('sample_y_col')
                sample_z_col = cfg.get('sample_z_col')
                sample_value_col = cfg.get('sample_value_col')
                sample_domain_col = cfg.get('sample_domain_col')
                wants_domain = bool(
                    cfg.get('algorithm') in ('string_theory', 'net_connector')
                    and str(cfg.get('string_theory_params', {}).get('interpolate_target', 'value')).strip().lower() == 'domain'
                )
                wants_ant_domain = bool(
                    cfg.get('algorithm') == 'ant_colony'
                    and str(cfg.get('ant_colony_interpolate_target', 'value')).strip().lower() == 'domain'
                )
                wants_domain_any = wants_domain or wants_ant_domain
                sample_filters = get_configured_sample_filters(cfg)
                block_filters = get_configured_block_filters(cfg)
                needs_sample_domains = should_resolve_sample_domains_for_interpolation(
                    wants_domain_any,
                    blocks_file=cfg.get('blocks_file'),
                    block_domain_col=cfg.get('block_domain_col'),
                )

                df, _, explicit_sample_map = load_samples_dataframe(
                    samples_file,
                    samples_delimiter=samples_delimiter,
                    samples_header_line=samples_header_line,
                    sample_x_col=sample_x_col,
                    sample_y_col=sample_y_col,
                    sample_z_col=sample_z_col,
                    sample_value_col=sample_value_col,
                    sample_domain_col=sample_domain_col,
                    sample_filters=sample_filters,
                    progress_label='Reading sample file',
                    extra_columns=[cfg.get('sample_weight_col')],
                )

                if needs_sample_domains and explicit_sample_map:
                    if samples_header_line and samples_header_line != 1 and samples_delimiter:
                        df, _ = read_csv_with_selected_header(
                            samples_file,
                            samples_delimiter,
                            samples_header_line,
                            expected_min_cols=4,
                            progress_label='Reading sample file',
                        )
                    else:
                        df = read_autodetect_csv(
                            samples_file,
                            forced_delimiter=samples_delimiter,
                            progress_label='Reading sample file',
                        )
                    if sample_filters:
                        df, _ = apply_sample_filters(df, sample_filters=sample_filters)
                    explicit_sample_map = None

                if needs_sample_domains:
                    df = normalize_selected_sample_domain_column(df, sample_domain_col=sample_domain_col)
                df = normalize_selected_sample_weight_column(
                    df,
                    sample_weight_col=cfg.get('sample_weight_col'),
                    sample_value_col=sample_value_col,
                )

                if explicit_sample_map:
                    pass
                elif sample_x_col and sample_y_col and sample_z_col and sample_value_col:
                    rename_map = {sample_x_col: 'x', sample_y_col: 'y', sample_z_col: 'z', sample_value_col: 'Value'}
                    df = df.rename(columns=rename_map)
                
                df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
                nan_before = int(df['Value'].isna().sum())

                if needs_sample_domains:
                    df, domain_resolution = ensure_sample_domains_for_domain_operations(
                        df,
                        sample_domain_col=sample_domain_col,
                        blank_domain_behavior=cfg.get('blank_sample_domain_behavior', 'skip'),
                        x_col='x',
                        y_col='y',
                        z_col='z',
                        blocks_file=cfg.get('blocks_file'),
                        blocks_delimiter=cfg.get('blocks_delimiter'),
                        blocks_header_line=cfg.get('blocks_header_line', 1),
                        block_x_col=cfg.get('block_x_col'),
                        block_y_col=cfg.get('block_y_col'),
                        block_z_col=cfg.get('block_z_col'),
                        block_domain_col=cfg.get('block_domain_col'),
                        block_size=cfg.get('block_size'),
                        block_filters=block_filters,
                    )
                    df['Domain'] = df['Domain'].astype(str).str.strip()
                    blank_domain = df['Domain'].isna() | (df['Domain'].str.strip() == '') | (df['Domain'].str.lower() == 'nan')
                    blank_count = int(blank_domain.sum())
                    if blank_count:
                        print(f"Detected {blank_count} sample rows with blank Domain; these will be excluded from domain interpolation.")
                        df = df.loc[~blank_domain].copy()

                    if nan_before:
                        print(f"Detected {nan_before} sample rows with blank/non-numeric Value; filling with 0.0 for domain interpolation.")
                        df['Value'] = df['Value'].fillna(0.0)
                else:
                    df = df.dropna(subset=['Value'])

                sample_weights = None
                selected_weight_column = str(cfg.get('sample_weight_col') or '').strip()
                if selected_weight_column and selected_weight_column != '(None)':
                    df['Weight'] = pd.to_numeric(df['Weight'], errors='coerce')
                    invalid_weight_mask = (~np.isfinite(df['Weight'])) | (df['Weight'] <= 0.0)
                    invalid_weight_count = int(invalid_weight_mask.sum())
                    if invalid_weight_count:
                        print(f"Detected {invalid_weight_count} sample rows with blank/non-numeric/non-positive weights; these will be excluded from samples.")
                        df = df.loc[~invalid_weight_mask].copy()
                    if len(df) == 0:
                        raise ValueError("After removing invalid sample weights, no samples remain. Check the configured 'sample_weight_col'.")
                    sample_weights = df['Weight'].to_numpy(dtype=float, copy=False)
                
                points = df[['x','y','z']].values
                values = df['Value'].values
                sample_domains = df['Domain'].values if 'Domain' in df.columns else None
                print(f"Loaded {len(values)} samples from {samples_file}.")
                
                # Create blocks with samples
                blocks = create_blocks(
                    points,
                    values,
                    block_size=cfg['block_size'],
                    verbose=cfg['verbose'],
                    range_size=cfg['range_size'],
                    max_pheromone=cfg['max_pheromone'],
                    ants_per_sample=cfg['ants_per_sample'],
                    blocks_file=cfg['blocks_file'],
                    background_value=cfg['background_value'],
                    background_distance=cfg['background_distance'],
                    average_with_blocks=cfg['average_with_blocks'],
                    blocks_delimiter=cfg.get('blocks_delimiter'),
                    avoid_visited_threshold_enabled=cfg.get('avoid_visited_threshold_enabled', False),
                    avoid_visited_threshold=cfg.get('avoid_visited_threshold', 100),
                    blocks_header_line=cfg.get('blocks_header_line', 1),
                    block_x_col=cfg.get('block_x_col'),
                    block_y_col=cfg.get('block_y_col'),
                    block_z_col=cfg.get('block_z_col'),
                    block_domain_col=cfg.get('block_domain_col'),
                    config=cfg,
                    sample_domains=sample_domains,
                    sample_weights=sample_weights,
                    build_visual_blocks=False,
                )
                
                # Run interpolation (handle both single and sequential domain processing)
                dims = tuple(blocks._block_info['dims'])
                iterations = cfg['iterations']
                
                # Check if we have multiple interpolators (sequential domain processing)
                if hasattr(blocks, '_interpolators'):
                    print(f"Running sequential domain interpolation for {len(blocks._interpolators)} domains...")
                    
                    for domain_idx, (domain, interpolator_list) in enumerate(blocks._interpolators.items(), 1):
                        # interpolator_list is [Pass1, (optional) Pass2]
                        
                        # --- Pass 1 ---
                        interp1 = interpolator_list[0]
                        algo_name1 = interp1.get_algorithm_name()
                        print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} - Pass 1 ({algo_name1}) ===")
                        seed_original_sample_provenance(interp1)
                        pass1_snapshot = snapshot_interpolator_state(interp1)

                        _run_interpolator_with_progress(
                            interp1,
                            dims,
                            iterations,
                            f"Domain {domain} - Pass 1",
                        )
                        finalize_phase_provenance(interp1, 'First Pass', algo_name1, pass1_snapshot)
                        
                        # Generate stats for Pass 1 if String Theory
                        output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                        _run_interpolator_statistics_with_retry(interp1, output_dir, f"{domain}_Pass1", parent=self)

                        last_interp = interp1

                        # --- Pass 2 ---
                        if len(interpolator_list) > 1:
                            interp2 = interpolator_list[1]
                            algo_name2 = interp2.get_algorithm_name()
                            print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} - Pass 2 ({algo_name2}) ===")
                            
                            print("  Transferring data from Pass 1 to Pass 2...")
                            pass1_values = interp1.get_interpolated_values()
                            print(f"  Pass 1 generated {len(pass1_values)} blocks.")
                            
                            # Re-initialize Pass 2 with merged samples
                            min_bounds = blocks._block_info['min_bounds']
                            block_size = blocks._block_info['block_size']
                            pass1_domain_mapping = {pos: domain for pos in pass1_values}
                            
                            # Determine if we should enforce domain mapping/grid restriction
                            use_mapping = False
                            if hasattr(interp2, 'allowed_grid_override') and interp2.allowed_grid_override is not None:
                                use_mapping = True
                            
                            interp2.initialize_blocks(
                                pass1_values,
                                dims,
                                min_bounds,
                                block_size,
                                use_domain_mapping=use_mapping,
                                sample_domain_mapping=pass1_domain_mapping,
                            )
                            copy_interpolator_provenance(interp1, interp2)
                            pass2_snapshot = snapshot_interpolator_state(interp2)
                            
                            if hasattr(interp2, 'create_ants'):
                                interp2.create_ants()

                            _run_interpolator_with_progress(
                                interp2,
                                dims,
                                iterations,
                                f"Domain {domain} - Pass 2",
                            )
                            finalize_phase_provenance(interp2, 'Second Pass', algo_name2, pass2_snapshot)
                            
                            # Generate stats for Pass 2 if String Theory
                            output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                            _run_interpolator_statistics_with_retry(interp2, output_dir, f"{domain}_Pass2", parent=self)
                            last_interp = interp2

                        post_process_mode = _get_domain_post_process_mode(cfg, domain)
                        post_process_snapshot = snapshot_interpolator_state(last_interp)
                        created, assigned = _run_domain_post_process(last_interp, dims, post_process_mode, domain_label=domain)
                        if created or assigned:
                            finalize_phase_provenance(last_interp, 'Post-process', 'Fill with Average', post_process_snapshot)
                            print(f"Applied post-process for {domain}: created={created}, assigned={assigned}")
            
                        # Print domain summary (of the last pass)
                        metadata = last_interp.get_metadata()
                        print(f"\n=== Domain {domain} Summary ===")
                        for key, value in metadata.items():
                            print(f"{key}: {value}")
                
                else:
                    # Single interpolator case
                    interpolator = blocks._ant_colony
                    algo_name = interpolator.get_algorithm_name()
                    seed_original_sample_provenance(interpolator)
                    first_pass_snapshot = snapshot_interpolator_state(interpolator)

                    _run_interpolator_with_progress(
                        interpolator,
                        dims,
                        iterations,
                        f"Interpolation ({algo_name})",
                    )
                    finalize_phase_provenance(interpolator, 'First Pass', algo_name, first_pass_snapshot)
                    
                    metadata = interpolator.get_metadata()
                    print(f"\n=== Summary ===")
                    for key, value in metadata.items():
                        print(f"{key}: {value}")
                        
                    # Generate stats if String Theory
                    output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                    _run_interpolator_statistics_with_retry(interpolator, output_dir, "Global", parent=self)
                
                # Export results (handles both single and multiple interpolators)
                interpolation_file = export_blocks_to_file(blocks, interpolation_file)
                cfg['interpolation_file'] = interpolation_file
                self.interp_edit.setText(interpolation_file)
                if block_evaluated_samples_file:
                    export_block_evaluated_samples_to_csv(blocks, block_evaluated_samples_file)
                print(f"Interpolation complete! Results saved to:\n  {interpolation_file}")
                if block_evaluated_samples_file:
                    print(f"Block-evaluated samples saved to:\n  {block_evaluated_samples_file}")
                print("=" * 60)
                self._prefer_interpolation_file_for_viewer = not _blocks_use_adaptive_leaf_cover(blocks)
                
                success_lines = [f"Interpolation complete!\nResults saved to:\n{interpolation_file}"]
                if block_evaluated_samples_file:
                    success_lines.append(f"Block-evaluated samples saved to:\n{block_evaluated_samples_file}")
                QtWidgets.QMessageBox.information(self, "Success", "\n\n".join(success_lines))

            except Exception as e:
                print(f"Error during interpolation: {e}")
                import traceback
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, "Error", f"An error occurred during interpolation:\n{str(e)}")

    # Create and run application
    app = QtWidgets.QApplication(sys.argv)
    dialog = ConfigDialog()
    app.aboutToQuit.connect(dialog._prepare_for_shutdown)
    atexit.register(dialog._prepare_for_shutdown)

    def _handle_sigint(*_args):
        if not dialog.isVisible():
            return
        QtCore.QTimer.singleShot(0, dialog.close)

    signal.signal(signal.SIGINT, _handle_sigint)
    sigint_pump = QtCore.QTimer()
    sigint_pump.setInterval(200)
    sigint_pump.timeout.connect(lambda: None)
    sigint_pump.start()

    dialog.show()

    try:
        exit_code = app.exec_()
    finally:
        if sigint_pump.isActive():
            sigint_pump.stop()
        dialog._prepare_for_shutdown()
    sys.exit(exit_code)