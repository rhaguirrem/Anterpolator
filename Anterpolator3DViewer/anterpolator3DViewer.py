import pandas as pd
import numpy as np
import csv
import atexit
from collections import Counter
import importlib.util
import io
import threading
import time
import tempfile
import subprocess
import traceback
import warnings
from tqdm import tqdm
import sys
import os
from datetime import datetime
import xml.etree.ElementTree as ET
from matplotlib.colors import ListedColormap
import json
from PyQt5 import QtWidgets, QtCore
sys.path.append("C:/Projects/Anterpolator")

pv = None
taichi_runtime_module = None
LARGE_BLOCK_FILE_THRESHOLD = 512 * 1024 * 1024
INITIAL_BLOCK_RENDER_THRESHOLD = 5000
INVALID_FILENAME_CHARS = str.maketrans({ch: '_' for ch in '<>:"/\\|?*'})


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
    try:
        with open(path, 'r', encoding='utf-8', errors='ignore') as f:
            for current, line in enumerate(f, start=1):
                if current == line_number:
                    raw = line.strip('\n\r')
                    tokens = [t.strip() for t in raw.split(delimiter)]
                    tokens = [t for t in tokens if t != '']
                    if not tokens:
                        raise ValueError(f"Parsed header line {line_number} in '{os.path.basename(path)}' produced no tokens.")
                    return tokens
        raise ValueError(f"Header line {line_number} exceeds total lines in file '{os.path.basename(path)}'.")
    except UnicodeDecodeError:
        raise ValueError(f"Could not decode file '{path}' with utf-8 encoding.")

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

def build_unique_column_names(headers):
    name_counts = {}
    final_names = []
    for header in headers:
        key = header if header else 'Unnamed'
        count = name_counts.get(key, 0)
        final_names.append(f"{key}_{count}" if count > 0 else key)
        name_counts[key] = count + 1
    return final_names

def read_selected_columns_with_header(path, delimiter, header_line, selected_columns, progress_label=None):
    headers = parse_header_line(path, delimiter, header_line)
    final_names = build_unique_column_names(headers)
    missing = [col for col in selected_columns if col not in final_names]
    if missing:
        raise ValueError(f"Selected columns not found in file '{os.path.basename(path)}': {missing}")
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=header_line,
        comment='#',
        usecols=selected_columns,
    )
    if progress_label:
        df = read_csv_with_progress(path, progress_label, **read_kwargs)
    else:
        df = pd.read_csv(path, **prepare_csv_read_kwargs(path, **read_kwargs))
    df = strip_leading_non_data_rows(df)
    df._detected_delimiter = delimiter
    return df, final_names

def load_samples_dataframe(samples_file, samples_delimiter=None, samples_header_line=1,
                          sample_x_col=None, sample_y_col=None, sample_z_col=None, sample_value_col=None,
                          progress_label=None):
    explicit_mapping = all([sample_x_col, sample_y_col, sample_z_col, sample_value_col])
    if explicit_mapping:
        delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
        selected_columns = [sample_x_col, sample_y_col, sample_z_col, sample_value_col]
        df, parsed_cols = read_selected_columns_with_header(
            samples_file,
            delimiter,
            samples_header_line or 1,
            selected_columns,
            progress_label=progress_label,
        )
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


def infer_sample_domains_from_blocks(sample_coords, blocks_file, block_size,
                                     blocks_delimiter=None, blocks_header_line=1,
                                     block_x_col=None, block_y_col=None, block_z_col=None,
                                     block_domain_col=None, progress_callback=None):
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
        config=None,
        progress_callback=progress_callback,
    )

    coords_for_mapping = coords
    if block_metadata.get('is_rotated') and len(coords_for_mapping) > 0:
        rotation_center = block_metadata['rotation_center']
        rotation_matrix = block_metadata['rotation_matrix']
        coords_for_mapping = (coords_for_mapping - rotation_center) @ rotation_matrix.T

    all_min_bounds = np.asarray(block_metadata['all_min_bounds'], dtype=float)
    unified_dims = np.asarray(block_metadata['unified_dims'], dtype=float)
    sample_block_indices = np.floor((coords_for_mapping - all_min_bounds) / unified_dims + 1e-6).astype(int)
    domain_mapping = block_metadata['domain_mapping']

    return np.array(
        [domain_mapping.get((int(idx[0]), int(idx[1]), int(idx[2])), '') for idx in sample_block_indices],
        dtype=object,
    )


def apply_blank_sample_domain_behavior(df, blank_domain_behavior='skip', domain_col='Domain',
                                       x_col='x', y_col='y', z_col='z',
                                       blocks_file=None, blocks_delimiter=None, blocks_header_line=1,
                                       block_x_col=None, block_y_col=None, block_z_col=None,
                                       block_domain_col=None, block_size=None, progress_callback=None):
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
                                                block_domain_col=None, block_size=None, progress_callback=None):
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

def build_domain_catalog_cache_signature(blocks_file, delimiter, header_line, domain_col):
    stats = os.stat(blocks_file)
    return {
        'path': os.path.abspath(blocks_file),
        'size': int(stats.st_size),
        'mtime_ns': int(getattr(stats, 'st_mtime_ns', int(stats.st_mtime * 1_000_000_000))),
        'delimiter': delimiter,
        'header_line': int(header_line or 1),
        'domain_col': str(domain_col or ''),
    }

def load_block_domain_catalog(blocks_file, delimiter, header_line, domain_col,
                              chunksize=250_000, progress_callback=None):
    headers = parse_header_line(blocks_file, delimiter, header_line)
    final_names = build_unique_column_names(headers)
    if domain_col not in final_names:
        raise ValueError(f'Domain column "{domain_col}" not found in blocks file.')

    domains = set()
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=header_line,
        comment='#',
        usecols=[domain_col],
        chunksize=chunksize,
    )

    for chunk in iterate_csv_with_progress(
        blocks_file,
        'Reading domain column',
        progress_callback=progress_callback,
        **read_kwargs,
    ):
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

def read_csv_with_selected_header(path, delimiter, header_line, expected_min_cols=1,
                                  progress_label=None, progress_callback=None):
    """Read CSV using a specific header line (1-based). Returns DataFrame.
    Uses manual header parsing to build names and then pandas read_csv with skiprows.
    """
    headers = parse_header_line(path, delimiter, header_line)
    final_names = build_unique_column_names(headers)
    read_kwargs = dict(
        delimiter=delimiter,
        header=None,
        names=final_names,
        skiprows=header_line-1,
        comment='#'
    )
    if progress_label:
        df = read_csv_with_progress(path, progress_label, progress_callback=progress_callback, **read_kwargs)
    else:
        df = pd.read_csv(path, **prepare_csv_read_kwargs(path, **read_kwargs))
    if df.shape[0] and all(str(df.iloc[0, i]).strip() == final_names[i] for i in range(min(len(final_names), df.shape[1]))):
        df = df.iloc[1:].reset_index(drop=True)
    df = strip_leading_non_data_rows(df)
    def is_all_empty(series):
        return series.isna().all() or (series.astype(str).str.strip() == '').all()
    empty_cols = [c for c in df.columns if is_all_empty(df[c])]
    if empty_cols:
        df = df.drop(columns=empty_cols)
    if df.shape[1] < expected_min_cols:
        raise ValueError(f"File '{path}' has fewer than {expected_min_cols} non-empty columns after cleanup.")
    df._detected_delimiter = delimiter
    return df, final_names

def detect_csv_delimiter(path):
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
                        progress_callback=None):
    """Read a CSV file detecting delimiter (comma, semicolon, tab, pipe) unless forced.
    Drops empty columns. Returns DataFrame."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")
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
    def is_all_empty(series):
        return series.isna().all() or (series.astype(str).str.strip() == '').all()
    empty_cols = [c for c in df.columns if is_all_empty(df[c])]
    if empty_cols:
        df = df.drop(columns=empty_cols)
    if df.shape[1] < min_cols:
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

    if len(points) > max_points:
        sample_indices = np.linspace(0, len(points) - 1, num=max_points, dtype=int)
        sampled = points[sample_indices]
    else:
        sampled = points

    z_keys = np.round(sampled[:, 2], decimals=6)
    _, inverse = np.unique(z_keys, return_inverse=True)
    vectors = []

    for level_idx in range(inverse.max() + 1 if len(inverse) else 0):
        level_points = sampled[inverse == level_idx]
        if len(level_points) < 2:
            continue

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
    
    best_y_angle = None
    
    for cand in candidates:
        cand_norm = cand % 360
        # Count vectors close to this angle
        diff = np.abs(angles_deg - cand_norm)
        diff = np.minimum(diff, 360 - diff)
        count = np.sum(diff < 5.0) # within 5 degrees
        
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
    
    # Check significance
    if abs(theta_deg % 360) < 1.0 or abs(theta_deg % 360) > 359.0:
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

def plan_block_file_columns(header_names, block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None):
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
    return selected_columns, rename_map, domain_copy_source, mapping_mode

def normalize_block_chunk(chunk, rename_map, domain_copy_source=None):
    if domain_copy_source and 'Domain' not in chunk.columns and domain_copy_source in chunk.columns:
        chunk['Domain'] = chunk[domain_copy_source]
    if rename_map:
        chunk = chunk.rename(columns=rename_map)
    keep_columns = [c for c in ['x', 'y', 'z', 'Domain'] if c in chunk.columns]
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
                               config=None, progress_callback=None):
    headers = parse_header_line(blocks_file, delimiter, header_line)
    final_names = build_unique_column_names(headers)
    selected_columns, rename_map, domain_copy_source, mapping_mode = plan_block_file_columns(
        final_names,
        block_x_col=block_x_col,
        block_y_col=block_y_col,
        block_z_col=block_z_col,
        block_domain_col=block_domain_col,
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
        skiprows=header_line,
        comment='#',
        usecols=selected_columns,
        chunksize=chunksize,
    )

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
        chunk, dropped = normalize_block_chunk(chunk, rename_map, domain_copy_source)
        if len(chunk) == 0:
            continue
        rotation_samples.append(chunk[['x', 'y', 'z']].to_numpy())
        if sum(len(arr) for arr in rotation_samples) >= rotation_sample_target:
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

    print("Building bounds...")
    all_min_bounds = None
    all_max_bounds = None

    full_read_kwargs = dict(base_read_kwargs)

    bounds_progress = _make_scaled_progress_callback(progress_callback, 25, 60, 'Reading grid file (bounds)')
    for chunk in iterate_csv_with_progress(
        blocks_file,
        "Reading grid file (bounds)",
        progress_callback=bounds_progress,
        **full_read_kwargs,
    ):
        chunk, dropped = normalize_block_chunk(chunk, rename_map, domain_copy_source)
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

    if all_min_bounds is None or all_max_bounds is None:
        raise ValueError('Could not determine grid bounds from blocks file.')

    _emit_progress(progress_callback, 62, 100, 'Calculating grid dimensions...')
    dims_grid = np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int)
    print('Calculated grid dimensions:', dims_grid)

    print("Building domain mapping...")
    grouped_domain_counts = {}
    skipped_count_total = 0
    resolved_rows = 0

    mapping_progress = _make_scaled_progress_callback(progress_callback, 65, 100, 'Reading grid file (domain mapping)')
    for chunk in iterate_csv_with_progress(
        blocks_file,
        "Reading grid file (domain mapping)",
        progress_callback=mapping_progress,
        **full_read_kwargs,
    ):
        chunk, dropped = normalize_block_chunk(chunk, rename_map, domain_copy_source)
        if len(chunk) == 0:
            continue

        coords = chunk[['x', 'y', 'z']].to_numpy(copy=False)
        if is_rotated:
            coords = (coords - rotation_center) @ rotation_matrix.T

        grid_indices = np.floor((coords - all_min_bounds) / unified_dims + 1e-6).astype(int)

        if 'Domain' in chunk.columns:
            domains = chunk['Domain'].fillna("Undomained").astype(str).str.strip().replace('', "Undomained")
        else:
            domains = pd.Series(["Undomained"] * len(chunk))

        if skipped_domains:
            keep_mask = ~domains.isin(skipped_domains)
            skipped_count_total += int((~keep_mask).sum())
            grid_indices = grid_indices[keep_mask.to_numpy()]
            domains = domains[keep_mask].reset_index(drop=True)

        if len(domains) == 0:
            continue

        resolved_rows += len(domains)
        grouped = pd.DataFrame(
            {
                'ix': grid_indices[:, 0],
                'iy': grid_indices[:, 1],
                'iz': grid_indices[:, 2],
                'Domain': domains.to_numpy(copy=False),
            }
        ).groupby(['ix', 'iy', 'iz', 'Domain'], sort=False).size()

        for (ix, iy, iz, domain), count in grouped.items():
            base_idx = (int(ix), int(iy), int(iz))
            domain_counts = grouped_domain_counts.setdefault(base_idx, Counter())
            domain_counts[str(domain)] += int(count)

    if skipped_domains:
        print(f"Skipping domains: {skipped_domains}")
        if skipped_count_total > 0:
            print(f"Skipped {skipped_count_total:,} sub-block rows due to domain overrides.")

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
        'min_quantized_idx': np.zeros(3, dtype=int),
    }

def create_blocks(points, values, block_size=10, verbose=False, range_size=10, max_pheromone=150,
                  ants_per_sample=3, blocks_file=None, background_value=0.0, background_distance=None, average_with_blocks=False,
                  blocks_delimiter=None,
                  avoid_visited_threshold_enabled=False,
                  avoid_visited_threshold=100,
                  blocks_header_line=1,
                  block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                  config=None,
                  sample_domains=None):
    pv = _require_pyvista()
    original_points_array = np.array(points, copy=True)
    original_values_array = np.array(values, copy=True)
    original_domains_array = np.array(sample_domains, copy=True) if sample_domains is not None else None
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
        use_streaming_blocks = os.path.getsize(blocks_file) >= LARGE_BLOCK_FILE_THRESHOLD
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
            load_pbar.update(4)
        elif blocks_header_line and blocks_header_line != 1 and blocks_delimiter:
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
                if len(cols) >= 4:
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

            all_min_bounds = np.min(centroids_all, axis=0)
            all_max_bounds = np.max(centroids_all, axis=0)
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

        # Group sample points using the same all_min_bounds (VECTORIZED)
        print("Assigning points to blocks...")
        # Compute all block indices at once
        points_array = np.array(points)
        values_array = np.array(values)
        block_indices = np.floor((points_array - all_min_bounds) / unified_dims + 1e-6).astype(int)
        allowed_mask = np.array([tuple(idx) in allowed_grid for idx in block_indices], dtype=bool)
        assigned_mask = compute_domain_sensitive_assignment_mask(
            block_indices,
            allowed_grid,
            domain_mapping=domain_mapping,
            sample_domains=sample_domains,
        )
        domain_mismatch_count = int(np.count_nonzero(allowed_mask & ~assigned_mask))
        if domain_mismatch_count:
            print(f"Rejected {domain_mismatch_count:,} samples whose domain does not match their target block domain.")
        
        # Create lookup for allowed blocks
        sample_blocks_dict = {}
        block_values = {}
        
        # Group by block index
        for i in tqdm(range(len(points_array)), desc="Assigning points to blocks"):
            block_idx = tuple(block_indices[i])
            if assigned_mask[i]:
                if block_idx not in sample_blocks_dict:
                    sample_blocks_dict[block_idx] = []
                    block_values[block_idx] = []
                sample_blocks_dict[block_idx].append(points_array[i])
                block_values[block_idx].append(values_array[i])
        
        # Count sample blocks per domain
        sample_domain_counts = {}
        for idx in sample_blocks_dict.keys():
            domain = domain_mapping.get(idx, "Undomained")
            sample_domain_counts[domain] = sample_domain_counts.get(domain, 0) + 1
        
        print(f"\nSample block distribution (samples assigned to blocks):")
        for domain, count in sorted(sample_domain_counts.items()):
            print(f"  {domain}: {count} sample blocks")
        print(f"  Total sample blocks: {len(sample_blocks_dict)}")
        load_pbar.update(1)
        load_pbar.set_postfix_str("creating blocks")

        block_data = []
        for idx in tqdm(sample_blocks_dict.keys(), desc="Creating blocks"):
            corner = all_min_bounds + np.array(idx) * unified_dims
            cell = pv.Box(bounds=(
                corner[0], corner[0] + unified_dims[0],
                corner[1], corner[1] + unified_dims[1],
                corner[2], corner[2] + unified_dims[2]
            ))
            avg_value = np.mean(block_values[idx])
            cell.cell_data['Value'] = np.full(cell.n_cells, avg_value)
            cell.cell_data['Raw_Value'] = np.full(cell.n_cells, avg_value)
            cell.cell_data['Is_Sample'] = np.full(cell.n_cells, True)
            cell.cell_data['Block_ID'] = np.full(cell.n_cells, 0)  # to be set later
            domain = domain_mapping.get(idx, "Undomained")
            cell.cell_data['Domain'] = np.full(cell.n_cells, domain)
            block_data.append(cell)
        block_info = {
            'min_bounds': all_min_bounds,
            'dims': np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int),
            'block_size': unified_dims.tolist(),
            'allowed_grid': list(allowed_grid),
            'rotation_matrix': rotation_matrix if is_rotated else None,
            'rotation_center': rotation_center if is_rotated else None,
            'domain_mapping': domain_mapping,
            'subblock_counts': subblock_counts,
            'mixed_domain_blocks': mixed_domain_blocks
        }
        multiblock = pv.MultiBlock(block_data)
        # Store metadata on multiblock with private-style names to avoid PyVista attribute restrictions
        multiblock._block_info = block_info
        multiblock._sample_blocks = {idx: np.mean(vals) for idx, vals in block_values.items()}
        multiblock._sample_assignment_data = {
            'points': original_points_array,
            'values': original_values_array,
            'domains': original_domains_array,
            'block_indices': block_indices,
            'assigned_mask': assigned_mask,
        }
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
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping)
                
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
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping_2)
                    
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
                    elif algo_type == 'string_theory':
                        interp1.allowed_grid_override = allowed_grid
                        interp1.domain_mapping = domain_mapping
                    use_mapping_1 = True
                else:
                    interp1.allowed_grid_override = None
                    interp1.domain_mapping = None
                    use_mapping_1 = False
                
                interp1.initialize_blocks(multiblock._sample_blocks, tuple(block_info['dims']),
                                         all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping_1)
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
                    elif algo2 == 'string_theory':
                        interp2.allowed_grid_override = allowed_grid
                        interp2.domain_mapping = domain_mapping
                    use_mapping_2 = True
                else:
                    interp2.allowed_grid_override = None
                    interp2.domain_mapping = None
                    use_mapping_2 = False
                
                interp2.initialize_blocks(multiblock._sample_blocks, tuple(block_info['dims']),
                                         all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping_2)
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
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=use_mapping)
                
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
        block_domains = {}
        block_info = {
            'min_bounds': min_bounds,
            'dims': dims,
            'block_size': block_size
        }
        # Vectorized block assignment
        points_array = np.array(points)
        values_array = np.array(values)
        domains_array = None
        if sample_domains is not None:
            domains_array = np.array(sample_domains)
        block_indices = ((points_array - min_bounds) // np.array(block_size)).astype(int)
        assigned_mask = np.ones(len(points_array), dtype=bool)
        
        for i in tqdm(range(len(points_array)), desc="Assigning points to blocks"):
            block_idx = tuple(block_indices[i])
            if block_idx not in blocks:
                blocks[block_idx] = []
                block_values[block_idx] = []
                if domains_array is not None:
                    block_domains[block_idx] = []
            blocks[block_idx].append(points_array[i])
            block_values[block_idx].append(values_array[i])
            if domains_array is not None:
                block_domains[block_idx].append(domains_array[i])
        block_data = []
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

        for idx in tqdm(blocks.keys(), desc="Creating blocks"):
            corner = min_bounds + np.array(idx) * np.array(block_size)
            cell = pv.Box(bounds=(
                corner[0], corner[0] + block_size[0],
                corner[1], corner[1] + block_size[1],
                corner[2], corner[2] + block_size[2]
            ))
            avg_value = np.mean(block_values[idx])
            cell.cell_data['Value'] = np.full(cell.n_cells, avg_value)
            cell.cell_data['Is_Sample'] = np.full(cell.n_cells, True)
            cell.cell_data['Block_ID'] = np.full(cell.n_cells, next_block_id)
            next_block_id += 1
            if is_st_domain_interpolation or is_ant_domain_interpolation:
                dom = sample_block_domain_mapping.get(idx, "")
                cell.cell_data['Domain'] = np.full(cell.n_cells, dom)
            block_data.append(cell)
        sample_blocks = {idx: np.mean(vals) for idx, vals in block_values.items()}
        multiblock = pv.MultiBlock(block_data)
        multiblock._block_info = block_info
        multiblock._sample_blocks = sample_blocks
        multiblock._sample_assignment_data = {
            'points': original_points_array,
            'values': original_values_array,
            'domains': original_domains_array,
            'block_indices': block_indices,
            'assigned_mask': assigned_mask,
        }

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
                                       min_bounds, block_size, use_domain_mapping=use_mapping)
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
                                       min_bounds, block_size, use_domain_mapping=use_mapping)
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
                                           min_bounds, block_size, use_domain_mapping=use_mapping)
            
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

        if getattr(plotter, '_fill_unvisited_domainwise', False):
            try:
                if hasattr(interp, 'fill_unvisited_blocks_domainwise'):
                    created, assigned = interp.fill_unvisited_blocks_domainwise(dims)
                    if created or assigned:
                        changes_made = True
            except Exception as e:
                print(f"Domain-wise fill error: {e}")

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


def _sanitize_filename_fragment(value, fallback='Domain'):
    text = str(value or '').strip()
    if not text or text == '(None)':
        text = fallback
    return text.translate(INVALID_FILENAME_CHARS)


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


def load_full_samples_dataframe(samples_file, samples_delimiter=None, samples_header_line=1,
                                progress_label=None, progress_callback=None):
    delimiter = samples_delimiter or detect_csv_delimiter(samples_file)
    if samples_header_line and samples_header_line != 1:
        df, _ = read_csv_with_selected_header(
            samples_file,
            delimiter,
            samples_header_line,
            expected_min_cols=3,
            progress_label=progress_label,
            progress_callback=progress_callback,
        )
        return df, delimiter

    df = read_autodetect_csv(
        samples_file,
        forced_delimiter=delimiter,
        progress_label=progress_label,
        progress_callback=progress_callback,
    )
    return df, delimiter


def load_full_blocks_dataframe(blocks_file, blocks_delimiter=None, blocks_header_line=1,
                               progress_label=None, progress_callback=None):
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
        return df, delimiter

    df = read_autodetect_csv(
        blocks_file,
        forced_delimiter=delimiter,
        progress_label=progress_label,
        progress_callback=progress_callback,
    )
    return df, delimiter


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


def apply_sample_filters(df_samples, sample_filters=None):
    if not sample_filters:
        return df_samples.copy(), []

    filtered_df = df_samples.copy()
    applied_filters = []

    for raw_filter in sample_filters:
        filter_spec = dict(raw_filter or {})
        field = str(filter_spec.get('field', '')).strip()
        filter_type = str(filter_spec.get('type', '')).strip().lower()

        if not field:
            raise ValueError('Each sample filter must define a field name.')
        if field not in filtered_df.columns:
            raise ValueError(f'Sample filter field not found in samples file: {field}')

        series = filtered_df[field]
        if filter_type == 'categorical':
            raw_values = filter_spec.get('values', [])
            allowed_values = [str(value) for value in raw_values]
            if not allowed_values:
                raise ValueError(f'Categorical filter for {field} must include at least one value.')
            mask = series.fillna('').astype(str).isin(allowed_values)
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
            raise ValueError(f'Unsupported sample filter type for {field}: {filter_type}')

        filtered_df = filtered_df.loc[mask].copy()
        applied_filters.append({
            'field': field,
            'type': filter_type,
            'summary': summarize_sample_filter_spec(filter_spec),
        })

    return filtered_df, applied_filters


def _compute_point_to_set_distance_stats(query_points, reference_points, query_chunk_size=1024,
                                         reference_chunk_size=4096, progress_callback=None,
                                         progress_label='Computing distance statistics',
                                         return_nearest_index=False):
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
            return nan_values.copy(), nan_values, np.full(len(query_points), -1, dtype=int)
        return nan_values.copy(), nan_values

    nearest = np.empty(len(query_points), dtype=float)
    average = np.empty(len(query_points), dtype=float)
    nearest_indices = np.full(len(query_points), -1, dtype=int) if return_nearest_index else None

    for query_start in range(0, len(query_points), query_chunk_size):
        query_end = min(query_start + query_chunk_size, len(query_points))
        query_chunk = query_points[query_start:query_end]
        chunk_nearest = np.full(len(query_chunk), np.inf, dtype=float)
        chunk_distance_sum = np.zeros(len(query_chunk), dtype=float)
        chunk_count = 0
        chunk_nearest_indices = np.full(len(query_chunk), -1, dtype=int) if return_nearest_index else None

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

        nearest[query_start:query_end] = chunk_nearest
        average[query_start:query_end] = chunk_distance_sum / max(chunk_count, 1)
        if return_nearest_index:
            nearest_indices[query_start:query_end] = chunk_nearest_indices
        if progress_callback:
            progress_callback(query_end, len(query_points), progress_label)

    if return_nearest_index:
        return nearest, average, nearest_indices
    return nearest, average


def build_concatenated_sample_ids(sample_df, id_columns):
    selected_columns = []
    for column_name in id_columns or []:
        normalized = str(column_name or '').strip()
        if not normalized or normalized == '(None)' or normalized in selected_columns:
            continue
        if normalized not in sample_df.columns:
            raise ValueError(f'Selected closest-sample ID column not found in samples file: {normalized}')
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


def export_block_domain_sample_metrics(samples_file, blocks_file, output_file=None,
                                       samples_delimiter=None, blocks_delimiter=None,
                                       samples_header_line=1, blocks_header_line=1,
                                       sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                       sample_domain_col=None,
                                       closest_sample_id_cols=None,
                                       block_x_col=None, block_y_col=None, block_z_col=None,
                                       block_domain_col=None, block_size=None,
                                       sample_filters=None, progress_callback=None,
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

    _emit_progress(progress_callback, 41 if use_explicit_sample_domains else 56, 100, 'Applying sample filters...')
    filtered_samples_df, applied_filters = apply_sample_filters(df_samples, sample_filters=sample_filters)

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
    candidate_sample_ids, selected_id_columns = build_concatenated_sample_ids(filtered_samples_df, closest_sample_id_cols)
    if candidate_sample_ids is not None:
        candidate_sample_ids = np.asarray(candidate_sample_ids, dtype=object)[valid_sample_mask.to_numpy()]

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

        all_min_bounds = np.asarray(block_metadata['all_min_bounds'], dtype=float)
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

    samples_by_domain = {}
    for domain in np.unique(matched_sample_domains):
        domain_mask = matched_sample_domains == domain
        samples_by_domain[str(domain)] = {
            'coords': matched_sample_coords[domain_mask],
            'ids': None if matched_sample_ids is None else matched_sample_ids[domain_mask],
        }

    nn_column = f'{domain_column_name}_NN_Distance'
    avg_column = f'{domain_column_name}_Avg_Distance'
    closest_id_column = f'{domain_column_name}_Closest_Sample_ID' if selected_id_columns else None

    output_df = df_blocks.copy()
    output_df[nn_column] = np.nan
    output_df[avg_column] = np.nan
    if closest_id_column:
        output_df[closest_id_column] = ''

    processed_block_count = 0
    populated_block_count = 0
    domains_to_process = [value for value in sorted(block_domain_series.unique()) if value]
    total_domain_blocks = int(sum(int((valid_block_mask & (block_domain_series == domain)).sum()) for domain in domains_to_process))
    completed_domain_blocks = 0
    compute_started_at = time.perf_counter()
    next_terminal_report_at = compute_started_at
    _emit_progress(progress_callback, 80 if use_explicit_sample_domains else 85, 100, 'Computing domain distance statistics...')

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
        distance_stats = _compute_point_to_set_distance_stats(
            query_points,
            domain_samples['coords'],
            progress_callback=(
                None if total_domain_blocks <= 0 else
                lambda current, total, label, completed=completed_domain_blocks: _emit_progress(
                    progress_callback,
                    (80 if use_explicit_sample_domains else 85) + int(round(((completed + current) / total_domain_blocks) * (18 if use_explicit_sample_domains else 13))),
                    100,
                    f'{label}... ({domain})',
                )
            ),
            return_nearest_index=bool(closest_id_column),
        )
        if closest_id_column:
            nearest_distances, average_distances, nearest_sample_indices = distance_stats
        else:
            nearest_distances, average_distances = distance_stats
        output_df.loc[domain_block_indices, nn_column] = nearest_distances
        output_df.loc[domain_block_indices, avg_column] = average_distances
        if closest_id_column:
            domain_ids = domain_samples['ids']
            closest_ids = [domain_ids[index] if index >= 0 else '' for index in nearest_sample_indices]
            output_df.loc[domain_block_indices, closest_id_column] = closest_ids
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

    output_dir = os.path.dirname(output_file) or '.'
    os.makedirs(output_dir, exist_ok=True)
    _emit_progress(progress_callback, 99, 100, 'Writing metrics export...')
    output_df.to_csv(output_file, index=False, sep=output_delimiter)
    _emit_progress(progress_callback, 100, 100, 'Block-domain sample metrics export complete.')

    return {
        'output_file': output_file,
        'domain_column': domain_column_name,
        'nearest_distance_column': nn_column,
        'average_distance_column': avg_column,
        'closest_sample_id_column': closest_id_column,
        'closest_sample_id_source_columns': list(selected_id_columns),
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


def export_domained_samples_from_blocks(samples_file, blocks_file, output_file=None,
                                        samples_delimiter=None, blocks_delimiter=None,
                                        samples_header_line=1, blocks_header_line=1,
                                        sample_x_col=None, sample_y_col=None, sample_z_col=None,
                                        block_x_col=None, block_y_col=None, block_z_col=None,
                                        block_domain_col=None, block_size=None,
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
        config=None,
        progress_callback=_make_scaled_progress_callback(progress_callback, 0, 45, 'Loading block domain mapping...'),
    )

    print(f"Loading samples from {samples_file}...")
    df_samples, sample_delimiter = load_full_samples_dataframe(
        samples_file,
        samples_delimiter=samples_delimiter,
        samples_header_line=samples_header_line,
        progress_label='Reading sample file',
        progress_callback=_make_scaled_progress_callback(progress_callback, 45, 80, 'Reading sample file...'),
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

    if block_metadata.get('is_rotated') and len(valid_coords) > 0:
        rotation_center = block_metadata['rotation_center']
        rotation_matrix = block_metadata['rotation_matrix']
        valid_coords = (valid_coords - rotation_center) @ rotation_matrix.T

    all_min_bounds = np.asarray(block_metadata['all_min_bounds'], dtype=float)
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


def export_blocks_to_csv(blocks, filepath):
    output_dir = os.path.dirname(filepath) or '.'
    os.makedirs(output_dir, exist_ok=True)

    print(f"Exporting blocks to {filepath}...")
    data = _collect_export_block_data(blocks)
    
    # Export to CSV
    df = pd.DataFrame([{k: v for k, v in row.items() if k != '_Grid_Index'} for row in data])
    df.to_csv(filepath, index=False)
    print(f"Exported {len(data)} blocks to {filepath}")


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
    name = str(algo_name or '').strip().lower()
    if 'ant colony' in name:
        return 'Anterpolator', 'ant_colony'
    if 'string theory' in name:
        return 'String Theory', 'string_theory'
    if 'molecular clock' in name or 'phylogeographic' in name or 'biochemical clock' in name:
        return 'Molecular Clock', 'molecular_clock'
    return 'Unknown', 'unknown'


def _add_interpolator_blocks_to_data(interpolator, min_bounds, block_size, data, rotation_matrix=None, rotation_center=None, domain_mapping=None, original_samples=None, pass_count=1, forced_domain=None, first_pass_algorithm_name=None, final_algorithm_name=None):
    """Process blocks from an interpolator and add to data list"""
    final_algorithm_label, final_algo_type = _normalize_export_algorithm_name(
        final_algorithm_name or interpolator.get_algorithm_name()
    )
    first_pass_algorithm_label, first_pass_algo_type = _normalize_export_algorithm_name(
        first_pass_algorithm_name or final_algorithm_name or interpolator.get_algorithm_name()
    )
    
    for pos, block in tqdm(interpolator.blocks.items(), desc="Processing blocks"):
        # Calculate block centroid - grid indices are calculated relative to centroids in min_bounds,
        # so we don't add block_size/2 (that would shift by half a block)
        centroid = min_bounds + np.array(pos) * np.array(block_size)
        
        # Apply inverse rotation if needed
        if rotation_matrix is not None and rotation_center is not None:
            # P_orig = P_aligned @ R + Center
            centroid = centroid @ rotation_matrix + rotation_center
        
        # Determine Source
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
            'Is_Feeder': None
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
                'Is_Feeder': is_feeder
            })
            
        data.append(row)

def silent_interpolation(plotter, iterations, interpolation_file):
    blocks = plotter._blocks_data
    dims = tuple(blocks._block_info['dims'])
    block_evaluated_samples_file = getattr(plotter, '_block_evaluated_samples_file', None)
    
    # Check if we have multiple interpolators (sequential domain processing)
    if hasattr(blocks, '_interpolators'):
        print(f"Running sequential domain interpolation for {len(blocks._interpolators)} domains...")
        
        for domain_idx, (domain, interpolator_list) in enumerate(blocks._interpolators.items(), 1):
            # interpolator_list is [Pass1, (optional) Pass2]
            
            # --- Pass 1 ---
            interp1 = interpolator_list[0]
            algo_name1 = interp1.get_algorithm_name()
            print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} - Pass 1 ({algo_name1}) ===")
            
            # Force verbose for first iteration if it's an AntColony
            if hasattr(interp1, 'verbose'):
                original_verbose = interp1.verbose
                interp1.verbose = True
            
            pbar = tqdm(range(iterations), desc=f"Domain {domain} - Pass 1")
            for i in pbar:
                should_continue = interp1.run_iteration(dims)
                
                if i == 0 and hasattr(interp1, 'verbose'):
                    interp1.verbose = original_verbose
                
                if not should_continue or interp1.is_converged():
                    pbar.set_description(f"Domain {domain} - Pass 1 (converged)")
                    print(f"Pass 1 converged at iteration {i+1}")
                    break
            
            # Generate stats for Pass 1 if String Theory
            output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
            if hasattr(interp1, 'generate_statistics'):
                interp1.generate_statistics(output_dir, domain_name=f"{domain}_Pass1")

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
                     
                     interp2.initialize_blocks(pass1_values, dims, min_bounds, block_size, use_domain_mapping=use_mapping)
                     
                     if hasattr(interp2, 'create_ants'):
                         interp2.create_ants()
                    
                     # Run Pass 2
                     if hasattr(interp2, 'verbose'):
                         original_verbose = interp2.verbose
                         interp2.verbose = True
                    
                     pbar = tqdm(range(iterations), desc=f"Domain {domain} - Pass 2")
                     for i in pbar:
                         should_continue = interp2.run_iteration(dims)
                         
                         if i == 0 and hasattr(interp2, 'verbose'):
                             interp2.verbose = original_verbose
                         
                         if not should_continue or interp2.is_converged():
                             pbar.set_description(f"Domain {domain} - Pass 2 (converged)")
                             print(f"Pass 2 converged at iteration {i+1}")
                             break
                
                # Generate stats for Pass 2 if String Theory
                output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                if hasattr(interp2, 'generate_statistics'):
                    interp2.generate_statistics(output_dir, domain_name=f"{domain}_Pass2")
            
            # Print domain summary (of the last pass)
            last_interp = interpolator_list[-1]
            metadata = last_interp.get_metadata()
            print(f"\n=== Domain {domain} Summary ===")
            for key, value in metadata.items():
                print(f"{key}: {value}")
    else:
        # Single interpolator
        interpolator = blocks._ant_colony
        algo_name = interpolator.get_algorithm_name()
        print(f"Running {algo_name} for {iterations} iterations...")
        
        pbar = tqdm(range(iterations), desc=f"Interpolation ({algo_name})")
        for i in pbar:
            should_continue = interpolator.run_iteration(dims)
            if not should_continue or interpolator.is_converged():
                print(f"Converged at iteration {i+1}")
                break
        
        metadata = interpolator.get_metadata()
        print(f"\n=== Summary ===")
        for key, value in metadata.items():
            print(f"{key}: {value}")
            
        # Generate stats if String Theory
        output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
        if hasattr(interpolator, 'generate_statistics'):
            interpolator.generate_statistics(output_dir, domain_name="Global")
    
    # Export results (handles both single and multiple interpolators)
    export_blocks_to_csv(blocks, interpolation_file)
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
        progress_label='Reading sample file',
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
        explicit_sample_map = None

    if needs_sample_domains:
        df = normalize_selected_sample_domain_column(df, sample_domain_col=config.get('sample_domain_col'))

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
            progress_label='Reading sample file',
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
            explicit_sample_map = None

        if wants_domain_any:
            df = normalize_selected_sample_domain_column(df, sample_domain_col=config.get('sample_domain_col') if config else None)

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
            sample_domains=sample_domains
        )

        # Store settings in plotter
        plotter._blocks_data = blocks
        plotter._verbose = verbose
        plotter._value_filter = value_filter  # Store value_filter
        plotter._fill_unvisited_domainwise = fill_unvisited_domainwise
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
            
            self.domain_configs = {}  # domain -> {'algorithm': str, 'second_pass_algorithm': str, 'skip': bool}
            
            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)
            
            # Info label
            info = QtWidgets.QLabel("Configure which algorithm to use for each domain.\n"
                                  "You can configure a second pass to run after the first one completes.\n"
                                  "The second pass uses the output of the first pass as input.")
            info.setWordWrap(True)
            layout.addWidget(info)
            
            # Table for domain mappings
            self.table = QtWidgets.QTableWidget()
            self.table.setColumnCount(3)
            self.table.setHorizontalHeaderLabels(['Domain', 'First Pass Algorithm', 'Second Pass Algorithm'])
            self.table.horizontalHeader().setStretchLastSection(True)
            self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
            layout.addWidget(self.table)
            
            self.populate_domains(domains)
            
            # Buttons
            btn_layout = QtWidgets.QHBoxLayout()
            self.apply_all_btn = QtWidgets.QPushButton('Apply First Pass to All')
            self.apply_all_btn.clicked.connect(self.apply_to_all)
            btn_layout.addWidget(self.apply_all_btn)
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
                algo1_combo.addItems(['(use default)', 'ant_colony', 'molecular_clock', 'string_theory', 'skip'])
                algo1_combo.setCurrentText('(use default)')
                self.table.setCellWidget(i, 1, algo1_combo)

                algo2_combo = QtWidgets.QComboBox()
                algo2_combo.addItems(['skip', 'ant_colony', 'molecular_clock', 'string_theory'])
                algo2_combo.setCurrentText('skip')
                self.table.setCellWidget(i, 2, algo2_combo)

                algo1_combo.currentTextChanged.connect(
                    lambda text, row=i: self.on_first_pass_changed(row, text)
                )
                algo2_combo.currentTextChanged.connect(
                    lambda text, row=i: self.on_second_pass_changed(row, text)
                )
        
        def on_first_pass_changed(self, row, text):
            """Handle changes to first pass algorithm"""
            algo2_combo = self.table.cellWidget(row, 2)
            if not algo2_combo:
                return
                
            if text == 'skip':
                # If first pass is skip, second pass must be skip and disabled
                algo2_combo.setCurrentText('skip')
                algo2_combo.setEnabled(False)
            else:
                algo2_combo.setEnabled(True)
                # Ensure second pass is not same as first (unless skip)
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
            algorithms = ['(use default)', 'ant_colony', 'molecular_clock', 'string_theory', 'skip']
            algo, ok = QtWidgets.QInputDialog.getItem(
                self, 'Apply to All', 'Select first pass algorithm for all domains:', 
                algorithms, 0, False
            )
            if ok:
                for i in range(self.table.rowCount()):
                    combo = self.table.cellWidget(i, 1)
                    if combo:
                        combo.setCurrentText(algo)
        
        def get_domain_configs(self):
            """Get domain algorithm configurations"""
            configs = {}
            for i in range(self.table.rowCount()):
                domain_item = self.table.item(i, 0)
                algo1_combo = self.table.cellWidget(i, 1)
                algo2_combo = self.table.cellWidget(i, 2)
                
                if domain_item and algo1_combo and algo2_combo:
                    domain = domain_item.text()
                    algo1 = algo1_combo.currentText()
                    algo2 = algo2_combo.currentText()
                    
                    config = {}
                    
                    # First Pass
                    if algo1 == 'skip':
                        config['skip'] = True
                    elif algo1 != '(use default)':
                        config['algorithm'] = algo1
                    
                    # Second Pass
                    if algo2 != 'skip':
                        config['second_pass_algorithm'] = algo2
                    
                    if config:
                        configs[domain] = config
            
            return configs
        
        def set_domain_configs(self, configs):
            """Set domain algorithm configurations from loaded config"""
            for i in range(self.table.rowCount()):
                domain_item = self.table.item(i, 0)
                algo1_combo = self.table.cellWidget(i, 1)
                algo2_combo = self.table.cellWidget(i, 2)
                
                if domain_item and algo1_combo and algo2_combo:
                    domain = domain_item.text()
                    if domain in configs:
                        config = configs[domain]
                        
                        # Set First Pass
                        if config.get('skip', False):
                            algo1_combo.setCurrentText('skip')
                        elif 'algorithm' in config:
                            algo1_combo.setCurrentText(config['algorithm'])
                        
                        # Set Second Pass
                        if 'second_pass_algorithm' in config:
                            algo2_combo.setCurrentText(config['second_pass_algorithm'])
                        else:
                            algo2_combo.setCurrentText('skip')

        def accept(self):
            super().accept()
        
        def reject(self):
            super().reject()

    class SampleFilterEditDialog(QtWidgets.QDialog):
        def __init__(self, df_samples, filter_spec=None, parent=None):
            super().__init__(parent)
            self.df_samples = df_samples
            self.setWindowTitle('Sample Filter')
            self.resize(520, 420)

            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)

            form = QtWidgets.QFormLayout()
            layout.addLayout(form)

            self.field_combo = QtWidgets.QComboBox()
            self.field_combo.addItems([str(col) for col in df_samples.columns])
            form.addRow('Field', self.field_combo)

            self.type_combo = QtWidgets.QComboBox()
            self.type_combo.addItems(['categorical', 'numeric'])
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
            series = self.df_samples[field] if field in self.df_samples.columns else pd.Series(dtype=object)

            if filter_type == 'categorical':
                self.criteria_stack.setCurrentIndex(0)
                unique_values = sorted({str(value) for value in series.dropna().astype(str)})
                truncated = False
                if len(unique_values) > 1000:
                    unique_values = unique_values[:1000]
                    truncated = True
                self.value_list.clear()
                self.value_list.addItems(unique_values)
                self.value_hint.setText(
                    f'Loaded {len(unique_values):,} distinct values.' +
                    (' Showing the first 1,000 values.' if truncated else '')
                )
            else:
                self.criteria_stack.setCurrentIndex(1)
                numeric_values = pd.to_numeric(series, errors='coerce').dropna()
                if len(numeric_values) == 0:
                    self.numeric_hint.setText('No numeric values detected in this field.')
                else:
                    self.numeric_hint.setText(
                        f'Available numeric range: {numeric_values.min():g} to {numeric_values.max():g}'
                    )

        def get_filter_spec(self):
            field = self.field_combo.currentText().strip()
            filter_type = self.type_combo.currentText().strip().lower()
            if not field:
                raise ValueError('Please select a field.')

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
        def __init__(self, df_samples, filters=None, parent=None):
            super().__init__(parent)
            self.df_samples = df_samples
            self.filters = [dict(entry) for entry in (filters or [])]
            self.setWindowTitle('Sample Filters')
            self.resize(760, 420)

            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)

            info = QtWidgets.QLabel(
                'Add one or more sample filters. Categorical filters keep selected values. '
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
            dialog = SampleFilterEditDialog(self.df_samples, parent=self)
            if dialog.exec_() == QtWidgets.QDialog.Accepted:
                self.filters.append(dialog.get_filter_spec())
                self._refresh_table()

        def edit_selected_filter(self):
            row = self.table.currentRow()
            if row < 0 or row >= len(self.filters):
                QtWidgets.QMessageBox.information(self, 'Edit Filter', 'Select a filter to edit.')
                return
            dialog = SampleFilterEditDialog(self.df_samples, filter_spec=self.filters[row], parent=self)
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
            self.setWindowTitle("Anterpolator 3D Viewer Configuration")
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
            operations_form = QtWidgets.QFormLayout()
            operations_tab.setLayout(operations_form)
            tabs.addTab(operations_tab, "Operations")
            
            # Tab 2: Ant Colony Parameters
            ant_tab = QtWidgets.QWidget()
            ant_form = QtWidgets.QFormLayout()
            ant_tab.setLayout(ant_form)
            tabs.addTab(ant_tab, "Ant Colony")
            
            # Tab 3: Molecular Clock Parameters
            mc_tab = QtWidgets.QWidget()
            mc_form = QtWidgets.QFormLayout()
            mc_tab.setLayout(mc_form)
            tabs.addTab(mc_tab, "Molecular Clock")

            st_tab = QtWidgets.QWidget()
            st_form = QtWidgets.QFormLayout()
            st_tab.setLayout(st_form)
            tabs.addTab(st_tab, "String Theory")

            display_tab = QtWidgets.QWidget()
            display_form = QtWidgets.QFormLayout()
            display_tab.setLayout(display_form)
            tabs.addTab(display_tab, "Display")

            # Tab 6: Advanced Options
            advanced_tab = QtWidgets.QWidget()
            advanced_form = QtWidgets.QFormLayout()
            advanced_tab.setLayout(advanced_form)
            tabs.addTab(advanced_tab, "Advanced")

            # === FILES & DATA TAB ===
            self.samples_edit = QtWidgets.QLineEdit('Data/ANT-Samples.csv')
            self.blocks_edit = QtWidgets.QLineEdit('Data/ANT-Domains.csv')
            self.color_edit = QtWidgets.QLineEdit('Data/Value.lfc')
            self.interp_edit = QtWidgets.QLineEdit('')

            def add_file_row(label, line_edit, filter_str, form_layout):
                h = QtWidgets.QHBoxLayout()
                h.addWidget(line_edit)
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
            
            add_file_row('Samples File', self.samples_edit, 'CSV Files (*.csv)', files_form)
            add_file_row('Blocks File', self.blocks_edit, 'CSV Files (*.csv)', files_form)
            add_file_row('Color File', self.color_edit, 'LFC Files (*.lfc);;All Files (*.*)', files_form)
            add_file_row('Interpolation File', self.interp_edit, 'CSV Files (*.csv)', files_form)

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

            # Algorithm selection
            self.algorithm_combo = QtWidgets.QComboBox()
            self.algorithm_combo.addItems(['ant_colony', 'molecular_clock', 'string_theory'])
            
            self.second_pass_combo = QtWidgets.QComboBox()
            self.second_pass_combo.addItems(['skip', 'ant_colony', 'molecular_clock', 'string_theory'])
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
            
            algo_layout = QtWidgets.QHBoxLayout()
            algo_layout.addWidget(self.algorithm_combo)
            algo_layout.addWidget(QtWidgets.QLabel("Second Pass:"))
            algo_layout.addWidget(self.second_pass_combo)
            
            files_form.addRow('First Pass', algo_layout)

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
            self.sample_x_col = QtWidgets.QComboBox(); self.sample_y_col = QtWidgets.QComboBox(); self.sample_z_col = QtWidgets.QComboBox(); self.sample_value_col = QtWidgets.QComboBox(); self.sample_domain_col = QtWidgets.QComboBox()
            self.block_domain_metrics_id_cols = [QtWidgets.QComboBox(), QtWidgets.QComboBox(), QtWidgets.QComboBox()]
            for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col, self.sample_domain_col]:
                configure_column_combo(cb)
            for cb in self.block_domain_metrics_id_cols:
                configure_column_combo(cb)
                cb.addItem('(None)')
            self.sample_domain_col.addItem('(None)')
            sample_map_layout = QtWidgets.QHBoxLayout(); sample_map_layout.addWidget(self.sample_x_col); sample_map_layout.addWidget(self.sample_y_col); sample_map_layout.addWidget(self.sample_z_col); sample_map_layout.addWidget(self.sample_value_col); sample_map_layout.addWidget(self.sample_domain_col)
            for index in range(5):
                sample_map_layout.setStretch(index, 1)
            files_form.addRow('Samples Columns (X Y Z Value Domain)', sample_map_layout)

            # Column mapping combo boxes for blocks
            self.block_x_col = QtWidgets.QComboBox(); self.block_y_col = QtWidgets.QComboBox(); self.block_z_col = QtWidgets.QComboBox(); self.block_domain_col = QtWidgets.QComboBox()
            for cb in [self.block_x_col, self.block_y_col, self.block_z_col, self.block_domain_col]:
                configure_column_combo(cb)
            self.block_domain_col.addItem('(None)')
            block_map_layout = QtWidgets.QHBoxLayout(); block_map_layout.addWidget(self.block_x_col); block_map_layout.addWidget(self.block_y_col); block_map_layout.addWidget(self.block_z_col); block_map_layout.addWidget(self.block_domain_col)
            for index in range(4):
                block_map_layout.setStretch(index, 1)
            files_form.addRow('Blocks Columns (X Y Z Domain)', block_map_layout)

            # Block Size
            self.block_x = QtWidgets.QSpinBox(); self.block_x.setRange(1, 10000); self.block_x.setValue(10)
            self.block_y = QtWidgets.QSpinBox(); self.block_y.setRange(1, 10000); self.block_y.setValue(10)
            self.block_z = QtWidgets.QSpinBox(); self.block_z.setRange(1, 10000); self.block_z.setValue(10)
            bx_layout = QtWidgets.QHBoxLayout(); bx_layout.addWidget(self.block_x); bx_layout.addWidget(self.block_y); bx_layout.addWidget(self.block_z)
            files_form.addRow('Block Size (x,y,z)', bx_layout)

            def refresh_sample_columns():
                path = self.samples_edit.text().strip()
                delim = self.samples_delim.currentText()
                header_line = self.samples_header_line.value()
                current_domain = self.sample_domain_col.currentText()
                current_metrics_id_cols = [cb.currentText() for cb in self.block_domain_metrics_id_cols]
                for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col]:
                    cb.clear()
                self.sample_domain_col.clear(); self.sample_domain_col.addItem('(None)')
                for cb in self.block_domain_metrics_id_cols:
                    cb.clear(); cb.addItem('(None)')
                if not os.path.isfile(path):
                    return
                try:
                    cols = parse_header_line(path, delim, header_line)
                    for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col, self.sample_domain_col]:
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
                    if current_domain and current_domain != '(None)':
                        idx = self.sample_domain_col.findText(current_domain)
                        if idx >= 0:
                            self.sample_domain_col.setCurrentIndex(idx)
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
                self.block_domain_col.clear(); self.block_domain_col.addItem('(None)')
                if not os.path.isfile(path):
                    return
                try:
                    cols = parse_header_line(path, delim, header_line)
                    for cb in [self.block_x_col, self.block_y_col, self.block_z_col, self.block_domain_col]:
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
                    # Restore domain selection if possible
                    if current_domain and current_domain != '(None)':
                        idx = self.block_domain_col.findText(current_domain)
                        if idx >= 0: self.block_domain_col.setCurrentIndex(idx)
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

            def suggested_block_domain_metrics_path():
                return resolve_block_domain_metrics_export_path(
                    None,
                    blocks_file=self.blocks_edit.text().strip(),
                    domain_col=self.block_domain_col.currentText(),
                )

            self.domain_samples_output_edit = QtWidgets.QLineEdit('')
            self.domain_samples_browse = QtWidgets.QPushButton('Browse')
            self.start_domaining_btn = QtWidgets.QPushButton('Start Domaining')
            self.block_domain_metrics_output_edit = QtWidgets.QLineEdit('')
            self.block_domain_metrics_browse = QtWidgets.QPushButton('Browse')
            self.configure_block_domain_metrics_filters_btn = QtWidgets.QPushButton('Configure Filters...')
            self.start_block_domain_metrics_btn = QtWidgets.QPushButton('Export Metrics')
            self.block_domain_sample_filters = []
            self.block_domain_metrics_filters_summary = QtWidgets.QLabel('No sample filters configured. All samples will be used.')
            self.block_domain_metrics_filters_summary.setWordWrap(True)
            self._domain_samples_auto_path = ''
            self._block_domain_metrics_auto_path = ''

            def refresh_domain_samples_output_path(force=False):
                suggested = suggested_domain_samples_path()
                current = self.domain_samples_output_edit.text().strip()
                if force or not current or current == self._domain_samples_auto_path:
                    self.domain_samples_output_edit.setText(suggested)
                self._domain_samples_auto_path = suggested

            def refresh_block_domain_metrics_output_path(force=False):
                suggested = suggested_block_domain_metrics_path()
                current = self.block_domain_metrics_output_edit.text().strip()
                if force or not current or current == self._block_domain_metrics_auto_path:
                    self.block_domain_metrics_output_edit.setText(suggested)
                self._block_domain_metrics_auto_path = suggested

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

            domain_samples_group = QtWidgets.QGroupBox('Domain Samples')
            domain_samples_form = QtWidgets.QFormLayout()
            domain_samples_group.setLayout(domain_samples_form)
            domain_output_layout = QtWidgets.QHBoxLayout()
            domain_output_layout.addWidget(self.domain_samples_output_edit)
            domain_output_layout.addWidget(self.domain_samples_browse)
            domain_samples_form.addRow('Output File', domain_output_layout)
            domain_samples_form.addRow('', self.start_domaining_btn)
            operations_form.addRow(domain_samples_group)

            block_metrics_group = QtWidgets.QGroupBox('Block Domain Sample Metrics')
            block_metrics_form = QtWidgets.QFormLayout()
            block_metrics_group.setLayout(block_metrics_form)
            block_metrics_output_layout = QtWidgets.QHBoxLayout()
            block_metrics_output_layout.addWidget(self.block_domain_metrics_output_edit)
            block_metrics_output_layout.addWidget(self.block_domain_metrics_browse)
            block_metrics_form.addRow('Output File', block_metrics_output_layout)
            block_metrics_id_layout = QtWidgets.QHBoxLayout()
            for index, cb in enumerate(self.block_domain_metrics_id_cols):
                block_metrics_id_layout.addWidget(cb)
                block_metrics_id_layout.setStretch(index, 1)
            block_metrics_form.addRow('Closest Sample ID Columns', block_metrics_id_layout)
            block_metrics_form.addRow('Sample Filters', self.configure_block_domain_metrics_filters_btn)
            block_metrics_form.addRow('', self.block_domain_metrics_filters_summary)
            block_metrics_form.addRow('', self.start_block_domain_metrics_btn)
            operations_form.addRow(block_metrics_group)

            self.domain_samples_browse.clicked.connect(browse_domain_samples_output)
            self.start_domaining_btn.clicked.connect(self.run_domain_samples_only)
            self.block_domain_metrics_browse.clicked.connect(browse_block_domain_metrics_output)
            self.configure_block_domain_metrics_filters_btn.clicked.connect(self.configure_block_domain_metrics_filters)
            self.start_block_domain_metrics_btn.clicked.connect(self.run_block_domain_sample_metrics_only)
            self.samples_edit.textChanged.connect(lambda _: refresh_domain_samples_output_path())
            self.blocks_edit.textChanged.connect(lambda _: refresh_block_domain_metrics_output_path())
            self.block_domain_col.currentTextChanged.connect(lambda _: refresh_domain_samples_output_path())
            self.block_domain_col.currentTextChanged.connect(lambda _: refresh_block_domain_metrics_output_path())
            refresh_domain_samples_output_path(force=True)
            refresh_block_domain_metrics_output_path(force=True)
            self._update_block_domain_metrics_filters_summary()

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
            self.average_with_blocks = QtWidgets.QCheckBox(); self.average_with_blocks.setChecked(True)
            self.verbose = QtWidgets.QCheckBox(); self.verbose.setChecked(False)
            self.fill_unvisited_domainwise = QtWidgets.QCheckBox(); self.fill_unvisited_domainwise.setChecked(False)
            self.process_domains_sequentially = QtWidgets.QCheckBox(); self.process_domains_sequentially.setChecked(True)
            self.blank_sample_domain_behavior = QtWidgets.QComboBox(); self.blank_sample_domain_behavior.addItems(['Skip', 'Infer From Blocks'])
            self.blank_sample_domain_behavior.setCurrentText('Skip')
            self.blank_sample_domain_behavior.setToolTip('How to handle blank sample domains during domain-based interpolation or metrics.\nSkip: exclude blank-domain rows.\nInfer From Blocks: infer missing sample domains from the blocks model when possible.')

            # === DISPLAY TAB ===
            self.taichi_sample_diameter = dbl_spin(1.0, 0.001, 1e6, 0.1)
            self.taichi_sample_diameter.setToolTip('Sample diameter in model units for the Taichi mesh viewer. Default is 1 unit.')
            display_form.addRow('Sample Diameter', self.taichi_sample_diameter)

            advanced_form.addRow('Average With Blocks', self.average_with_blocks)
            advanced_form.addRow('Fill Unvisited (Domain-wise)', self.fill_unvisited_domainwise)
            advanced_form.addRow('Process Domains Sequentially', self.process_domains_sequentially)
            advanced_form.addRow('Blank Sample Domains', self.blank_sample_domain_behavior)
            advanced_form.addRow('Verbose', self.verbose)

            # Domain Algorithm Mapping
            self.domain_overrides = {}  # Store domain -> algorithm mapping
            self.domain_mapping_btn = QtWidgets.QPushButton('Configure Domain Algorithms...')
            self.domain_mapping_btn.clicked.connect(self.open_domain_mapping)
            advanced_form.addRow('Domain-Specific Algorithms', self.domain_mapping_btn)

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
            self._terminate_viewer_process()
            super().closeEvent(event)

        def _invalidate_domain_catalog_cache(self):
            self._domain_catalog_cache = None

        def to_dict(self, include_runtime_state=False):
            config = {
                'samples_file': self.samples_edit.text(),
                'blocks_file': self.blocks_edit.text(),
                'color_file': self.color_edit.text(),
                'interpolation_file': self.interp_edit.text(),
                'export_block_evaluated_samples': self.block_evaluated_samples_enabled.isChecked(),
                'block_evaluated_samples_file': self.block_evaluated_samples_edit.text(),
                'domain_samples_file': self.domain_samples_output_edit.text(),
                'block_domain_metrics_file': self.block_domain_metrics_output_edit.text(),
                'block_domain_metrics_closest_sample_id_cols': [cb.currentText() for cb in self.block_domain_metrics_id_cols],
                'block_domain_sample_filters': [dict(entry) for entry in self.block_domain_sample_filters],
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
                'block_x_col': self.block_x_col.currentText(),
                'block_y_col': self.block_y_col.currentText(),
                'block_z_col': self.block_z_col.currentText(),
                'block_domain_col': self.block_domain_col.currentText(),
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
                'fill_unvisited_domainwise': self.fill_unvisited_domainwise.isChecked(),
                'process_domains_sequentially': self.process_domains_sequentially.isChecked(),
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
            return config

        def from_dict(self, config):
            was_running = self._is_viewer_running()
            self._suspend_auto_viewer_refresh = True
            try:
                self.block_domain_sample_filters = []
                self._update_block_domain_metrics_filters_summary()
                if 'samples_file' in config: self.samples_edit.setText(config['samples_file'])
                if 'blocks_file' in config: self.blocks_edit.setText(config['blocks_file'])
                if 'color_file' in config: self.color_edit.setText(config['color_file'])
                if 'interpolation_file' in config: self.interp_edit.setText(config['interpolation_file'])
                if 'export_block_evaluated_samples' in config: self.block_evaluated_samples_enabled.setChecked(bool(config['export_block_evaluated_samples']))
                if 'block_evaluated_samples_file' in config: self.block_evaluated_samples_edit.setText(config['block_evaluated_samples_file'])
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
                if 'block_x_col' in config: self.block_x_col.setCurrentText(config['block_x_col'])
                if 'block_y_col' in config: self.block_y_col.setCurrentText(config['block_y_col'])
                if 'block_z_col' in config: self.block_z_col.setCurrentText(config['block_z_col'])
                if 'block_domain_col' in config: self.block_domain_col.setCurrentText(config['block_domain_col'])
                if 'domain_samples_file' in config: self.domain_samples_output_edit.setText(config['domain_samples_file'])
                if 'block_domain_metrics_file' in config: self.block_domain_metrics_output_edit.setText(config['block_domain_metrics_file'])
                if 'block_domain_metrics_closest_sample_id_cols' in config:
                    for cb, column_name in zip(self.block_domain_metrics_id_cols, config['block_domain_metrics_closest_sample_id_cols']):
                        cb.setCurrentText(column_name)
                if 'block_domain_sample_filters' in config:
                    self.block_domain_sample_filters = [dict(entry) for entry in config['block_domain_sample_filters']]
                    self._update_block_domain_metrics_filters_summary()
            
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
                if 'fill_unvisited_domainwise' in config: self.fill_unvisited_domainwise.setChecked(config['fill_unvisited_domainwise'])
                if 'process_domains_sequentially' in config: self.process_domains_sequentially.setChecked(config['process_domains_sequentially'])
                if 'blank_sample_domain_behavior' in config:
                    behavior = str(config['blank_sample_domain_behavior']).strip().lower()
                    self.blank_sample_domain_behavior.setCurrentText('Infer From Blocks' if behavior == 'infer_from_blocks' else 'Skip')
                if 'verbose' in config: self.verbose.setChecked(config['verbose'])
                if 'viewer_backend' in config: self.viewer_backend = _normalize_viewer_backend(config['viewer_backend'])
                if 'taichi_sample_diameter' in config: self.taichi_sample_diameter.setValue(config['taichi_sample_diameter'])
                if 'domain_algorithm_overrides' in config: self.domain_overrides = config['domain_algorithm_overrides']
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
                dialog.set_domain_configs(self.domain_overrides)
            
            if dialog.exec_() == QtWidgets.QDialog.Accepted:
                self.domain_overrides = dialog.get_domain_configs()
                count = len(self.domain_overrides)
                if count > 0:
                    self.domain_mapping_btn.setText(f'Configure Domain Algorithms... ({count} configured)')
                else:
                    self.domain_mapping_btn.setText('Configure Domain Algorithms...')

        def _load_samples_dataframe_for_operations(self):
            cfg = self.to_dict()
            samples_file = cfg.get('samples_file')
            if not samples_file or not os.path.isfile(samples_file):
                raise ValueError('Please select a valid samples file first.')
            df_samples, _ = load_full_samples_dataframe(
                samples_file,
                samples_delimiter=cfg.get('samples_delimiter'),
                samples_header_line=cfg.get('samples_header_line', 1),
                progress_label='Reading sample file',
            )
            return df_samples

        def _set_operation_buttons_enabled(self, enabled):
            self.start_domaining_btn.setEnabled(enabled)
            self.start_block_domain_metrics_btn.setEnabled(enabled)

        def _run_operation_with_progress(self, title, initial_message, operation, kwargs,
                                         success_handler, error_context):
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
                    progress.setLabelText(message)

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
                QtWidgets.QMessageBox.critical(self, 'Error', f'{error_context}:\n{error_text}')

            worker.progress.connect(handle_progress, QtCore.Qt.QueuedConnection)
            worker.finished.connect(handle_finished, QtCore.Qt.QueuedConnection)
            worker.failed.connect(handle_failed, QtCore.Qt.QueuedConnection)
            thread.started.connect(worker.run)
            thread.start()

        def _update_block_domain_metrics_filters_summary(self):
            count = len(self.block_domain_sample_filters)
            if count == 0:
                self.block_domain_metrics_filters_summary.setText('No sample filters configured. All samples will be used.')
                self.configure_block_domain_metrics_filters_btn.setText('Configure Filters...')
                return

            summaries = [summarize_sample_filter_spec(spec) for spec in self.block_domain_sample_filters]
            self.block_domain_metrics_filters_summary.setText('\n'.join(summaries[:4]))
            self.configure_block_domain_metrics_filters_btn.setText(f'Configure Filters... ({count} active)')

        def configure_block_domain_metrics_filters(self):
            cursor_set = False
            try:
                QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
                cursor_set = True
                df_samples = self._load_samples_dataframe_for_operations()
                QtWidgets.QApplication.restoreOverrideCursor()
                cursor_set = False

                dialog = SampleFiltersDialog(df_samples, filters=self.block_domain_sample_filters, parent=self)
                if dialog.exec_() == QtWidgets.QDialog.Accepted:
                    self.block_domain_sample_filters = dialog.get_filters()
                    self._update_block_domain_metrics_filters_summary()
            except Exception as exc:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Could not configure sample filters:\n{exc}')
            finally:
                if cursor_set:
                    QtWidgets.QApplication.restoreOverrideCursor()

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
                    },
                    handle_success,
                    'An error occurred during sample domaining',
                )
            except Exception as e:
                print(f"Error during sample domaining: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during sample domaining:\n{str(e)}')

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
                    closest_id_text = 'None'
                    if result.get('closest_sample_id_source_columns'):
                        closest_id_text = ', '.join(result['closest_sample_id_source_columns'])

                    QtWidgets.QMessageBox.information(
                        self,
                        'Success',
                        (
                            f"Block metrics export complete!\nResults saved to:\n{result['output_file']}\n\n"
                            f"Domain column: {result['domain_column']}\n"
                            f"Filtered samples: {result['filtered_samples']:,} of {result['input_samples']:,}\n"
                            f"Matched samples: {result['matched_samples']:,}\n"
                            f"Processed blocks: {result['processed_blocks']:,}\n"
                            f"Blocks with samples in domain: {result['blocks_with_samples_in_domain']:,}\n"
                            f"Invalid sample coordinates: {result['invalid_coordinate_samples']:,}\n"
                            f"Closest sample ID columns: {closest_id_text}\n\n"
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
                        'closest_sample_id_cols': cfg.get('block_domain_metrics_closest_sample_id_cols'),
                        'block_x_col': cfg.get('block_x_col'),
                        'block_y_col': cfg.get('block_y_col'),
                        'block_z_col': cfg.get('block_z_col'),
                        'block_domain_col': cfg.get('block_domain_col'),
                        'block_size': cfg.get('block_size'),
                        'sample_filters': cfg.get('block_domain_sample_filters'),
                        'blank_sample_domain_behavior': cfg.get('blank_sample_domain_behavior', 'skip'),
                    },
                    handle_success,
                    'An error occurred during block domain metrics export',
                )
            except Exception as e:
                print(f"Error during block domain metrics export: {e}")
                traceback.print_exc()
                QtWidgets.QMessageBox.critical(self, 'Error', f'An error occurred during block domain metrics export:\n{str(e)}')
        
        def run_interpolation_only(self):
            """Run interpolation without visualization"""
            try:
                cfg = self.to_dict(include_runtime_state=True)
                interpolation_file = cfg['interpolation_file']
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
                    progress_label='Reading sample file',
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
                    explicit_sample_map = None

                if needs_sample_domains:
                    df = normalize_selected_sample_domain_column(df, sample_domain_col=sample_domain_col)

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
                    sample_domains=sample_domains
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
                        
                        print(f"Running {algo_name1} for {iterations} iterations...")
                        pbar = tqdm(range(iterations), desc=f"Domain {domain} - Pass 1")
                        for i in pbar:
                            should_continue = interp1.run_iteration(dims)
                            if not should_continue or interp1.is_converged():
                                print(f"Pass 1 converged at iteration {i+1}")
                                break
                        
                        # Generate stats for Pass 1 if String Theory
                        output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                        if hasattr(interp1, 'generate_statistics'):
                            interp1.generate_statistics(output_dir, domain_name=f"{domain}_Pass1")

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
                            
                            # Determine if we should enforce domain mapping/grid restriction
                            use_mapping = False
                            if hasattr(interp2, 'allowed_grid_override') and interp2.allowed_grid_override is not None:
                                use_mapping = True
                            
                            interp2.initialize_blocks(pass1_values, dims, min_bounds, block_size, use_domain_mapping=use_mapping)
                            
                            if hasattr(interp2, 'create_ants'):
                                interp2.create_ants()
                            
                            print(f"Running {algo_name2} for {iterations} iterations...")
                            pbar = tqdm(range(iterations), desc=f"Domain {domain} - Pass 2")
                            for i in pbar:
                                should_continue = interp2.run_iteration(dims)
                                if not should_continue or interp2.is_converged():
                                    print(f"Pass 2 converged at iteration {i+1}")
                                    break
                            
                            # Generate stats for Pass 2 if String Theory
                            output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                            if hasattr(interp2, 'generate_statistics'):
                                interp2.generate_statistics(output_dir, domain_name=f"{domain}_Pass2")
            
                        # Print domain summary (of the last pass)
                        last_interp = interpolator_list[-1]
                        metadata = last_interp.get_metadata()
                        print(f"\n=== Domain {domain} Summary ===")
                        for key, value in metadata.items():
                            print(f"{key}: {value}")
                
                else:
                    # Single interpolator case
                    interpolator = blocks._ant_colony
                    algo_name = interpolator.get_algorithm_name()
                    print(f"Running {algo_name} for {iterations} iterations...")
                    
                    pbar = tqdm(range(iterations), desc=f"Interpolation ({algo_name})")
                    for i in pbar:
                        should_continue = interpolator.run_iteration(dims)
                        if not should_continue or interpolator.is_converged():
                            print(f"Converged at iteration {i+1}")
                            break
                    
                    metadata = interpolator.get_metadata()
                    print(f"\n=== Summary ===")
                    for key, value in metadata.items():
                        print(f"{key}: {value}")
                        
                    # Generate stats if String Theory
                    output_dir = os.path.dirname(interpolation_file) if interpolation_file else "."
                    if hasattr(interpolator, 'generate_statistics'):
                        interpolator.generate_statistics(output_dir, domain_name="Global")
                
                # Export results (handles both single and multiple interpolators)
                export_blocks_to_csv(blocks, interpolation_file)
                if block_evaluated_samples_file:
                    export_block_evaluated_samples_to_csv(blocks, block_evaluated_samples_file)
                print(f"Interpolation complete! Results saved to:\n  {interpolation_file}")
                if block_evaluated_samples_file:
                    print(f"Block-evaluated samples saved to:\n  {block_evaluated_samples_file}")
                print("=" * 60)
                
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
    app.aboutToQuit.connect(dialog._terminate_viewer_process)
    atexit.register(dialog._terminate_viewer_process)
    dialog.show()

    sys.exit(app.exec_())