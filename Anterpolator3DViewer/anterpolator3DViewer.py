import pandas as pd
import pyvista as pv
import numpy as np
import csv
from tqdm import tqdm
import sys
import os
from datetime import datetime
import xml.etree.ElementTree as ET
from matplotlib.colors import ListedColormap
import json
from PyQt5 import QtWidgets, QtCore
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
sys.path.append("C:/Projects/Anterpolator")
from ant_colony import AntColonyInterpolator
from molecular_clock_interpolator import MolecularClockInterpolator
from interpolator_base import InterpolatorBase

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
    
    # Check for domain-specific override
    if domain and 'domain_algorithm_overrides' in config:
        domain_config = config['domain_algorithm_overrides'].get(domain, {})
        if 'algorithm' in domain_config:
            algo_type = domain_config['algorithm']
            # Merge domain-specific parameters
            if algo_type == 'molecular_clock' and 'molecular_clock_params' in config:
                mc_params = {**config['molecular_clock_params'], **domain_config}
                config = {**config, 'molecular_clock_params': mc_params}
    
    if algo_type == 'ant_colony':
        return AntColonyInterpolator(
            range_size=config.get('range_size', 10),
            max_pheromone=config.get('max_pheromone', 150),
            ants_per_sample=config.get('ants_per_sample', 3),
            verbose=config.get('verbose', False),
            background_value=config.get('background_value', 0.0),
            background_distance=config.get('background_distance', None),
            average_with_blocks=config.get('average_with_blocks', False),
            avoid_visited_threshold_enabled=config.get('avoid_visited_threshold_enabled', False),
            avoid_visited_threshold=config.get('avoid_visited_threshold', 100)
        )
    
    elif algo_type == 'molecular_clock':
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

def read_csv_with_selected_header(path, delimiter, header_line, expected_min_cols=1):
    """Read CSV using a specific header line (1-based). Returns DataFrame.
    Uses manual header parsing to build names and then pandas read_csv with skiprows.
    """
    headers = parse_header_line(path, delimiter, header_line)
    # Build a name list; allow duplicate names by enumerating duplicates
    name_counts = {}
    final_names = []
    for h in headers:
        key = h if h else 'Unnamed'
        count = name_counts.get(key, 0)
        if count > 0:
            new_name = f"{key}_{count}"
        else:
            new_name = key
        name_counts[key] = count + 1
        final_names.append(new_name)
    df = pd.read_csv(
        path,
        delimiter=delimiter,
        engine='python',
        header=None,
        names=final_names,
        skiprows=header_line-1,
        comment='#'
    )
    # Drop the first data row if it actually contained the header again (defensive)
    if df.shape[0] and all(str(df.iloc[0, i]).strip() == final_names[i] for i in range(min(len(final_names), df.shape[1]))):
        df = df.iloc[1:].reset_index(drop=True)
    # Remove fully empty columns
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
                # Collect non-comment, non-empty for delimiter sensing
                if line.strip() and not line.lstrip().startswith(('#','//')):
                    lines.append(line)
                if len(lines) >= 3:
                    break
            sample = ''.join(lines)
    except StopIteration:
        sample = ''
    counts = {d: sample.count(d) for d in [',',';','\t','|']}
    # Prefer delimiter with highest count; fallback comma
    delim = max(counts, key=counts.get)
    return delim if counts[delim] > 0 else ','

def read_autodetect_csv(path, min_cols=1, forced_delimiter=None):
    """Read a CSV file detecting delimiter (comma, semicolon, tab, pipe) unless forced.
    Drops empty columns. Returns DataFrame."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")
    delim = forced_delimiter if forced_delimiter else detect_csv_delimiter(path)
    # Read, skipping comment-like lines
    def base_read(delimiter):
        return pd.read_csv(path, delimiter=delimiter, engine='python', comment='#')
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
    # If single column but delimiter chars appear inside, attempt multi-try
    if df.shape[1] == 1:
        text_sample = str(df.iloc[0,0]) if len(df) else ''
        candidates = ['; ',',','\t','|']
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

def middle_click_callback(plotter, event_id):
    if hasattr(plotter, 'picked_point'):
        point = plotter.picked_point
        if point is not None:
            # Update camera focus
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

    # Use a contiguous chunk to ensure we have neighbors
    chunk_size = min(len(points), 10000)
    chunk = points[:chunk_size]
    
    # 2. Find nearest neighbors in 2D (XY plane) to speed up
    # We only care about rotation around Z, so we look at XY vectors
    nbrs = NearestNeighbors(n_neighbors=8, algorithm='auto').fit(chunk[:, :2])
    distances, indices = nbrs.kneighbors(chunk[:, :2])
    
    # 3. Collect candidate vectors
    vectors = []
    for i in range(len(chunk)):
        # Skip self (index 0)
        for j in range(1, indices.shape[1]):
            neighbor_idx = indices[i][j]
            
            # Check if on same Z level
            if not np.isclose(chunk[i, 2], chunk[neighbor_idx, 2]):
                continue
                
            vec = chunk[neighbor_idx] - chunk[i]
            dist = distances[i][j]
            
            if dist > 0:
                vectors.append(vec)
    
    if not vectors:
        print("No vectors found between neighbors on same Z level.")
        return np.eye(3), np.zeros(3), False
        
    vectors = np.array(vectors)
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

def create_blocks(points, values, block_size=10, verbose=False, range_size=10, max_pheromone=150,
                  ants_per_sample=3, blocks_file=None, background_value=0.0, background_distance=None, average_with_blocks=False,
                  blocks_delimiter=None,
                  avoid_visited_threshold_enabled=False,
                  avoid_visited_threshold=100,
                  blocks_header_line=1,
                  block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                  config=None):
    if blocks_file is not None:
        print("Loading predefined cells from blocks_file...")
        if blocks_header_line and blocks_header_line != 1 and blocks_delimiter:
            # Use custom header line parsing
            df_blocks, parsed_cols = read_csv_with_selected_header(blocks_file, blocks_delimiter, blocks_header_line, expected_min_cols=3)
            print(f"Blocks file (custom header line {blocks_header_line}) parsed columns: {parsed_cols}")
        else:
            df_blocks = read_autodetect_csv(blocks_file, forced_delimiter=blocks_delimiter)
            print(f"Blocks file delimiter used: '{df_blocks._detected_delimiter}'")
        # If only one column but delimiter appears inside, attempt alternate reparse
        if df_blocks.shape[1] == 1:
            sole_col = df_blocks.columns[0]
            sample_val = str(df_blocks.iloc[0,0]) if len(df_blocks) else ''
            if (';' in sample_val and df_blocks._detected_delimiter != ';') or (',' in sample_val and df_blocks._detected_delimiter != ','):
                alt = ';' if df_blocks._detected_delimiter == ',' else ','
                try:
                    df_blocks = read_autodetect_csv(blocks_file, forced_delimiter=alt)
                    print(f"Reparsed blocks file with alternate delimiter '{alt}' -> columns: {list(df_blocks.columns)}")
                except Exception as e:
                    print(f"Alternate delimiter parse failed: {e}")
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
                    # Domain column could coincide; keep as Domain copy
                    if block_domain_col in df_blocks.columns:
                        df_blocks['Domain'] = df_blocks[block_domain_col]
                else:
                    if block_domain_col not in df_blocks.columns:
                        raise ValueError(f"Selected domain column '{block_domain_col}' not found in blocks file.")
                    rename_map[block_domain_col] = 'Domain'
            df_blocks = df_blocks.rename(columns=rename_map)
            print(f"Applied user block column mapping: {rename_map}")
        else:
            # Legacy auto-mapping fallback
            needed = ['x','y','z']
            if not all(n in df_blocks.columns for n in needed):
                cols = [c for c in df_blocks.columns if str(c).strip() != '']
                if len(cols) >= 3:
                    mapping = {}
                    for new, old in zip(['x','y','z'], cols[:3]):
                        mapping[old] = new
                    if len(cols) >= 4:
                        mapping[cols[3]] = 'Domain'
                    df_blocks = df_blocks.rename(columns=mapping)
                    print(f"Blocks header auto-mapped: {mapping}")
            # Generic positional mapping (first 4 columns -> x,y,z,Domain) if not already set
            cols = df_blocks.columns.tolist()
            if len(cols) >= 4:
                positional_targets = ['x', 'y', 'z', 'Domain']
                if any(cols[i] != positional_targets[i] for i in range(4)):
                    pos_map = {cols[i]: positional_targets[i] for i in range(4)}
                    df_blocks = df_blocks.rename(columns=pos_map)
                    print(f"Applied generic positional mapping to first four columns: {pos_map}")

        print(f"Blocks final headers: {list(df_blocks.columns)}")
        # Ensure coordinate columns are numeric
        coord_before = len(df_blocks)
        for c in ['x','y','z']:
            if c in df_blocks.columns:
                df_blocks[c] = pd.to_numeric(df_blocks[c], errors='coerce')
        # Drop rows with non-numeric coordinates
        df_blocks = df_blocks.dropna(subset=['x','y','z'])
        dropped = coord_before - len(df_blocks)
        if dropped > 0:
            print(f"Dropped {dropped} block rows with non-numeric coordinates.")
        if len(df_blocks) == 0:
            raise ValueError("All block rows have non-numeric coordinates after conversion.")
        missing_coords = [c for c in ['x','y','z'] if c not in df_blocks.columns]
        if missing_coords:
            raise ValueError(f"Blocks file missing required coordinate columns after mapping: {missing_coords}. Parsed headers: {list(df_blocks.columns)}")
        # Use ALL rows from df_blocks for boundaries
        centroids_all = df_blocks[['x','y','z']].values
        
        # Detect and apply rotation if needed
        rotation_matrix, rotation_center, is_rotated = detect_grid_rotation(centroids_all, block_size_hint=block_size)
        
        # Calculate azimuth of the Grid Y-axis (Row 1)
        # theta is angle of Grid Y-axis relative to World X-axis (East)
        theta_rad = np.arctan2(rotation_matrix[1, 1], rotation_matrix[1, 0])
        theta_deg = np.degrees(theta_rad)
        # Azimuth is clockwise from North (World Y)
        # Math angle is CCW from East (World X)
        # Azimuth = 90 - theta
        azimuth = (90 - theta_deg) % 360
        
        if is_rotated:
            print(f"Detected rotated block model (Azimuth: {azimuth:.2f}°). Applying rotation correction to align with axes.")
            print(f"Rotation center: {rotation_center}")
            print(f"Rotation matrix:\n{rotation_matrix}")
            
            # Rotate blocks
            # P_aligned = (P - Center) @ R.T
            # We use R.T because rotation_matrix rows are the new axes
            centroids_all = (centroids_all - rotation_center) @ rotation_matrix.T
            df_blocks['x'] = centroids_all[:, 0]
            df_blocks['y'] = centroids_all[:, 1]
            df_blocks['z'] = centroids_all[:, 2]
            
            # Rotate samples
            # Use in-place update so the caller sees the rotated points (for visualization consistency)
            points[:] = (points - rotation_center) @ rotation_matrix.T
            print("Applied rotation to sample points.")
            
        else:
            print(f"No significant rotation detected (Azimuth: {azimuth:.2f}°).")
            
        all_min_bounds = np.min(centroids_all, axis=0)
        all_max_bounds = np.max(centroids_all, axis=0)
        
        # Determine unified block dimensions
        # Prefer user-provided block_size if available, as auto-detection on rotated/float coordinates is unreliable
        if block_size is not None:
            if isinstance(block_size, (list, tuple, np.ndarray)):
                unified_dims = np.array(block_size, dtype=float)
            else:
                unified_dims = np.array([block_size, block_size, block_size], dtype=float)
            print(f"Using configured block size: {unified_dims}")
        else:
            # Fallback to auto-detection (legacy behavior)
            dims_all = []
            for col in ['x','y','z']:
                uniq = np.sort(df_blocks[col].unique())
                diffs = np.diff(uniq)
                positive_diffs = diffs[diffs > 0]
                dim = positive_diffs.min() if len(positive_diffs) > 0 else 10
                dims_all.append(dim)
            unified_dims = np.array(dims_all)
            print(f"Auto-detected block dimensions: {unified_dims}")

        dims_grid = np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int)
        print("Calculated grid dimensions:", dims_grid)
        # Build mapping from grid index to Domain using ALL rows from df_blocks (VECTORIZED)
        print("Building domain mapping...")
        centroids = df_blocks[['x', 'y', 'z']].values
        grid_indices = np.floor((centroids - all_min_bounds) / unified_dims + 1e-6).astype(int)
        
        # Handle domain column with default "Undomained"
        if 'Domain' in df_blocks.columns:
            domains = df_blocks['Domain'].fillna("Undomained").astype(str).str.strip()
            domains = domains.replace("", "Undomained")
        else:
            domains = pd.Series(["Undomained"] * len(df_blocks))
        
        # Create mapping dictionary
        domain_mapping = {tuple(idx): domain for idx, domain in zip(grid_indices, domains)}
        
        # Filter out skipped domains if config is provided
        if config and 'domain_algorithm_overrides' in config:
            skipped_domains = {domain for domain, cfg in config['domain_algorithm_overrides'].items() 
                             if cfg.get('skip', False)}
            if skipped_domains:
                print(f"Skipping domains: {skipped_domains}")
                domain_mapping = {idx: domain for idx, domain in domain_mapping.items() 
                                if domain not in skipped_domains}
        
        allowed_grid = set(domain_mapping.keys())
        
        # Count blocks per domain for debugging
        domain_block_counts = {}
        for idx, domain in domain_mapping.items():
            domain_block_counts[domain] = domain_block_counts.get(domain, 0) + 1
        
        print(f"\nDomain block distribution (from blocks file):")
        for domain, count in sorted(domain_block_counts.items()):
            print(f"  {domain}: {count} blocks")
        print(f"  Total allowed blocks: {len(allowed_grid)}")
        
        # Group sample points using the same all_min_bounds (VECTORIZED)
        print("Assigning points to blocks...")
        # Compute all block indices at once
        points_array = np.array(points)
        values_array = np.array(values)
        block_indices = np.floor((points_array - all_min_bounds) / unified_dims + 1e-6).astype(int)
        
        # Create lookup for allowed blocks
        sample_blocks_dict = {}
        block_values = {}
        
        # Group by block index
        for i in tqdm(range(len(points_array)), desc="Assigning points to blocks"):
            block_idx = tuple(block_indices[i])
            if block_idx in allowed_grid:
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
            'rotation_center': rotation_center if is_rotated else None
        }
        multiblock = pv.MultiBlock(block_data)
        # Store metadata on multiblock with private-style names to avoid PyVista attribute restrictions
        multiblock._block_info = block_info
        multiblock._sample_blocks = {idx: np.mean(vals) for idx, vals in block_values.items()}
        
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
            
            # Store multiple interpolators
            multiblock._interpolators = {}
            
            for domain, domain_samples in samples_by_domain.items():
                # Determine which algorithm for this domain
                domain_config = config['domain_algorithm_overrides'].get(domain, {})
                if domain_config.get('skip', False):
                    print(f"  Skipping domain: {domain}")
                    continue
                
                algo = domain_config.get('algorithm', config.get('algorithm', 'ant_colony'))
                print(f"  Domain '{domain}': {len(domain_samples)} samples, algorithm={algo}")
                
                # Create config for this domain
                domain_cfg = config.copy()
                domain_cfg['algorithm'] = algo
                
                # Create interpolator for this domain
                interpolator = create_interpolator(domain_cfg, domain=domain)
                
                # Filter allowed_grid to this domain only
                domain_allowed_grid = {pos for pos in allowed_grid if domain_mapping.get(pos) == domain}
                
                # Attach domain-specific attributes
                interpolator.allowed_grid_override = domain_allowed_grid
                interpolator.domain_mapping = {pos: domain for pos in domain_allowed_grid}
                
                # Initialize with only this domain's samples
                interpolator.initialize_blocks(domain_samples, tuple(block_info['dims']),
                                             all_min_bounds, unified_dims.tolist(), use_domain_mapping=True)
                
                if hasattr(interpolator, 'create_ants'):
                    interpolator.create_ants()
                
                multiblock._interpolators[domain] = interpolator
            
            # Keep reference to primary interpolator (for compatibility)
            if multiblock._interpolators:
                multiblock._ant_colony = list(multiblock._interpolators.values())[0]
            
        else:
            # Original single interpolator approach
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
            
            interpolator = create_interpolator(config)
            
            # Attach domain-specific attributes BEFORE initialize_blocks
            # For molecular clock, don't filter by allowed_grid (let it create blocks anywhere)
            # For ant colony, use allowed_grid to restrict movement
            algo_type = config.get('algorithm', 'ant_colony')
            if algo_type == 'ant_colony':
                interpolator.allowed_grid_override = allowed_grid
                interpolator.domain_mapping = domain_mapping
            elif algo_type == 'molecular_clock':
                # Enable domain sensitivity for molecular clock
                interpolator.allowed_grid_override = allowed_grid
                interpolator.domain_mapping = domain_mapping
            
            # Use use_domain_mapping=True to respect allowed_grid_override (geometry)
            interpolator.initialize_blocks(multiblock._sample_blocks, tuple(block_info['dims']),
                                         all_min_bounds, unified_dims.tolist(), use_domain_mapping=True)
            
            if hasattr(interpolator, 'create_ants'):
                interpolator.create_ants()
            
            multiblock._ant_colony = interpolator
        
        return multiblock
    else:
        print("Assigning points to blocks...")
        min_bounds = np.min(points, axis=0)
        max_bounds = np.max(points, axis=0)
        dims = np.ceil((max_bounds - min_bounds) / np.array(block_size)).astype(int)
        blocks = {}
        block_values = {}
        block_info = {
            'min_bounds': min_bounds,
            'dims': dims,
            'block_size': block_size
        }
        # Vectorized block assignment
        points_array = np.array(points)
        values_array = np.array(values)
        block_indices = ((points_array - min_bounds) // np.array(block_size)).astype(int)
        
        for i in tqdm(range(len(points_array)), desc="Assigning points to blocks"):
            block_idx = tuple(block_indices[i])
            if block_idx not in blocks:
                blocks[block_idx] = []
                block_values[block_idx] = []
            blocks[block_idx].append(points_array[i])
            block_values[block_idx].append(values_array[i])
        block_data = []
        next_block_id = 1
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
            block_data.append(cell)
        sample_blocks = {idx: np.mean(vals) for idx, vals in block_values.items()}
        multiblock = pv.MultiBlock(block_data)
        multiblock._block_info = block_info
        multiblock._sample_blocks = sample_blocks
        
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
        
        interpolator = create_interpolator(config)
        interpolator.initialize_blocks(sample_blocks, tuple(block_info['dims']),
                                       min_bounds, block_size, use_domain_mapping=False)
        
        if hasattr(interpolator, 'create_ants'):
            interpolator.create_ants()
        
        multiblock._ant_colony = interpolator
        return multiblock

def toggle_blocks(plotter):
    if hasattr(plotter, '_blocks_actor'):
        is_visible = not plotter._blocks_actor.GetVisibility()
        plotter._blocks_actor.SetVisibility(is_visible)
        plotter.render()

def update_interpolation(plotter):
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
    
    # Helper for classification of continuous to discrete indices
    def classify(val_array, bins, n_colors):
        bins_np = np.array(bins)
        threshold_style = (n_colors == len(bins_np) + 1)
        if threshold_style:
            idx = np.searchsorted(bins_np, val_array, side='right')
        else:
            if len(bins_np) == n_colors + 1:
                idx = np.digitize(val_array, bins_np, right=False) - 1
            else:
                idx = np.digitize(val_array, bins_np, right=False) - 1
        return np.clip(idx, 0, n_colors-1)

    # Move ants and track if changes were made
    print("Moving ants...")
    changes_made = False
    
    # Handle multiple interpolators (Scenario 3) or single (Scenario 1 & 2)
    interpolators = []
    if hasattr(blocks, '_interpolators') and blocks._interpolators:
        interpolators = list(blocks._interpolators.values())
    else:
        interpolators = [interpolator]
        
    for interp in interpolators:
        if interp.run_iteration(dims):
            changes_made = True
        
        # Optional domain-wise fill of unvisited blocks (ant colony specific)
        if getattr(plotter, '_fill_unvisited_domainwise', False):
            try:
                if hasattr(interp, 'fill_unvisited_blocks_domainwise'):
                    created, assigned = interp.fill_unvisited_blocks_domainwise(dims)
                    if (created or assigned):
                        changes_made = True
            except Exception as e:
                print(f"Domain-wise fill error: {e}")
    
    if changes_made:
        colony_blocks = {}
        # Map positions to their source interpolator to use correct next_block_id
        pos_to_interpolator = {}
        
        for interp in interpolators:
            vals = interp.get_interpolated_values()
            colony_blocks.update(vals)
            for pos in vals:
                pos_to_interpolator[pos] = interp
        
        # Track blocks by both position and ID
        block_mapping = {}
        id_mapping = {}
        
        # First pass: map existing blocks and track IDs - no value filtering here
        for i, block in enumerate(blocks):
            corner = block.bounds[::2]
            pos = tuple(np.floor((corner - min_bounds) / np.array(block_size)).astype(int))
            if 'Block_ID' in block.cell_data:
                block_id = block.cell_data['Block_ID'][0]
            else:
                block_id = -1
            
            if pos in blocks._sample_blocks:
                block_mapping[pos] = (i, block)
                id_mapping[block_id] = pos
            else:
                # Keep track of all non-sample blocks regardless of value
                block_mapping[pos] = (i, block)
        
        # Process blocks
        new_blocks = []
        modified = False
        
        # Update blocks - create all blocks, filter only for display
        for pos, value in colony_blocks.items():
            try:
                if pos in blocks._sample_blocks:
                    continue
                
                if pos not in block_mapping:
                    # Create new interpolated block regardless of value
                    corner = min_bounds + np.array(pos) * np.array(block_size)
                    half_size = np.array(block_size) / 2
                    center = corner + half_size
                    
                    new_block = pv.Box(bounds=(
                        center[0] - half_size[0]/2, center[0] + half_size[0]/2,
                        center[1] - half_size[1]/2, center[1] + half_size[1]/2,
                        center[2] - half_size[2]/2, center[2] + half_size[2]/2
                    ))
                    if getattr(plotter, '_value_is_indexed', False) and hasattr(plotter, '_lfc_bins'):
                        raw_val = value
                        if isinstance(plotter._colormap, ListedColormap):
                            n_colors = len(plotter._colormap.colors)
                            idx_val = classify(np.array([raw_val]), plotter._lfc_bins, n_colors)[0]
                            new_block.cell_data['Raw_Value'] = np.full(new_block.n_cells, raw_val)
                            new_block.cell_data['Value'] = np.full(new_block.n_cells, idx_val)
                        else:
                            new_block.cell_data['Value'] = np.full(new_block.n_cells, raw_val)
                    else:
                        new_block.cell_data['Value'] = np.full(new_block.n_cells, value)
                    new_block.cell_data['Is_Sample'] = np.full(new_block.n_cells, False)
                    
                    target_interp = pos_to_interpolator.get(pos, interpolator)
                    new_block.cell_data['Block_ID'] = np.full(new_block.n_cells, target_interp.next_block_id)
                    target_interp.next_block_id += 1
                    
                    # Set domain if available
                    if hasattr(target_interp, 'domain_mapping'):
                        domain = target_interp.domain_mapping.get(pos, "Undomained")
                        new_block.cell_data['Domain'] = np.full(new_block.n_cells, domain)
                    
                    new_blocks.append(new_block)
                elif pos in block_mapping:
                    _, block = block_mapping[pos]
                    base_val = block.cell_data['Raw_Value'][0] if 'Raw_Value' in block.cell_data else block.cell_data['Value'][0]
                    if abs(base_val - value) > 0.0001:
                        if getattr(plotter, '_value_is_indexed', False) and hasattr(plotter, '_lfc_bins') and isinstance(plotter._colormap, ListedColormap):
                            n_colors = len(plotter._colormap.colors)
                            idx_val = classify(np.array([value]), plotter._lfc_bins, n_colors)[0]
                            block.cell_data['Raw_Value'][:] = value
                            block.cell_data['Value'][:] = idx_val
                        else:
                            block.cell_data['Value'][:] = value
                        modified = True
            except Exception as e:
                print(f"Error processing position {pos}: {str(e)}")
                continue

        # Add all new blocks to the multiblock
        if new_blocks:
            print(f"Adding {len(new_blocks)} new blocks...")
            for block in new_blocks:
                blocks.append(block)
        
        # Filter blocks for display only
        visible_blocks = pv.MultiBlock()
        for block in blocks:
            if block.cell_data['Is_Sample'][0]:  # Always show sample blocks
                visible_blocks.append(block)
            else:
                raw_val = block.cell_data['Raw_Value'][0] if 'Raw_Value' in block.cell_data else block.cell_data['Value'][0]
                if raw_val >= plotter._value_filter:
                    visible_blocks.append(block)
        
        # Use the colormap stored in the plotter
        colormap = getattr(plotter, '_colormap', 'rainbow')

        # Create fresh mesh actor with filtered blocks
        plotter.remove_actor(plotter._blocks_actor)
        if getattr(plotter, '_value_is_indexed', False) and hasattr(plotter, '_lfc_bins') and isinstance(colormap, ListedColormap):
            annotations = {i: plotter._lfc_tick_labels[i] for i in range(len(plotter._lfc_tick_labels))} if hasattr(plotter, '_lfc_tick_labels') else None
            clim = [-0.5, len(colormap.colors)-0.5]
            plotter._blocks_actor = plotter.add_mesh(
                visible_blocks,
                style='surface',
                scalars='Value',
                opacity=0.5,
                show_edges=True,
                cmap=colormap,
                categories=True,
                annotations=annotations,
                clim=clim,
            )
        else:
            clim = [0, max(colony_blocks.values())]
            plotter._blocks_actor = plotter.add_mesh(
                visible_blocks,
                style='surface',
                scalars='Value',
                opacity=0.5,
                show_edges=True,
                cmap=colormap,
                clim=clim,
            )
        plotter.render()
        print("Visualization updated")
    else:
        print("No blocks to update")

def export_blocks_to_csv(blocks, filepath):
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    print(f"Exporting blocks to {filepath}...")
    
    # Prepare data for export
    data = []
    min_bounds = blocks._block_info['min_bounds']
    block_size = blocks._block_info['block_size']
    rotation_matrix = blocks._block_info.get('rotation_matrix')
    rotation_center = blocks._block_info.get('rotation_center')
    
    # Check if we have multiple interpolators (sequential domain processing)
    if hasattr(blocks, '_interpolators'):
        print(f"  Processing {len(blocks._interpolators)} domain interpolators...")
        for domain, interpolator in blocks._interpolators.items():
            print(f"    Exporting domain: {domain} ({len(interpolator.blocks)} blocks)")
            _add_interpolator_blocks_to_data(interpolator, min_bounds, block_size, data, rotation_matrix, rotation_center)
    else:
        # Single interpolator
        interpolator = blocks._ant_colony
        # Pass domain_mapping from block_info to restore original domains
        domain_mapping = blocks._block_info.get('domain_mapping')
        _add_interpolator_blocks_to_data(interpolator, min_bounds, block_size, data, rotation_matrix, rotation_center, domain_mapping=domain_mapping)
    
    # Export to CSV
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Exported {len(data)} blocks to {filepath}")

def _add_interpolator_blocks_to_data(interpolator, min_bounds, block_size, data, rotation_matrix=None, rotation_center=None, domain_mapping=None):
    """Process blocks from an interpolator and add to data list"""
    # Determine algorithm type
    algo_name = interpolator.get_algorithm_name()
    if "Ant Colony" in algo_name:
        algo_type = "ant_colony"
    elif "Biochemical Clock" in algo_name or "Phylogeographic" in algo_name:
        algo_type = "bio_clock"
    else:
        algo_type = "unknown"
    
    for pos, block in tqdm(interpolator.blocks.items(), desc="Processing blocks"):
        # Calculate block centroid - grid indices are calculated relative to centroids in min_bounds,
        # so we don't add block_size/2 (that would shift by half a block)
        centroid = min_bounds + np.array(pos) * np.array(block_size)
        
        # Apply inverse rotation if needed
        if rotation_matrix is not None and rotation_center is not None:
            # P_orig = P_aligned @ R + Center
            centroid = centroid @ rotation_matrix + rotation_center
        
        # Initialize common fields with None/NaN for all possible columns
        row = {
            'x': centroid[0],
            'y': centroid[1],
            'z': centroid[2],
            'Algo_Type': algo_type,
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
            
            # Determine domain: use mapping if provided (for single interpolator mode), else use block's domain
            domain = block.domain
            if domain_mapping and pos in domain_mapping:
                domain = domain_mapping[pos]
            
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
            domain = block.get('domain', 'Undomained')
            if domain_mapping and pos in domain_mapping:
                domain = domain_mapping[pos]
            
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
    
    # Check if we have multiple interpolators (sequential domain processing)
    if hasattr(blocks, '_interpolators'):
        print(f"Running sequential domain interpolation for {len(blocks._interpolators)} domains...")
        
        for domain_idx, (domain, interpolator) in enumerate(blocks._interpolators.items(), 1):
            algo_name = interpolator.get_algorithm_name()
            print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} ({algo_name}) ===")
            
            # Force verbose for first iteration if it's an AntColony
            if hasattr(interpolator, 'verbose'):
                original_verbose = interpolator.verbose
                interpolator.verbose = True
            
            pbar = tqdm(range(iterations), desc=f"Domain {domain} ({algo_name})")
            for i in pbar:
                # Use the base class interface
                should_continue = interpolator.run_iteration(dims)
                
                # Restore verbose after first iteration
                if i == 0 and hasattr(interpolator, 'verbose'):
                    interpolator.verbose = original_verbose
                    # Print block count after first iteration
                    metadata = interpolator.get_metadata()
                    print(f"\nAfter iteration 1:")
                    print(f"  Total blocks: {metadata.get('total_blocks', 0)}")
                    print(f"  Sample blocks: {metadata.get('sample_blocks', 0)}")
                    print(f"  Interpolated blocks: {metadata.get('interpolated_blocks', 0)}")
                    print(f"  Remaining ants: {metadata.get('remaining_ants', 0)}")
                
                # Check for early termination
                if not should_continue or interpolator.is_converged():
                    pbar.set_description(f"Domain {domain} (converged)")
                    print(f"Early stop: algorithm converged. Iterations run: {i+1}/{iterations}")
                    break
            
            # Optionally fill domain-wise after silent run (ant colony specific)
            if getattr(plotter, '_fill_unvisited_domainwise', False):
                try:
                    if hasattr(interpolator, 'fill_unvisited_blocks_domainwise'):
                        interpolator.fill_unvisited_blocks_domainwise(dims)
                except Exception as e:
                    print(f"Domain-wise fill error: {e}")
            
            # Print domain metadata
            metadata = interpolator.get_metadata()
            print(f"\n=== Domain {domain} Summary ===")
            for key, value in metadata.items():
                print(f"{key}: {value}")
    else:
        # Single interpolator
        interpolator = blocks._ant_colony
        algo_name = interpolator.get_algorithm_name()
        print(f"Running {algo_name} interpolation...")
        
        # Force verbose for first iteration if it's an AntColony
        if hasattr(interpolator, 'verbose'):
            original_verbose = interpolator.verbose
            interpolator.verbose = True
        
        pbar = tqdm(range(iterations), desc=f"Interpolating ({algo_name})")
        for i in pbar:
            # Use the base class interface
            should_continue = interpolator.run_iteration(dims)
            
            # Restore verbose after first iteration
            if i == 0 and hasattr(interpolator, 'verbose'):
                interpolator.verbose = original_verbose
                # Print block count after first iteration
                metadata = interpolator.get_metadata()
                print(f"\nAfter iteration 1:")
                print(f"  Total blocks: {metadata.get('total_blocks', 0)}")
                print(f"  Sample blocks: {metadata.get('sample_blocks', 0)}")
                print(f"  Interpolated blocks: {metadata.get('interpolated_blocks', 0)}")
                print(f"  Remaining ants: {metadata.get('remaining_ants', 0)}")
            
            # Check for early termination
            if not should_continue or interpolator.is_converged():
                pbar.set_description(f"Interpolating (converged)")
                print(f"Early stop: algorithm converged. Iterations run: {i+1}/{iterations}")
                break
        
        # Optionally fill domain-wise after silent run (ant colony specific)
        if getattr(interpolator, '_fill_unvisited_domainwise', False):
            try:
                if hasattr(interpolator, 'fill_unvisited_blocks_domainwise'):
                    interpolator.fill_unvisited_blocks_domainwise(dims)
            except Exception as e:
                print(f"Domain-wise fill error: {e}")
        
        # Print metadata
        metadata = interpolator.get_metadata()
        print(f"\n=== Interpolation Summary ===")
        for key, value in metadata.items():
            print(f"{key}: {value}")
    
    # Export results (handles both single and multiple interpolators)
    export_blocks_to_csv(blocks, interpolation_file)

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

def load_and_visualize_samples(samples_file, block_size=10, value_filter=60, verbose=False, iterations=100, range_size=10, max_pheromone=150, ants_per_sample=3, blocks_file=None, color_file=None, background_value=0.0, background_distance=None, average_with_blocks=False,
                               samples_delimiter=None, blocks_delimiter=None, fill_unvisited_domainwise=False,
                               avoid_visited_threshold_enabled=False,
                               avoid_visited_threshold=100,
                               samples_header_line=1,
                               sample_x_col=None, sample_y_col=None, sample_z_col=None, sample_value_col=None,
                               blocks_header_line=1,
                               block_x_col=None, block_y_col=None, block_z_col=None, block_domain_col=None,
                               config=None):
    try:
        # Display message about loading the sample file
        print(f"Loading sample file from {samples_file}...")
        
        # Read CSV
        if samples_header_line and samples_header_line != 1 and samples_delimiter:
            df, parsed_cols = read_csv_with_selected_header(samples_file, samples_delimiter, samples_header_line, expected_min_cols=4)
            print(f"Samples file (custom header line {samples_header_line}) parsed columns: {parsed_cols}")
        else:
            df = read_autodetect_csv(samples_file, forced_delimiter=samples_delimiter)
            print(f"Samples file delimiter used: '{df._detected_delimiter}'")

        # Apply explicit mapping if provided
        if sample_x_col and sample_y_col and sample_z_col and sample_value_col:
            for chosen in [sample_x_col, sample_y_col, sample_z_col, sample_value_col]:
                if chosen not in df.columns:
                    raise ValueError(f"Selected samples column '{chosen}' not present in file.")
            rename_map = {sample_x_col: 'x', sample_y_col: 'y', sample_z_col: 'z', sample_value_col: 'Value'}
            df = df.rename(columns=rename_map)
            print(f"Applied user sample column mapping: {rename_map}")
        else:
            # Legacy fallback auto mapping
            expected = ['x','y','z','Value']
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
        print(f"Samples final headers: {list(df.columns)}")
        # --- Ensure the selected Value column is strictly honored and blanks handled ---
        # Coerce Value to numeric; blank strings or invalid entries become NaN
        try:
            df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
        except Exception as e:
            print(f"Warning: could not coerce 'Value' column to numeric directly: {e}")
        # Count and drop rows with NaN Value (true blanks)
        nan_before = df['Value'].isna().sum()
        if nan_before:
            print(f"Detected {nan_before} sample rows with blank/non-numeric values in the selected value column; these will be excluded from samples.")
            df = df.dropna(subset=['Value'])
        if len(df) == 0:
            raise ValueError("After removing blank/non-numeric sample values, no samples remain. Check the configured 'sample_value_col'.")
        # Explicitly report the effective value range post-cleaning
        print(f"Cleaned sample count: {len(df)} (removed {nan_before} blank value rows)")
        points = df[['x','y','z']].values
        values = df['Value'].values
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
            config=config
        )

        # Store settings in plotter
        plotter._blocks_data = blocks
        plotter._verbose = verbose
        plotter._value_filter = value_filter  # Store value_filter
        plotter._fill_unvisited_domainwise = fill_unvisited_domainwise
        plotter._avoid_visited_threshold_enabled = avoid_visited_threshold_enabled
        plotter._avoid_visited_threshold = avoid_visited_threshold

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

        # Add blocks without scalar bar
        if getattr(plotter, '_value_is_indexed', False):
            annotations = {i: plotter._lfc_tick_labels[i] for i in range(len(plotter._lfc_tick_labels))}
            plotter._blocks_actor = plotter.add_mesh(
                blocks,
                style='surface',
                scalars='Value',
                opacity=0.5,
                show_edges=True,
                cmap=plotter._colormap,
                categories=True,
                annotations=annotations,
                clim=[-0.5, len(plotter._colormap.colors)-0.5],
                show_scalar_bar=False,
            )
        else:
            plotter._blocks_actor = plotter.add_mesh(
                blocks,
                style='surface',
                scalars='Value',
                opacity=0.5,
                show_edges=True,
                cmap=plotter._colormap,
                clim=[data_min, data_max],
                show_scalar_bar=False,
            )
        # Persist range on plotter for later updates
        plotter._scalar_range = (data_min, data_max)

        # If discrete mode active, classify existing sample blocks too
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
                    if len(bins_np) == n_colors + 1:
                        idx_val = np.digitize([raw_val], bins_np, right=False)[0] - 1
                    else:
                        idx_val = np.digitize([raw_val], bins_np, right=False)[0] - 1
                idx_val = max(0, min(idx_val, n_colors-1))
                block.cell_data['Value'][:] = idx_val

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
    class DomainAlgorithmDialog(QtWidgets.QDialog):
        """Dialog for configuring algorithm per domain"""
        def __init__(self, blocks_file, blocks_delimiter, blocks_header_line, block_domain_col, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Domain Algorithm Mapping")
            self.resize(700, 500)
            
            self.domain_configs = {}  # domain -> {'algorithm': str, 'skip': bool, 'params': dict}
            
            layout = QtWidgets.QVBoxLayout()
            self.setLayout(layout)
            
            # Info label
            info = QtWidgets.QLabel("Configure which algorithm to use for each domain.\nYou can also skip domains to exclude them from interpolation.")
            info.setWordWrap(True)
            layout.addWidget(info)
            
            # Table for domain mappings
            self.table = QtWidgets.QTableWidget()
            self.table.setColumnCount(3)
            self.table.setHorizontalHeaderLabels(['Domain', 'Algorithm', 'Skip Domain'])
            self.table.horizontalHeader().setStretchLastSection(False)
            self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
            self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
            self.table.horizontalHeader().setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeToContents)
            layout.addWidget(self.table)
            
            # Load domains from blocks file
            self.load_domains(blocks_file, blocks_delimiter, blocks_header_line, block_domain_col)
            
            # Buttons
            btn_layout = QtWidgets.QHBoxLayout()
            self.apply_all_btn = QtWidgets.QPushButton('Apply Algorithm to All')
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
        
        def load_domains(self, blocks_file, delimiter, header_line, domain_col):
            """Load unique domains from blocks file"""
            domains = set()
            
            if not blocks_file or not os.path.isfile(blocks_file):
                QtWidgets.QMessageBox.warning(self, 'Warning', 
                    'Blocks file not found. Please select a valid blocks file first.')
                return
            
            if not domain_col or domain_col == '(None)':
                QtWidgets.QMessageBox.warning(self, 'Warning',
                    'No domain column selected. Please select a domain column in the main dialog first.')
                return
            
            try:
                # Read domains from file
                df, _ = read_csv_with_selected_header(blocks_file, delimiter, header_line, expected_min_cols=1)
                
                if domain_col not in df.columns:
                    QtWidgets.QMessageBox.warning(self, 'Warning',
                        f'Domain column "{domain_col}" not found in blocks file.')
                    return
                
                # Get unique domains
                domains = set(df[domain_col].dropna().unique())
                domains = {str(d).strip() for d in domains if str(d).strip() and str(d).strip().lower() != 'nan'}
                
                if not domains:
                    QtWidgets.QMessageBox.warning(self, 'Warning', 'No domains found in blocks file.')
                    return
                
                # Populate table
                self.table.setRowCount(len(domains))
                for i, domain in enumerate(sorted(domains)):
                    # Domain name (read-only)
                    domain_item = QtWidgets.QTableWidgetItem(domain)
                    domain_item.setFlags(domain_item.flags() & ~QtCore.Qt.ItemIsEditable)
                    self.table.setItem(i, 0, domain_item)
                    
                    # Algorithm selector
                    algo_combo = QtWidgets.QComboBox()
                    algo_combo.addItems(['(use default)', 'ant_colony', 'molecular_clock', 'skip'])
                    algo_combo.setCurrentText('(use default)')
                    self.table.setCellWidget(i, 1, algo_combo)
                    
                    # Skip checkbox
                    skip_check = QtWidgets.QCheckBox()
                    skip_check.setEnabled(False)  # Controlled by algorithm selection
                    skip_widget = QtWidgets.QWidget()
                    skip_layout = QtWidgets.QHBoxLayout(skip_widget)
                    skip_layout.addWidget(skip_check)
                    skip_layout.setAlignment(QtCore.Qt.AlignCenter)
                    skip_layout.setContentsMargins(0, 0, 0, 0)
                    self.table.setCellWidget(i, 2, skip_widget)
                    
                    # Connect algorithm change to update skip checkbox
                    algo_combo.currentTextChanged.connect(
                        lambda text, row=i: self.update_skip_state(row, text)
                    )
                    
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Failed to load domains: {e}')
        
        def update_skip_state(self, row, algorithm):
            """Update skip checkbox based on algorithm selection"""
            skip_widget = self.table.cellWidget(row, 2)
            if skip_widget:
                skip_check = skip_widget.findChild(QtWidgets.QCheckBox)
                if skip_check:
                    if algorithm == 'skip':
                        skip_check.setChecked(True)
                        skip_check.setEnabled(False)
                    else:
                        skip_check.setChecked(False)
                        skip_check.setEnabled(False)
        
        def apply_to_all(self):
            """Apply same algorithm to all domains"""
            algorithms = ['(use default)', 'ant_colony', 'molecular_clock', 'skip']
            algo, ok = QtWidgets.QInputDialog.getItem(
                self, 'Apply to All', 'Select algorithm for all domains:', 
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
                algo_combo = self.table.cellWidget(i, 1)
                
                if domain_item and algo_combo:
                    domain = domain_item.text()
                    algorithm = algo_combo.currentText()
                    
                    if algorithm != '(use default)':
                        if algorithm == 'skip':
                            configs[domain] = {'skip': True}
                        else:
                            configs[domain] = {'algorithm': algorithm}
            
            return configs
        
        def set_domain_configs(self, configs):
            """Set domain algorithm configurations from loaded config"""
            for i in range(self.table.rowCount()):
                domain_item = self.table.item(i, 0)
                algo_combo = self.table.cellWidget(i, 1)
                
                if domain_item and algo_combo:
                    domain = domain_item.text()
                    if domain in configs:
                        config = configs[domain]
                        if config.get('skip', False):
                            algo_combo.setCurrentText('skip')
                        elif 'algorithm' in config:
                            algo_combo.setCurrentText(config['algorithm'])

    class ConfigDialog(QtWidgets.QDialog):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("Anterpolator 3D Viewer Configuration")
            self.resize(700, 600)
            
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
            
            # Tab 4: Advanced Options
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

            # Algorithm selection
            self.algorithm_combo = QtWidgets.QComboBox()
            self.algorithm_combo.addItems(['ant_colony', 'molecular_clock'])
            files_form.addRow('Algorithm', self.algorithm_combo)

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

            # Column mapping combo boxes for samples
            self.sample_x_col = QtWidgets.QComboBox(); self.sample_y_col = QtWidgets.QComboBox(); self.sample_z_col = QtWidgets.QComboBox(); self.sample_value_col = QtWidgets.QComboBox()
            for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col]:
                cb.setEditable(False)
            sample_map_layout = QtWidgets.QHBoxLayout(); sample_map_layout.addWidget(self.sample_x_col); sample_map_layout.addWidget(self.sample_y_col); sample_map_layout.addWidget(self.sample_z_col); sample_map_layout.addWidget(self.sample_value_col)
            files_form.addRow('Samples Columns (X Y Z Value)', sample_map_layout)

            # Column mapping combo boxes for blocks
            self.block_x_col = QtWidgets.QComboBox(); self.block_y_col = QtWidgets.QComboBox(); self.block_z_col = QtWidgets.QComboBox(); self.block_domain_col = QtWidgets.QComboBox()
            for cb in [self.block_x_col, self.block_y_col, self.block_z_col, self.block_domain_col]:
                cb.setEditable(False)
            self.block_domain_col.addItem('(None)')
            block_map_layout = QtWidgets.QHBoxLayout(); block_map_layout.addWidget(self.block_x_col); block_map_layout.addWidget(self.block_y_col); block_map_layout.addWidget(self.block_z_col); block_map_layout.addWidget(self.block_domain_col)
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
                for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col]:
                    cb.clear()
                if not os.path.isfile(path):
                    return
                try:
                    cols = parse_header_line(path, delim, header_line)
                    for cb in [self.sample_x_col, self.sample_y_col, self.sample_z_col, self.sample_value_col]:
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

            # Initial refresh (silent if files missing)
            refresh_sample_columns()
            refresh_block_columns()

            # === ANT COLONY TAB ===
            def dbl_spin(default, minv, maxv, step=0.1):
                s = QtWidgets.QDoubleSpinBox(); s.setRange(minv, maxv); s.setSingleStep(step); s.setValue(default); return s
            def int_spin(default, minv, maxv, step=1):
                s = QtWidgets.QSpinBox(); s.setRange(minv, maxv); s.setSingleStep(step); s.setValue(default); return s

            self.range_size = dbl_spin(0.2, 0.0001, 1e6, 0.1)
            self.max_pheromone = int_spin(1000, 1, 10_000_000, 10)
            self.ants_per_sample = int_spin(16, 1, 10000, 1)
            self.iterations = int_spin(500, 1, 10_000_000, 50)
            self.background_value = dbl_spin(0.0, -1e9, 1e9, 0.1)
            self.background_distance = dbl_spin(32.0, 0.0, 1e9, 1.0)
            self.value_filter = dbl_spin(0.0, -1e9, 1e9, 1.0)
            self.avoid_visited_enabled = QtWidgets.QCheckBox(); self.avoid_visited_enabled.setChecked(False)
            self.avoid_visited_threshold = QtWidgets.QSpinBox(); self.avoid_visited_threshold.setRange(1, 10_000_000); self.avoid_visited_threshold.setValue(100)

            ant_form.addRow('Range Size', self.range_size)
            ant_form.addRow('Max Pheromone', self.max_pheromone)
            ant_form.addRow('Ants per Sample', self.ants_per_sample)
            ant_form.addRow('Iterations (silent)', self.iterations)
            ant_form.addRow('Background Value', self.background_value)
            ant_form.addRow('Background Distance', self.background_distance)
            ant_form.addRow('Value Filter', self.value_filter)
            ant_form.addRow('Avoid Heavily-Visited', self.avoid_visited_enabled)
            ant_form.addRow('Visited Threshold', self.avoid_visited_threshold)

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

            # === ADVANCED TAB ===
            self.average_with_blocks = QtWidgets.QCheckBox(); self.average_with_blocks.setChecked(True)
            self.verbose = QtWidgets.QCheckBox(); self.verbose.setChecked(False)
            self.fill_unvisited_domainwise = QtWidgets.QCheckBox(); self.fill_unvisited_domainwise.setChecked(False)
            self.process_domains_sequentially = QtWidgets.QCheckBox(); self.process_domains_sequentially.setChecked(True)
            advanced_form.addRow('Average With Blocks', self.average_with_blocks)
            advanced_form.addRow('Fill Unvisited (Domain-wise)', self.fill_unvisited_domainwise)
            advanced_form.addRow('Process Domains Sequentially', self.process_domains_sequentially)
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
            self.start_btn = QtWidgets.QPushButton('Start with Visualization')
            self.cancel_btn = QtWidgets.QPushButton('Cancel')
            btn_box.addWidget(self.load_btn); btn_box.addWidget(self.save_btn)
            btn_box.addStretch(1)
            btn_box.addWidget(self.run_only_btn); btn_box.addWidget(self.start_btn); btn_box.addWidget(self.cancel_btn)
            main_layout.addLayout(btn_box)

            self.load_btn.clicked.connect(self.load_config)
            self.save_btn.clicked.connect(self.save_config)
            self.run_only_btn.clicked.connect(self.run_interpolation_only)
            self.start_btn.clicked.connect(self.accept)
            self.cancel_btn.clicked.connect(self.reject)

        def open_domain_mapping(self):
            """Open dialog to configure domain-specific algorithms"""
            blocks_file = self.blocks_edit.text().strip()
            blocks_delimiter = self.blocks_delim.currentText()
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
            
            dialog = DomainAlgorithmDialog(
                blocks_file, blocks_delimiter, blocks_header_line, 
                block_domain_col, self
            )
            
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
        
        def run_interpolation_only(self):
            """Run interpolation without visualization and close dialog"""
            try:
                cfg = self.to_dict()
                interpolation_file = cfg['interpolation_file']
                
                # Run interpolation directly without visualization
                print("=" * 60)
                print("Running interpolation without visualization...")
                print("=" * 60)
                
                # Load samples first
                samples_file = cfg['samples_file']
                samples_delimiter = cfg.get('samples_delimiter')
                samples_header_line = cfg.get('samples_header_line', 1)
                
                print(f"Loading sample file from {samples_file}...")
                if samples_header_line and samples_header_line != 1 and samples_delimiter:
                    df, _ = read_csv_with_selected_header(samples_file, samples_delimiter, samples_header_line, expected_min_cols=4)
                else:
                    df = read_autodetect_csv(samples_file, forced_delimiter=samples_delimiter)
                
                # Apply column mapping
                sample_x_col = cfg.get('sample_x_col')
                sample_y_col = cfg.get('sample_y_col')
                sample_z_col = cfg.get('sample_z_col')
                sample_value_col = cfg.get('sample_value_col')
                
                if sample_x_col and sample_y_col and sample_z_col and sample_value_col:
                    rename_map = {sample_x_col: 'x', sample_y_col: 'y', sample_z_col: 'z', sample_value_col: 'Value'}
                    df = df.rename(columns=rename_map)
                
                df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
                df = df.dropna(subset=['Value'])
                
                points = df[['x','y','z']].values
                values = df['Value'].values
                print(f"Loaded {len(values)} samples from {samples_file}.")
                
                # Create blocks with samples
                blocks = create_blocks(
                    points, values,
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
                    config=cfg
                )
                
                # Run interpolation (handle both single and sequential domain processing)
                dims = tuple(blocks._block_info['dims'])
                iterations = cfg['iterations']
                
                # Check if we have multiple interpolators (sequential domain processing)
                if hasattr(blocks, '_interpolators'):
                    print(f"Running sequential domain interpolation for {len(blocks._interpolators)} domains...")
                    
                    for domain_idx, (domain, interpolator) in enumerate(blocks._interpolators.items(), 1):
                        algo_name = interpolator.get_algorithm_name()
                        print(f"\n=== Domain {domain_idx}/{len(blocks._interpolators)}: {domain} ({algo_name}) ===")
                        
                        print(f"Running {algo_name} for {iterations} iterations...")
                        pbar = tqdm(range(iterations), desc=f"Domain {domain}")
                        for i in pbar:
                            should_continue = interpolator.run_iteration(dims)
                            if not should_continue or interpolator.is_converged():
                                print(f"Converged at iteration {i+1}")
                                break
                        
                        # Print domain summary
                        metadata = interpolator.get_metadata()
                        print(f"\n=== Domain {domain} Summary ===")
                        for key, value in metadata.items():
                            if key in ['intrusion_trees', 'ancestor_nodes']:
                                print(f"{key}: <{len(value)} items>")
                            else:
                                print(f"{key}: {value}")
                else:
                    # Single interpolator
                    interpolator = blocks._ant_colony
                    algo_name = interpolator.get_algorithm_name()
                    print(f"Running {algo_name} for {iterations} iterations...")
                    pbar = tqdm(range(iterations), desc="Interpolating")
                    for i in pbar:
                        should_continue = interpolator.run_iteration(dims)
                        if not should_continue or interpolator.is_converged():
                            print(f"Converged at iteration {i+1}")
                            break
                    
                    # Print summary
                    metadata = interpolator.get_metadata()
                    print(f"\n=== Interpolation Summary ===")
                    for key, value in metadata.items():
                        if key in ['intrusion_trees', 'ancestor_nodes']:
                            print(f"{key}: <{len(value)} items>")
                        else:
                            print(f"{key}: {value}")
                
                # Export results (handles both single and multiple interpolators)
                export_blocks_to_csv(blocks, interpolation_file)
                
                print("=" * 60)
                print(f"Interpolation complete! Results saved to:")
                print(f"  {interpolation_file}")
                print("=" * 60)
                
                QtWidgets.QMessageBox.information(self, 'Complete', 
                    f'Interpolation complete!\n\nResults saved to:\n{interpolation_file}')
                
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Interpolation failed: {e}')
                import traceback
                traceback.print_exc()

        def to_dict(self):
            samples_file = self.samples_edit.text().strip()
            if not os.path.isfile(samples_file):
                raise ValueError('Samples file does not exist.')
            blocks_file = self.blocks_edit.text().strip() or None
            if blocks_file and not os.path.isfile(blocks_file):
                raise ValueError('Blocks file does not exist.')
            color_file = self.color_edit.text().strip() or None
            if color_file and not os.path.isfile(color_file):
                color_file = None
            interp = self.interp_edit.text().strip()
            if not interp:
                base = os.path.splitext(os.path.basename(samples_file))[0]
                interp = os.path.join(os.path.dirname(samples_file), f"{base}_anterpolation.csv")
            os.makedirs(os.path.dirname(interp), exist_ok=True)
            return {
                'samples_file': samples_file,
                'blocks_file': blocks_file,
                'color_file': color_file,
                'interpolation_file': interp,
                'block_size': (self.block_x.value(), self.block_y.value(), self.block_z.value()),
                'range_size': self.range_size.value(),
                'max_pheromone': self.max_pheromone.value(),
                'ants_per_sample': self.ants_per_sample.value(),
                'iterations': self.iterations.value(),
                'background_value': self.background_value.value(),
                'background_distance': self.background_distance.value(),
                'average_with_blocks': self.average_with_blocks.isChecked(),
                'value_filter': self.value_filter.value(),
                'verbose': self.verbose.isChecked(),
                'fill_unvisited_domainwise': self.fill_unvisited_domainwise.isChecked(),
                'avoid_visited_threshold_enabled': self.avoid_visited_enabled.isChecked(),
                'avoid_visited_threshold': self.avoid_visited_threshold.value()
                ,'samples_delimiter': self.samples_delim.currentText()
                ,'blocks_delimiter': self.blocks_delim.currentText()
                ,'samples_header_line': self.samples_header_line.value()
                ,'blocks_header_line': self.blocks_header_line.value()
                ,'sample_x_col': self.sample_x_col.currentText() if self.sample_x_col.count() else None
                ,'sample_y_col': self.sample_y_col.currentText() if self.sample_y_col.count() else None
                ,'sample_z_col': self.sample_z_col.currentText() if self.sample_z_col.count() else None
                ,'sample_value_col': self.sample_value_col.currentText() if self.sample_value_col.count() else None
                ,'block_x_col': self.block_x_col.currentText() if self.block_x_col.count() else None
                ,'block_y_col': self.block_y_col.currentText() if self.block_y_col.count() else None
                ,'block_z_col': self.block_z_col.currentText() if self.block_z_col.count() else None
                ,'block_domain_col': (self.block_domain_col.currentText() if (self.block_domain_col.count() and self.block_domain_col.currentText() != '(None)') else None)
                ,'process_domains_sequentially': self.process_domains_sequentially.isChecked()
                ,'algorithm': self.algorithm_combo.currentText()
                ,'molecular_clock_params': {
                    'spatial_weight': self.mc_spatial_weight.value(),
                    'attr_weight': self.mc_attr_weight.value(),
                    'ancestor_depth_offset': self.mc_ancestor_depth_offset.value(),
                    'branch_threshold': self.mc_branch_threshold.value(),
                    'min_samples_per_event': self.mc_min_samples.value(),
                    'max_samples_per_event': self.mc_max_samples.value(),
                    'detect_multiple_events': self.mc_detect_multiple.isChecked(),
                    'interpolation_method': self.mc_interp_method.currentText(),
                    'fill_background': False,
                    'background_value': self.background_value.value()
                }
                ,'domain_algorithm_overrides': self.domain_overrides
            }

        def load_config(self):
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load Config', '.', 'JSON Files (*.json)')
            if not path:
                return
            try:
                with open(path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self.samples_edit.setText(data.get('samples_file', self.samples_edit.text()))
                self.blocks_edit.setText(data.get('blocks_file', self.blocks_edit.text() or ''))
                self.color_edit.setText(data.get('color_file', self.color_edit.text() or ''))
                self.interp_edit.setText(data.get('interpolation_file', self.interp_edit.text() or ''))
                # Delimiters & header lines
                if 'samples_delimiter' in data: self.samples_delim.setCurrentText(data['samples_delimiter'])
                if 'blocks_delimiter' in data: self.blocks_delim.setCurrentText(data['blocks_delimiter'])
                self.samples_header_line.setValue(int(data.get('samples_header_line', 1)))
                self.blocks_header_line.setValue(int(data.get('blocks_header_line', 1)))
                # Refresh column choices after header changes
                # (Signals will auto refresh; we proceed to set mapping selections)
                def set_if_exists(cb, value):
                    if value is None: return
                    idx = cb.findText(value)
                    if idx >= 0: cb.setCurrentIndex(idx)
                set_if_exists(self.sample_x_col, data.get('sample_x_col'))
                set_if_exists(self.sample_y_col, data.get('sample_y_col'))
                set_if_exists(self.sample_z_col, data.get('sample_z_col'))
                set_if_exists(self.sample_value_col, data.get('sample_value_col'))
                set_if_exists(self.block_x_col, data.get('block_x_col'))
                set_if_exists(self.block_y_col, data.get('block_y_col'))
                set_if_exists(self.block_z_col, data.get('block_z_col'))
                if data.get('block_domain_col'):
                    set_if_exists(self.block_domain_col, data.get('block_domain_col'))
                self.block_domain_col.setCurrentText(data.get('block_domain_col') or '(None)')
                # Numeric sizes
                bs = data.get('block_size', (10,10,10))
                if isinstance(bs, (list, tuple)) and len(bs) == 3:
                    self.block_x.setValue(int(bs[0])); self.block_y.setValue(int(bs[1])); self.block_z.setValue(int(bs[2]))
                self.range_size.setValue(float(data.get('range_size', self.range_size.value())))
                self.max_pheromone.setValue(int(data.get('max_pheromone', self.max_pheromone.value())))
                self.ants_per_sample.setValue(int(data.get('ants_per_sample', self.ants_per_sample.value())))
                self.iterations.setValue(int(data.get('iterations', self.iterations.value())))
                self.background_value.setValue(float(data.get('background_value', self.background_value.value())))
                self.background_distance.setValue(float(data.get('background_distance', self.background_distance.value())))
                self.average_with_blocks.setChecked(bool(data.get('average_with_blocks', self.average_with_blocks.isChecked())))
                self.value_filter.setValue(float(data.get('value_filter', self.value_filter.value())))
                self.verbose.setChecked(bool(data.get('verbose', self.verbose.isChecked())))
                self.fill_unvisited_domainwise.setChecked(bool(data.get('fill_unvisited_domainwise', self.fill_unvisited_domainwise.isChecked())))
                self.process_domains_sequentially.setChecked(bool(data.get('process_domains_sequentially', True)))
                self.avoid_visited_enabled.setChecked(bool(data.get('avoid_visited_threshold_enabled', False)))
                try:
                    self.avoid_visited_threshold.setValue(int(data.get('avoid_visited_threshold', 100)))
                except Exception:
                    self.avoid_visited_threshold.setValue(100)
                
                # Algorithm selection
                if 'algorithm' in data:
                    algo = data['algorithm']
                    if algo == 'biochemical_clock':
                        algo = 'molecular_clock'
                    idx = self.algorithm_combo.findText(algo)
                    if idx >= 0:
                        self.algorithm_combo.setCurrentIndex(idx)
                
                # Molecular clock parameters
                mc_params = data.get('molecular_clock_params', data.get('biochemical_clock_params', {}))
                self.mc_spatial_weight.setValue(float(mc_params.get('spatial_weight', 1.0)))
                self.mc_attr_weight.setValue(float(mc_params.get('attr_weight', 1.0)))
                self.mc_ancestor_depth_offset.setValue(float(mc_params.get('ancestor_depth_offset', 1.0)))
                self.mc_branch_threshold.setValue(float(mc_params.get('branch_threshold', 2.0)))
                self.mc_min_samples.setValue(int(mc_params.get('min_samples_per_event', 3)))
                self.mc_max_samples.setValue(int(mc_params.get('max_samples_per_event', 1000)))
                self.mc_detect_multiple.setChecked(bool(mc_params.get('detect_multiple_events', True)))
                interp_method = mc_params.get('interpolation_method', 'linear')
                idx = self.mc_interp_method.findText(interp_method)
                if idx >= 0:
                    self.mc_interp_method.setCurrentIndex(idx)
                
                # Domain algorithm overrides
                if 'domain_algorithm_overrides' in data:
                    self.domain_overrides = data['domain_algorithm_overrides']
                    count = len(self.domain_overrides)
                    if count > 0:
                        self.domain_mapping_btn.setText(f'Configure Domain Algorithms... ({count} configured)')
                
                QtWidgets.QMessageBox.information(self, 'Config', 'Configuration loaded.')
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Failed to load config: {e}')

        def save_config(self):
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Config', 'config.json', 'JSON Files (*.json)')
            if not path:
                return
            try:
                with open(path, 'w', encoding='utf-8') as f:
                    json.dump(self.to_dict(), f, indent=2)
                QtWidgets.QMessageBox.information(self, 'Config', 'Configuration saved.')
            except Exception as e:
                QtWidgets.QMessageBox.critical(self, 'Error', f'Failed to save config: {e}')

    app = QtWidgets.QApplication(sys.argv)
    dialog = ConfigDialog()
    if dialog.exec_() == QtWidgets.QDialog.Accepted:
        config = dialog.to_dict()
        load_and_visualize_samples(
            samples_file=config['samples_file'],
            block_size=config['block_size'],
            value_filter=config['value_filter'],
            verbose=config['verbose'],
            iterations=config['iterations'],
            range_size=config['range_size'],
            max_pheromone=config['max_pheromone'],
            ants_per_sample=config['ants_per_sample'],
            blocks_file=config['blocks_file'],
            color_file=config['color_file'],
            background_value=config['background_value'],
            background_distance=config['background_distance'],
            average_with_blocks=config['average_with_blocks'],
            samples_delimiter=config['samples_delimiter'],
            blocks_delimiter=config['blocks_delimiter'],
            fill_unvisited_domainwise=config['fill_unvisited_domainwise'],
            avoid_visited_threshold_enabled=config['avoid_visited_threshold_enabled'],
            avoid_visited_threshold=config['avoid_visited_threshold'],
            samples_header_line=config['samples_header_line'],
            sample_x_col=config['sample_x_col'],
            sample_y_col=config['sample_y_col'],
            sample_z_col=config['sample_z_col'],
            sample_value_col=config['sample_value_col'],
            blocks_header_line=config['blocks_header_line'],
            block_x_col=config['block_x_col'],
            block_y_col=config['block_y_col'],
            block_z_col=config['block_z_col'],
            block_domain_col=config['block_domain_col'],
            config=config
        )
    sys.exit()