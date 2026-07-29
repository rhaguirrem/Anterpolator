import numpy as np
from scipy.spatial import cKDTree
from typing import Dict, Tuple, List, Set, Any
from interpolator_base import InterpolatorBase

class StringTheoryInterpolator(InterpolatorBase):
    def __init__(
        self,
        distance_threshold: float = 10.0,
        grade_difference: float = 1.0,
        connect_to_all: bool = True,
        max_connections: int = 1,
        min_connections: int = 1,
        collision_policy: str = 'average',
        processing_order: str = 'ascending',
        filter_by_frequency: bool = False,
        min_azimuth_freq: float = 10.0,
        min_dip_freq: float = 10.0,
        interpolation_target: str = 'value',
        verbose: bool = False,
    ):
        super().__init__(verbose)
        self.distance_threshold = distance_threshold
        self.grade_difference = grade_difference
        self.connect_to_all = connect_to_all
        self.max_connections = max_connections
        self.min_connections = min_connections
        self.collision_policy = collision_policy
        self.processing_order = processing_order
        self.filter_by_frequency = filter_by_frequency
        self.min_azimuth_freq = min_azimuth_freq
        self.min_dip_freq = min_dip_freq

        target = (interpolation_target or 'value').strip().lower()
        if target in ('numeric', 'number', 'grade', 'value'):
            target = 'value'
        elif target in ('domain', 'domains', 'category', 'categorical'):
            target = 'domain'
        else:
            target = 'value'
        self.interpolation_target = target
        
        # State
        self.sorted_samples: List[Tuple[Tuple[int, int, int], float]] = []
        self.sample_coords: np.ndarray = None
        self.sample_values: np.ndarray = None
        self.kdtree: cKDTree = None
        self.current_index = 0
        self.sample_locations: Set[Tuple[int, int, int]] = set()
        self.progress_callback = None
        self._interpolated_block_count = 0
        
        # Statistics
        self.stats = {
            'neighbors_found': 0,
            'grade_rejected': 0,
            'adjacent_rejected': 0,
            'domain_rejected': 0,
            'obstacle_rejected': 0,
            'connections_made': 0,
            'no_candidates': 0
        }
        
        self.paths = [] # List of (p1, p2, path_points)
        self.kept_paths_indices = None # Set of indices of paths that passed filtering
        
        # Result storage: (x, y, z) -> dict
        # We initialize with sample values.
        self.blocks: Dict[Tuple[int, int, int], dict] = {}

        # Domain interpolation state
        self.sample_domain_mapping: Dict[Tuple[int, int, int], str] = {}
        # In domain mode, we also use self.domain_mapping as the evolving (pos -> domain) output mapping.

    def get_algorithm_name(self) -> str:
        return "String Theory"

    def initialize_blocks(self, sample_blocks: Dict[Tuple[int, int, int], float], 
                         dims: Tuple[int, int, int], min_bounds, block_size, **kwargs):
        self.dims = dims
        self.min_bounds = min_bounds
        self.block_size = block_size

        # Numeric-mode domain mapping is a constraint (pos -> domain), typically coming from blocks file.
        # Domain-mode domain mapping is the *target* being interpolated.
        self.use_domain_mapping = kwargs.get('use_domain_mapping', False)
        self.sample_domain_mapping = kwargs.get('sample_domain_mapping', None) or {}

        self.sample_locations = set(sample_blocks.keys())
        self.blocks = {}

        if self.interpolation_target == 'domain':
            # For domain interpolation, every sample block must have an assigned domain.
            # Store domains both on blocks and on a mapping for quick access.
            self.domain_mapping = {}
            for pos, val in sample_blocks.items():
                domain = self.sample_domain_mapping.get(pos)
                if domain is None or str(domain).strip() == '':
                    # Leave it unmapped; it will be skipped later.
                    continue
                domain = str(domain).strip()
                self.domain_mapping[pos] = domain
                self.blocks[pos] = {'value': float(val) if val is not None else 0.0, 'is_sample': True, 'count': 1, 'domain': domain}
            # Domain interpolation should not use domain-mapping-as-constraint.
            self.use_domain_mapping = False
        else:
            for pos, val in sample_blocks.items():
                self.blocks[pos] = {'value': val, 'is_sample': True, 'count': 1}
        
        # Prepare samples for processing
        reverse_order = (self.processing_order == 'descending')
        if self.processing_order == 'random':
            import random
            items = list(sample_blocks.items())
            random.shuffle(items)
            self.sorted_samples = items
            if self.verbose:
                print("String Theory: Processing samples in RANDOM order.")
        else:
            if self.interpolation_target == 'domain':
                # No numeric value ordering in domain mode; keep deterministic coordinate ordering.
                self.sorted_samples = sorted(sample_blocks.items(), key=lambda x: x[0], reverse=reverse_order)
                if self.verbose:
                    direction = 'DESCENDING' if reverse_order else 'ASCENDING'
                    print(f"String Theory: Processing samples in {direction} coordinate order (domain mode).")
            else:
                # Sort by value ascending/descending.
                self.sorted_samples = sorted(sample_blocks.items(), key=lambda x: x[1], reverse=reverse_order)
                if self.verbose:
                    direction = 'DESCENDING' if reverse_order else 'ASCENDING'
                    print(f"String Theory: Processing samples in {direction} value order.")
        
        # Prepare KDTree for fast spatial queries
        coords = np.array([list(k) for k in sample_blocks.keys()])
        self.sample_coords = coords
        # We need a mapping from coord index back to value/pos if we use the tree directly,
        # or we can just use the tree to find indices in the coords array.
        self.sample_values = np.array([sample_blocks[tuple(c)] for c in coords])
        
        if len(coords) > 0:
            self.kdtree = cKDTree(coords)
        
        self.current_index = 0
        self._interpolated_block_count = sum(1 for block in self.blocks.values() if not block.get('is_sample'))
        if self.verbose:
            print(f"Initialized StringTheory with {len(self.sorted_samples)} samples.")
            if self.interpolation_target == 'domain':
                print(f"Distance Threshold: {self.distance_threshold} (domain mode)")
            else:
                print(f"Distance Threshold: {self.distance_threshold}, Grade Diff: {self.grade_difference}")

    def _emit_progress(self, processed_samples: int, total_samples: int):
        callback = getattr(self, 'progress_callback', None)
        if callback is None:
            return
        callback(processed_samples, total_samples, self._interpolated_block_count)

    def _bresenham_3d(self, p1: Tuple[int, int, int], p2: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        """
        Generate points on a line between p1 and p2 (inclusive) using Bresenham's algorithm.
        """
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        points = []
        
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        dz = abs(z2 - z1)
        
        xs = 1 if x2 > x1 else -1
        ys = 1 if y2 > y1 else -1
        zs = 1 if z2 > z1 else -1
        
        # Driving axis
        if dx >= dy and dx >= dz:
            p1_err = 2 * dy - dx
            p2_err = 2 * dz - dx
            while x1 != x2:
                points.append((x1, y1, z1))
                x1 += xs
                if p1_err >= 0:
                    y1 += ys
                    p1_err -= 2 * dx
                if p2_err >= 0:
                    z1 += zs
                    p2_err -= 2 * dx
                p1_err += 2 * dy
                p2_err += 2 * dz
            points.append((x1, y1, z1))
            
        elif dy >= dx and dy >= dz:
            p1_err = 2 * dx - dy
            p2_err = 2 * dz - dy
            while y1 != y2:
                points.append((x1, y1, z1))
                y1 += ys
                if p1_err >= 0:
                    x1 += xs
                    p1_err -= 2 * dy
                if p2_err >= 0:
                    z1 += zs
                    p2_err -= 2 * dy
                p1_err += 2 * dx
                p2_err += 2 * dz
            points.append((x1, y1, z1))
            
        else:
            p1_err = 2 * dy - dz
            p2_err = 2 * dx - dz
            while z1 != z2:
                points.append((x1, y1, z1))
                z1 += zs
                if p1_err >= 0:
                    y1 += ys
                    p1_err -= 2 * dz
                if p2_err >= 0:
                    x1 += xs
                    p2_err -= 2 * dz
                p1_err += 2 * dy
                p2_err += 2 * dx
            points.append((x1, y1, z1))
            
        return points

    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        # Process ALL remaining samples in one go
        # String Theory is not iterative in the sense of convergence; it's a single pass algorithm.

        if self.interpolation_target == 'domain':
            return self._run_iteration_domain(dims)
        
        if self.kdtree is None:
            return False

        total_samples = len(self.sorted_samples)
        processed_count = 0
        
        # Debug counters
        # self.stats is initialized in __init__
        
        print(f"String Theory Debug: use_domain_mapping={self.use_domain_mapping}")
        if self.use_domain_mapping:
             print(f"String Theory Debug: domain_mapping size={len(self.domain_mapping) if getattr(self, 'domain_mapping', None) else 'None'}")

        while self.current_index < total_samples:
            # Process one sample
            s1_pos, s1_val = self.sorted_samples[self.current_index]
            self.current_index += 1
            processed_count += 1
            
            if self.verbose and processed_count % 1000 == 0:
                print(f"String Theory: Processed {self.current_index}/{total_samples} samples...", end='\r')

            # Find neighbors within distance threshold
            # query_ball_point returns indices of points in self.sample_coords
            neighbor_indices = self.kdtree.query_ball_point(s1_pos, self.distance_threshold)
            
            candidates = []
            for idx in neighbor_indices:
                s2_pos = tuple(self.sample_coords[idx])
                # Skip self
                if s2_pos == s1_pos:
                    continue
                    
                s2_val = self.sample_values[idx]
                
                # Check grade difference and direction (equal or greater)
                # We only connect to samples with equal or higher grade to avoid
                # redundant paths (A->B and B->A) and follow the gradient.
                # NOTE: We use <= grade_difference to include the boundary case
                if s2_val >= s1_val and abs(s2_val - s1_val) <= self.grade_difference:
                    dist = np.linalg.norm(np.array(s1_pos) - np.array(s2_pos))
                    candidates.append((dist, s2_pos, s2_val))
                else:
                    self.stats['grade_rejected'] += 1
            
            self.stats['neighbors_found'] += len(candidates)
            if not candidates:
                self.stats['no_candidates'] += 1
            
            # Sort candidates by grade difference (primary) and distance (secondary)
            # "search the equal or greater sample that makes the module of the difference smaller"
            candidates.sort(key=lambda x: (abs(x[2] - s1_val), x[0]))
            
            valid_connections_buffer = []
            
            for dist, s2_pos, s2_val in candidates:
                # Get path
                path = self._bresenham_3d(s1_pos, s2_pos)
                
                # Skip if samples are adjacent (path length <= 2)
                if len(path) <= 2:
                    self.stats['adjacent_rejected'] += 1
                    continue
                
                # Check domain constraint
                if self.use_domain_mapping and getattr(self, 'domain_mapping', None) is not None:
                    d1 = self.domain_mapping.get(s1_pos)
                    d2 = self.domain_mapping.get(s2_pos)
                    
                    if d1 is not None and d2 is not None and d1 == d2:
                        # Check if path leaves domain
                        path_leaves_domain = False
                        for pt in path:
                            # Note: path includes s1_pos and s2_pos, which we know are in d1/d2.
                            # We check intermediate points too.
                            if self.domain_mapping.get(pt) != d1:
                                path_leaves_domain = True
                                break
                        if path_leaves_domain:
                            self.stats['domain_rejected'] += 1
                            continue # Skip this candidate
                    elif d1 is None or d2 is None:
                         # If one of the points is not in the domain map (and we are using mapping), reject?
                         # Or if they are in different domains?
                         # Current logic only enforces if d1 == d2. 
                         # If d1 != d2, we might want to reject too? 
                         # Assuming we only connect within same domain.
                         if d1 != d2:
                             self.stats['domain_rejected'] += 1
                             continue

                # Check for obstacles (samples in the middle)
                # Exclude start and end points from check
                obstacle_found = False
                # path includes s1_pos and s2_pos. 
                # We check points strictly between them.
                # If path length is 2, there are no points in between.
                
                points_to_fill = []
                if len(path) > 2:
                    for pt in path[1:-1]:
                        if pt in self.sample_locations:
                            obstacle_found = True
                            break
                        points_to_fill.append(pt)
                
                if obstacle_found:
                    self.stats['obstacle_rejected'] += 1
                else:
                    # Valid path found!
                    # Buffer it
                    valid_connections_buffer.append((s1_pos, s1_val, s2_pos, s2_val, points_to_fill, path))
                    
                    # If not connecting to all, stop after max_connections
                    if not self.connect_to_all and len(valid_connections_buffer) >= self.max_connections:
                        break
            
            # Apply connections if minimum requirement is met
            if len(valid_connections_buffer) >= self.min_connections:
                for s1_p, s1_v, s2_p, s2_v, pts, pth in valid_connections_buffer:
                    self.stats['connections_made'] += 1
                    
                    # Store path for statistics
                    self.paths.append((s1_p, s2_p, pth))
                    
                    # Fill blocks if NOT filtering by frequency (otherwise we do it later)
                    if not self.filter_by_frequency:
                        self._fill_path_blocks(s1_p, s1_v, s2_p, s2_v, pts)

            if processed_count == 1 or processed_count == total_samples or processed_count % 100 == 0:
                self._emit_progress(processed_count, total_samples)
        
        # Apply filtering if enabled
        if self.filter_by_frequency:
            self._filter_and_fill_paths()

        self._interpolated_block_count = sum(1 for block in self.blocks.values() if not block.get('is_sample'))
        self._emit_progress(total_samples, total_samples)
        
        print(f"String Theory: Finished processing {total_samples} samples.        ")
        print("String Theory Stats:")
        for k, v in self.stats.items():
            print(f"  {k}: {v}")
            
        return False

    def _run_iteration_domain(self, dims: Tuple[int, int, int]) -> bool:
        if self.kdtree is None:
            return False

        total_samples = len(self.sorted_samples)
        processed_count = 0

        if self.verbose:
            print("String Theory: Running DOMAIN interpolation.")

        while self.current_index < total_samples:
            s1_pos, s1_val = self.sorted_samples[self.current_index]
            self.current_index += 1
            processed_count += 1

            d1 = None
            if getattr(self, 'domain_mapping', None) is not None:
                d1 = self.domain_mapping.get(s1_pos)
            if d1 is None:
                # Sample block has no domain; skip.
                continue

            neighbor_indices = self.kdtree.query_ball_point(s1_pos, self.distance_threshold)

            candidates = []
            for idx in neighbor_indices:
                s2_pos = tuple(self.sample_coords[idx])
                if s2_pos == s1_pos:
                    continue

                d2 = self.domain_mapping.get(s2_pos) if getattr(self, 'domain_mapping', None) is not None else None
                if d2 != d1:
                    self.stats['domain_rejected'] += 1
                    continue

                dist = float(np.linalg.norm(np.array(s1_pos) - np.array(s2_pos)))
                candidates.append((dist, s2_pos))

            self.stats['neighbors_found'] += len(candidates)
            if not candidates:
                self.stats['no_candidates'] += 1

            candidates.sort(key=lambda x: x[0])

            valid_connections_buffer = []

            for dist, s2_pos in candidates:
                path = self._bresenham_3d(s1_pos, s2_pos)

                if len(path) <= 2:
                    self.stats['adjacent_rejected'] += 1
                    continue

                obstacle_found = False
                points_to_fill = []
                if len(path) > 2:
                    for pt in path[1:-1]:
                        if pt in self.sample_locations:
                            obstacle_found = True
                            break
                        points_to_fill.append(pt)

                if obstacle_found:
                    self.stats['obstacle_rejected'] += 1
                    continue

                valid_connections_buffer.append((s1_pos, s2_pos, d1, points_to_fill, path))
                if not self.connect_to_all and len(valid_connections_buffer) >= self.max_connections:
                    break

            if len(valid_connections_buffer) >= self.min_connections:
                for s1_p, s2_p, domain, pts, pth in valid_connections_buffer:
                    self.stats['connections_made'] += 1
                    self.paths.append((s1_p, s2_p, pth))
                    if not self.filter_by_frequency:
                        self._fill_path_blocks_domain(domain, pts)

            if processed_count == 1 or processed_count == total_samples or processed_count % 100 == 0:
                self._emit_progress(processed_count, total_samples)

        if self.filter_by_frequency:
            self._filter_and_fill_paths_domain()

        self._interpolated_block_count = sum(1 for block in self.blocks.values() if not block.get('is_sample'))
        self._emit_progress(total_samples, total_samples)

        if self.verbose:
            print(f"String Theory: Finished processing {total_samples} samples (domain mode).")
            print("String Theory Stats:")
            for k, v in self.stats.items():
                print(f"  {k}: {v}")

        return False

    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        return {pos: b['value'] for pos, b in self.blocks.items()}
    
    def get_metadata(self) -> Dict[str, Any]:
        metadata = dict(self.stats)
        metadata['interpolated_blocks'] = int(self._interpolated_block_count)
        return metadata

    def _collect_sample_connection_rows(self) -> List[Tuple[Tuple[int, int, int], float, int]]:
        connection_counts = {
            pos: 0
            for pos, block in self.blocks.items()
            if block.get('is_sample')
        }

        for source_pos, _, _ in self.paths:
            if source_pos in connection_counts:
                connection_counts[source_pos] += 1

        rows = []
        for pos, connection_count in connection_counts.items():
            value = self.blocks.get(pos, {}).get('value')
            if value is None:
                continue
            rows.append((pos, float(value), int(connection_count)))

        rows.sort(key=lambda row: (row[1], row[0]))
        return rows
    
    def generate_statistics(self, output_dir: str, domain_name: str = "Global"):
        import matplotlib.pyplot as plt
        import os
        import csv
        
        if not self.paths:
            print(f"No paths generated for statistics in domain {domain_name}.")
            return

        lengths = []
        azimuths = []
        dips = []
        
        for p1, p2, path in self.paths:
            # Calculate length in BLOCKS (grid units)
            length = np.linalg.norm(np.array(p2) - np.array(p1))
            lengths.append(length)

            # Calculate orientation using REAL coordinates if available
            if hasattr(self, 'block_size') and self.block_size is not None:
                bs = np.array(self.block_size)
                p1_real = np.array(p1) * bs
                p2_real = np.array(p2) * bs
                
                # Ensure p_start has higher Z (or equal)
                if p1_real[2] < p2_real[2]:
                    p_start = p2_real
                    p_end = p1_real
                else:
                    p_start = p1_real
                    p_end = p2_real
                
                dx = p_end[0] - p_start[0]
                dy = p_end[1] - p_start[1]
                dz = p_end[2] - p_start[2]
            else:
                # Ensure p_start has higher Z (or equal)
                if p1[2] < p2[2]:
                    p_start = p2
                    p_end = p1
                else:
                    p_start = p1
                    p_end = p2

                dx = p_end[0] - p_start[0]
                dy = p_end[1] - p_start[1]
                dz = p_end[2] - p_start[2]
            
            # Azimuth (0-360)
            # arctan2(dx, dy) gives angle from North (Y axis) clockwise if we consider:
            # Y is North, X is East.
            az = np.degrees(np.arctan2(dx, dy))
            if az < 0:
                az += 360
            azimuths.append(az)
            
            # Dip (0 to 90) - positive downwards
            h_dist = np.sqrt(dx**2 + dy**2)
            if h_dist == 0:
                dip = 90.0 # Vertical
            else:
                # dz is <= 0 because we swapped to make start higher than end.
                # We want positive dip [0, 90].
                dip = np.degrees(np.arctan(abs(dz) / h_dist))
            dips.append(dip)
            
        # Create output directory
        stats_dir = os.path.join(output_dir, "StringTheory_Stats")
        os.makedirs(stats_dir, exist_ok=True)
        
        # Sanitize domain name for filename
        safe_domain = "".join([c for c in domain_name if c.isalnum() or c in (' ', '_', '-')]).strip()
        sample_connection_rows = self._collect_sample_connection_rows()
        
        # 1. Path Length Histogram
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, color='skyblue', edgecolor='black')
        plt.title(f'Path Length Distribution - {domain_name}')
        plt.xlabel('Length (blocks)')
        plt.ylabel('Frequency')
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(stats_dir, f'path_length_histogram_{safe_domain}.png'))
        plt.close()
        
        # 2. Rose Diagram (Azimuth)
        plt.figure(figsize=(8, 8))
        ax = plt.subplot(111, polar=True)
        
        # Convert Azimuth (deg, CW from North) to Radians for plotting
        # Matplotlib polar plot with set_theta_zero_location('N') and set_theta_direction(-1)
        # expects angles in radians, 0 at North, increasing Clockwise.
        
        az_rad = np.radians(azimuths)
        bins = np.linspace(0.0, 2 * np.pi, 37) # 10 degree bins
        hist, _ = np.histogram(az_rad, bins=bins)
        width = 2 * np.pi / 36
        
        bars = ax.bar(bins[:-1], hist, width=width, bottom=0.0, color='coral', edgecolor='black')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1) # Clockwise
        plt.title(f'Path Azimuth Distribution - {domain_name}')
        plt.savefig(os.path.join(stats_dir, f'azimuth_rose_diagram_{safe_domain}.png'))
        plt.close()
        
        # 3. Dip Histogram
        plt.figure(figsize=(10, 6))
        plt.hist(dips, bins=36, range=(0, 90), color='lightgreen', edgecolor='black')
        plt.title(f'Path Dip Distribution - {domain_name}')
        plt.xlabel('Dip (degrees)')
        plt.ylabel('Frequency')
        plt.grid(True, alpha=0.3)
        plt.savefig(os.path.join(stats_dir, f'dip_histogram_{safe_domain}.png'))
        plt.close()

        # 4. Connection-count histogram per sample block
        if sample_connection_rows:
            connection_counts = [row[2] for row in sample_connection_rows]
            unique_connection_counts = np.unique(connection_counts)
            plt.figure(figsize=(10, 6))
            if len(unique_connection_counts) == 1:
                connection_count = int(unique_connection_counts[0])
                plt.bar([connection_count], [len(connection_counts)], width=0.8, color='mediumpurple', edgecolor='black')
            else:
                min_count = int(np.min(connection_counts))
                max_count = int(np.max(connection_counts))
                bins = np.arange(min_count - 0.5, max_count + 1.5, 1.0)
                plt.hist(connection_counts, bins=bins, color='mediumpurple', edgecolor='black')
                if max_count - min_count <= 25:
                    plt.xticks(range(min_count, max_count + 1))
            plt.title(f'Sample Block Connection Count Distribution - {domain_name}')
            plt.xlabel('Number of Connections')
            plt.ylabel('Sample Block Count')
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(stats_dir, f'connection_count_histogram_{safe_domain}.png'))
            plt.close()

            with open(os.path.join(stats_dir, f'sample_connection_stats_{safe_domain}.csv'), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sample_x', 'sample_y', 'sample_z', 'sample_value', 'connection_count'])
                for pos, sample_value, connection_count in sample_connection_rows:
                    writer.writerow([pos[0], pos[1], pos[2], sample_value, connection_count])
        
        # Save raw stats to CSV
        with open(os.path.join(stats_dir, f'path_stats_{safe_domain}.csv'), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Length', 'Azimuth', 'Dip', 'mid_x', 'mid_y', 'mid_z', 'Polarity', 'Category', 'Accepted'])
            
            for i, (l, a, d) in enumerate(zip(lengths, azimuths, dips)):
                p1, p2, _ = self.paths[i]
                
                # Determine acceptance status
                is_accepted = True
                if self.filter_by_frequency and self.kept_paths_indices is not None:
                    is_accepted = i in self.kept_paths_indices
                
                # Calculate midpoint
                p1_arr = np.array(p1)
                p2_arr = np.array(p2)
                mid_block = (p1_arr + p2_arr) / 2.0
                
                if hasattr(self, 'block_size') and self.block_size is not None and hasattr(self, 'min_bounds') and self.min_bounds is not None:
                    bs = np.array(self.block_size)
                    mb = np.array(self.min_bounds)
                    mid_real = mb + mid_block * bs
                    mid_x, mid_y, mid_z = mid_real
                else:
                    mid_x, mid_y, mid_z = mid_block
                
                writer.writerow([l, a, d, mid_x, mid_y, mid_z, 1, domain_name, is_accepted])
                
        print(f"String Theory statistics for {domain_name} saved to {stats_dir}")
    
    def _calculate_path_orientation(self, p1: Tuple[int, int, int], p2: Tuple[int, int, int]) -> Tuple[float, float]:
        """
        Calculate Azimuth and Dip for a path between p1 and p2.
        Returns (azimuth, dip) in degrees.
        Dip is always positive [0, 90] (downwards from horizontal).
        """
        # Calculate orientation using REAL coordinates if available
        if hasattr(self, 'block_size') and self.block_size is not None:
            bs = np.array(self.block_size)
            p1_real = np.array(p1) * bs
            p2_real = np.array(p2) * bs
            
            # Ensure p_start has higher Z (or equal)
            if p1_real[2] < p2_real[2]:
                p_start = p2_real
                p_end = p1_real
            else:
                p_start = p1_real
                p_end = p2_real
            
            dx = p_end[0] - p_start[0]
            dy = p_end[1] - p_start[1]
            dz = p_end[2] - p_start[2]
        else:
            # Ensure p_start has higher Z (or equal)
            if p1[2] < p2[2]:
                p_start = p2
                p_end = p1
            else:
                p_start = p1
                p_end = p2

            dx = p_end[0] - p_start[0]
            dy = p_end[1] - p_start[1]
            dz = p_end[2] - p_start[2]
        
        # Azimuth (0-360)
        az = np.degrees(np.arctan2(dx, dy))
        if az < 0:
            az += 360
        
        # Dip (0 to 90)
        h_dist = np.sqrt(dx**2 + dy**2)
        if h_dist == 0:
            dip = 90.0 # Vertical
        else:
            dip = np.degrees(np.arctan(abs(dz) / h_dist))
            
        return az, dip

    def _fill_path_blocks(self, s1_pos, s1_val, s2_pos, s2_val, points_to_fill):
        """Helper to fill blocks along a path."""
        for pt in points_to_fill:
            # Calculate weighted average
            d1 = np.linalg.norm(np.array(pt) - np.array(s1_pos))
            d2 = np.linalg.norm(np.array(pt) - np.array(s2_pos))
            total_d = d1 + d2
            if total_d == 0:
                val = (s1_val + s2_val) / 2
            else:
                val = (s1_val * d2 + s2_val * d1) / total_d
            
            if pt in self.blocks:
                # Handle collision
                existing = self.blocks[pt]
                if existing.get('is_sample', False):
                    continue # Don't overwrite samples
                
                old_val = existing['value']
                count = existing.get('count', 1)
                
                if self.collision_policy == 'average':
                    # Running average
                    new_val = (old_val * count + val) / (count + 1)
                    self.blocks[pt]['value'] = new_val
                    self.blocks[pt]['count'] = count + 1
                elif self.collision_policy == 'min':
                    self.blocks[pt]['value'] = min(old_val, val)
                    self.blocks[pt]['count'] = count + 1
                elif self.collision_policy == 'max':
                    self.blocks[pt]['value'] = max(old_val, val)
                    self.blocks[pt]['count'] = count + 1
                else: # 'overwrite' or unknown
                    self.blocks[pt]['value'] = val
                    self.blocks[pt]['count'] = count + 1
            else:
                self.blocks[pt] = {'value': val, 'is_sample': False, 'count': 1}
                self._interpolated_block_count += 1

    def _fill_path_blocks_domain(self, domain: str, points_to_fill: List[Tuple[int, int, int]]):
        domain = str(domain).strip()
        for pt in points_to_fill:
            if pt in self.blocks and self.blocks[pt].get('is_sample', False):
                continue

            existing = self.blocks.get(pt)
            if existing is None:
                self.blocks[pt] = {'value': 0.0, 'is_sample': False, 'count': 1, 'domain': domain}
                if getattr(self, 'domain_mapping', None) is None:
                    self.domain_mapping = {}
                self.domain_mapping[pt] = domain
                self._interpolated_block_count += 1
                continue

            existing_domain = existing.get('domain')
            if existing_domain == domain:
                continue

            # Collision between domains: only 'overwrite' is meaningful here.
            if self.collision_policy == 'overwrite':
                existing['domain'] = domain
                existing['value'] = 0.0
                if getattr(self, 'domain_mapping', None) is None:
                    self.domain_mapping = {}
                self.domain_mapping[pt] = domain

    def _filter_and_fill_paths_domain(self):
        if not self.paths:
            return

        azimuths = []
        dips = []
        for p1, p2, _ in self.paths:
            az, dip = self._calculate_path_orientation(p1, p2)
            azimuths.append(az)
            dips.append(dip)

        az_bins = np.linspace(0, 360, 37)
        az_hist, az_bin_edges = np.histogram(azimuths, bins=az_bins)
        max_az_freq = np.max(az_hist) if len(az_hist) > 0 else 0

        dip_bins = np.linspace(0, 90, 37)
        dip_hist, dip_bin_edges = np.histogram(dips, bins=dip_bins)
        max_dip_freq = np.max(dip_hist) if len(dip_hist) > 0 else 0

        az_threshold_count = max_az_freq * (self.min_azimuth_freq / 100.0)
        dip_threshold_count = max_dip_freq * (self.min_dip_freq / 100.0)

        kept_count = 0
        rejected_count = 0
        self.kept_paths_indices = set()

        for i, (p1, p2, path) in enumerate(self.paths):
            az, dip = azimuths[i], dips[i]

            az_idx = np.digitize([az], az_bin_edges)[0] - 1
            if az_idx >= len(az_hist):
                az_idx = len(az_hist) - 1
            if az_idx < 0:
                az_idx = 0

            dip_idx = np.digitize([dip], dip_bin_edges)[0] - 1
            if dip_idx >= len(dip_hist):
                dip_idx = len(dip_hist) - 1
            if dip_idx < 0:
                dip_idx = 0

            current_az_freq = az_hist[az_idx]
            current_dip_freq = dip_hist[dip_idx]

            if current_az_freq >= az_threshold_count and current_dip_freq >= dip_threshold_count:
                kept_count += 1
                self.kept_paths_indices.add(i)
                domain = None
                if getattr(self, 'domain_mapping', None) is not None:
                    domain = self.domain_mapping.get(p1)
                if domain is None:
                    continue
                if len(path) > 2:
                    points_to_fill = path[1:-1]
                    self._fill_path_blocks_domain(domain, points_to_fill)
            else:
                rejected_count += 1

        if self.verbose:
            print(f"String Theory Filtering (domain): Kept {kept_count} paths, Rejected {rejected_count} paths.")

    def _filter_and_fill_paths(self):
        """
        Analyze all generated paths, calculate Azimuth/Dip histograms,
        and filter out paths that fall into bins with frequency lower than
        the specified percentage of the maximum frequency.
        """
        if not self.paths:
            return

        azimuths = []
        dips = []
        
        # 1. Calculate orientations for all paths
        for p1, p2, _ in self.paths:
            az, dip = self._calculate_path_orientation(p1, p2)
            azimuths.append(az)
            dips.append(dip)
            
        # 2. Calculate Histograms
        
        # Azimuth bins (10 degrees)
        az_bins = np.linspace(0, 360, 37)
        az_hist, az_bin_edges = np.histogram(azimuths, bins=az_bins)
        max_az_freq = np.max(az_hist) if len(az_hist) > 0 else 0
        
        # Dip bins (5 degrees)
        dip_bins = np.linspace(0, 90, 37)
        dip_hist, dip_bin_edges = np.histogram(dips, bins=dip_bins)
        max_dip_freq = np.max(dip_hist) if len(dip_hist) > 0 else 0
        
        # Calculate thresholds (counts)
        az_threshold_count = max_az_freq * (self.min_azimuth_freq / 100.0)
        dip_threshold_count = max_dip_freq * (self.min_dip_freq / 100.0)
        
        if self.verbose:
            print(f"String Theory Filtering: Max Azimuth Freq={max_az_freq}, Threshold={az_threshold_count:.1f} ({self.min_azimuth_freq}%)")
            print(f"String Theory Filtering: Max Dip Freq={max_dip_freq}, Threshold={dip_threshold_count:.1f} ({self.min_dip_freq}%)")
            
        # 3. Filter and Fill
        kept_count = 0
        rejected_count = 0
        self.kept_paths_indices = set()
        
        for i, (p1, p2, path) in enumerate(self.paths):
            az, dip = azimuths[i], dips[i]
            
            # Find bin indices
            # np.digitize returns 1-based index. 
            az_idx = np.digitize([az], az_bin_edges)[0] - 1
            if az_idx >= len(az_hist): az_idx = len(az_hist) - 1
            if az_idx < 0: az_idx = 0 # Should not happen given range 0-360
            
            dip_idx = np.digitize([dip], dip_bin_edges)[0] - 1
            if dip_idx >= len(dip_hist): dip_idx = len(dip_hist) - 1
            if dip_idx < 0: dip_idx = 0
            
            current_az_freq = az_hist[az_idx]
            current_dip_freq = dip_hist[dip_idx]
            
            if current_az_freq >= az_threshold_count and current_dip_freq >= dip_threshold_count:
                # Keep path
                kept_count += 1
                self.kept_paths_indices.add(i)
                
                # Re-construct points_to_fill (excluding start/end which are samples)
                # path contains [p1, ..., p2]
                if len(path) > 2:
                    points_to_fill = path[1:-1]
                    # Get values
                    v1 = self.blocks[p1]['value']
                    v2 = self.blocks[p2]['value']
                    
                    self._fill_path_blocks(p1, v1, p2, v2, points_to_fill)
            else:
                rejected_count += 1
                
        if self.verbose:
            print(f"String Theory Filtering: Kept {kept_count} paths, Rejected {rejected_count} paths.")
