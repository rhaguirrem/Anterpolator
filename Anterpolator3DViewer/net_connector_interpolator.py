import numpy as np
from scipy.spatial import cKDTree
from typing import Dict, Tuple, List, Set
from interpolator_base import InterpolatorBase

class NetConnectorInterpolator(InterpolatorBase):
    def __init__(self, distance_threshold: float = 10.0, grade_difference: float = 1.0, verbose=False):
        super().__init__(verbose)
        self.distance_threshold = distance_threshold
        self.grade_difference = grade_difference
        
        # State
        self.sorted_samples: List[Tuple[Tuple[int, int, int], float]] = []
        self.sample_coords: np.ndarray = None
        self.sample_values: np.ndarray = None
        self.kdtree: cKDTree = None
        self.current_index = 0
        self.sample_locations: Set[Tuple[int, int, int]] = set()
        
        # Result storage: (x, y, z) -> dict
        # We initialize with sample values.
        self.blocks: Dict[Tuple[int, int, int], dict] = {}

    def get_algorithm_name(self) -> str:
        return "Net Connector"

    def initialize_blocks(self, sample_blocks: Dict[Tuple[int, int, int], float], 
                         dims: Tuple[int, int, int], min_bounds, block_size, **kwargs):
        self.dims = dims
        self.use_domain_mapping = kwargs.get('use_domain_mapping', False)
        self.sample_locations = set(sample_blocks.keys())
        self.blocks = {}
        for pos, val in sample_blocks.items():
            self.blocks[pos] = {'value': val, 'is_sample': True}
        
        # Prepare samples for processing
        # Sort by grade (value) ascending
        self.sorted_samples = sorted(sample_blocks.items(), key=lambda x: x[1])
        
        # Prepare KDTree for fast spatial queries
        coords = np.array([list(k) for k in sample_blocks.keys()])
        self.sample_coords = coords
        # We need a mapping from coord index back to value/pos if we use the tree directly,
        # or we can just use the tree to find indices in the coords array.
        self.sample_values = np.array([sample_blocks[tuple(c)] for c in coords])
        
        if len(coords) > 0:
            self.kdtree = cKDTree(coords)
        
        self.current_index = 0
        if self.verbose:
            print(f"Initialized NetConnector with {len(self.sorted_samples)} samples.")
            print(f"Distance Threshold: {self.distance_threshold}, Grade Diff: {self.grade_difference}")

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
        if self.current_index >= len(self.sorted_samples):
            return False
            
        # Process one sample
        s1_pos, s1_val = self.sorted_samples[self.current_index]
        self.current_index += 1
        
        if self.kdtree is None:
            return True

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
            
            # Check grade difference
            if abs(s2_val - s1_val) < self.grade_difference:
                dist = np.linalg.norm(np.array(s1_pos) - np.array(s2_pos))
                candidates.append((dist, s2_pos, s2_val))
        
        # Sort candidates by distance (nearest first)
        candidates.sort(key=lambda x: x[0])
        
        path_found = False
        for dist, s2_pos, s2_val in candidates:
            # Get path
            path = self._bresenham_3d(s1_pos, s2_pos)
            
            # Check domain constraint
            if self.use_domain_mapping and hasattr(self, 'domain_mapping'):
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
                        continue # Skip this candidate

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
            
            if not obstacle_found:
                # Valid path found!
                # Fill blocks
                for pt in points_to_fill:
                    # Calculate weighted average
                    d1 = np.linalg.norm(np.array(pt) - np.array(s1_pos))
                    d2 = np.linalg.norm(np.array(pt) - np.array(s2_pos))
                    total_d = d1 + d2
                    if total_d == 0:
                        val = (s1_val + s2_val) / 2
                    else:
                        val = (s1_val * d2 + s2_val * d1) / total_d
                    
                    self.blocks[pt] = {'value': val, 'is_sample': False}
                
                path_found = True
                break # Stop after finding the first valid connection for this sample
        
        return True

    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        return {pos: b['value'] for pos, b in self.blocks.items()}
