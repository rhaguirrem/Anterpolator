"""
Molecular Clock Interpolator - Phylogeographic approach for intrusion modeling.

This interpolator uses spatial parsimony principles to find "common ancestor" points
(feeders) and builds minimum spanning trees connecting samples through branches.
Ideal for modeling magmatic intrusions, dikes, and vein systems.
"""
import numpy as np
from typing import Dict, Tuple, List, Any
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from scipy.ndimage import label
from interpolator_base import InterpolatorBase

try:
    from sklearn.cluster import DBSCAN
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    print("Warning: sklearn not available. Multi-event detection disabled.")


class MolecularClockInterpolator(InterpolatorBase):
    """
    Phylogeographic interpolator inspired by molecular clocks.
    
    Finds feeder sources (common ancestors) and builds branching tree structures
    connecting samples. Uses minimum spanning trees and spatial parsimony.
    
    Perfect for modeling:
    - Magmatic intrusions
    - Hydrothermal vein systems
    - Dike/sill networks
    - Any system with point-source dispersion
    """
    
    def __init__(self, 
                 spatial_weight: float = 1.0,
                 attr_weight: float = 1.0,
                 detect_multiple_events: bool = True,
                 branch_threshold: float = 2.0,
                 min_samples_per_event: int = 3,
                 max_samples_per_event: int = 1000,
                 interpolation_method: str = "linear",
                 fill_background: bool = True,
                 background_value: float = 0.0,
                 ancestor_depth_offset: float = 1.0,
                 verbose: bool = False):
        """
        Parameters:
        -----------
        spatial_weight : float
            Weight for spatial distance in composite distance calculation
        attr_weight : float
            Weight for attribute (value) distance
        detect_multiple_events : bool
            If True, attempts to identify multiple intrusion events
        branch_threshold : float
            Multiplier for clustering threshold (higher = fewer clusters)
        min_samples_per_event : int
            Minimum samples required to define an intrusion event
        max_samples_per_event : int
            Maximum samples allowed per event. Events exceeding this will be randomly
            subsampled to prevent performance issues. Set to 0 to disable limit.
        interpolation_method : str
            Method for interpolating along branches: "linear", "inverse_distance"
        fill_background : bool
            If True, fills blocks not on branches with background value
        background_value : float
            Value for background blocks
        ancestor_depth_offset : float
            Minimum additional depth (in grid units) used when placing ancestor nodes
        verbose : bool
            Print detailed progress information
        """
        super().__init__(verbose)
        self.spatial_weight = spatial_weight
        self.attr_weight = attr_weight
        self.detect_multiple_events = detect_multiple_events and HAS_SKLEARN
        self.branch_threshold = branch_threshold
        self.depth_weight = 0.0  # Deprecated
        self.min_samples_per_event = min_samples_per_event
        self.max_samples_per_event = max_samples_per_event
        self.interpolation_method = interpolation_method
        self.fill_background = fill_background
        self.background_value = background_value
        self.ancestor_depth_offset = ancestor_depth_offset
        self.min_branch_angle = 0.0  # Deprecated/Disabled by default
        
        # Data storage
        self.sample_coords = None
        self.sample_values = None
        self.sample_domains = None
        self.intrusion_trees: List[Dict[str, Any]] = []
        self.ancestor_nodes: List[Dict[str, Any]] = []
        self.min_bounds = None
        self.block_size = None
        self.converged = False
        self.mutation_rate = 1.0  # Will be inferred from data

    def _calculate_mutation_rate(self):
        """
        Infer mutation rate (grade change per unit distance) from data.
        Uses median gradient between nearest neighbors.
        """
        if self.sample_coords is None or len(self.sample_coords) < 2:
            return 1.0
            
        # Convert to world coordinates for accurate distance
        world_coords = self._grid_to_world(self.sample_coords)
        
        # Find nearest neighbor for each point
        from scipy.spatial import cKDTree
        tree = cKDTree(world_coords)
        dists, indices = tree.query(world_coords, k=2)  # k=2 because closest is itself
        
        # dists[:, 1] is distance to nearest neighbor
        # indices[:, 1] is index of nearest neighbor
        neighbor_dists = dists[:, 1]
        neighbor_indices = indices[:, 1]
        
        # Calculate grade differences
        grade_diffs = np.abs(self.sample_values - self.sample_values[neighbor_indices])
        
        # Calculate gradients (avoid division by zero)
        valid_mask = neighbor_dists > 1e-6
        if np.sum(valid_mask) == 0:
            return 1.0
            
        gradients = grade_diffs[valid_mask] / neighbor_dists[valid_mask]
        
        # Use 75th percentile gradient as the "natural" mutation rate
        # This allows for steeper gradients (higher variance) without penalizing connection too much
        rate = np.percentile(gradients, 75)
        
        # Fallback if rate is too small (e.g. constant grade)
        if rate < 1e-6:
            rate = 1.0
            
        if self.verbose:
            print(f"  Inferred mutation rate: {rate:.4f} grade_units/distance_unit (75th percentile)")
            
        return rate
        
    def initialize_blocks(self, sample_blocks: Dict[Tuple[int, int, int], float],
                         dims: Tuple[int, int, int], min_bounds, block_size, 
                         use_domain_mapping: bool = False, **kwargs):
        """Initialize with sample data"""
        self.dims = np.array(dims, dtype=int)
        self.min_bounds = np.array(min_bounds)
        self.block_size = np.array(block_size)
        
        # Extract sample coordinates and values
        coords_list = []
        values_list = []
        domains_list = []
        
        for pos, value in sample_blocks.items():
            coords_list.append(pos)
            values_list.append(value)
            # Extract domain if available
            domain = ''
            if hasattr(self, 'domain_mapping') and self.domain_mapping:
                domain = self.domain_mapping.get(pos, 'Undomained')
            else:
                domain = kwargs.get('block_domains', {}).get(pos, '')
            domains_list.append(domain)
        
        self.sample_coords = np.array(coords_list)
        self.sample_values = np.array(values_list)
        self.sample_domains = np.array(domains_list)
        
        # Initialize sample blocks
        for pos, value in sample_blocks.items():
            # Find domain for this position
            domain = 'Undomained'
            if hasattr(self, 'domain_mapping') and self.domain_mapping:
                domain = self.domain_mapping.get(pos, 'Undomained')
            else:
                domain = kwargs.get('block_domains', {}).get(pos, 'Undomained')
                
            self.blocks[pos] = {
                'value': value,
                'is_sample': True,
                'event_id': -1,
                'distance_to_feeder': 0,
                'branch_id': -1,
                'is_feeder': False,
                'visit_count': 1,
                'domain': domain
            }
        
        if self.verbose:
            print(f"Initialized {len(sample_blocks)} sample blocks")
            print(f"Grid dimensions: {dims}")
            print(f"Spatial weight: {self.spatial_weight}, Attribute weight: {self.attr_weight}")
            
        # Infer mutation rate
        self.mutation_rate = self._calculate_mutation_rate()
    
    def _grid_to_world(self, coords: np.ndarray) -> np.ndarray:
        if self.min_bounds is None or self.block_size is None:
            return np.array(coords, dtype=float)
        return self.min_bounds + np.array(coords, dtype=float) * self.block_size

    def _world_to_grid(self, coords: np.ndarray) -> Tuple[int, int, int]:
        if self.min_bounds is None or self.block_size is None:
            return tuple(np.round(coords).astype(int))
        relative = (np.array(coords, dtype=float) - self.min_bounds) / self.block_size
        relative = np.clip(relative, np.zeros(3), self.dims - 1)
        return tuple(np.round(relative).astype(int))

    def detect_intrusion_events(self, sample_indices: np.ndarray = None) -> List[np.ndarray]:
        """
        Identify multiple intrusion events using DBSCAN clustering.
        Returns list of sample indices for each event.
        """
        if sample_indices is None:
            sample_indices = np.arange(len(self.sample_coords))
            
        n_samples = len(sample_indices)
        coords = self.sample_coords[sample_indices]
        values = self.sample_values[sample_indices]
        
        if not self.detect_multiple_events or n_samples < self.min_samples_per_event:
            return [sample_indices]
        
        # Compute composite distance matrix using Evolutionary Distance (Divergence)
        
        # 1. Spatial Distance (World Coordinates)
        world_coords = self._grid_to_world(coords)
        spatial_dist = cdist(world_coords, world_coords, 'euclidean')
        
        # 2. Grade Distance (Normalized by Mutation Rate)
        value_diffs = cdist(values.reshape(-1, 1), 
                           values.reshape(-1, 1), 'cityblock') # Absolute difference
        grade_dist = value_diffs / self.mutation_rate
        
        # 3. Total Divergence
        # sqrt(Spatial^2 + (Weight * Grade_Spatial)^2)
        composite_dist = np.sqrt(
            self.spatial_weight * spatial_dist**2 + 
            (self.attr_weight * grade_dist)**2
        )
        
        # DBSCAN clustering with more conservative eps
        # Use minimum distance to k-th nearest neighbor as baseline
        # Sort distances for each point and take the k-th closest
        # Dynamic k based on dataset size to capture broader structure (inter-hole spacing)
        # sqrt(N) is a standard heuristic for density estimation neighborhood
        target_k = int(np.sqrt(n_samples))
        target_k = np.clip(target_k, 5, 100)  # Clamp between 5 and 100
        
        k = min(target_k, n_samples - 1)
        sorted_dists = np.sort(composite_dist, axis=1)
        kth_neighbor_dists = sorted_dists[:, k]  # k-th nearest neighbor distance for each point
        
        # Use a percentile of k-nearest distances, scaled by branch_threshold
        # Use a higher percentile (90th) to better capture inter-hole spacing in drillhole data
        # With larger k, this becomes a very robust estimate of "regional" connectivity
        eps = np.percentile(kth_neighbor_dists, 90) * self.branch_threshold
        
        # Enforce minimum epsilon based on project extent (to bridge large gaps in sparse data)
        if self.min_bounds is not None and self.block_size is not None and self.dims is not None:
            # Calculate project diagonal in world units
            extent = self.dims * self.block_size
            diagonal = np.linalg.norm(extent)
            
            # Minimum spatial reach: 5% of project diagonal
            min_spatial_reach = diagonal * 0.05
            
            # Convert to composite distance units
            # composite = sqrt(spatial_weight * spatial^2 + ...)
            # We assume grade contribution is zero for this baseline reach
            min_eps = np.sqrt(self.spatial_weight * min_spatial_reach**2)
            
            if eps < min_eps:
                if self.verbose:
                    print(f"  Boosting eps from {eps:.2f} to {min_eps:.2f} (5% of project extent: {min_spatial_reach:.1f}m)")
                eps = min_eps
        
        if self.verbose:
            print(f"  DBSCAN eps: {eps:.2f} (90th percentile of {k}-nearest neighbor distances × {self.branch_threshold})")
        
        clustering = DBSCAN(
            eps=eps,
            min_samples=self.min_samples_per_event,
            metric='precomputed'
        ).fit(composite_dist)
        
        # Group samples by cluster
        events = []
        unique_labels = set(clustering.labels_)
        
        # Temporary storage for clusters
        clusters = []
        
        for label in sorted(unique_labels):
            if label == -1:  # Noise points - treat as individual events if enough samples
                noise_indices = np.where(clustering.labels_ == -1)[0]
                # For noise, we might want to try merging them into nearest clusters later
                # But for now, treat as potential events
                if len(noise_indices) >= self.min_samples_per_event:
                    clusters.append(noise_indices)
                continue
            
            mask = clustering.labels_ == label
            indices = np.where(mask)[0]
            if len(indices) >= self.min_samples_per_event:
                clusters.append(indices)
        
        # Post-processing: Merge clusters that are "close enough"
        # DBSCAN is strict (density-connected). We want to allow "hops" between dense clusters.
        # If the minimum distance between two clusters is < eps * 2.0, merge them.
        
        if len(clusters) > 1:
            if self.verbose:
                print(f"  Found {len(clusters)} initial clusters. Checking for merges...")
            
            # Iteratively merge closest clusters
            while True:
                min_dist = float('inf')
                merge_pair = None
                
                # Check all pairs
                for i in range(len(clusters)):
                    for j in range(i + 1, len(clusters)):
                        # Calculate min distance between cluster i and cluster j
                        # We use the precomputed composite_dist matrix
                        # Extract submatrix for these two clusters
                        indices_i = clusters[i]
                        indices_j = clusters[j]
                        
                        # Use min() of the submatrix
                        # This is the "Single Linkage" distance between the two sets
                        dists = composite_dist[np.ix_(indices_i, indices_j)]
                        d = np.min(dists)
                        
                        if d < min_dist:
                            min_dist = d
                            merge_pair = (i, j)
                
                # Threshold for merging: 2.0x the DBSCAN epsilon
                # This allows for a "gap" twice as large as the internal density requirement
                merge_threshold = eps * 2.0
                
                if merge_pair and min_dist < merge_threshold:
                    i, j = merge_pair
                    if self.verbose:
                        print(f"    Merging cluster {i} and {j} (dist: {min_dist:.2f} < {merge_threshold:.2f})")
                    
                    # Merge j into i
                    clusters[i] = np.concatenate([clusters[i], clusters[j]])
                    clusters.pop(j)
                else:
                    break
        
        # Map local cluster indices back to global sample indices
        events = [sample_indices[cluster] for cluster in clusters]
        
        # If no valid events found, treat all as single event
        if len(events) == 0:
            events = [sample_indices]
        
        if self.verbose:
            print(f"Detected {len(events)} intrusion event(s)")
            for i, event in enumerate(events):
                print(f"  Event {i+1}: {len(event)} samples")
        
        return events
    
    def build_hierarchical_tree(self, sample_indices: np.ndarray) -> Dict[str, Any]:
        """
        Build hierarchical phylogenetic tree using agglomerative clustering.
        Finds most parsimonious common ancestors for samples, preferring deeper ancestors.
        
        WARNING: Computational complexity is O(n²) per merge iteration.
        For large datasets (>200 samples), enable detect_multiple_events to split into smaller groups.
        
        Returns:
        --------
        tree_info : Dict
            Contains:
            - 'nodes': List of all nodes (samples and ancestors)
            - 'edges': List of (parent_idx, child_idx) tuples
            - 'root_idx': Index of the root ancestor node
            - 'node_data': Dict mapping node_idx to {coords, value, is_sample, depth}
        """
        coords = self.sample_coords[sample_indices].astype(float)
        values = self.sample_values[sample_indices]
        n = len(coords)
        
        # Warn about performance for large datasets
        if n > 200:
            if self.verbose:
                print(f"  WARNING: Building tree for {n} samples - this may take several minutes!")
                print(f"  Recommendation: Enable 'detect_multiple_events' and reduce 'branch_threshold' to split into smaller groups")
        
        if n == 1:
            # Single sample: no ancestors needed, it's the root
            tree_info = {
                'nodes': [0],
                'edges': [],
                'root_idx': 0,
                'node_data': {
                    0: {
                        'coords': coords[0],
                        'value': values[0],
                        'is_sample': True,
                        'depth': coords[0][2],
                        'sample_indices': [sample_indices[0]]
                    }
                }
            }
            return tree_info
        
        # Initialize: each sample is its own cluster
        clusters = []
        for i in range(n):
            clusters.append({
                'members': [i],  # Original sample indices (0 to n-1)
                'coords': coords[i].copy(),
                'value': values[i],
                'is_sample': True,
                'node_idx': i  # Node index in final tree
            })
        
        # Track all nodes and edges
        all_nodes = list(range(n))  # Start with sample nodes
        all_edges = []
        node_data = {}
        
        # Initialize node data for samples
        for i in range(n):
            node_data[i] = {
                'coords': coords[i].copy(),
                'value': values[i],
                'is_sample': True,
                'depth': coords[i][2],
                'sample_indices': [sample_indices[i]]  # Original sample index
            }
        
        next_node_idx = n  # Next available node index for ancestors
        
        # Agglomerative clustering: repeatedly merge closest clusters
        merge_iteration = 0
        if self.verbose and n > 10:
            print(f"  Building hierarchical tree for {n} samples...")
        
        while len(clusters) > 1:
            merge_iteration += 1
            
            # Progress reporting for large datasets
            if self.verbose and n > 50 and merge_iteration % 10 == 0:
                print(f"    Merge iteration {merge_iteration}/{n-1}: {len(clusters)} clusters remaining")
            
            # Find the two closest clusters that satisfy angle constraint (VECTORIZED)
            nc = len(clusters)
            
            # Early optimization: for small cluster counts, use original loop
            if nc <= 10:
                min_dist = float('inf')
                merge_i, merge_j = 0, 1
                best_found = False
                
                for i in range(nc):
                    for j in range(i+1, nc):
                        if not self._check_branch_angle_valid(clusters[i], clusters[j]):
                            continue
                        dist = self._compute_parsimony_distance(clusters[i], clusters[j], coords, values)
                        if dist < min_dist:
                            min_dist = dist
                            merge_i, merge_j = i, j
                            best_found = True
            else:
                # Vectorized approach for large cluster counts
                cluster_coords = np.array([c['coords'] for c in clusters])
                cluster_values = np.array([c['value'] for c in clusters])
                
                # Compute all pairwise distances at once
                # Note: This vectorized calculation is an approximation of _compute_parsimony_distance
                # It assumes grid coordinates are proportional to world coordinates, which is true if block_size is uniform
                # For exact results, we should convert to world coords, but that's expensive inside the loop
                
                # Spatial distances: ||coords_i - coords_j||
                # We use grid coordinates here for speed, assuming they are roughly proportional
                diff = cluster_coords[:, np.newaxis, :] - cluster_coords[np.newaxis, :, :]
                spatial_dists = np.linalg.norm(diff, axis=2)
                
                # If block_size is available, scale spatial_dists to approximate world units
                if self.block_size is not None:
                    avg_block_size = np.mean(self.block_size)
                    spatial_dists_world = spatial_dists * avg_block_size
                else:
                    spatial_dists_world = spatial_dists
                
                # Value differences
                value_diffs = np.abs(cluster_values[:, np.newaxis] - cluster_values[np.newaxis, :])
                
                # Grade distance in spatial units
                grade_dists_world = value_diffs / self.mutation_rate
                
                # Composite distances (Divergence)
                composite_dists = np.sqrt(
                    self.spatial_weight * spatial_dists_world**2 + 
                    (self.attr_weight * grade_dists_world)**2
                )
                
                # Mask out invalid pairs (i >= j)
                mask = np.ones((nc, nc), dtype=bool)
                np.fill_diagonal(mask, False)
                mask[np.tril_indices(nc, -1)] = False  # Keep only upper triangle
                
                composite_dists[~mask] = np.inf
                
                # Find minimum valid distance
                min_dist = float('inf')
                merge_i, merge_j = 0, 1
                best_found = False
                
                # Sort pairs by distance and check angle constraint
                valid_pairs = np.where(mask)
                if len(valid_pairs[0]) > 0:
                    pair_dists = composite_dists[valid_pairs]
                    sorted_indices = np.argsort(pair_dists)
                    
                    # Check top candidates until we find a valid one
                    for idx in sorted_indices[:min(100, len(sorted_indices))]:  # Check top 100 candidates
                        i, j = valid_pairs[0][idx], valid_pairs[1][idx]
                        if self._check_branch_angle_valid(clusters[i], clusters[j]):
                            min_dist = pair_dists[idx]
                            merge_i, merge_j = i, j
                            best_found = True
                            break
            
            # If no valid merge found (all angles too shallow), break
            if not best_found:
                if self.verbose:
                    print(f"  Warning: No valid merges found satisfying {self.min_branch_angle}° angle constraint")
                    print(f"  Stopping with {len(clusters)} separate trees")
                break
            
            # Merge the two closest clusters
            cluster_i = clusters[merge_i]
            cluster_j = clusters[merge_j]
            
            # Compute optimal ancestor position and grade
            ancestor_coords, ancestor_value = self._compute_optimal_ancestor(
                cluster_i, cluster_j, coords, values
            )
            
            # Create new ancestor node
            ancestor_node_idx = next_node_idx
            next_node_idx += 1
            
            # Store node data
            node_data[ancestor_node_idx] = {
                'coords': ancestor_coords,
                'value': ancestor_value,
                'is_sample': False,
                'depth': ancestor_coords[2],
                'sample_indices': cluster_i['members'] + cluster_j['members']
            }
            
            # Add edges from ancestor to children
            all_edges.append((ancestor_node_idx, cluster_i['node_idx']))
            all_edges.append((ancestor_node_idx, cluster_j['node_idx']))
            all_nodes.append(ancestor_node_idx)
            
            # Create merged cluster
            merged_cluster = {
                'members': cluster_i['members'] + cluster_j['members'],
                'coords': ancestor_coords,
                'value': ancestor_value,
                'is_sample': False,
                'node_idx': ancestor_node_idx
            }
            
            # Remove old clusters and add merged one
            # Remove in reverse order to maintain indices
            if merge_j > merge_i:
                clusters.pop(merge_j)
                clusters.pop(merge_i)
            else:
                clusters.pop(merge_i)
                clusters.pop(merge_j)
            clusters.append(merged_cluster)
        
        # Final cluster is the root
        root_idx = clusters[0]['node_idx']
        
        if self.verbose:
            print(f"  Built hierarchical tree: {n} samples, {len(all_nodes)} total nodes, {len(all_edges)} edges")
            print(f"  Root ancestor at: {node_data[root_idx]['coords']} (depth: {node_data[root_idx]['depth']:.2f})")
        
        tree_info = {
            'nodes': all_nodes,
            'edges': all_edges,
            'root_idx': root_idx,
            'node_data': node_data
        }
        
        return tree_info
    
    def _compute_parsimony_distance(self, cluster_i: Dict, cluster_j: Dict, 
                                     coords: np.ndarray, values: np.ndarray) -> float:
        """
        Compute parsimony distance between two clusters using inferred mutation rate.
        Distance represents the "Evolutionary Distance" (divergence).
        """
        # Spatial distance between cluster centers (Horizontal + Vertical)
        # Note: We use full 3D distance as the baseline spatial separation
        coords_world_i = self._grid_to_world(cluster_i['coords'])
        coords_world_j = self._grid_to_world(cluster_j['coords'])
        spatial_dist = np.linalg.norm(coords_world_i - coords_world_j)
        
        # Attribute distance converted to spatial units
        value_diff = abs(cluster_i['value'] - cluster_j['value'])
        grade_dist = value_diff / self.mutation_rate
        
        # Apply user weights
        # We treat grade_dist as an orthogonal dimension to spatial_dist
        # Total Divergence = sqrt(Spatial^2 + (Weight * Grade_Spatial)^2)
        
        divergence = np.sqrt(
            self.spatial_weight * spatial_dist**2 + 
            (self.attr_weight * grade_dist)**2
        )
        
        return divergence
    
    def _check_branch_angle_valid(self, cluster_i: Dict, cluster_j: Dict) -> bool:
        """
        Check if branches from ancestor to both children would satisfy minimum angle constraint.
        With divergence-based depth, this is less critical but still useful to prevent
        extremely flat branches if grade difference is small.
        """
        if self.min_branch_angle <= 0:
            return True  # No constraint
        
        # Compute potential ancestor position in world coordinates
        coords_world_i = self._grid_to_world(cluster_i['coords'])
        coords_world_j = self._grid_to_world(cluster_j['coords'])
        
        # Use the new divergence-based ancestor logic to check the angle
        ancestor_coords_grid, _ = self._compute_optimal_ancestor(
            cluster_i, cluster_j, None, None
        )
        ancestor_world = self._grid_to_world(ancestor_coords_grid)
        
        # Check angle from ancestor to each child
        for child_world in [coords_world_i, coords_world_j]:
            # Vector from ancestor to child
            branch_vector = child_world - ancestor_world
            
            # Horizontal distance (XY plane)
            horizontal_dist = np.linalg.norm(branch_vector[:2])
            
            # Vertical distance (Z axis)
            vertical_dist = abs(branch_vector[2])
            
            if horizontal_dist < 1e-6:
                continue
                
            angle_degrees = np.degrees(np.arctan(vertical_dist / horizontal_dist))
            
            if angle_degrees < self.min_branch_angle:
                return False
        
        return True
    
    def _compute_optimal_ancestor(self, cluster_i: Dict, cluster_j: Dict,
                                   coords: np.ndarray, values: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Compute optimal ancestor position and grade for two clusters.
        Uses Divergence-Based Depth:
        - Horizontal: Midpoint
        - Depth: Derived from total divergence (Pythagorean theorem)
        """
        coords_i = cluster_i['coords']
        coords_j = cluster_j['coords']
        value_i = cluster_i['value']
        value_j = cluster_j['value']
        
        # 1. Horizontal Position: Midpoint
        centroid = (coords_i + coords_j) / 2.0
        avg_value = (value_i + value_j) / 2.0
        
        # Convert to world coords for depth calculation
        coords_world_i = self._grid_to_world(coords_i)
        coords_world_j = self._grid_to_world(coords_j)
        
        # 2. Calculate Divergence (Target Branch Length)
        # Spatial separation (Horizontal only for the base of the triangle)
        # We project the problem into a 2D plane: Horizontal Separation vs Depth
        horizontal_sep = np.linalg.norm(coords_world_i[:2] - coords_world_j[:2])
        
        # Grade separation converted to spatial units
        grade_diff = abs(value_i - value_j)
        grade_dist = (grade_diff / self.mutation_rate) * self.attr_weight
        
        # Total "Evolutionary Path" required (Divergence)
        # We assume the path length is driven by the grade difference + spatial separation
        # If we assume a clock, the path from Ancestor to Child is half the total divergence
        # But we need to account for the existing Z separation too
        z_sep = abs(coords_world_i[2] - coords_world_j[2])
        
        # Total divergence between children
        total_divergence = np.sqrt(horizontal_sep**2 + z_sep**2 + grade_dist**2)
        
        # Target branch length (Ancestor -> Child)
        target_branch_len = total_divergence / 2.0
        
        # 3. Calculate Required Vertical Drop
        # We place ancestor at horizontal midpoint.
        # Horizontal distance from midpoint to child
        dist_to_midpoint_horiz = horizontal_sep / 2.0
        
        # We need to find V (vertical drop from the *deepest* child)
        # Actually, let's place it relative to the Z-midpoint first
        z_midpoint = (coords_world_i[2] + coords_world_j[2]) / 2.0
        
        # The ancestor is at (Mid_X, Mid_Y, Z_anc)
        # Distance to child i: sqrt( (dist_to_midpoint_horiz)^2 + (Z_anc - Z_i)^2 )
        # We want this to equal target_branch_len (ignoring grade dimension for the physical placement)
        # Wait, the grade dimension IS the reason we push it deep.
        # The "Physical Path" length in 3D space should explain the divergence?
        # No, the "Evolutionary Path" explains the divergence.
        # We want PhysicalPath to be such that it fits?
        
        # Let's use the "Excess Divergence" logic:
        # The grade difference implies a certain "Time".
        # That time must be spent traveling.
        # If horizontal travel isn't enough, we must travel vertically.
        
        # Minimum physical distance required to explain the grade change?
        # min_physical_dist = grade_dist (if rate is 1:1)
        
        # Let's simply use the Pythagorean approach on the "Missing Leg":
        # We have a triangle with Hypotenuse = target_branch_len
        # Base = dist_to_midpoint_horiz
        # Height = V (vertical distance from Z-midpoint to Ancestor)
        
        # V^2 + Base^2 = Hypotenuse^2
        # V = sqrt( Hypotenuse^2 - Base^2 )
        
        term = target_branch_len**2 - dist_to_midpoint_horiz**2
        if term > 0:
            vertical_drop_from_midpoint = np.sqrt(term)
        else:
            vertical_drop_from_midpoint = 0
            
        # Calculate target depth
        # We want the ancestor to be deeper than the midpoint by V
        # But we also want to respect the "flowed upwards" constraint (deeper than BOTH if possible)
        # And we have a minimum offset
        
        min_offset_world = self.ancestor_depth_offset * abs(self.block_size[2] if self.block_size is not None else 1.0)
        
        # Base depth is the deepest child (to ensure we go deeper)
        min_child_z = min(coords_world_i[2], coords_world_j[2])
        
        # Calculate proposed Z
        # If we use Z-midpoint - V, we might end up shallower than the deepest child if V is small
        # So let's force it to be at least min_offset below the deepest child
        
        proposed_z = z_midpoint - vertical_drop_from_midpoint
        
        # Enforce minimum depth constraint (Geological constraint)
        max_allowed_z = min_child_z - min_offset_world
        
        final_z = min(proposed_z, max_allowed_z)
        
        # Do NOT clamp to grid bottom - allow ancestors to be "virtual" deep sources
        # min_world = self.min_bounds if self.min_bounds is not None else np.zeros(3)
        # final_z = max(final_z, min_world[2])
        
        # Construct result
        ancestor_world = self._grid_to_world(centroid)
        ancestor_world[2] = final_z
        
        # Convert back to grid coordinates manually to avoid Z-clipping
        if self.min_bounds is not None and self.block_size is not None:
            rel = (ancestor_world - self.min_bounds) / self.block_size
            # Clip X/Y to grid, but allow Z to be negative (deep)
            if self.dims is not None:
                rel[:2] = np.clip(rel[:2], 0, self.dims[:2] - 1)
            ancestor_coords = np.round(rel).astype(int)
        else:
            ancestor_coords = np.round(ancestor_world).astype(int)
            
        # Option 2: Snap ancestor to nearest valid domain block
        if hasattr(self, 'domain_kdtree') and self.domain_kdtree is not None:
            # Check if current ancestor is already in domain
            # (Convert to tuple for set lookup)
            anc_tuple = tuple(ancestor_coords)
            if hasattr(self, 'allowed_grid_override') and self.allowed_grid_override:
                if anc_tuple not in self.allowed_grid_override:
                    # Not in domain, find nearest valid block
                    # Note: self.domain_kdtree is built from the current sub-volume's blocks
                    # so this ensures we snap to the same closed volume as the samples.
                    dist, idx = self.domain_kdtree.query(ancestor_coords)
                    snapped_coords = self.domain_grid_points[idx]
                    
                    if self.verbose:
                        print(f"      Snapping ancestor from {ancestor_coords} to {snapped_coords} (dist: {dist:.2f}) to stay in sub-volume")
                    
                    ancestor_coords = snapped_coords
                    # Note: We keep the calculated avg_value, assuming the "virtual" ancestor
                    # projects its grade to this nearest physical point.
        
        return ancestor_coords, avg_value

    
    def interpolate_hierarchical_tree(self, tree_info: Dict[str, Any],
                                       sample_indices: np.ndarray,
                                       event_id: int,
                                       domain: str = 'Undomained'):
        """
        Interpolate values along hierarchical tree from ancestors to samples.
        Traverses parent→child edges and creates blocks along each path.
        """
        node_data = tree_info['node_data']
        edges = tree_info['edges']
        root_idx = tree_info['root_idx']
        
        n_samples = sum(1 for nd in node_data.values() if nd['is_sample'])
        print(f"  Interpolating hierarchical tree: {n_samples} samples, {len(node_data)} nodes, {len(edges)} edges...")
        
        # Get allowed positions if filtering is enabled
        allowed_positions = getattr(self, 'allowed_grid_override', None)
        if allowed_positions:
            print(f"  Domain filtering: ENABLED ({len(allowed_positions)} allowed positions)")
        else:
            print(f"  Domain filtering: DISABLED (blocks can be created anywhere)")
        
        total_paths = 0
        total_blocks_created = 0
        total_blocks_filtered = 0
        
        # Create or update the root ancestor block
        root_data = node_data[root_idx]
        root_coords = root_data['coords']
        root_value = root_data['value']
        root_pos = tuple(root_coords.astype(int))
        
        if root_pos not in self.blocks:
            self.blocks[root_pos] = {
                'value': root_value,
                'is_sample': False,
                'event_id': event_id,
                'distance_to_feeder': 0.0,
                'branch_id': 0,
                'is_feeder': True,
                'visit_count': 1,
                'domain': domain
            }
        else:
            # Update existing block (might be a sample)
            block = self.blocks[root_pos]
            if not block.get('is_sample', False):
                block['value'] = root_value
            block['event_id'] = event_id
            block['distance_to_feeder'] = 0.0
            block['branch_id'] = 0
            block['is_feeder'] = True
            block['domain'] = domain
            # Don't reset visit_count, just update properties
        
        # Build adjacency list: parent -> children
        children_map = {}
        for parent_idx, child_idx in edges:
            if parent_idx not in children_map:
                children_map[parent_idx] = []
            children_map[parent_idx].append(child_idx)
        
        # Traverse tree using BFS from root
        visited = set()
        queue = [(root_idx, 0.0, 0)]  # (node_idx, cumulative_distance, branch_id)
        next_branch_id = 1
        
        while queue:
            current_idx, parent_distance, current_branch_id = queue.pop(0)
            
            if current_idx in visited:
                continue
            visited.add(current_idx)
            
            current_data = node_data[current_idx]
            current_coords = current_data['coords']
            current_value = current_data['value']
            current_is_sample = current_data['is_sample']
            
            # Get children of current node
            children = children_map.get(current_idx, [])
            
            # Process each child edge
            for i, child_idx in enumerate(children):
                child_data = node_data[child_idx]
                child_coords = child_data['coords']
                child_value = child_data['value']
                
                # Determine branch ID: new branch for each child after the first
                if i == 0:
                    edge_branch_id = current_branch_id
                else:
                    edge_branch_id = next_branch_id
                    next_branch_id += 1
                
                # Create path from parent to child
                path_blocks = self._create_path_blocks(
                    current_coords, child_coords,
                    current_value, child_value,
                    event_id, edge_branch_id
                )
                
                if path_blocks:
                    total_paths += 1
                
                blocks_created = 0
                blocks_filtered = 0
                
                for block_pos, block_value, dist_to_parent in path_blocks:
                    # Filter out blocks outside the grid (e.g. deep virtual paths)
                    if self.dims is not None:
                        if not (0 <= block_pos[0] < self.dims[0] and 
                                0 <= block_pos[1] < self.dims[1] and 
                                0 <= block_pos[2] < self.dims[2]):
                            blocks_filtered += 1
                            continue

                    if allowed_positions is not None and block_pos not in allowed_positions:
                        blocks_filtered += 1
                        continue
                    
                    # Calculate cumulative distance from root
                    cumulative_dist = parent_distance + dist_to_parent
                    
                    # Check if block already exists
                    if block_pos not in self.blocks:
                        self.blocks[block_pos] = {
                            'value': block_value,
                            'is_sample': False,
                            'event_id': event_id,
                            'distance_to_feeder': cumulative_dist,
                            'branch_id': edge_branch_id,
                            'is_feeder': False,
                            'visit_count': 1,
                            'domain': domain
                        }
                        blocks_created += 1
                    else:
                        # Update existing block if not a sample
                        existing = self.blocks[block_pos]
                        if not existing.get('is_sample', False):
                            # Magma Mixing: Average the values
                            old_count = existing.get('visit_count', 1)
                            new_count = old_count + 1
                            
                            # Running average
                            existing['value'] = (existing['value'] * old_count + block_value) / new_count
                            existing['visit_count'] = new_count
                            existing['domain'] = domain
                            
                            # Keep the shortest path to feeder (dominant flow)
                            if cumulative_dist < existing['distance_to_feeder']:
                                existing['distance_to_feeder'] = cumulative_dist
                                existing['event_id'] = event_id
                                existing['branch_id'] = edge_branch_id
                
                total_blocks_created += blocks_created
                total_blocks_filtered += blocks_filtered
                
                # Calculate edge distance for next iteration
                edge_distance = np.linalg.norm(child_coords - current_coords)
                
                # Add child to queue
                queue.append((child_idx, parent_distance + edge_distance, edge_branch_id))
        
        print(f"  Created {total_blocks_created} blocks across {total_paths} paths")
        if total_blocks_filtered > 0:
            print(f"  Filtered {total_blocks_filtered} blocks (outside allowed domains)")
    
        """Create or tag the synthetic ancestor node for an event."""
        coords = self.sample_coords[sample_indices].astype(float)
        values = self.sample_values[sample_indices]
        coords_world = self._grid_to_world(coords)
        if self.dims is not None:
            grid_dims = self.dims
        else:
            grid_dims = np.ceil(coords.max(axis=0) + 1).astype(int)
        if len(coords) == 0:
            raise ValueError("Cannot compute ancestor without samples")
        ancestor_world = coords_world.mean(axis=0)
        depth_offset_world = max(self.ancestor_depth_offset, 0.0) * abs(self.block_size[2] if self.block_size is not None else 1.0)
        min_world = self.min_bounds if self.min_bounds is not None else np.zeros(3)
        max_world = self._grid_to_world(grid_dims - 1)
        deepest = coords_world[:, 2].min()
        target_depth = max(deepest - depth_offset_world, min_world[2])
        ancestor_world[2] = target_depth
        ancestor_world = np.clip(ancestor_world, min_world, max_world)
        ancestor_idx = self._world_to_grid(ancestor_world)
        ancestor_center = np.array(ancestor_idx, dtype=float)
        dists = np.linalg.norm(coords_world - ancestor_world, axis=1)
        weights = 1.0 / (dists + 1e-6)
        ancestor_value = float(np.average(values, weights=weights))
        matched_sample = False
        existing_block = self.blocks.get(ancestor_idx)
        if existing_block is not None and existing_block.get('is_sample', False):
            ancestor_value = existing_block['value']
            matched_sample = True
        else:
            if ancestor_idx not in self.blocks:
                self.blocks[ancestor_idx] = {
                    'value': ancestor_value,
                    'is_sample': False,
                    'event_id': event_id,
                    'distance_to_feeder': 0,
                    'branch_id': 0,
                    'is_feeder': True
                }
            else:
                block = self.blocks[ancestor_idx]
                block.update({
                    'value': ancestor_value,
                    'is_sample': False,
                    'event_id': event_id,
                    'distance_to_feeder': 0,
                    'branch_id': 0,
                    'is_feeder': True
                })
        ancestor_block = self.blocks[ancestor_idx]
        ancestor_block['event_id'] = event_id
        ancestor_block['branch_id'] = 0
        ancestor_block['distance_to_feeder'] = 0
        ancestor_block['is_feeder'] = True
        ancestor_info = {
            'event_id': event_id,
            'grid_coords': ancestor_center,
            'block_idx': ancestor_idx,
            'value': ancestor_value,
            'matched_sample': matched_sample,
            'world_coords': self._grid_to_world(ancestor_center)
        }
        self.ancestor_nodes = [node for node in self.ancestor_nodes if node.get('event_id') != event_id]
        self.ancestor_nodes.append(ancestor_info)
        return ancestor_info
    
    def interpolate_along_branches(self, feeder_coords: np.ndarray, 
                                   tree_structure: np.ndarray, 
                                   sample_indices: np.ndarray,
                                   event_id: int,
                                   ancestor_info: Dict[str, Any]):
        """
        Interpolate values along tree branches from feeder to samples.
        Creates blocks along each path in the MST.
        """
        coords = self.sample_coords[sample_indices]
        values = self.sample_values[sample_indices]
        n = len(coords)
        
        print(f"  Interpolating branches for {n} samples...")
        
        # Get allowed positions if filtering is enabled
        allowed_positions = getattr(self, 'allowed_grid_override', None)
        if allowed_positions:
            print(f"  Domain filtering: ENABLED ({len(allowed_positions)} allowed positions)")
        else:
            print(f"  Domain filtering: DISABLED (blocks can be created anywhere)")
        
        total_paths = 0
        total_blocks_created = 0
        total_blocks_filtered = 0
        
        # Mark feeder
        feeder_idx = np.argmin(np.linalg.norm(coords - feeder_coords, axis=1))
        feeder_pos = tuple(coords[feeder_idx].astype(int))
        if feeder_pos in self.blocks:
            self.blocks[feeder_pos]['is_feeder'] = True
            self.blocks[feeder_pos]['event_id'] = event_id
            self.blocks[feeder_pos]['branch_id'] = 0
        ancestor_coords = np.array(ancestor_info['grid_coords'], dtype=float)
        ancestor_value = ancestor_info['value']
        ancestor_block_idx = ancestor_info['block_idx']
        if ancestor_block_idx in self.blocks:
            self.blocks[ancestor_block_idx]['event_id'] = event_id
            self.blocks[ancestor_block_idx]['is_feeder'] = True
            self.blocks[ancestor_block_idx]['branch_id'] = 0
            self.blocks[ancestor_block_idx]['distance_to_feeder'] = 0
        trunk_blocks = self._create_path_blocks(
            ancestor_coords, coords[feeder_idx].astype(float),
            ancestor_value, values[feeder_idx],
            event_id, 0
        )
        if trunk_blocks:
            total_paths += 1
        blocks_created = 0
        blocks_filtered = 0
        parent_distance = 0.0
        for block_pos, block_value, dist_to_feeder in trunk_blocks:
            if allowed_positions is not None and block_pos not in allowed_positions:
                blocks_filtered += 1
                continue
            if block_pos not in self.blocks:
                self.blocks[block_pos] = {
                    'value': block_value,
                    'is_sample': False,
                    'event_id': event_id,
                    'distance_to_feeder': parent_distance + dist_to_feeder,
                    'branch_id': 0,
                    'is_feeder': False
                }
                blocks_created += 1
        total_blocks_created += blocks_created
        total_blocks_filtered += blocks_filtered
        if feeder_pos in self.blocks:
            feeder_distance = np.linalg.norm(coords[feeder_idx].astype(float) - ancestor_coords)
            self.blocks[feeder_pos]['distance_to_feeder'] = feeder_distance
        
        # Build adjacency list from MST (VECTORIZED)
        adjacency = {i: [] for i in range(n)}
        # Find all edges in tree_structure
        edges_i, edges_j = np.where(tree_structure > 0)
        for i, j in zip(edges_i, edges_j):
            adjacency[i].append(j)
        
        # BFS from feeder to trace all branches
        visited = set()
        queue = [(feeder_idx, None, 0)]  # (node_index, parent_index, branch_id)
        next_branch_id = 1
        
        while queue:
            current_idx, parent_idx, current_branch_id = queue.pop(0)
            
            if current_idx in visited:
                continue
            visited.add(current_idx)
            
            current_pos = coords[current_idx]
            current_value = values[current_idx]
            
            # Mark sample with event and branch
            sample_pos = tuple(current_pos.astype(int))
            if sample_pos in self.blocks:
                self.blocks[sample_pos]['event_id'] = event_id
                self.blocks[sample_pos]['branch_id'] = (current_branch_id
                                                        if current_branch_id is not None else -1)
                if parent_idx is None:
                    # Root sample gets distance from ancestor trunk
                    existing_distance = self.blocks[sample_pos].get('distance_to_feeder', 0.0)
                    self.blocks[sample_pos]['distance_to_feeder'] = existing_distance
            
            # Create interpolated blocks from parent to current
            if parent_idx is not None:
                total_paths += 1
                parent_pos = coords[parent_idx]
                parent_value = values[parent_idx]
                parent_block_idx = tuple(parent_pos.astype(int))
                parent_distance = self.blocks.get(parent_block_idx, {}).get('distance_to_feeder', 0.0)
                
                path_blocks = self._create_path_blocks(
                    parent_pos, current_pos,
                    parent_value, current_value,
                    event_id, current_branch_id if current_branch_id is not None else -1
                )
                
                blocks_created = 0
                blocks_filtered = 0
                for block_pos, block_value, dist_to_feeder in path_blocks:
                    # Check if position is in allowed grid (if filtering enabled)
                    if allowed_positions is not None and block_pos not in allowed_positions:
                        blocks_filtered += 1
                        continue
                    
                    if block_pos not in self.blocks:
                        self.blocks[block_pos] = {
                            'value': block_value,
                            'is_sample': False,
                            'event_id': event_id,
                            'distance_to_feeder': parent_distance + dist_to_feeder,
                            'branch_id': current_branch_id if current_branch_id is not None else -1,
                            'is_feeder': False
                        }
                        blocks_created += 1
                
                total_blocks_created += blocks_created
                total_blocks_filtered += blocks_filtered
                
                if self.verbose and (blocks_created > 0 or blocks_filtered > 0):
                    print(f"    Branch: created {blocks_created} blocks, filtered {blocks_filtered} (not in domain)")

                # Update sample distance using parent's cumulative distance
                if sample_pos in self.blocks:
                    edge_length = np.linalg.norm(current_pos - parent_pos)
                    self.blocks[sample_pos]['distance_to_feeder'] = parent_distance + edge_length

            
            # Add children to queue
            for neighbor_idx in adjacency[current_idx]:
                if neighbor_idx not in visited:
                    neighbor_branch_id = next_branch_id
                    next_branch_id += 1
                    queue.append((neighbor_idx, current_idx, neighbor_branch_id))
        
        print(f"  Summary: {total_paths} paths processed, {total_blocks_created} blocks created, {total_blocks_filtered} blocks filtered")
    
    def _create_path_blocks(self, start: np.ndarray, end: np.ndarray,
                           start_value: float, end_value: float,
                           event_id: int, branch_id: int) -> List[Tuple[Tuple[int, int, int], float, float]]:
        """
        Create interpolated blocks along 3D line from start to end.
        
        Returns:
        --------
        List of (block_position, interpolated_value, distance_to_feeder)
        """
        start = np.array(start)
        end = np.array(end)
        
        total_dist = np.linalg.norm(end - start)
        
        if self.verbose:
            print(f"    Path from {start} to {end}: distance = {total_dist:.2f}")
        
        if total_dist < 0.1:  # Same block - no interpolation needed
            if self.verbose:
                print(f"      Skipped (same block)")
            return []
        
        # Number of steps based on distance (ensure at least 1 step)
        num_steps = max(2, int(np.ceil(total_dist * 2)))  # 2 steps per grid unit
        path_blocks = []
        seen_positions = set()
        
        if self.verbose:
            print(f"      Creating {num_steps-1} intermediate blocks...")
        
        for step in range(1, num_steps):  # Skip start (0) and end (num_steps)
            t = step / num_steps
            
            # Linear interpolation of position
            pos = start + t * (end - start)
            block_pos = tuple(np.round(pos).astype(int))
            
            # Skip if we've already seen this position
            if block_pos in seen_positions:
                continue
            seen_positions.add(block_pos)
            
            # Skip if this is the start or end position
            if block_pos == tuple(start.astype(int)) or block_pos == tuple(end.astype(int)):
                continue
            
            # Interpolate value based on method
            if self.interpolation_method == "linear":
                interpolated_value = start_value + t * (end_value - start_value)
            elif self.interpolation_method == "inverse_distance":
                # Inverse distance weighting
                dist_to_start = t * total_dist
                dist_to_end = (1 - t) * total_dist
                w_start = 1 / (dist_to_start + 1e-6)
                w_end = 1 / (dist_to_end + 1e-6)
                interpolated_value = (w_start * start_value + w_end * end_value) / (w_start + w_end)
            else:
                interpolated_value = start_value + t * (end_value - start_value)
            
            # Distance from feeder (approximate)
            dist_to_feeder = t * total_dist
            
            path_blocks.append((block_pos, interpolated_value, dist_to_feeder))
        
        if self.verbose:
            print(f"      Created {len(path_blocks)} unique blocks")
        
        return path_blocks
    
    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        """
        Execute the molecular clock algorithm.
        This runs once as it's a deterministic algorithm.
        
        Returns False after completion to signal no more iterations needed.
        """
        if self.converged:
            return False
        
        if self.verbose:
            print("\n=== Running Molecular Clock Interpolation ===")
        
        # Get unique domains
        unique_domains = np.unique(self.sample_domains)
        if len(unique_domains) == 0:
            unique_domains = ['']
            
        print(f"Processing {len(unique_domains)} domains: {unique_domains}")
        
        # Store original allowed_grid_override
        original_allowed = getattr(self, 'allowed_grid_override', None)
        
        for domain in unique_domains:
            if self.verbose:
                print(f"\n=== Processing Domain: {domain} ===")
            
            # Filter samples for this domain
            domain_mask = (self.sample_domains == domain)
            domain_indices = np.where(domain_mask)[0]
            
            if len(domain_indices) == 0:
                if self.verbose:
                    print(f"  No samples in domain {domain}, skipping.")
                continue
                
            # Determine allowed blocks for this domain
            domain_allowed = None
            if original_allowed is not None:
                if hasattr(self, 'domain_mapping') and self.domain_mapping:
                    domain_allowed = {pos for pos in original_allowed if self.domain_mapping.get(pos, 'Undomained') == domain}
                else:
                    domain_allowed = original_allowed
            
            # Identify disconnected sub-volumes within the domain
            sub_volumes = []
            
            if domain_allowed and len(domain_allowed) > 0 and self.dims is not None:
                # Check if grid is too massive for dense array
                total_cells = np.prod(self.dims)
                if total_cells > 1e8: # 100 million cells limit
                    if self.verbose:
                        print(f"  Grid too large ({total_cells} cells) for volume analysis. Treating domain as single volume.")
                    sub_volumes = [(domain_allowed, domain_indices)]
                else:
                    try:
                        grid_mask = np.zeros(self.dims, dtype=bool)
                        for pos in domain_allowed:
                            if 0 <= pos[0] < self.dims[0] and 0 <= pos[1] < self.dims[1] and 0 <= pos[2] < self.dims[2]:
                                grid_mask[pos] = True
                        
                        labeled_array, num_features = label(grid_mask)
                        
                        if self.verbose:
                            print(f"  Domain '{domain}' has {num_features} disconnected sub-volumes.")
                        
                        if num_features <= 1:
                             sub_volumes = [(domain_allowed, domain_indices)]
                        else:
                            # Group blocks by volume
                            volume_blocks_map = {i: set() for i in range(1, num_features + 1)}
                            for pos in domain_allowed:
                                if 0 <= pos[0] < self.dims[0] and 0 <= pos[1] < self.dims[1] and 0 <= pos[2] < self.dims[2]:
                                    lbl = labeled_array[pos]
                                    if lbl > 0:
                                        volume_blocks_map[lbl].add(pos)
                            
                            # Assign samples to volumes
                            volume_samples_map = {i: [] for i in range(1, num_features + 1)}
                            unassigned_samples = []
                            
                            for idx in domain_indices:
                                sample_pos = tuple(self.sample_coords[idx])
                                if 0 <= sample_pos[0] < self.dims[0] and 0 <= sample_pos[1] < self.dims[1] and 0 <= sample_pos[2] < self.dims[2]:
                                    lbl = labeled_array[sample_pos]
                                    if lbl > 0:
                                        volume_samples_map[lbl].append(idx)
                                    else:
                                        unassigned_samples.append(idx)
                                else:
                                    unassigned_samples.append(idx)
                            
                            # Create sub-volume list
                            for i in range(1, num_features + 1):
                                if len(volume_samples_map[i]) > 0:
                                    sub_volumes.append((volume_blocks_map[i], np.array(volume_samples_map[i])))
                            
                            # Handle unassigned samples
                            if unassigned_samples:
                                if self.verbose:
                                    print(f"  {len(unassigned_samples)} samples fell outside defined volumes. Processing separately.")
                                sub_volumes.append((domain_allowed, np.array(unassigned_samples)))

                    except Exception as e:
                        print(f"  Error in volume analysis: {e}. Treating domain as single volume.")
                        sub_volumes = [(domain_allowed, domain_indices)]
            else:
                sub_volumes = [(domain_allowed, domain_indices)]
            
            # Process each sub-volume
            for vol_idx, (vol_allowed, vol_indices) in enumerate(sub_volumes):
                if self.verbose and len(sub_volumes) > 1:
                    print(f"  Processing Sub-volume {vol_idx + 1}/{len(sub_volumes)}: {len(vol_indices)} samples")
                
                # Set allowed grid for this sub-volume
                self.allowed_grid_override = vol_allowed
                
                # Build KDTree for domain snapping (Option 2: Force ancestors into domain)
                if vol_allowed:
                    # Convert set of tuples to array for KDTree
                    # This represents the valid blocks for the CURRENT SUB-VOLUME only
                    self.domain_grid_points = np.array(list(vol_allowed))
                    self.domain_kdtree = cKDTree(self.domain_grid_points)
                else:
                    self.domain_kdtree = None
                    self.domain_grid_points = None
                
                # Detect intrusion events for this sub-volume
                events = self.detect_intrusion_events(sample_indices=vol_indices)
                
                # Build hierarchical tree for each event
                for event_id, sample_indices in enumerate(events):
                    num_samples = len(sample_indices)
                    
                    # Check if event is too large and needs subsampling
                    if self.max_samples_per_event > 0 and num_samples > self.max_samples_per_event:
                        if self.verbose:
                            print(f"\nEvent {event_id + 1}/{len(events)}: {num_samples} samples (exceeds limit {self.max_samples_per_event})")
                            print(f"  Randomly subsampling to {self.max_samples_per_event} samples...")
                        # Randomly subsample (truly random, no fixed seed)
                        sample_indices = np.random.choice(sample_indices, size=self.max_samples_per_event, replace=False)
                        num_samples = len(sample_indices)
                    
                    if self.verbose:
                        print(f"\nProcessing Event {event_id + 1}/{len(events)}: {num_samples} samples")
                    
                    # Build hierarchical phylogenetic tree
                    tree_info = self.build_hierarchical_tree(sample_indices)
                    
                    # Store tree structure
                    self.intrusion_trees.append({
                        'event_id': event_id,
                        'domain': domain,
                        'volume_id': vol_idx,
                        'tree_info': tree_info,
                        'sample_indices': sample_indices,
                        'num_samples': len(sample_indices)
                    })
                    
                    # Interpolate along hierarchical branches
                    self.interpolate_hierarchical_tree(tree_info, sample_indices, event_id, domain=domain)
        
        # Restore allowed_grid_override
        if original_allowed is not None:
            self.allowed_grid_override = original_allowed
        
        # Fill background blocks if requested
        if self.fill_background:
            self._fill_background_blocks()
        
        self.converged = True
        
        if self.verbose:
            print(f"\n=== Interpolation Complete ===")
            print(f"Total blocks created: {len(self.blocks)}")
            print(f"Intrusion events: {len(self.intrusion_trees)}")
        
        return False  # Signal completion
    
    def _fill_background_blocks(self):
        """Fill remaining grid positions with configured background value."""
        allowed_positions = getattr(self, 'allowed_grid_override', None)
        if allowed_positions is not None:
            candidate_positions = allowed_positions
        elif self.dims is not None:
            candidate_positions = np.ndindex(tuple(self.dims))
        else:
            if self.verbose:
                print("    Skipped background fill (unknown grid dimensions)")
            return
        filled = 0
        for pos in candidate_positions:
            if pos in self.blocks:
                continue
            
            domain = 'Undomained'
            if hasattr(self, 'domain_mapping') and self.domain_mapping:
                domain = self.domain_mapping.get(pos, 'Undomained')
                
            self.blocks[pos] = {
                'value': self.background_value,
                'is_sample': False,
                'event_id': -1,
                'distance_to_feeder': np.inf,
                'branch_id': -1,
                'is_feeder': False,
                'domain': domain
            }
            filled += 1
        if self.verbose:
            print(f"    Background fill added {filled} blocks with value {self.background_value}")
    
    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        """Return interpolated block values"""
        return {pos: block['value'] for pos, block in self.blocks.items()}
    
    def get_algorithm_name(self) -> str:
        return "Molecular Clock (Phylogeographic MST)"
    
    def is_converged(self) -> bool:
        """Check if algorithm has finished"""
        return self.converged
    
    def get_metadata(self) -> Dict[str, Any]:
        """Return algorithm-specific metadata including tree structures"""
        return {
            'algorithm': self.get_algorithm_name(),
            'num_intrusion_events': len(self.intrusion_trees),
            'intrusion_trees': self.intrusion_trees,
            'spatial_weight': self.spatial_weight,
            'attribute_weight': self.attr_weight,
            'total_blocks': len(self.blocks),
            'sample_blocks': sum(1 for b in self.blocks.values() if b['is_sample']),
            'interpolated_blocks': sum(1 for b in self.blocks.values() if not b['is_sample']),
            'ancestor_nodes': self.ancestor_nodes
        }
    
    def get_intrusion_trees(self) -> List[Dict[str, Any]]:
        """Return tree structures for visualization/export"""
        return self.intrusion_trees
    
    def export_tree_structure(self, filepath: str):
        """
        Export hierarchical tree structures to CSV for visualization in Leapfrog.
        Creates a CSV with all nodes (samples and ancestors) and their connections.
        """
        import pandas as pd
        
        branch_data = []
        
        for tree in self.intrusion_trees:
            event_id = tree['event_id']
            tree_info = tree['tree_info']
            node_data = tree_info['node_data']
            edges = tree_info['edges']
            root_idx = tree_info['root_idx']
            
            # Export all nodes
            for node_idx, data in node_data.items():
                coords = data['coords']
                world_coords = self._grid_to_world(coords)
                
                # Determine point type
                if data['is_sample']:
                    point_type = 'Sample'
                elif node_idx == root_idx:
                    point_type = 'Root_Ancestor'
                else:
                    point_type = 'Intermediate_Ancestor'
                
                branch_data.append({
                    'Event_ID': event_id,
                    'Node_ID': node_idx,
                    'Point_Type': point_type,
                    'X': world_coords[0],
                    'Y': world_coords[1],
                    'Z': world_coords[2],
                    'Depth': data['depth'],
                    'Value': data['value'],
                    'From_Node_ID': None,
                    'From_X': None,
                    'From_Y': None,
                    'From_Z': None
                })
            
            # Export edges (parent->child connections)
            for parent_idx, child_idx in edges:
                parent_data = node_data[parent_idx]
                child_data = node_data[child_idx]
                
                parent_coords = self._grid_to_world(parent_data['coords'])
                child_coords = self._grid_to_world(child_data['coords'])
                
                branch_data.append({
                    'Event_ID': event_id,
                    'Node_ID': child_idx,
                    'Point_Type': 'Edge',
                    'X': child_coords[0],
                    'Y': child_coords[1],
                    'Z': child_coords[2],
                    'Depth': child_data['depth'],
                    'Value': child_data['value'],
                    'From_Node_ID': parent_idx,
                    'From_X': parent_coords[0],
                    'From_Y': parent_coords[1],
                    'From_Z': parent_coords[2]
                })
        
        df = pd.DataFrame(branch_data)
        df.to_csv(filepath, index=False)
        
        if self.verbose:
            print(f"Exported {len(branch_data)} tree points to {filepath}")

