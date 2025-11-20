# Hierarchical Phylogenetic Tree Implementation

## Summary

Replaced the Minimum Spanning Tree (MST) approach with a true hierarchical phylogenetic tree construction algorithm for the BiochemicalClockInterpolator. The new implementation correctly models geological intrusion systems as evolutionary trees, inferring common ancestors at branch points using parsimony minimization.

## Changes Made

### 1. Core Algorithm Replacement

**Old Approach (MST):**
- Connected samples directly to each other using minimum spanning tree
- Created only one "feeder" node per event
- Violated phylogenetic principle: samples shouldn't connect to samples

**New Approach (Hierarchical Clustering):**
- Builds tree bottom-up using agglomerative clustering
- Infers optimal ancestor positions/grades at each merge point
- Creates internal ancestor nodes at all branch points
- Correctly models ancestor→descendant relationships

### 2. New Methods

#### `build_hierarchical_tree(sample_indices)`
- Implements UPGMA-style agglomerative clustering
- Iteratively merges closest cluster pairs based on parsimony distance
- Computes optimal ancestor position/grade at each merge:
  - Weighted centroid of child positions/grades
  - Biased toward deeper (lower Z) positions
  - Respects grid bounds and world coordinates
- Returns full tree structure with:
  - `nodes`: List of all node indices (samples + ancestors)
  - `edges`: Parent→child connections
  - `root_idx`: Root ancestor node
  - `node_data`: Full metadata for each node

#### `_compute_parsimony_distance(cluster_i, cluster_j)`
- Calculates composite distance between clusters
- Combines spatial distance + attribute (grade) difference
- Applies depth bias to favor merges creating deeper ancestors
- Formula: `√(spatial_weight × spatial² + attr_weight × grade²) - depth_bonus`

#### `_compute_optimal_ancestor(cluster_i, cluster_j)`
- Computes ancestor position/grade for merged cluster
- Uses weighted centroid with depth bias
- Pushes ancestor deeper than both children by `ancestor_depth_offset`
- Handles grid/world coordinate conversion correctly

#### `interpolate_hierarchical_tree(tree_info, sample_indices, event_id)`
- Replaces `interpolate_along_branches`
- Traverses tree using BFS from root ancestor
- Follows parent→child edges (not sample→sample)
- Creates interpolated blocks along each edge
- Assigns proper branch IDs at each fork point

### 3. Updated Methods

#### `run_iteration()`
- Now calls `build_hierarchical_tree()` instead of `build_mst_tree()`
- Stores full tree structure in `intrusion_trees`
- Calls `interpolate_hierarchical_tree()` for traversal

#### `export_tree_structure()`
- Exports all nodes (samples + all ancestors)
- Includes node types: Sample, Intermediate_Ancestor, Root_Ancestor, Edge
- Shows parent→child connections with From_Node_ID
- Enables full tree visualization in Leapfrog

### 4. Test Updates

#### `test_biochemical_clock.py`
- Updated `run_depth_bias_case()`:
  - Verifies hierarchical structure (n samples → n-1 ancestors)
  - Checks root is deeper than all samples
  - Validates correct number of edges
  
- Updated `run_branch_and_background_case()`:
  - Verifies 3 samples create 2 internal ancestors
  - Checks branch IDs are assigned correctly
  - Tests background fill still works

#### `test_hierarchical_visual.py` (new)
- Visual test with clustered samples
- Prints full tree structure with node positions/grades
- Verifies tree properties (node counts, depth bias)
- Demonstrates intermediate ancestors at merge points

## Algorithm Properties

### Correctness
✅ **Phylogenetic principle**: Samples connect to ancestors, not to each other  
✅ **Parsimony**: Minimizes total spatial + grade "mutation" distance  
✅ **Depth bias**: Prefers deeper ancestors when multiple solutions exist  
✅ **Multi-tree**: DBSCAN clustering creates separate events  
✅ **Sample preservation**: Never overwrites original sample values  

### Tree Structure
- Binary tree with n samples has n-1 internal ancestor nodes
- Root ancestor is deepest node in the tree
- Intermediate ancestors appear at merge points
- Each merge creates optimal ancestor via weighted centroid

### Performance
- Time complexity: O(n² log n) for n samples per event
- Space complexity: O(n) for tree storage
- Deterministic: same input always produces same tree

## Coordinate Convention

**Z is always elevation** (hardcoded, not parametric):
- Larger Z = shallower (closer to surface)
- Smaller Z = deeper (further from surface)
- Ancestors biased toward smaller Z values

## Export Format

Tree CSV contains:
- **Nodes**: All samples and ancestors with positions/grades
- **Edges**: Parent→child connections
- **Point_Type**: Sample / Intermediate_Ancestor / Root_Ancestor / Edge
- **Node_ID**: Unique identifier for each node
- **From_Node_ID**: Parent node (for edges)
- **Depth**: Z-coordinate (elevation)

## Testing Results

All tests pass:
```
✓ Depth bias test passed: root ancestor is deeper than samples
✓ Branch and background test passed: 3 samples, 2 ancestors, 2 branch IDs
Biochemical clock regression suite passed
```

Visual test with 5 samples:
```
Nodes: 9 (5 samples + 4 ancestors)
Edges: 8 (parent→child connections)
Root ancestor depth: 0.0 (deeper than deepest sample at 20.0)
Depth bias satisfied: True
```

## Documentation Updated

- `README_biochemical_clock.md`: 
  - Added "Phylogenetic Tree Inference" section
  - Documented hierarchical clustering algorithm
  - Updated algorithm steps
  - Added tree export format
  - Clarified coordinate convention

## Backward Compatibility

⚠️ **Breaking change**: 
- Old `intrusion_trees` format replaced with new structure
- `build_mst_tree()` no longer used (kept for reference, could be removed)
- `_ensure_ancestor_node()` no longer used (kept for potential future use)
- Tree export CSV format changed (new columns, different node types)

Projects using the old format will need to regenerate interpolations.

## Future Enhancements

Possible improvements:
1. Add cluster size weighting when computing ancestor centroids
2. Implement non-binary tree support (more than 2 children per node)
3. Add tree pruning for very small branches
4. Support for asymmetric depth bias (different weights per axis)
5. Interactive tree visualization in the GUI

## Files Modified

1. `biochemical_clock_interpolator.py`:
   - Added `build_hierarchical_tree()` (223 lines)
   - Added `_compute_parsimony_distance()` (43 lines)
   - Added `_compute_optimal_ancestor()` (49 lines)
   - Added `interpolate_hierarchical_tree()` (130 lines)
   - Updated `run_iteration()` (9 lines)
   - Updated `export_tree_structure()` (64 lines)

2. `README_biochemical_clock.md`:
   - Updated overview and key concepts (47 lines)
   - Added algorithm steps section (15 lines)
   - Updated output documentation (30 lines)

3. `tests/test_biochemical_clock.py`:
   - Updated `run_depth_bias_case()` (34 lines)
   - Updated `run_branch_and_background_case()` (35 lines)

4. `tests/test_hierarchical_visual.py` (new):
   - Created visual test (106 lines)

**Total changes**: ~750 lines modified/added across 4 files
