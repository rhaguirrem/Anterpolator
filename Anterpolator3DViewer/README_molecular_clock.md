# Molecular Clock Interpolator

## Overview

The **Molecular Clock Interpolator** is a phylogeographic approach for spatial interpolation, particularly suited for modeling **intrusion systems** such as:
- Magmatic intrusions and feeder zones
- Hydrothermal vein systems
- Dike and sill networks
- Any geological feature with point-source dispersion

This algorithm treats spatial dispersion like evolutionary divergence, inferring hierarchical "common ancestor" points (feeders) and building branching tree structures connecting samples via parsimony minimization.

## Key Concepts

### Phylogenetic Tree Inference
- Builds hierarchical trees bottom-up using agglomerative clustering
- Each sample has a "genome" defined by [grade, x, y, z] coordinates
- Finds most parsimonious common ancestors that minimize total spatial + geochemical "mutation" distance
- Creates both leaf nodes (samples) and internal nodes (ancestors) at branch points

### Parsimony-Based Distance
- Composite distance combines spatial distance + attribute (grade) difference
- Depth bias: prefers deeper ancestors when multiple equally-parsimonious solutions exist
- Distance metric:
  ```
  distance = √(spatial_weight × spatial² + attr_weight × grade_diff²) - depth_bonus
  ```

### Multi-Event Detection
- Automatically identifies multiple intrusion events using DBSCAN clustering
- Each event has its own phylogenetic tree with root ancestor
- Useful for distinguishing different magmatic pulses or phases

### Depth-Biased Ancestors
- Ancestor nodes are inferred at merge points during bottom-up clustering
- Ancestors positioned deeper than both child nodes (configurable offset)
- Root ancestor is the deepest common source for all samples in an event
- Ancestor coordinates and grades computed via weighted centroid with depth bias

### Branch Interpolation
- Creates interpolated blocks along parent→child edges in the hierarchical tree
- Values interpolated linearly from ancestor to descendant
- Models realistic conduit geometries from deep feeders to surface samples

## Algorithm Steps

1. **Event Detection**: DBSCAN clusters samples into separate intrusion events
2. **Hierarchical Clustering**: For each event:
   - Initialize: each sample is its own cluster
   - Iteratively merge the two closest clusters:
     - Compute optimal ancestor position/grade at merge point
     - Create internal ancestor node with depth bias
     - Record parent→child edges
   - Continue until single root remains
3. **Tree Traversal**: BFS from root ancestor, interpolating blocks along each edge
4. **Export**: Output all nodes (samples + ancestors) and edges for visualization

## Configuration

### Algorithm Selection

In `config.json`, set the algorithm type:

```json
{
  "algorithm": "biochemical_clock",
  ...
}
```

### Parameters

Add biochemical clock parameters to your config:

```json
{
  "biochemical_clock_params": {
    "spatial_weight": 1.0,
    "attr_weight": 1.0,
    "detect_multiple_events": true,
    "branch_threshold": 2.0,
    "depth_weight": 1.0,
    "min_samples_per_event": 3,
    "interpolation_method": "linear",
    "fill_background": false,
    "background_value": 0.0,
    "ancestor_depth_offset": 1.0
  }
}
```

#### Parameter Descriptions

- **spatial_weight**: Weight for spatial distance in composite distance calculation (default: 1.0)
- **attr_weight**: Weight for attribute (geochemical) distance (default: 1.0)
- **detect_multiple_events**: If true, attempts to identify multiple intrusion events via DBSCAN clustering (default: true). Set to false to treat all samples as a single event.
- **branch_threshold**: DBSCAN clustering sensitivity multiplier (default: 2.0)
  - **Lower values (0.5-1.5)**: More sensitive, creates more separate events (good for widely-spaced intrusions)
  - **Higher values (2.0-5.0)**: Less sensitive, groups more samples together (good for closely-spaced samples)
  - If you see one giant tree connecting distant samples, **decrease this value**
  - If samples are being split into too many tiny events, **increase this value**
- **depth_weight**: Weight favoring deeper feeders; positive values prefer deeper Z coordinates (default: 1.0)
- **min_samples_per_event**: Minimum samples required to define an intrusion event (default: 3)
- **interpolation_method**: Method for branch interpolation - "linear" or "inverse_distance" (default: "linear")
- **fill_background**: Whether to fill blocks not on branches with background value (default: false)
- **background_value**: Value for background blocks if fill_background is true (default: 0.0)
- **ancestor_depth_offset**: Depth offset (in grid units) to push ancestor nodes below the deepest child at each merge. Higher values create deeper feeder systems (default: 1.0)
- **min_branch_angle**: Minimum dip angle in degrees for branches (default: 30.0)
  - Prevents near-horizontal connections that don't make geological sense
  - Angle measured from horizontal plane: 0° = horizontal, 90° = vertical
  - If a potential merge would create branches with angles below this threshold, the merge is rejected
  - **Typical values**: 20-45° for typical intrusions, 0° to disable constraint
  - **Effect**: Higher values force steeper, more vertical tree structures

**Note on Coordinate Convention**: Z coordinate is always treated as elevation. Larger Z = shallower, smaller Z = deeper. This is hardcoded and not configurable.

### Domain-Specific Algorithms

You can use different algorithms for different domains:

#### Using the UI

1. In the main configuration dialog, click **"Configure Domain Algorithms..."**
2. A dialog will open showing all domains from your Blocks file
3. For each domain, select:
   - **(use default)** - Use the main algorithm selected above
   - **ant_colony** - Force ant colony for this domain
   - **biochemical_clock** - Force biochemical clock for this domain
   - **skip** - Exclude this domain from interpolation entirely
4. You can use **"Apply Algorithm to All"** to quickly set all domains
5. Click **OK** to apply the configuration

The button will show how many domains are configured (e.g., "Configure Domain Algorithms... (3 configured)")

#### Using config.json (Manual)

Alternatively, edit the config.json file directly:

```json
{
  "algorithm": "ant_colony",
  "domain_algorithm_overrides": {
    "Intrusion": {
      "algorithm": "biochemical_clock",
      "spatial_weight": 1.0,
      "attr_weight": 0.5,
      "detect_multiple_events": true,
      "branch_threshold": 2.5,
      "depth_weight": 2.0
    },
    "Waste_Rock": {
      "skip": true
    },
    "Sedimentary": {
      "algorithm": "ant_colony"
    }
  }
}
```

#### Skipping Domains

Domains marked as **"skip"** will be:
- Excluded from the interpolation grid
- Not assigned to any blocks
- Not processed by any algorithm
- Not included in the output CSV

This is useful for:
- Excluding waste rock or barren zones
- Focusing interpolation on ore zones only
- Reducing computation time
- Avoiding interpolation in geologically inappropriate areas

## Output

The molecular clock exports additional fields in the CSV:

- **Event_ID**: Intrusion event identifier (-1 for unassigned)
- **Distance_To_Feeder**: Distance from block to root ancestor
- **Branch_ID**: Branch identifier within the tree structure
- **Is_Feeder**: Boolean indicating if block is the root ancestor
- **Node_ID**: Unique identifier for each node in the tree (samples and ancestors)
- **Point_Type**: "Sample", "Intermediate_Ancestor", "Root_Ancestor", or "Edge" in tree export

### Tree Export CSV

The `export_tree_structure()` method creates a detailed CSV showing the full hierarchical tree:

- **Nodes**: All samples and inferred ancestors with their positions and grades
- **Edges**: Parent→child connections showing the phylogenetic structure
- **Depth**: Z-coordinate (elevation) for each node
- **From_Node_ID**: Parent node for edge entries, enabling tree reconstruction

This CSV can be used to:
1. Visualize the inferred phylogenetic tree structure
2. Verify ancestor positions and grades
3. Analyze branching patterns and dispersion paths
4. Compare multiple intrusion events

## Visualization in Leapfrog

The exported CSV can be imported into Leapfrog Geo for visualization:

1. Import the interpolation CSV as a points file
2. Use **Event_ID** to color-code different intrusion events
3. Use **Distance_To_Feeder** to show dispersion gradients from root
4. Use **Branch_ID** to visualize individual conduits
5. Filter **Is_Feeder** = True to show root ancestor locations
6. Import tree export CSV and use **From_X/Y/Z** to draw parent→child connections
7. Filter by **Point_Type** to highlight ancestors vs samples

## Algorithm Comparison

| Feature | Ant Colony | Molecular Clock |
|---------|-----------|-------------------|
| **Best for** | Continuous interpolation, complex geometries | Intrusions, feeders, vein systems |
| **Iteration** | Multiple iterations needed | Single deterministic run |
| **Output** | Smooth gradients | Branching structures |
| **Feeders** | Not explicitly identified | Automatically located |
| **Multiple events** | Not detected | Can identify multiple sources |
| **Computation** | Slower (iterative) | Faster (deterministic) |

## Tips for Intrusion Modeling

1. **Increase depth_weight** (e.g., 2.0-5.0) to strongly prefer deeper feeders
2. **Lower branch_threshold** (e.g., 1.5) to detect more intrusion events
3. **Adjust spatial_weight vs attr_weight** based on whether spatial proximity or geochemical similarity is more important
4. Use **detect_multiple_events=true** when you expect multiple intrusive phases
5. Tune **ancestor_depth_offset** if ancestors should be significantly deeper (or shallower) than observed samples
6. Ensure your coordinate system uses Z as elevation (larger Z = shallower). If it does not, transform the input data before running the molecular clock.
7. Set **min_samples_per_event** to avoid spurious single-sample events

## Dependencies

The molecular clock requires additional Python packages:

```bash
pip install scipy scikit-learn
```

These are included in the updated `requirements.txt`.

## References

The algorithm is inspired by:
- Phylogeographic parsimony methods
- Minimum spanning trees in spatial analysis  
- Molecular clock concepts from molecular biology
- Geometric median optimization

---

**Note**: This is an experimental algorithm. Results should be validated against geological knowledge and other interpolation methods.
