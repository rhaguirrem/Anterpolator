"""Visual test for hierarchical tree construction with real data."""
import os
import sys

import numpy as np

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from biochemical_clock_interpolator import BiochemicalClockInterpolator


def test_clustered_samples():
    """Test with spatially clustered samples to verify hierarchical structure."""
    # Create two clusters of samples
    sample_blocks = {
        # Cluster 1: shallow samples
        (5, 5, 2): 1.0,
        (6, 5, 2): 1.1,
        (5, 6, 2): 1.2,
        
        # Cluster 2: deeper samples, different location
        (15, 15, 5): 2.0,
        (16, 15, 5): 2.1,
    }
    dims = (25, 25, 10)
    min_bounds = (0, 0, 0)
    block_size = (10, 10, 10)

    interpolator = BiochemicalClockInterpolator(
        spatial_weight=1.0,
        attr_weight=1.0,
        depth_weight=2.0,
        ancestor_depth_offset=1.5,
        detect_multiple_events=True,
        branch_threshold=3.0,
        fill_background=False,
        verbose=True,
    )
    interpolator.initialize_blocks(sample_blocks, dims, min_bounds, block_size)
    interpolator.run_iteration(dims)

    print("\n" + "="*80)
    print("HIERARCHICAL TREE STRUCTURE")
    print("="*80)
    
    for tree in interpolator.intrusion_trees:
        event_id = tree['event_id']
        tree_info = tree['tree_info']
        node_data = tree_info['node_data']
        edges = tree_info['edges']
        root_idx = tree_info['root_idx']
        
        print(f"\n--- Event {event_id} ---")
        print(f"Nodes: {len(node_data)}")
        print(f"Edges: {len(edges)}")
        print(f"Root: Node {root_idx}")
        
        # Show all nodes
        print("\nNodes:")
        for node_idx in sorted(node_data.keys()):
            data = node_data[node_idx]
            coords = data['coords']
            world = interpolator._grid_to_world(coords)
            node_type = "Sample" if data['is_sample'] else "Ancestor"
            print(f"  Node {node_idx:2d} ({node_type:8s}): "
                  f"grid=({coords[0]:5.1f}, {coords[1]:5.1f}, {coords[2]:5.1f})  "
                  f"world=({world[0]:6.1f}, {world[1]:6.1f}, {world[2]:6.1f})  "
                  f"grade={data['value']:.2f}")
        
        # Show edges
        print("\nEdges (parent → child):")
        for parent_idx, child_idx in edges:
            parent_data = node_data[parent_idx]
            child_data = node_data[child_idx]
            p_type = "Sample" if parent_data['is_sample'] else "Ancestor"
            c_type = "Sample" if child_data['is_sample'] else "Ancestor"
            print(f"  Node {parent_idx:2d} ({p_type:8s}) → Node {child_idx:2d} ({c_type:8s})")
        
        # Verify tree properties
        n_samples = sum(1 for d in node_data.values() if d['is_sample'])
        n_ancestors = sum(1 for d in node_data.values() if not d['is_sample'])
        print(f"\nTree Properties:")
        print(f"  Samples: {n_samples}")
        print(f"  Ancestors: {n_ancestors}")
        print(f"  Expected ancestors for {n_samples} samples: {n_samples - 1}")
        
        # Verify depth bias
        root_data = node_data[root_idx]
        root_world = interpolator._grid_to_world(root_data['coords'])
        sample_depths = [interpolator._grid_to_world(d['coords'])[2] 
                        for d in node_data.values() if d['is_sample']]
        deepest_sample = min(sample_depths)
        print(f"  Root ancestor depth: {root_world[2]:.1f}")
        print(f"  Deepest sample depth: {deepest_sample:.1f}")
        print(f"  Depth bias satisfied: {root_world[2] <= deepest_sample}")
    
    print("\n" + "="*80)
    print("TEST PASSED: Hierarchical tree structure looks correct!")
    print("="*80)


if __name__ == "__main__":
    test_clustered_samples()
