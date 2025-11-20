"""Lightweight regression checks for BiochemicalClockInterpolator."""
import os
import sys

import numpy as np

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from biochemical_clock_interpolator import BiochemicalClockInterpolator


def run_depth_bias_case():
    """Test that ancestors are placed deeper than samples."""
    sample_blocks = {
        (0, 0, 0): 0.2,
        (0, 0, 4): 1.0,
    }
    dims = (1, 1, 7)
    min_bounds = (0, 0, 0)
    block_size = (1, 1, 1)

    interpolator = BiochemicalClockInterpolator(
        depth_weight=5.0,
        fill_background=False,
        detect_multiple_events=False,
        verbose=False,
    )
    interpolator.initialize_blocks(sample_blocks, dims, min_bounds, block_size)
    interpolator.run_iteration(dims)

    # Check hierarchical tree structure
    assert len(interpolator.intrusion_trees) == 1, "Should have one event"
    tree = interpolator.intrusion_trees[0]
    tree_info = tree['tree_info']
    node_data = tree_info['node_data']
    
    # Should have 2 samples + 1 ancestor = 3 nodes total
    assert len(node_data) == 3, f"Expected 3 nodes (2 samples + 1 ancestor), got {len(node_data)}"
    
    # Find root ancestor
    root_idx = tree_info['root_idx']
    root_data = node_data[root_idx]
    assert not root_data['is_sample'], "Root should be an ancestor, not a sample"
    
    # Check that root is deeper (lower Z) than samples
    sample_world = [interpolator._grid_to_world(np.array(pos)) for pos in sample_blocks.keys()]
    deepest_sample = min(world[2] for world in sample_world)
    root_world = interpolator._grid_to_world(root_data['coords'])
    assert root_world[2] <= deepest_sample, f"Root ancestor Z={root_world[2]} should be <= deepest sample Z={deepest_sample}"
    
    # Check edges: should have 2 edges (root to each sample)
    edges = tree_info['edges']
    assert len(edges) == 2, f"Expected 2 edges, got {len(edges)}"
    
    print("  ✓ Depth bias test passed: root ancestor is deeper than samples")


def run_branch_and_background_case():
    """Test branch IDs and background fill with 3 samples."""
    sample_blocks = {
        (0, 0, 0): 0.2,
        (0, 0, 2): 0.5,
        (0, 0, 4): 1.0,
    }
    dims = (1, 1, 7)
    min_bounds = (0, 0, 0)
    block_size = (1, 1, 1)

    interpolator = BiochemicalClockInterpolator(
        depth_weight=5.0,
        fill_background=True,
        background_value=-99.0,
        detect_multiple_events=False,
        verbose=False,
    )
    interpolator.initialize_blocks(sample_blocks, dims, min_bounds, block_size)
    interpolator.run_iteration(dims)

    # Check hierarchical tree structure
    assert len(interpolator.intrusion_trees) == 1, "Should have one event"
    tree = interpolator.intrusion_trees[0]
    tree_info = tree['tree_info']
    node_data = tree_info['node_data']
    
    # Should have 3 samples + 2 intermediate ancestors + 1 root = 6 nodes total
    # (binary tree with 3 leaves has 2 internal nodes)
    n_samples = sum(1 for nd in node_data.values() if nd['is_sample'])
    n_ancestors = sum(1 for nd in node_data.values() if not nd['is_sample'])
    assert n_samples == 3, f"Expected 3 samples, got {n_samples}"
    assert n_ancestors == 2, f"Expected 2 ancestors for 3 samples, got {n_ancestors}"
    
    # Check branch IDs are assigned
    branch_ids = sorted(
        {
            data["branch_id"]
            for data in interpolator.blocks.values()
            if data.get("event_id") == 0 and not data.get("is_sample")
        }
    )
    assert branch_ids and branch_ids[0] == 0, "Branch IDs not populated as expected"
    
    # Check background fill
    expected_total = int(np.prod(dims))
    assert len(interpolator.blocks) == expected_total, (
        f"Background fill missing cells: {len(interpolator.blocks)} != {expected_total}"
    )
    
    print(f"  ✓ Branch and background test passed: {n_samples} samples, {n_ancestors} ancestors, {len(branch_ids)} branch IDs")
def run_regression_suite():
    run_depth_bias_case()
    run_branch_and_background_case()
    print("Biochemical clock regression suite passed")


if __name__ == "__main__":
    run_regression_suite()
