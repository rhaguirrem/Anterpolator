"""Regression checks for the bottom-up adaptive octree interpolator."""
import os
import sys

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from octree_domain_interpolator import OctreeDomainInterpolator


def run_dense_cover_case():
    interpolator = OctreeDomainInterpolator(output_mode="dense_blocks_cover")
    interpolator.allowed_grid_override = {(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)}
    interpolator.domain_mapping = {
        (0, 0, 0): "A",
        (1, 0, 0): "A",
        (2, 0, 0): "A",
        (3, 0, 0): "A",
    }
    interpolator.initialize_blocks(
        {
            (0, 0, 0): 10.0,
            (1, 0, 0): 14.0,
        },
        (4, 1, 1),
        (0, 0, 0),
        (1, 1, 1),
        use_domain_mapping=True,
        sample_domain_mapping={(0, 0, 0): "A", (1, 0, 0): "A"},
    )

    assert interpolator.run_iteration((4, 1, 1)) is False
    values = interpolator.get_interpolated_values()

    assert len(values) == 4, "Dense cover should assign every allowed finest cell"
    assert values[(0, 0, 0)] == 10.0
    assert values[(1, 0, 0)] == 14.0
    assert abs(values[(2, 0, 0)] - 12.0) < 1e-9, "Unsampled sibling should inherit parent aggregate"
    assert abs(values[(3, 0, 0)] - 12.0) < 1e-9, "Unsampled sibling should inherit parent aggregate"
    assert interpolator.blocks[(2, 0, 0)]["adaptive_level"] == 1


def run_adaptive_leaf_cover_case():
    interpolator = OctreeDomainInterpolator(output_mode="adaptive_leaf_cover")
    interpolator.allowed_grid_override = {(0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)}
    interpolator.domain_mapping = {
        (0, 0, 0): "A",
        (1, 0, 0): "A",
        (2, 0, 0): "A",
        (3, 0, 0): "A",
    }
    interpolator.initialize_blocks(
        {
            (0, 0, 0): 10.0,
            (1, 0, 0): 14.0,
        },
        (4, 1, 1),
        (0, 0, 0),
        (1, 1, 1),
        use_domain_mapping=True,
        sample_domain_mapping={(0, 0, 0): "A", (1, 0, 0): "A"},
    )

    interpolator.run_iteration((4, 1, 1))
    values = interpolator.get_interpolated_values()

    assert len(values) == 3, "Adaptive cover should keep sampled leaves and one inherited coarse sibling"
    assert (0, 0, 0) in values
    assert (1, 0, 0) in values
    assert (2, 0, 0) in values, "Coarse sibling should be emitted at its finest-grid origin"
    inherited_candidates = [
        block for block in interpolator.blocks.values() if block.get("is_inherited")
    ]
    assert inherited_candidates, "Adaptive cover should include inherited unsampled sibling blocks"
    inherited_block = inherited_candidates[0]
    assert inherited_block["relative_size"] == (2, 2, 2)
    assert abs(inherited_block["value"] - 12.0) < 1e-9


def run_domain_separation_case():
    interpolator = OctreeDomainInterpolator(output_mode="dense_blocks_cover")
    interpolator.allowed_grid_override = {(0, 0, 0), (1, 0, 0)}
    interpolator.domain_mapping = {
        (0, 0, 0): "A",
        (1, 0, 0): "B",
    }
    interpolator.initialize_blocks(
        {(0, 0, 0): 5.0, (1, 0, 0): 25.0},
        (2, 1, 1),
        (0, 0, 0),
        (1, 1, 1),
        use_domain_mapping=True,
        sample_domain_mapping={(0, 0, 0): "A", (1, 0, 0): "B"},
    )

    interpolator.run_iteration((2, 1, 1))
    values = interpolator.get_interpolated_values()

    assert values[(0, 0, 0)] == 5.0
    assert values[(1, 0, 0)] == 25.0


def run_support_density_alpha_case():
    allowed_positions = {
        (x, y, z)
        for x in (0, 1)
        for y in (0, 1)
        for z in (0, 1)
    }
    allowed_positions.update({(2, 0, 0), (0, 2, 0)})
    domain_mapping = {pos: "A" for pos in allowed_positions}
    sample_blocks = {
        (0, 0, 0): 0.0,
        (2, 0, 0): 10.0,
    }
    sample_domains = {pos: "A" for pos in sample_blocks}

    baseline = OctreeDomainInterpolator(output_mode="dense_blocks_cover", support_density_alpha=0.0)
    baseline.allowed_grid_override = set(allowed_positions)
    baseline.domain_mapping = dict(domain_mapping)
    baseline.initialize_blocks(
        sample_blocks,
        (4, 4, 4),
        (0, 0, 0),
        (1, 1, 1),
        use_domain_mapping=True,
        sample_domain_mapping=sample_domains,
    )
    baseline.run_iteration((4, 4, 4))

    density_weighted = OctreeDomainInterpolator(output_mode="dense_blocks_cover", support_density_alpha=1.0)
    density_weighted.allowed_grid_override = set(allowed_positions)
    density_weighted.domain_mapping = dict(domain_mapping)
    density_weighted.initialize_blocks(
        sample_blocks,
        (4, 4, 4),
        (0, 0, 0),
        (1, 1, 1),
        use_domain_mapping=True,
        sample_domain_mapping=sample_domains,
    )
    density_weighted.run_iteration((4, 4, 4))

    inherited_pos = (0, 2, 0)
    baseline_value = baseline.get_interpolated_values()[inherited_pos]
    density_weighted_value = density_weighted.get_interpolated_values()[inherited_pos]

    assert abs(baseline_value - 5.0) < 1e-9
    assert abs(density_weighted_value - (80.0 / 9.0)) < 1e-9
    assert density_weighted_value > baseline_value
    assert density_weighted.get_metadata()["support_density_alpha"] == 1.0


def run_regression_suite():
    run_dense_cover_case()
    run_adaptive_leaf_cover_case()
    run_domain_separation_case()
    run_support_density_alpha_case()
    print("Adaptive octree regression suite passed")


if __name__ == "__main__":
    run_regression_suite()