"""Regression checks for GaussianKernelInterpolator."""
import os
import sys

CURRENT_DIR = os.path.dirname(__file__)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from gaussian_kernel_interpolator import GaussianKernelInterpolator


def run_single_pass_case():
    sample_blocks = {
        (0, 0, 0): 0.0,
        (2, 0, 0): 10.0,
    }
    interpolator = GaussianKernelInterpolator(bandwidth=1.0, cutoff_sigma=3.0)
    interpolator.initialize_blocks(sample_blocks, (3, 1, 1), (0, 0, 0), (1, 1, 1))

    should_continue = interpolator.run_iteration((3, 1, 1))
    values = interpolator.get_interpolated_values()

    assert should_continue is False, "Gaussian kernel should finish in one pass"
    assert values[(0, 0, 0)] == 0.0, "Sample blocks must be preserved"
    assert values[(2, 0, 0)] == 10.0, "Sample blocks must be preserved"
    assert (1, 0, 0) in values, "Interpolator should fill the interior block"
    assert 4.0 < values[(1, 0, 0)] < 6.0, f"Unexpected smoothed midpoint: {values[(1, 0, 0)]}"


def run_domain_constraint_case():
    sample_blocks = {
        (0, 0, 0): 0.0,
        (2, 0, 0): 100.0,
    }
    interpolator = GaussianKernelInterpolator(bandwidth=1.0, cutoff_sigma=3.0, use_nearest_fallback=False)
    interpolator.allowed_grid_override = {(0, 0, 0), (1, 0, 0), (2, 0, 0)}
    interpolator.domain_mapping = {
        (0, 0, 0): 'A',
        (1, 0, 0): 'A',
        (2, 0, 0): 'B',
    }
    interpolator.initialize_blocks(
        sample_blocks,
        (3, 1, 1),
        (0, 0, 0),
        (1, 1, 1),
        use_domain_mapping=True,
        sample_domain_mapping={(0, 0, 0): 'A', (2, 0, 0): 'B'},
    )
    interpolator.run_iteration((3, 1, 1))
    values = interpolator.get_interpolated_values()

    assert abs(values[(1, 0, 0)] - 0.0) < 1e-9, "Domain filtering should prevent smoothing across domains"


def run_regression_suite():
    run_single_pass_case()
    run_domain_constraint_case()
    print("Gaussian kernel regression suite passed")


if __name__ == "__main__":
    run_regression_suite()