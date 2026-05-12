import math
from typing import Any, Dict, Tuple

import numpy as np
from scipy.spatial import cKDTree

from interpolator_base import InterpolatorBase


class GaussianKernelInterpolator(InterpolatorBase):
    """Single-pass Gaussian kernel smoother for scalar block interpolation."""

    def __init__(
        self,
        bandwidth: float = 3.0,
        cutoff_sigma: float = 3.0,
        use_nearest_fallback: bool = True,
        background_value: float | None = None,
        fill_background: bool = False,
        verbose: bool = False,
    ):
        super().__init__(verbose)
        self.bandwidth = max(float(bandwidth), 1e-6)
        self.cutoff_sigma = max(float(cutoff_sigma), 0.5)
        self.use_nearest_fallback = bool(use_nearest_fallback)
        self.background_value = None if background_value is None else float(background_value)
        self.fill_background = bool(fill_background)

        self.min_bounds = None
        self.block_size = None
        self.sample_blocks: Dict[Tuple[int, int, int], float] = {}
        self.sample_positions = []
        self.sample_values = None
        self.sample_domains: Dict[Tuple[int, int, int], str] = {}
        self.sample_tree: cKDTree | None = None
        self._ran = False

    def get_algorithm_name(self) -> str:
        return "Gaussian Kernel"

    def initialize_blocks(
        self,
        sample_blocks: Dict[Tuple[int, int, int], float],
        dims: Tuple[int, int, int],
        min_bounds,
        block_size,
        use_domain_mapping: bool = False,
        sample_domain_mapping: Dict[Tuple[int, int, int], str] | None = None,
        **kwargs,
    ):
        self.dims = tuple(int(axis) for axis in dims)
        self.min_bounds = None if min_bounds is None else np.asarray(min_bounds, dtype=float)
        self.block_size = None if block_size is None else np.asarray(block_size, dtype=float)
        self.use_domain_mapping = bool(use_domain_mapping)
        self.sample_blocks = {tuple(map(int, pos)): float(value) for pos, value in sample_blocks.items()}
        self.sample_domains = {
            tuple(map(int, pos)): str(domain).strip()
            for pos, domain in (sample_domain_mapping or {}).items()
            if str(domain).strip()
        }
        self.blocks = {}

        self.sample_positions = sorted(self.sample_blocks.keys())
        if self.sample_positions:
            sample_coords = np.asarray(self.sample_positions, dtype=float)
            self.sample_values = np.asarray([self.sample_blocks[pos] for pos in self.sample_positions], dtype=float)
            self.sample_tree = cKDTree(sample_coords)
        else:
            self.sample_values = np.empty(0, dtype=float)
            self.sample_tree = None

        for pos, value in self.sample_blocks.items():
            block_data: Dict[str, Any] = {
                'value': float(value),
                'is_sample': True,
                'visit_count': 1,
            }
            if self.sample_domains.get(pos):
                block_data['domain'] = self.sample_domains[pos]
            self.blocks[pos] = block_data

        self._ran = False

    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        if self._ran:
            return False
        self._ran = True

        if not self.sample_positions:
            return False

        allowed_positions = getattr(self, 'allowed_grid_override', None)
        if allowed_positions is None:
            positions = (
                (x, y, z)
                for x in range(self.dims[0])
                for y in range(self.dims[1])
                for z in range(self.dims[2])
            )
        else:
            positions = sorted(tuple(map(int, pos)) for pos in allowed_positions)

        search_radius = self.bandwidth * self.cutoff_sigma
        domain_mapping = getattr(self, 'domain_mapping', None) or {}

        for pos in positions:
            if pos in self.blocks:
                continue

            target_domain = None
            if self.use_domain_mapping:
                target_domain = str(domain_mapping.get(pos, '')).strip() or None
                if target_domain is None and not self.fill_background:
                    continue

            value = self._interpolate_position(pos, search_radius, target_domain)
            if value is None:
                if not self.fill_background or self.background_value is None:
                    continue
                value = self.background_value

            block_data: Dict[str, Any] = {
                'value': float(value),
                'is_sample': False,
                'visit_count': 1,
            }
            if target_domain is not None:
                block_data['domain'] = target_domain
            self.blocks[pos] = block_data

        return False

    def _interpolate_position(self, pos: Tuple[int, int, int], search_radius: float, target_domain: str | None) -> float | None:
        if self.sample_tree is None:
            return None

        target = np.asarray(pos, dtype=float)
        candidate_indices = self.sample_tree.query_ball_point(target, r=search_radius)

        if target_domain is not None:
            candidate_indices = [
                index for index in candidate_indices
                if self.sample_domains.get(self.sample_positions[index]) == target_domain
            ]

        if not candidate_indices:
            if not self.use_nearest_fallback:
                return None
            _, nearest_index = self.sample_tree.query(target, k=1)
            nearest_index = int(nearest_index)
            if target_domain is not None and self.sample_domains.get(self.sample_positions[nearest_index]) != target_domain:
                matching_indices = [
                    index for index, sample_pos in enumerate(self.sample_positions)
                    if self.sample_domains.get(sample_pos) == target_domain
                ]
                if not matching_indices:
                    return None
                candidate_indices = matching_indices
            else:
                candidate_indices = [nearest_index]

        neighbor_coords = np.asarray([self.sample_positions[index] for index in candidate_indices], dtype=float)
        deltas = neighbor_coords - target
        distances_sq = np.einsum('ij,ij->i', deltas, deltas)
        weights = np.exp(-0.5 * distances_sq / (self.bandwidth ** 2))
        total_weight = float(np.sum(weights))

        if not math.isfinite(total_weight) or total_weight <= 0.0:
            return None

        values = self.sample_values[candidate_indices]
        return float(np.dot(weights, values) / total_weight)

    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        return {pos: float(data['value']) for pos, data in self.blocks.items()}

    def is_converged(self) -> bool:
        return self._ran

    def get_metadata(self) -> Dict[str, Any]:
        return {
            'algorithm': self.get_algorithm_name(),
            'bandwidth': self.bandwidth,
            'cutoff_sigma': self.cutoff_sigma,
            'fill_background': self.fill_background,
            'total_blocks': len(self.blocks),
            'sample_blocks': len(self.sample_blocks),
            'interpolated_blocks': max(len(self.blocks) - len(self.sample_blocks), 0),
        }