from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Tuple

import numpy as np

from interpolator_base import InterpolatorBase


GridIndex = Tuple[int, int, int]
NodeKey = Tuple[int, int, int, int]


@dataclass
class _LeafAssignment:
    leaf_id: int
    level: int
    size: Tuple[int, int, int]
    value: float
    is_inherited: bool


class OctreeDomainInterpolator(InterpolatorBase):
    """Bottom-up domain-aware octree aggregation over the configured finest grid."""

    def __init__(
        self,
        output_mode: str = "dense_blocks_cover",
        max_levels: int = 0,
        support_density_alpha: float = 0.0,
        include_dense_provenance: bool = True,
        verbose: bool = False,
    ):
        super().__init__(verbose)
        normalized_mode = str(output_mode or "dense_blocks_cover").strip().lower()
        if normalized_mode not in {"dense_blocks_cover", "adaptive_leaf_cover"}:
            raise ValueError(
                "output_mode must be 'dense_blocks_cover' or 'adaptive_leaf_cover'"
            )

        self.output_mode = normalized_mode
        self.max_levels = max(int(max_levels or 0), 0)
        alpha_value = float(support_density_alpha or 0.0)
        if not math.isfinite(alpha_value) or alpha_value < 0.0 or alpha_value > 1.0:
            raise ValueError("support_density_alpha must be a finite value between 0 and 1")
        self.support_density_alpha = alpha_value
        self.include_dense_provenance = bool(include_dense_provenance)

        self.min_bounds = None
        self.block_size = None
        self.sample_blocks: Dict[GridIndex, float] = {}
        self.sample_domains: Dict[GridIndex, str] = {}
        self.use_domain_mapping = False
        self.allowed_positions: set[GridIndex] = set()
        self.allowed_positions_by_domain: Dict[str, set[GridIndex]] = {}
        self.domain_roots: Dict[str, list[NodeKey]] = {}
        self.valid_counts_by_domain_level: Dict[str, Dict[int, Dict[GridIndex, int]]] = {}
        self.sample_stats_by_domain_level: Dict[str, Dict[int, Dict[GridIndex, Dict[str, float]]]] = {}
        self.output_blocks_by_domain: Dict[str, Dict[GridIndex, Dict[str, Any]]] = {}
        self.output_nodes_by_domain: Dict[str, Dict[NodeKey, Dict[str, Any]]] = {}
        self._ran = False
        self._metadata: Dict[str, Any] = {}

    def get_algorithm_name(self) -> str:
        return "Adaptive Octree"

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
        self.output_blocks_by_domain = {}
        self.output_nodes_by_domain = {}
        self._metadata = {}
        self._ran = False
        sample_block_weight_sums = {
            tuple(map(int, pos)): float(weight)
            for pos, weight in (kwargs.get("sample_block_weight_sums") or {}).items()
            if float(weight) > 0.0
        }
        sample_block_counts = {
            tuple(map(int, pos)): int(count)
            for pos, count in (kwargs.get("sample_block_counts") or {}).items()
            if int(count) > 0
        }

        allowed_positions = getattr(self, "allowed_grid_override", None)
        if allowed_positions is None:
            self.allowed_positions = {
                (x, y, z)
                for x in range(self.dims[0])
                for y in range(self.dims[1])
                for z in range(self.dims[2])
            }
        else:
            self.allowed_positions = {tuple(map(int, pos)) for pos in allowed_positions}

        domain_mapping = getattr(self, "domain_mapping", None) or {}
        self.allowed_positions_by_domain = defaultdict(set)
        if self.use_domain_mapping:
            for pos in self.allowed_positions:
                domain_name = str(domain_mapping.get(pos, "")).strip()
                if domain_name:
                    self.allowed_positions_by_domain[domain_name].add(pos)
        else:
            default_domain = "Undomained"
            self.allowed_positions_by_domain[default_domain] = set(self.allowed_positions)

        for pos, value in self.sample_blocks.items():
            domain_name = self._resolve_sample_domain(pos)
            if not domain_name:
                continue
            if pos not in self.allowed_positions_by_domain.get(domain_name, set()):
                continue
            support_weight = float(sample_block_weight_sums.get(pos, 1.0))
            sample_count = int(sample_block_counts.get(pos, 1))
            self.blocks[pos] = {
                "value": float(value),
                "is_sample": True,
                "domain": domain_name,
                "level": 0,
                "relative_size": (1, 1, 1),
                "sample_count": sample_count,
                "support_weight": support_weight,
                "leaf_id": 0,
                "is_inherited": False,
            }

    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        if self._ran:
            return False
        self._ran = True

        self.output_blocks_by_domain = {}
        self.output_nodes_by_domain = {}
        all_output_blocks: Dict[GridIndex, Dict[str, Any]] = {}
        total_valid = 0
        total_output = 0
        total_leaf_nodes = 0
        max_depth_reached = 0

        for domain_name, allowed_positions in sorted(self.allowed_positions_by_domain.items()):
            if not allowed_positions:
                continue

            valid_levels = self._build_valid_levels(allowed_positions)
            sample_levels = self._build_sample_levels(domain_name, valid_levels)
            self.valid_counts_by_domain_level = self.valid_counts_by_domain_level or {}
            self.sample_stats_by_domain_level = self.sample_stats_by_domain_level or {}
            self.valid_counts_by_domain_level[domain_name] = valid_levels
            self.sample_stats_by_domain_level[domain_name] = sample_levels

            domain_roots, level_count = self._build_domain_roots(valid_levels)
            self.domain_roots[domain_name] = domain_roots

            domain_output_nodes: Dict[NodeKey, Dict[str, Any]] = {}
            leaf_assignments: Dict[GridIndex, _LeafAssignment] = {}
            next_leaf_id = 1

            for root_key in domain_roots:
                next_leaf_id = self._emit_cover_for_node(
                    domain_name,
                    root_key,
                    inherited_value=None,
                    leaf_assignments=leaf_assignments,
                    output_nodes=domain_output_nodes,
                    next_leaf_id=next_leaf_id,
                )

            if self.output_mode == "dense_blocks_cover":
                dense_blocks = self._build_dense_blocks(domain_name, leaf_assignments)
                self.output_blocks_by_domain[domain_name] = dense_blocks
                all_output_blocks.update(dense_blocks)
                total_output += len(dense_blocks)
            else:
                adaptive_blocks = self._build_adaptive_blocks(domain_output_nodes)
                self.output_blocks_by_domain[domain_name] = adaptive_blocks
                all_output_blocks.update(adaptive_blocks)
                total_output += len(adaptive_blocks)

            self.output_nodes_by_domain[domain_name] = domain_output_nodes
            total_valid += len(allowed_positions)
            total_leaf_nodes += len(domain_output_nodes)
            max_depth_reached = max(max_depth_reached, level_count)

        self.blocks = all_output_blocks
        self._metadata = {
            "algorithm": self.get_algorithm_name(),
            "output_mode": self.output_mode,
            "support_density_alpha": self.support_density_alpha,
            "include_dense_provenance": self.include_dense_provenance,
            "domains": len(self.output_blocks_by_domain),
            "total_blocks": len(self.blocks),
            "sample_blocks": len(self.sample_blocks),
            "interpolated_blocks": max(len(self.blocks) - len(self.sample_blocks), 0),
            "valid_finest_blocks": total_valid,
            "leaf_nodes": total_leaf_nodes,
            "max_levels_used": max_depth_reached,
        }
        return False

    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        return {pos: float(block["value"]) for pos, block in self.blocks.items()}

    def is_converged(self) -> bool:
        return self._ran

    def get_metadata(self) -> Dict[str, Any]:
        return dict(self._metadata)

    def _resolve_sample_domain(self, pos: GridIndex) -> str:
        if self.sample_domains.get(pos):
            return self.sample_domains[pos]
        if not self.use_domain_mapping:
            return "Undomained"
        domain_mapping = getattr(self, "domain_mapping", None) or {}
        return str(domain_mapping.get(pos, "")).strip()

    def _build_valid_levels(self, allowed_positions: Iterable[GridIndex]) -> Dict[int, Dict[GridIndex, int]]:
        levels: Dict[int, Dict[GridIndex, int]] = {0: {pos: 1 for pos in allowed_positions}}
        current = levels[0]
        level = 0
        max_levels = self.max_levels if self.max_levels > 0 else None
        while current and len(current) > 1 and (max_levels is None or level < max_levels):
            parent_counts: Dict[GridIndex, int] = defaultdict(int)
            for pos, count in current.items():
                parent_counts[self._parent_pos(pos)] += int(count)
            level += 1
            levels[level] = dict(parent_counts)
            current = levels[level]
        return levels

    def _build_sample_levels(
        self,
        domain_name: str,
        valid_levels: Dict[int, Dict[GridIndex, int]],
    ) -> Dict[int, Dict[GridIndex, Dict[str, float]]]:
        level_zero: Dict[GridIndex, Dict[str, float]] = {}
        for pos, value in self.sample_blocks.items():
            if self._resolve_sample_domain(pos) != domain_name:
                continue
            base_block = self.blocks.get(pos, {})
            raw_support = float(base_block.get("support_weight", 1.0))
            sample_count = float(base_block.get("sample_count", 1.0))
            effective_weight = self._compute_effective_support(
                raw_support,
                int(valid_levels.get(0, {}).get(pos, 1)),
            )
            level_zero[pos] = {
                "value": float(value),
                "raw_support": raw_support,
                "effective_weight": effective_weight,
                "sample_count": sample_count,
            }

        levels: Dict[int, Dict[GridIndex, Dict[str, float]]] = {0: level_zero}
        current = level_zero
        level = 0
        max_levels = self.max_levels if self.max_levels > 0 else None
        while current and len(current) > 1 and (max_levels is None or level < max_levels):
            parent_stats: Dict[GridIndex, Dict[str, float]] = {}
            for pos, stats in current.items():
                parent_pos = self._parent_pos(pos)
                parent = parent_stats.setdefault(
                    parent_pos,
                    {
                        "value_weight_sum": 0.0,
                        "effective_weight_sum": 0.0,
                        "raw_support": 0.0,
                        "sample_count": 0.0,
                    },
                )
                child_effective_weight = float(stats.get("effective_weight", 0.0))
                parent["value_weight_sum"] += float(stats.get("value", 0.0)) * child_effective_weight
                parent["effective_weight_sum"] += child_effective_weight
                parent["raw_support"] += float(stats.get("raw_support", 0.0))
                parent["sample_count"] += float(stats.get("sample_count", 0.0))
            level += 1
            finalized_parent_stats: Dict[GridIndex, Dict[str, float]] = {}
            for parent_pos, stats in parent_stats.items():
                effective_weight_sum = float(stats.get("effective_weight_sum", 0.0))
                if effective_weight_sum <= 0.0:
                    continue
                raw_support = float(stats.get("raw_support", 0.0))
                finalized_parent_stats[parent_pos] = {
                    "value": float(stats["value_weight_sum"] / effective_weight_sum),
                    "raw_support": raw_support,
                    "effective_weight": self._compute_effective_support(
                        raw_support,
                        int(valid_levels.get(level, {}).get(parent_pos, 0)),
                    ),
                    "sample_count": float(stats.get("sample_count", 0.0)),
                }
            levels[level] = finalized_parent_stats
            current = finalized_parent_stats
        return levels

    def _compute_effective_support(self, raw_support: float, represented_volume: int) -> float:
        raw_support_value = float(raw_support)
        if raw_support_value <= 0.0:
            return 0.0
        volume = max(int(represented_volume), 1)
        return raw_support_value / (float(volume) ** self.support_density_alpha)

    def _build_domain_roots(self, valid_levels: Dict[int, Dict[GridIndex, int]]) -> Tuple[list[NodeKey], int]:
        max_level = max(valid_levels.keys()) if valid_levels else 0
        roots = [
            (max_level, int(pos[0]), int(pos[1]), int(pos[2]))
            for pos in sorted(valid_levels.get(max_level, {}).keys())
        ]
        return roots, max_level

    def _emit_cover_for_node(
        self,
        domain_name: str,
        node_key: NodeKey,
        inherited_value: float | None,
        leaf_assignments: Dict[GridIndex, _LeafAssignment],
        output_nodes: Dict[NodeKey, Dict[str, Any]],
        next_leaf_id: int,
    ) -> int:
        level, ix, iy, iz = node_key
        pos = (ix, iy, iz)
        node_stats = self.sample_stats_by_domain_level.get(domain_name, {}).get(level, {}).get(pos)
        node_valid_count = int(self.valid_counts_by_domain_level.get(domain_name, {}).get(level, {}).get(pos, 0))
        if node_valid_count <= 0:
            return next_leaf_id

        node_value = inherited_value
        support_weight = 0.0
        sample_count = 0
        has_samples = False
        if node_stats and float(node_stats.get("effective_weight", 0.0)) > 0.0:
            support_weight = float(node_stats["effective_weight"])
            node_value = float(node_stats["value"])
            sample_count = int(round(float(node_stats.get("sample_count", 0.0))))
            has_samples = True

        children = []
        if level > 0:
            for child_pos in self._child_positions(pos):
                child_valid_count = int(
                    self.valid_counts_by_domain_level.get(domain_name, {}).get(level - 1, {}).get(child_pos, 0)
                )
                if child_valid_count > 0:
                    children.append((level - 1, int(child_pos[0]), int(child_pos[1]), int(child_pos[2])))

        if not children:
            if node_value is None:
                return next_leaf_id
            return self._record_leaf_node(
                domain_name,
                node_key,
                node_value,
                bool(not has_samples and inherited_value is not None),
                sample_count,
                support_weight,
                node_valid_count,
                leaf_assignments,
                output_nodes,
                next_leaf_id,
            )

        children_with_samples = [
            child_key
            for child_key in children
            if float(
                self.sample_stats_by_domain_level.get(domain_name, {})
                .get(child_key[0], {})
                .get((child_key[1], child_key[2], child_key[3]), {})
                .get("effective_weight", 0.0)
            ) > 0.0
        ]

        if not children_with_samples:
            if node_value is None:
                return next_leaf_id
            return self._record_leaf_node(
                domain_name,
                node_key,
                node_value,
                bool(inherited_value is not None and not has_samples),
                sample_count,
                support_weight,
                node_valid_count,
                leaf_assignments,
                output_nodes,
                next_leaf_id,
            )

        if node_value is None and has_samples is False:
            sibling_values = []
            for child_key in children_with_samples:
                child_stats = self.sample_stats_by_domain_level.get(domain_name, {}).get(child_key[0], {}).get(
                    (child_key[1], child_key[2], child_key[3]),
                    {},
                )
                if float(child_stats.get("effective_weight", 0.0)) > 0.0:
                    sibling_values.append(float(child_stats["value"]))
            if sibling_values:
                node_value = float(np.mean(sibling_values))

        for child_key in children:
            child_pos = (child_key[1], child_key[2], child_key[3])
            child_has_samples = float(
                self.sample_stats_by_domain_level.get(domain_name, {})
                .get(child_key[0], {})
                .get(child_pos, {})
                .get("effective_weight", 0.0)
            ) > 0.0
            if child_has_samples:
                next_leaf_id = self._emit_cover_for_node(
                    domain_name,
                    child_key,
                    inherited_value=node_value,
                    leaf_assignments=leaf_assignments,
                    output_nodes=output_nodes,
                    next_leaf_id=next_leaf_id,
                )
            else:
                if node_value is None:
                    continue
                next_leaf_id = self._record_leaf_node(
                    domain_name,
                    child_key,
                    node_value,
                    True,
                    0,
                    0.0,
                    int(
                        self.valid_counts_by_domain_level.get(domain_name, {}).get(child_key[0], {}).get(child_pos, 0)
                    ),
                    leaf_assignments,
                    output_nodes,
                    next_leaf_id,
                )

        return next_leaf_id

    def _record_leaf_node(
        self,
        domain_name: str,
        node_key: NodeKey,
        value: float,
        is_inherited: bool,
        sample_count: int,
        support_weight: float,
        valid_count: int,
        leaf_assignments: Dict[GridIndex, _LeafAssignment],
        output_nodes: Dict[NodeKey, Dict[str, Any]],
        next_leaf_id: int,
    ) -> int:
        level, ix, iy, iz = node_key
        relative_size = 2 ** int(level)
        node_record = {
            "value": float(value),
            "is_sample": level == 0 and sample_count > 0,
            "domain": domain_name,
            "level": int(level),
            "relative_size": (relative_size, relative_size, relative_size),
            "sample_count": int(sample_count),
            "support_weight": float(support_weight),
            "leaf_id": int(next_leaf_id),
            "is_inherited": bool(is_inherited),
            "valid_descendant_count": int(valid_count),
        }
        output_nodes[node_key] = node_record

        fine_cells = self._iter_descendant_finest_cells((ix, iy, iz), level)
        for cell_pos in fine_cells:
            if cell_pos not in self.allowed_positions_by_domain.get(domain_name, set()):
                continue
            leaf_assignments[cell_pos] = _LeafAssignment(
                leaf_id=int(next_leaf_id),
                level=int(level),
                size=(relative_size, relative_size, relative_size),
                value=float(value),
                is_inherited=bool(is_inherited),
            )
        return next_leaf_id + 1

    def _build_dense_blocks(
        self,
        domain_name: str,
        leaf_assignments: Dict[GridIndex, _LeafAssignment],
    ) -> Dict[GridIndex, Dict[str, Any]]:
        dense_blocks: Dict[GridIndex, Dict[str, Any]] = {}
        for pos in sorted(self.allowed_positions_by_domain.get(domain_name, set())):
            assignment = leaf_assignments.get(pos)
            if assignment is None:
                continue
            block_record: Dict[str, Any] = {
                "value": float(assignment.value),
                "is_sample": pos in self.sample_blocks,
                "domain": domain_name,
                "level": 0,
                "relative_size": (1, 1, 1),
                "sample_count": 1 if pos in self.sample_blocks else 0,
                "support_weight": 1.0 if pos in self.sample_blocks else 0.0,
                "is_inherited": bool(assignment.is_inherited),
            }
            if self.include_dense_provenance:
                block_record.update({
                    "adaptive_leaf_id": int(assignment.leaf_id),
                    "adaptive_level": int(assignment.level),
                    "adaptive_relative_size": tuple(int(v) for v in assignment.size),
                })
            dense_blocks[pos] = block_record
        return dense_blocks

    def _build_adaptive_blocks(self, output_nodes: Dict[NodeKey, Dict[str, Any]]) -> Dict[GridIndex, Dict[str, Any]]:
        adaptive_blocks: Dict[GridIndex, Dict[str, Any]] = {}
        for node_key, block_record in sorted(output_nodes.items()):
            _, ix, iy, iz = node_key
            relative_size = tuple(int(v) for v in block_record.get("relative_size", (1, 1, 1)))
            fine_origin = (
                int(ix) * relative_size[0],
                int(iy) * relative_size[1],
                int(iz) * relative_size[2],
            )
            adaptive_blocks[fine_origin] = {
                **block_record,
                "is_sample": bool(block_record.get("level", 0) == 0 and block_record.get("sample_count", 0) > 0),
            }
        return adaptive_blocks

    @staticmethod
    def _parent_pos(pos: GridIndex) -> GridIndex:
        return (int(pos[0]) // 2, int(pos[1]) // 2, int(pos[2]) // 2)

    @staticmethod
    def _child_positions(pos: GridIndex) -> list[GridIndex]:
        base_x = int(pos[0]) * 2
        base_y = int(pos[1]) * 2
        base_z = int(pos[2]) * 2
        return [
            (base_x + dx, base_y + dy, base_z + dz)
            for dx in (0, 1)
            for dy in (0, 1)
            for dz in (0, 1)
        ]

    @staticmethod
    def _iter_descendant_finest_cells(pos: GridIndex, level: int) -> Iterable[GridIndex]:
        scale = 2 ** int(level)
        base_x = int(pos[0]) * scale
        base_y = int(pos[1]) * scale
        base_z = int(pos[2]) * scale
        for dx in range(scale):
            for dy in range(scale):
                for dz in range(scale):
                    yield (base_x + dx, base_y + dy, base_z + dz)
