"""Rectified Flow implementation removed.
This file is intentionally disabled.
"""
"""
import numpy as np
from typing import Dict, Tuple, Any, List, Optional

from interpolator_base import InterpolatorBase
from molecular_clock_interpolator import MolecularClockInterpolator


class RectifiedFlowInterpolator(InterpolatorBase):
    """Rectified Flow (intrusion) interpolator.

    Two-stage algorithm:
    1) Feeder roots are inferred using the Molecular Clock algorithm (multi-event clustering + LUCA/root).
    2) Each sample is connected to its closest feeder root using rectified-flow-style straight paths.

    Notes:
    - This is a deterministic, single-pass interpolator (run_iteration executes once).
    - "Rectified" here is implemented as a straight path to a per-sample sub-voxel target inside a
      Gaussian blob around the feeder root, to diversify endpoints without inflating the feeder.
    """

    def __init__(
        self,
        molecular_clock_params: Optional[Dict[str, Any]] = None,
        sigma_blocks: float = 0.5,
        z_sigma_multiplier: float = 1.0,
        seed: int = 0,
        steps_per_block: float = 2.0,
        collision_policy: str = 'average',
        processing_order: str = 'ascending',
        path_mode: str = 'straight',
        curve_strength: float = 0.0,
        interpolation_target: str = 'value',
        verbose: bool = False,
    ):
        super().__init__(verbose)
        self.molecular_clock_params = molecular_clock_params or {}

        self.sigma_blocks = float(sigma_blocks)
        self.z_sigma_multiplier = float(z_sigma_multiplier)
        self.seed = int(seed)
        self.steps_per_block = float(steps_per_block)

        self.collision_policy = (collision_policy or 'average').strip().lower()
        self.processing_order = (processing_order or 'ascending').strip().lower()
        self.path_mode = (path_mode or 'straight').strip().lower()
        self.curve_strength = float(curve_strength)

        target = (interpolation_target or 'value').strip().lower()
        if target in ('numeric', 'number', 'grade', 'value'):
            target = 'value'
        elif target in ('domain', 'domains', 'category', 'categorical'):
            target = 'domain'
        else:
            target = 'value'
        self.interpolation_target = target

        self.converged = False

        # Populated during initialize_blocks
        self.dims = None
        self.min_bounds = None
        self.block_size = None
        self.sample_blocks: Dict[Tuple[int, int, int], float] = {}
        self.sample_domain_mapping: Dict[Tuple[int, int, int], str] = {}

        # Output domain mapping when interpolating categorical domains
        self.domain_mapping: Dict[Tuple[int, int, int], str] = {}

        # Feeder roots and event membership (computed in run_iteration)
        self.feeder_roots: List[Dict[str, Any]] = []
        self.sample_to_event: Dict[Tuple[int, int, int], int] = {}

    def get_algorithm_name(self) -> str:
        return 'Rectified Flow (Feeder Paths)'

    def initialize_blocks(
        self,
        sample_blocks: Dict[Tuple[int, int, int], float],
        dims: Tuple[int, int, int],
        min_bounds,
        block_size,
        **kwargs,
    ):
        self.dims = tuple(int(x) for x in dims)
        self.min_bounds = np.array(min_bounds, dtype=float) if min_bounds is not None else None
        self.block_size = np.array(block_size, dtype=float) if block_size is not None else None

        self.sample_blocks = dict(sample_blocks)
        self.blocks = {}

        self.use_domain_mapping = kwargs.get('use_domain_mapping', False)
        self.sample_domain_mapping = kwargs.get('sample_domain_mapping', None) or {}

        if self.interpolation_target == 'domain':
            # For domain interpolation, initialize samples with their domain.
            self.domain_mapping = {}
            for pos, val in self.sample_blocks.items():
                domain = self.sample_domain_mapping.get(pos)
                if domain is None or str(domain).strip() == '':
                    continue
                domain = str(domain).strip()
                self.domain_mapping[pos] = domain
                self.blocks[pos] = {
                    'value': float(val) if val is not None else 0.0,
                    'is_sample': True,
                    'visit_count': 1,
                    'domain': domain,
                }
            # Do not use domain constraint mapping when the target is the domain itself.
            self.use_domain_mapping = False
        else:
            for pos, val in self.sample_blocks.items():
                self.blocks[pos] = {
                    'value': float(val),
                    'is_sample': True,
                    'visit_count': 1,
                }

        self.converged = False

    def _grid_to_world(self, coords: np.ndarray) -> np.ndarray:
        if self.min_bounds is None or self.block_size is None:
            return np.array(coords, dtype=float)
        return self.min_bounds + np.array(coords, dtype=float) * self.block_size

    def _world_to_grid_float(self, world: np.ndarray) -> np.ndarray:
        if self.min_bounds is None or self.block_size is None:
            return np.array(world, dtype=float)
        return (np.array(world, dtype=float) - self.min_bounds) / self.block_size

    def _iter_samples_in_order(self, event_samples: List[Tuple[Tuple[int, int, int], float]]):
        if self.interpolation_target == 'domain':
            # Deterministic ordering for categorical mode.
            return sorted(event_samples, key=lambda x: x[0])

        if self.processing_order == 'random':
            rng = np.random.default_rng(self.seed)
            idx = np.arange(len(event_samples))
            rng.shuffle(idx)
            return [event_samples[i] for i in idx]

        if self.processing_order == 'descending':
            return sorted(event_samples, key=lambda x: x[1], reverse=True)

        # Default: ascending
        return sorted(event_samples, key=lambda x: x[1])

    def _apply_collision_numeric(self, block: Dict[str, Any], new_value: float):
        policy = self.collision_policy
        if policy == 'overwrite':
            block['value'] = float(new_value)
            block['visit_count'] = int(block.get('visit_count', 1)) + 1
            return

        if policy == 'min':
            block['value'] = float(min(block.get('value', new_value), new_value))
            block['visit_count'] = int(block.get('visit_count', 1)) + 1
            return

        if policy == 'max':
            block['value'] = float(max(block.get('value', new_value), new_value))
            block['visit_count'] = int(block.get('visit_count', 1)) + 1
            return

        # Default: average
        old_count = int(block.get('visit_count', 1))
        new_count = old_count + 1
        block['value'] = (float(block.get('value', 0.0)) * old_count + float(new_value)) / new_count
        block['visit_count'] = new_count

    def _apply_collision_domain(self, block: Dict[str, Any], new_domain: str):
        policy = self.collision_policy
        new_domain = str(new_domain).strip()
        if new_domain == '':
            return

        if policy in ('overwrite', 'last'):
            block['domain'] = new_domain
            return

        # Default: vote/majority
        counts = block.get('domain_counts')
        if counts is None:
            counts = {}
            # Seed with existing domain if present
            existing = block.get('domain')
            if existing is not None and str(existing).strip() != '':
                counts[str(existing).strip()] = 1
            block['domain_counts'] = counts

        counts[new_domain] = counts.get(new_domain, 0) + 1
        # Pick winner by max count, tie-break lexicographically for determinism
        winner = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        block['domain'] = winner

    def _create_path_blocks(
        self,
        start_grid: np.ndarray,
        end_grid: np.ndarray,
        start_value: float,
        event_id: int,
        feeder_grid_int: Tuple[int, int, int],
        allowed_positions: Optional[set],
        sample_domain: Optional[str],
        rng: np.random.Generator,
    ):
        # Straight segment in grid coordinates.
        start = np.array(start_grid, dtype=float)
        end = np.array(end_grid, dtype=float)

        vec = end - start
        dist = float(np.linalg.norm(vec))
        if dist < 1e-6:
            return

        # Optional curved path (quadratic Bézier)
        use_curve = self.path_mode in ('curved', 'bezier') and self.curve_strength > 0.0
        control = None
        if use_curve:
            # Build a control point offset from the midpoint in a direction perpendicular to the path.
            dir_hat = vec / dist
            rand_vec = rng.normal(size=3)
            ortho = rand_vec - np.dot(rand_vec, dir_hat) * dir_hat
            ortho_norm = float(np.linalg.norm(ortho))
            if ortho_norm < 1e-6:
                # Fallback orthogonal vector
                ortho = np.array([dir_hat[1], -dir_hat[0], 0.0])
                ortho_norm = float(np.linalg.norm(ortho))
                if ortho_norm < 1e-6:
                    ortho = np.array([0.0, 0.0, 1.0])
                    ortho_norm = 1.0
            ortho = ortho / ortho_norm
            offset = ortho * (self.curve_strength * dist)
            control = 0.5 * (start + end) + offset

        steps = max(2, int(np.ceil(dist * max(self.steps_per_block, 0.5))))

        seen = set()
        for step in range(1, steps):
            t = step / steps
            if use_curve and control is not None:
                # Quadratic Bézier: (1-t)^2*start + 2(1-t)t*control + t^2*end
                pos = (1 - t) ** 2 * start + 2 * (1 - t) * t * control + t ** 2 * end
            else:
                pos = start + t * vec
            block_pos = tuple(np.round(pos).astype(int))

            # Skip start voxel; allow feeder voxel.
            if block_pos == tuple(np.round(start).astype(int)):
                continue

            if self.dims is not None:
                if not (0 <= block_pos[0] < self.dims[0] and 0 <= block_pos[1] < self.dims[1] and 0 <= block_pos[2] < self.dims[2]):
                    continue

            if allowed_positions is not None and block_pos not in allowed_positions:
                continue

            if block_pos in seen:
                continue
            seen.add(block_pos)

            is_feeder = (block_pos == feeder_grid_int)

            if self.interpolation_target == 'domain':
                if sample_domain is None or str(sample_domain).strip() == '':
                    continue
                if block_pos not in self.blocks:
                    self.blocks[block_pos] = {
                        'value': 0.0,
                        'is_sample': False,
                        'visit_count': 1,
                        'event_id': event_id,
                        'is_feeder': is_feeder,
                        'distance_to_feeder': float(dist * (1.0 - t)),
                        'domain': str(sample_domain).strip(),
                        'domain_counts': {str(sample_domain).strip(): 1},
                    }
                    self.domain_mapping[block_pos] = str(sample_domain).strip()
                else:
                    block = self.blocks[block_pos]
                    if not block.get('is_sample', False):
                        self._apply_collision_domain(block, sample_domain)
                        self.domain_mapping[block_pos] = block.get('domain', '')
                        # Prefer closer-to-feeder distance as dominant flow
                        cur = block.get('distance_to_feeder', float('inf'))
                        d2 = float(dist * (1.0 - t))
                        if d2 < cur:
                            block['distance_to_feeder'] = d2
                            block['event_id'] = event_id
                            block['is_feeder'] = is_feeder
                continue

            # Numeric mode: deposit start_value along the path.
            if block_pos not in self.blocks:
                self.blocks[block_pos] = {
                    'value': float(start_value),
                    'is_sample': False,
                    'visit_count': 1,
                    'event_id': event_id,
                    'is_feeder': is_feeder,
                    'distance_to_feeder': float(dist * (1.0 - t)),
                }
            else:
                block = self.blocks[block_pos]
                if not block.get('is_sample', False):
                    self._apply_collision_numeric(block, float(start_value))
                    cur = block.get('distance_to_feeder', float('inf'))
                    d2 = float(dist * (1.0 - t))
                    if d2 < cur:
                        block['distance_to_feeder'] = d2
                        block['event_id'] = event_id
                        block['is_feeder'] = is_feeder

    def _compute_feeders_with_molecular_clock(self) -> List[Dict[str, Any]]:
        """Run Molecular Clock just to obtain event splits + feeder roots."""
        mc_params = dict(self.molecular_clock_params or {})

        # Map UI naming (viewer uses min_samples/max_samples/detect_multiple/interp_method)
        spatial_weight = mc_params.get('spatial_weight', 1.0)
        attr_weight = mc_params.get('attr_weight', 1.0)
        ancestor_depth_offset = mc_params.get('ancestor_depth_offset', 1.0)
        branch_threshold = mc_params.get('branch_threshold', 2.0)
        min_samples = mc_params.get('min_samples', 3)
        max_samples = mc_params.get('max_samples', 1000)
        detect_multiple = mc_params.get('detect_multiple', True)
        interp_method = mc_params.get('interp_method', 'linear')

        # In categorical mode we must ignore grade differences.
        if self.interpolation_target == 'domain':
            attr_weight = 0.0

        mc = MolecularClockInterpolator(
            spatial_weight=float(spatial_weight),
            attr_weight=float(attr_weight),
            detect_multiple_events=bool(detect_multiple),
            branch_threshold=float(branch_threshold),
            min_samples_per_event=int(min_samples),
            max_samples_per_event=int(max_samples),
            interpolation_method=str(interp_method),
            fill_background=False,
            background_value=0.0,
            ancestor_depth_offset=float(ancestor_depth_offset),
            verbose=bool(self.verbose),
        )

        # Attach domain constraints if present
        mc.allowed_grid_override = getattr(self, 'allowed_grid_override', None)
        mc.domain_mapping = getattr(self, 'domain_mapping', None)

        # Initialize with our samples
        mc.initialize_blocks(self.sample_blocks, self.dims, self.min_bounds, self.block_size, block_domains={})

        # If we are constrained to a domain volume, Molecular Clock will also do its own sub-volume splitting.
        # For feeder roots, we just run detect_intrusion_events + build root per event at the top level.
        sample_indices_all = np.arange(len(mc.sample_coords))
        events = mc.detect_intrusion_events(sample_indices=sample_indices_all)

        feeder_roots: List[Dict[str, Any]] = []

        for event_id, sample_indices in enumerate(events):
            if len(sample_indices) == 0:
                continue

            # Subsample large events to keep performance bounded
            if max_samples and max_samples > 0 and len(sample_indices) > max_samples:
                rng = np.random.default_rng(self.seed + event_id)
                sample_indices = rng.choice(sample_indices, size=max_samples, replace=False)

            tree_info = mc.build_hierarchical_tree(sample_indices)
            root_idx = tree_info['root_idx']
            root_coords = tree_info['node_data'][root_idx]['coords']

            # Convert coords to grid int tuple
            feeder_grid = tuple(np.round(np.array(root_coords, dtype=float)).astype(int))

            feeder_roots.append({
                'event_id': event_id,
                'feeder_grid': feeder_grid,
                'root_node': root_idx,
                'num_samples': int(len(sample_indices)),
                'sample_indices': sample_indices,
            })

        return feeder_roots

    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        if self.converged:
            return False

        if self.verbose:
            print("\n=== Running Rectified Flow Interpolation ===")

        # Stage 1: feeder roots via Molecular Clock
        self.feeder_roots = self._compute_feeders_with_molecular_clock()

        if self.verbose:
            print(f"Rectified Flow: detected {len(self.feeder_roots)} feeder event(s)")

        # If no feeders found, fall back to a single "feeder" at the centroid of samples.
        if len(self.feeder_roots) == 0:
            if len(self.sample_blocks) == 0:
                self.converged = True
                return False
            coords = np.array(list(self.sample_blocks.keys()), dtype=float)
            centroid = np.mean(coords, axis=0)
            feeder_grid = tuple(np.round(centroid).astype(int))
            self.feeder_roots = [{
                'event_id': 0,
                'feeder_grid': feeder_grid,
                'root_node': None,
                'num_samples': len(self.sample_blocks),
                'sample_indices': np.arange(len(self.sample_blocks)),
            }]

        allowed_positions = getattr(self, 'allowed_grid_override', None) if self.use_domain_mapping else None

        # Build a mapping from sample position -> event_id based on Molecular Clock clustering
        # Molecular Clock returns indices into its internal sample list; we reconstruct positions.
        # We stored sample_indices in feeder_roots; use that to tag samples.
        all_sample_positions = list(self.sample_blocks.keys())
        for feeder in self.feeder_roots:
            eid = int(feeder['event_id'])
            for idx in feeder.get('sample_indices', []):
                if 0 <= int(idx) < len(all_sample_positions):
                    self.sample_to_event[all_sample_positions[int(idx)]] = eid

        # Stage 2: rectified-flow-style paths to feeder roots
        rng = np.random.default_rng(self.seed)

        # Group samples by event
        samples_by_event: Dict[int, List[Tuple[Tuple[int, int, int], float]]] = {}
        for pos, val in self.sample_blocks.items():
            eid = self.sample_to_event.get(pos, 0)
            samples_by_event.setdefault(eid, []).append((pos, float(val) if val is not None else 0.0))

        for feeder in self.feeder_roots:
            eid = int(feeder['event_id'])
            feeder_grid = feeder['feeder_grid']

            feeder_grid_arr = np.array(feeder_grid, dtype=float)

            event_samples = samples_by_event.get(eid, [])
            if not event_samples:
                continue

            ordered_samples = self._iter_samples_in_order(event_samples)

            for sample_pos, sample_value in ordered_samples:
                # In domain mode, value is irrelevant; we propagate sample's domain label.
                sample_domain = None
                if self.interpolation_target == 'domain':
                    sample_domain = self.sample_domain_mapping.get(sample_pos)
                    if sample_domain is None or str(sample_domain).strip() == '':
                        continue

                start_grid = np.array(sample_pos, dtype=float)

                # Sub-voxel Gaussian endpoint around feeder.
                sigma = max(self.sigma_blocks, 0.0)
                noise = rng.normal(0.0, sigma, size=3)
                noise[2] *= max(self.z_sigma_multiplier, 0.0)
                end_grid = feeder_grid_arr + noise

                self._create_path_blocks(
                    start_grid=start_grid,
                    end_grid=end_grid,
                    start_value=sample_value,
                    event_id=eid,
                    feeder_grid_int=tuple(int(x) for x in feeder_grid),
                    allowed_positions=allowed_positions,
                    sample_domain=sample_domain,
                    rng=rng,
                )

            # Ensure feeder cell exists and is tagged.
            if feeder_grid not in self.blocks:
                self.blocks[feeder_grid] = {
                    'value': float(np.mean([v for _, v in event_samples])) if self.interpolation_target != 'domain' else 0.0,
                    'is_sample': False,
                    'visit_count': 1,
                    'event_id': eid,
                    'is_feeder': True,
                    'distance_to_feeder': 0.0,
                }
            else:
                b = self.blocks[feeder_grid]
                b['event_id'] = eid
                b['is_feeder'] = True
                b['distance_to_feeder'] = 0.0
                if self.interpolation_target == 'domain':
                    # For domain mode, keep domain as current majority.
                    dom = b.get('domain', '')
                    if dom:
                        self.domain_mapping[feeder_grid] = dom

        self.converged = True
        return False

    def is_converged(self) -> bool:
        return self.converged

    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        if self.interpolation_target == 'domain':
            # Viewer expects numeric values here; domains are carried separately.
            return {pos: 0.0 for pos in self.blocks.keys()}
        return {pos: float(block.get('value', 0.0)) for pos, block in self.blocks.items()}

    def get_metadata(self) -> Dict[str, Any]:
        return {
            'algorithm': self.get_algorithm_name(),
            'num_feeders': len(self.feeder_roots),
            'feeder_roots': self.feeder_roots,
            'sigma_blocks': self.sigma_blocks,
            'seed': self.seed,
            'collision_policy': self.collision_policy,
            'processing_order': self.processing_order,
            'interpolation_target': self.interpolation_target,
            'total_blocks': len(self.blocks),
        }
"""
