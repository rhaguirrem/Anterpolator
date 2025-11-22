import numpy as np
from functools import lru_cache
from dataclasses import dataclass
from typing import Dict, List, Tuple
import random
from tqdm import tqdm
from interpolator_base import InterpolatorBase

@dataclass
class Block:
    value: float
    is_sample: bool
    mark_class: float
    block_id: int = 0
    pheromone: int = 0
    visited: bool = False
    visit_count: int = 0
    ant_count: int = 0
    distance_to_sample: int = 999999
    nearest_sample_value: float = 0.0
    domain: str = ""
    heavy_reached: bool = False

@dataclass
class Ant:
    origin_block: Tuple[int, int, int]
    current_pos: Tuple[int, int, int]
    mark_class: float
    steps: int = 0
    origin_sample_id: int = 0
    domain: str = ""
    steps_history: List[Tuple[int, int, int]] = None

    def __post_init__(self):
        self.steps_history = []

class AntColonyInterpolator(InterpolatorBase):
    def __init__(self, range_size=10, max_pheromone=150, ants_per_sample=3, verbose=False, background_value=0.0, background_distance=None, average_with_blocks=False,
                 avoid_visited_threshold_enabled: bool = False,
                 avoid_visited_threshold: int = 100):
        super().__init__(verbose)
        self.range_size = range_size
        self.max_pheromone = max_pheromone
        self.ants_per_sample = ants_per_sample
        # self.blocks inherited from InterpolatorBase
        self.ants: List[Ant] = []
        self.next_block_id = 1
        self.allowed_positions = set()
        self.background_value = background_value
        self.background_distance = background_distance
        self.average_with_blocks = average_with_blocks
        self.avoid_visited_threshold_enabled = avoid_visited_threshold_enabled
        self.avoid_visited_threshold = int(avoid_visited_threshold) if avoid_visited_threshold is not None else 100
        # Domain-level tracking for early stopping
        self.domain_totals: Dict[str, int] = {}
        self.domain_created_counts: Dict[str, int] = {}
        self.domain_heavy_counts: Dict[str, int] = {}
        self.domains_frozen: set[str] = set()

    def _maybe_freeze_domain(self, domain: str):
        if not self.avoid_visited_threshold_enabled:
            return
        total = self.domain_totals.get(domain)
        created = self.domain_created_counts.get(domain, 0)
        heavy = self.domain_heavy_counts.get(domain, 0)
        if total is None:
            return
        if created >= total and heavy >= total:
            self.domains_frozen.add(domain)

    def get_mark_class(self, value: float) -> float:
        decimals = len(str(self.range_size).split('.')[-1]) if '.' in str(self.range_size) else 0
        mark_class = value - (value % self.range_size) + self.range_size / 2
        return round(mark_class, decimals)

    def initialize_blocks(self, sample_blocks: Dict[Tuple[int, int, int], float], dims: Tuple[int, int, int], min_bounds, block_size, use_domain_mapping: bool = False):
        print(f"Total volume: {dims[0]*dims[1]*dims[2]} blocks")
        print(f"Initializing {len(sample_blocks)} sample blocks...")
        self.dims = dims
        if use_domain_mapping:
            if hasattr(self, 'allowed_grid_override'):
                self.allowed_positions = self.allowed_grid_override
            else:
                self.allowed_positions = set(sample_blocks.keys())
        else:
            self.allowed_positions = {(i, j, k) for i in range(dims[0]) for j in range(dims[1]) for k in range(dims[2])}
        # Build domain totals for early stopping
        self.domain_totals.clear(); self.domain_created_counts.clear(); self.domain_heavy_counts.clear(); self.domains_frozen.clear()
        if use_domain_mapping and hasattr(self, 'domain_mapping'):
            for pos in self.allowed_positions:
                d = self.domain_mapping.get(pos, "Undomained")
                self.domain_totals[d] = self.domain_totals.get(d, 0) + 1
        else:
            self.domain_totals["default"] = len(self.allowed_positions)
        for d in list(self.domain_totals.keys()):
            self.domain_created_counts[d] = 0
            self.domain_heavy_counts[d] = 0
        self.blocks = {}
        for pos, value in sample_blocks.items():
            mark_class = self.get_mark_class(value)
            domain = self.domain_mapping.get(pos, "Undomained") if use_domain_mapping and hasattr(self, 'domain_mapping') else "default"
            self.blocks[pos] = Block(
                value=value,
                is_sample=True,
                mark_class=mark_class,
                block_id=self.next_block_id,
                pheromone=self.max_pheromone,
                visited=True,
                visit_count=1,
                ant_count=0,
                distance_to_sample=0,
                nearest_sample_value=value,
                domain=domain
            )
            # Track domain created/heavy counters
            self.domain_created_counts[domain] = self.domain_created_counts.get(domain, 0) + 1
            if self.avoid_visited_threshold_enabled and self.blocks[pos].visit_count >= self.avoid_visited_threshold:
                self.blocks[pos].heavy_reached = True
                self.domain_heavy_counts[domain] = self.domain_heavy_counts.get(domain, 0) + 1
                self._maybe_freeze_domain(domain)
            self.next_block_id += 1
            # if self.verbose:
            #     print(f"Sample at {pos}: value={value:.2f}, mark_class={mark_class}, domain={domain}")

    @lru_cache(maxsize=1024)
    def get_neighbors(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        x, y, z = pos
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    if dx == dy == dz == 0:
                        continue
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if 0 <= nx < dims[0] and 0 <= ny < dims[1] and 0 <= nz < dims[2]:
                        neighbors.append((nx, ny, nz))
        return neighbors

    def create_ants(self):
        self.ants.clear()
        for pos, block in self.blocks.items():
            if block.is_sample:
                for _ in range(self.ants_per_sample):
                    ant = Ant(pos, pos, block.mark_class)
                    ant.origin_sample_id = block.block_id
                    ant.domain = block.domain
                    self.ants.append(ant)
        
        # Debug: count ants per domain
        if self.verbose:
            domain_ant_counts = {}
            for ant in self.ants:
                domain_ant_counts[ant.domain] = domain_ant_counts.get(ant.domain, 0) + 1
            print(f"\nAnts created per domain:")
            for domain, count in sorted(domain_ant_counts.items()):
                print(f"  {domain}: {count} ants")
            print(f"  Total ants: {len(self.ants)}")

    def calculate_block_value(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int], ant: Ant) -> float:
        neighbors = self.get_neighbors(pos, dims)
        values = []
        weights = []
        current_block = self.blocks.get(pos)
        current_mark_class = current_block.mark_class if current_block else None
        origin_block_value = self.blocks[ant.origin_block].value if ant and ant.origin_block in self.blocks else None
        for npos in neighbors:
            if npos in self.blocks:
                nblock = self.blocks[npos]
                values.append(nblock.value)
                weight = 1.0
                # Mark-class similarity
                if current_mark_class is not None:
                    mark_class_diff = abs(nblock.value - nblock.mark_class)
                    weight *= 1.0 / (1.0 + mark_class_diff / self.range_size)
                # Distance to nearest sample
                weight *= 1.0 / (1.0 + nblock.distance_to_sample)
                # Consistency with nearest sample value (favor closer to sample)
                if not nblock.is_sample:
                    sample_diff = abs(nblock.value - nblock.nearest_sample_value)
                    weight *= 1.0 / (1.0 + sample_diff / self.range_size)
                else:
                    weight *= 3.0
                # Proximity of neighbor's nearest sample to current cell's nearest sample (from backup)
                if current_block is not None:
                    nearest_sample_diff = abs(nblock.nearest_sample_value - current_block.nearest_sample_value)
                    weight *= 1.0 / (1.0 + nearest_sample_diff / self.range_size)
                # Value closeness to neighbor's nearest sample
                value_to_sample_diff = abs(nblock.value - nblock.nearest_sample_value)
                weight *= 1.0 / (1.0 + value_to_sample_diff / self.range_size)
                # Value closeness to ant's origin sample value
                if origin_block_value is not None:
                    value_to_origin_diff = abs(nblock.value - origin_block_value)
                    weight *= 1.0 / (1.0 + value_to_origin_diff / self.range_size)
                weights.append(weight)
        block = self.blocks.get(pos)
        if block and self.background_distance is not None and block.distance_to_sample >= self.background_distance:
            return float(self.background_value)
        if values:
            weights = np.array(weights, dtype=float)
            values = np.array(values, dtype=float)
            neighbor_average = float(np.average(values, weights=weights))
            if self.average_with_blocks and block and block.value is not None:
                return float(round((neighbor_average + float(block.value)) / 2, 2))
            return float(round(neighbor_average, 2))
        return float(self.background_value)

    def move_ants(self, dims: Tuple[int, int, int], quiet=False):
        if not quiet:
            print("Moving ants...")
        changes_made = False
        
        # Debug: Track why ants can't move (only for verbose mode and first call)
        debug_first_iteration = self.verbose and not hasattr(self, '_debug_done')
        if debug_first_iteration:
            self._debug_done = True
            debug_stats = {
                'total_ants_checked': 0,
                'neighbors_not_in_allowed': 0,
                'neighbors_domain_mismatch': 0,
                'neighbors_frozen_domain': 0,
                'blocks_created': 0,
                'no_valid_neighbors': 0
            }
        
        for block in self.blocks.values():
            block.ant_count = 0
        # Optionally drop ants belonging to frozen domains to speed up
        if self.avoid_visited_threshold_enabled and self.domains_frozen:
            self.ants = [a for a in self.ants if a.domain not in self.domains_frozen]
        for ant in self.ants:
            if ant.current_pos in self.blocks:
                self.blocks[ant.current_pos].ant_count += 1
        ant_iterator = self.ants if quiet else tqdm(self.ants)
        for ant in ant_iterator:
            if debug_first_iteration:
                debug_stats['total_ants_checked'] += 1
            
            if self.avoid_visited_threshold_enabled and ant.domain in self.domains_frozen:
                continue
            if ant.current_pos not in self.blocks:
                continue
            current_block = self.blocks[ant.current_pos]
            neighbors = self.get_neighbors(ant.current_pos, dims)
            random.shuffle(neighbors)
            next_pos = None
            
            # Track why neighbors are rejected (debug)
            if debug_first_iteration:
                found_valid = False
            
            # Priority 1: Forward expansion to uncreated or unvisited neighbor within domain
            for npos in neighbors:
                if npos not in self.blocks and npos in self.allowed_positions:
                    domain = "default"
                    if hasattr(self, 'domain_mapping'):
                        domain = self.domain_mapping.get(npos, "Undomained")
                        if self.avoid_visited_threshold_enabled and domain in self.domains_frozen:
                            if debug_first_iteration:
                                debug_stats['neighbors_frozen_domain'] += 1
                            continue
                        if domain != ant.domain:
                            if debug_first_iteration:
                                debug_stats['neighbors_domain_mismatch'] += 1
                            continue
                    self.blocks[npos] = Block(
                        value=None,
                        is_sample=False,
                        mark_class=current_block.mark_class,
                        block_id=self.next_block_id,
                        pheromone=max(0, current_block.pheromone - 1),
                        visited=False,
                        visit_count=0,
                        ant_count=0,
                        distance_to_sample=current_block.distance_to_sample + 1,
                        nearest_sample_value=current_block.nearest_sample_value,
                        domain=domain
                    )
                    # Track domain creation
                    self.domain_created_counts[domain] = self.domain_created_counts.get(domain, 0) + 1
                    self.next_block_id += 1
                    next_pos = npos
                    if debug_first_iteration:
                        debug_stats['blocks_created'] += 1
                        found_valid = True
                    break
                elif npos not in self.blocks and npos not in self.allowed_positions:
                    if debug_first_iteration:
                        debug_stats['neighbors_not_in_allowed'] += 1
                elif npos in self.blocks:
                    nblock = self.blocks[npos]
                    if not nblock.visited and nblock.domain == ant.domain and nblock.mark_class <= ant.mark_class:
                        if self.avoid_visited_threshold_enabled and nblock.domain in self.domains_frozen:
                            continue
                        next_pos = npos
                        break
            if next_pos is None:
                # Priority 2: same mark_class, different nearest sample, shortest distance-to-sample
                candidates = [
                    (npos, self.blocks[npos]) for npos in neighbors
                    if npos in self.blocks and self.blocks[npos].mark_class == ant.mark_class and
                    self.blocks[npos].nearest_sample_value != current_block.nearest_sample_value and
                    self.blocks[npos].domain == ant.domain and (not ant.steps_history or npos != ant.steps_history[-1]) and
                    (not self.avoid_visited_threshold_enabled or self.blocks[npos].domain not in self.domains_frozen)
                ]
                if candidates:
                    if self.avoid_visited_threshold_enabled:
                        under = [c for c in candidates if c[1].visit_count < self.avoid_visited_threshold]
                        if under:
                            next_pos = min(under, key=lambda x: x[1].distance_to_sample)[0]
                        # else: do not select here; allow later fallback to random
                    else:
                        next_pos = min(candidates, key=lambda x: x[1].distance_to_sample)[0]
            if next_pos is None:
                # Priority 3: value closer to ant's origin sample
                candidates = [
                    (npos, self.blocks[npos]) for npos in neighbors
                    if npos in self.blocks and self.blocks[npos].domain == ant.domain and self.blocks[npos].mark_class <= ant.mark_class and (not ant.steps_history or npos != ant.steps_history[-1]) and (not self.avoid_visited_threshold_enabled or self.blocks[npos].domain not in self.domains_frozen)
                ]
                if candidates:
                    if self.avoid_visited_threshold_enabled:
                        under = [c for c in candidates if c[1].visit_count < self.avoid_visited_threshold]
                        if under:
                            next_pos = min(under, key=lambda x: abs(x[1].value - current_block.nearest_sample_value))[0]
                        # else: defer to random fallback
                    else:
                        next_pos = min(candidates, key=lambda x: abs(x[1].value - current_block.nearest_sample_value))[0]
            if next_pos is None:
                # Priority 4: random valid neighbor
                valid_neighbors = [
                    npos for npos in neighbors
                    if npos in self.blocks and self.blocks[npos].domain == ant.domain and self.blocks[npos].mark_class <= ant.mark_class and (not ant.steps_history or npos != ant.steps_history[-1]) and (not self.avoid_visited_threshold_enabled or self.blocks[npos].domain not in self.domains_frozen)
                ]
                if valid_neighbors:
                    if self.avoid_visited_threshold_enabled:
                        under = [n for n in valid_neighbors if self.blocks[n].visit_count < self.avoid_visited_threshold]
                        if under:
                            next_pos = random.choice(under)
                        else:
                            # all neighbors are heavily visited, pick random neighbor as requested
                            next_pos = random.choice(valid_neighbors)
                    else:
                        next_pos = random.choice(valid_neighbors)
                else:
                    next_pos = ant.steps_history[-1] if ant.steps_history else ant.current_pos

            # Apply move
            if next_pos not in self.blocks:
                if next_pos in self.allowed_positions:
                    domain = "default"
                    if hasattr(self, 'domain_mapping'):
                        domain = self.domain_mapping.get(next_pos, "Undomained")
                        if self.avoid_visited_threshold_enabled and domain in self.domains_frozen:
                            ant.steps_history.append(ant.current_pos)
                            continue
                        if domain != ant.domain:
                            ant.steps_history.append(ant.current_pos)
                            continue
                    self.blocks[next_pos] = Block(
                        value=0.0,
                        is_sample=False,
                        mark_class=current_block.mark_class,
                        block_id=self.next_block_id,
                        pheromone=max(0, current_block.pheromone - 1),
                        visited=True,
                        visit_count=1,
                        ant_count=1,
                        distance_to_sample=current_block.distance_to_sample + 1,
                        nearest_sample_value=current_block.nearest_sample_value,
                        domain=domain
                    )
                    # Track domain creation and heavy
                    self.domain_created_counts[domain] = self.domain_created_counts.get(domain, 0) + 1
                    if self.avoid_visited_threshold_enabled and self.blocks[next_pos].visit_count >= self.avoid_visited_threshold and not self.blocks[next_pos].heavy_reached:
                        self.blocks[next_pos].heavy_reached = True
                        self.domain_heavy_counts[domain] = self.domain_heavy_counts.get(domain, 0) + 1
                        self._maybe_freeze_domain(domain)
                    self.next_block_id += 1
                    self.blocks[next_pos].value = self.calculate_block_value(next_pos, dims, ant)
                    self.blocks[next_pos].mark_class = self.get_mark_class(self.blocks[next_pos].value)
                    changes_made = True
                else:
                    ant.steps_history.append(ant.current_pos)
                    continue
            else:
                next_block = self.blocks[next_pos]
                if not next_block.is_sample:
                    if next_block.pheromone < current_block.pheromone:
                        next_block.pheromone = max(0, current_block.pheromone - 1)
                        next_block.mark_class = current_block.mark_class
                    elif current_block.pheromone < next_block.pheromone and not current_block.is_sample:
                        current_block.pheromone = max(0, next_block.pheromone - 1)
                        current_block.mark_class = next_block.mark_class
                    if current_block.distance_to_sample < next_block.distance_to_sample:
                        next_block.distance_to_sample = current_block.distance_to_sample + 1
                        next_block.nearest_sample_value = current_block.nearest_sample_value
                    elif next_block.distance_to_sample < current_block.distance_to_sample:
                        current_block.distance_to_sample = next_block.distance_to_sample + 1
                        current_block.nearest_sample_value = next_block.nearest_sample_value
                    # If threshold enabled and this block is heavily visited, optionally avoid updating its value
                    should_update_value = True
                    if self.avoid_visited_threshold_enabled and next_block.visit_count >= self.avoid_visited_threshold:
                        should_update_value = False
                    if should_update_value:
                        next_block.value = self.calculate_block_value(next_pos, dims, ant)
                        next_block.mark_class = self.get_mark_class(next_block.value)
                    next_block.visited = True
                    next_block.visit_count += 1
                    # Heavy tracking
                    if self.avoid_visited_threshold_enabled and not next_block.heavy_reached and next_block.visit_count >= self.avoid_visited_threshold:
                        next_block.heavy_reached = True
                        self.domain_heavy_counts[next_block.domain] = self.domain_heavy_counts.get(next_block.domain, 0) + 1
                        self._maybe_freeze_domain(next_block.domain)
                    next_block.ant_count += 1
                    changes_made = True
            if ant.current_pos in self.blocks:
                self.blocks[ant.current_pos].ant_count = max(0, self.blocks[ant.current_pos].ant_count - 1)
            ant.current_pos = next_pos
            ant.steps += 1
            ant.steps_history.append(next_pos)
            
            # Track ants that found no valid neighbors
            if debug_first_iteration and not found_valid:
                debug_stats['no_valid_neighbors'] += 1
        
        # Print debug summary
        if debug_first_iteration:
            print(f"\n=== First Iteration Debug Statistics ===")
            print(f"Total ants checked: {debug_stats['total_ants_checked']}")
            print(f"Blocks created: {debug_stats['blocks_created']}")
            print(f"Ants with no valid neighbors: {debug_stats['no_valid_neighbors']}")
            print(f"Neighbor rejections:")
            print(f"  - Not in allowed_positions: {debug_stats['neighbors_not_in_allowed']}")
            print(f"  - Domain mismatch: {debug_stats['neighbors_domain_mismatch']}")
            print(f"  - Frozen domain: {debug_stats['neighbors_frozen_domain']}")
            print(f"Changes made: {changes_made}\n")
        
        return changes_made

    def update_block_value(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int]):
        block = self.blocks[pos]
        if block.is_sample:
            return
        neighbors = self.get_neighbors(pos, dims)
        weights = []
        values = []
        for npos in neighbors:
            nblock = self.blocks.get(npos)
            if nblock is None:
                continue
            if nblock.is_sample:
                weights.append(self.max_pheromone * 3.0)
                values.append(nblock.value)
            elif nblock.visited:
                weights.append(1.0)
                values.append(nblock.value)
        if values:
            try:
                weights = np.array(weights, dtype=float)
                if weights.sum() > 0:
                    block.value = np.average(values, weights=weights)
                    block.mark_class = self.get_mark_class(block.value)
            except Exception as e:
                print(f"Error updating block {pos}: {str(e)}")
                print(f"Values: {values}")
                print(f"Weights: {weights}")

    def update_block(self, current_pos, previous_pos):
        current_block = self.blocks[current_pos]
        previous_block = self.blocks[previous_pos]
        if current_block.visit_count == 0:
            current_block.value = max(0, previous_block.value - 1)
            current_block.mark_class = previous_block.mark_class
        else:
            if current_block.pheromone > previous_block.pheromone:
                previous_block.value = max(0, current_block.value - 1)
                previous_block.mark_class = current_block.mark_class
            else:
                current_block.value = max(0, previous_block.value - 1)
                current_block.mark_class = previous_block.mark_class

    def decay_pheromone(self, quiet=False):
        if not quiet:
            print("Updating pheromones...")
        blocks_to_remove = []
        for pos in (tqdm(list(self.blocks.keys())) if not quiet else list(self.blocks.keys())):
            block = self.blocks[pos]
            if not block.is_sample:
                block.pheromone = max(0, block.pheromone - 1)
                if block.pheromone <= 0:
                    blocks_to_remove.append(pos)
        for pos in blocks_to_remove:
            if pos in self.blocks and not self.blocks[pos].is_sample:
                del self.blocks[pos]
        for ant in self.ants:
            if ant.current_pos not in self.blocks:
                ant.current_pos = ant.origin_block
                ant.steps = 0

    def get_blocks_data(self) -> Dict[Tuple[int, int, int], float]:
        return {pos: block.value for pos, block in self.blocks.items() if block.is_sample or block.visited}

    def fill_unvisited_blocks(self, dims: Tuple[int, int, int], quiet=False):
        if not quiet:
            print("Filling unvisited blocks with neighbor averages...")
        created = 0
        updated = 0
        allowed = self.allowed_positions if hasattr(self, 'allowed_positions') and self.allowed_positions else None
        if allowed is None or len(allowed) == 0:
            allowed = {(i, j, k) for i in range(dims[0]) for j in range(dims[1]) for k in range(dims[2])}
        for pos in (tqdm(list(allowed)) if not quiet else allowed):
            b = self.blocks.get(pos)
            domain = self.domain_mapping.get(pos, "Undomained") if hasattr(self, 'domain_mapping') else "default"
            if b is None:
                self.blocks[pos] = Block(
                    value=0.0,
                    is_sample=False,
                    mark_class=0,
                    block_id=self.next_block_id,
                    pheromone=0,
                    visited=False,
                    visit_count=0,
                    ant_count=0,
                    distance_to_sample=999999,
                    nearest_sample_value=0.0,
                    domain=domain,
                )
                self.next_block_id += 1
                neighbors = self.get_neighbors(pos, dims)
                nbr_samples = [self.blocks[n] for n in neighbors if n in self.blocks and self.blocks[n].is_sample]
                if nbr_samples:
                    self.blocks[pos].distance_to_sample = 1 + min(nb.distance_to_sample for nb in nbr_samples)
                    nearest = min(nbr_samples, key=lambda nb: nb.distance_to_sample)
                    self.blocks[pos].nearest_sample_value = nearest.value
                ant_origin = next((p for p, blk in self.blocks.items() if blk.is_sample and blk.domain == domain), None)
                if ant_origin is None:
                    ant_origin = next((p for p, blk in self.blocks.items() if blk.is_sample), None)
                if ant_origin is None:
                    self.blocks[pos].value = self.background_value
                    self.blocks[pos].mark_class = self.get_mark_class(self.blocks[pos].value)
                else:
                    fake_ant = Ant(origin_block=ant_origin, current_pos=pos, mark_class=self.blocks[ant_origin].mark_class)
                    fake_ant.domain = self.blocks[ant_origin].domain
                    val = self.calculate_block_value(pos, dims, fake_ant)
                    self.blocks[pos].value = val
                    self.blocks[pos].mark_class = self.get_mark_class(val)
                created += 1
            else:
                if not b.is_sample and not b.visited:
                    ant_origin = next((p for p, blk in self.blocks.items() if blk.is_sample and blk.domain == b.domain), None)
                    if ant_origin is None:
                        ant_origin = next((p for p, blk in self.blocks.items() if blk.is_sample), None)
                    if ant_origin is not None:
                        fake_ant = Ant(origin_block=ant_origin, current_pos=pos, mark_class=self.blocks[ant_origin].mark_class)
                        fake_ant.domain = self.blocks[ant_origin].domain
                        val = self.calculate_block_value(pos, dims, fake_ant)
                        b.value = val
                        b.mark_class = self.get_mark_class(val)
                        updated += 1
        if not quiet:
            print(f"Fill complete. Created: {created}, Updated: {updated}")
        return created, updated

    def fill_unvisited_blocks_domainwise(self, dims: Tuple[int, int, int]):
        print("Filling unvisited blocks (domain-wise propagation)...")
        # Determine domain sets
        if hasattr(self, 'allowed_positions') and self.allowed_positions:
            all_positions = set(self.allowed_positions)
        else:
            all_positions = {(i, j, k) for i in range(dims[0]) for j in range(dims[1]) for k in range(dims[2])}
        # Helper: get domain for a position
        def get_domain(pos):
            if hasattr(self, 'domain_mapping'):
                return self.domain_mapping.get(pos, "Undomained")
            return "default"
        # Group positions by domain
        domains = {}
        for pos in all_positions:
            d = get_domain(pos)
            domains.setdefault(d, set()).add(pos)
        total_created = 0
        total_assigned = 0
        # For each domain, propagate values via BFS from frontier cells
        domain_iter = tqdm(list(domains.items()), desc="Domains")
        for domain, positions in domain_iter:
            # Known valued positions within this domain
            valued = set()
            for pos in positions:
                b = self.blocks.get(pos)
                if b is not None and b.value is not None:
                    valued.add(pos)
            # Target unvisited positions lacking value
            targets = set()
            for pos in positions:
                b = self.blocks.get(pos)
                if b is None:
                    targets.add(pos)
                else:
                    if not b.is_sample and not b.visited and b.value is None:
                        targets.add(pos)
            if not targets:
                continue
            # Initialize frontier: unvisited cells that have at least one valued neighbor in same domain
            frontier = [pos for pos in targets if any((n in valued) for n in self.get_neighbors(pos, dims) if n in positions)]
            queued = set(frontier)
            # Progress bar for this domain
            pbar = None
            if len(targets) > 0:
                try:
                    pbar = tqdm(total=len(targets), desc=f"Filling {domain}", leave=False)
                except Exception:
                    pbar = None
            created_d = 0
            assigned_d = 0
            # BFS propagation
            while frontier:
                next_frontier = []
                for pos in frontier:
                    if pos not in positions:
                        continue
                    # Collect neighbor values from same domain
                    neighbor_vals = []
                    neighbor_info = []
                    for n in self.get_neighbors(pos, dims):
                        if n not in positions:
                            continue
                        nb = self.blocks.get(n)
                        if nb is not None and nb.value is not None:
                            neighbor_vals.append(float(nb.value))
                            neighbor_info.append((n, nb))
                    if not neighbor_vals:
                        continue
                    avg_val = float(sum(neighbor_vals) / len(neighbor_vals))
                    # Create or update block at pos
                    b = self.blocks.get(pos)
                    updated_now = False
                    if b is None:
                        # derive distance and nearest sample value from nearest neighbor by distance_to_sample
                        dist = 999999
                        ns_val = avg_val
                        if neighbor_info:
                            nearest_nb = min(neighbor_info, key=lambda t: t[1].distance_to_sample)[1]
                            dist = min(999999, nearest_nb.distance_to_sample + 1)
                            ns_val = nearest_nb.nearest_sample_value
                        self.blocks[pos] = Block(
                            value=avg_val,
                            is_sample=False,
                            mark_class=self.get_mark_class(avg_val),
                            block_id=self.next_block_id,
                            pheromone=0,
                            visited=False,
                            visit_count=0,
                            ant_count=0,
                            distance_to_sample=dist,
                            nearest_sample_value=ns_val,
                            domain=domain,
                        )
                        self.next_block_id += 1
                        total_created += 1
                        created_d += 1
                        updated_now = True
                    else:
                        if b.value is None:
                            b.value = avg_val
                            b.mark_class = self.get_mark_class(avg_val)
                            total_assigned += 1
                            assigned_d += 1
                            updated_now = True
                    valued.add(pos)
                    if pbar is not None and updated_now:
                        try:
                            pbar.update(1)
                            pbar.set_postfix(created=created_d, assigned=assigned_d)
                        except Exception:
                            pass
                    # Expand to unvisited same-domain neighbors
                    for nn in self.get_neighbors(pos, dims):
                        if nn in positions and nn not in valued and nn in targets and nn not in queued:
                            # Enqueue if it now has at least one valued neighbor
                            next_frontier.append(nn)
                            queued.add(nn)
                frontier = next_frontier
            if pbar is not None:
                try:
                    pbar.close()
                except Exception:
                    pass
        print(f"Domain-wise fill complete. Created: {total_created}, Newly assigned: {total_assigned}")
        return total_created, total_assigned
    
    # ========== InterpolatorBase abstract methods implementation ==========
    
    def run_iteration(self, dims: Tuple[int, int, int]) -> bool:
        """
        Execute one iteration of the ant colony algorithm.
        Returns True to continue iterating, False when finished.
        """
        self.move_ants(dims, quiet=True)
        self.decay_pheromone(quiet=True)
        
        # Check if should continue
        if self.avoid_visited_threshold_enabled:
            if len(self.domain_totals) > 0 and len(self.domains_frozen) >= len(self.domain_totals):
                return False  # All domains frozen
        
        if len(self.ants) == 0:
            return False  # No ants left
        
        return True  # Continue
    
    def get_interpolated_values(self) -> Dict[Tuple[int, int, int], float]:
        """Return interpolated values for all blocks"""
        return {pos: block.value for pos, block in self.blocks.items()}
    
    def get_algorithm_name(self) -> str:
        """Return algorithm name"""
        return "Ant Colony Optimization"
    
    def is_converged(self) -> bool:
        """Check if algorithm has converged"""
        if self.avoid_visited_threshold_enabled:
            if len(self.domain_totals) > 0 and len(self.domains_frozen) >= len(self.domain_totals):
                return True
        return len(self.ants) == 0
    
    def get_metadata(self) -> Dict:
        """Return algorithm-specific metadata"""
        return {
            'algorithm': self.get_algorithm_name(),
            'total_blocks': len(self.blocks),
            'sample_blocks': sum(1 for b in self.blocks.values() if b.is_sample),
            'interpolated_blocks': sum(1 for b in self.blocks.values() if not b.is_sample),
            'remaining_ants': len(self.ants),
            'range_size': self.range_size,
            'max_pheromone': self.max_pheromone,
            'ants_per_sample': self.ants_per_sample,
            'domains_frozen': list(self.domains_frozen) if self.avoid_visited_threshold_enabled else []
        }
