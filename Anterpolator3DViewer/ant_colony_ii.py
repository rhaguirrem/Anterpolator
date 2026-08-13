import random
from typing import Dict, Optional, Tuple

from tqdm import tqdm

from ant_colony import AntColonyInterpolator, Block


class AntColonyIIInterpolator(AntColonyInterpolator):
    def __init__(
        self,
        range_size=10,
        max_pheromone=1000,
        ants_per_sample=3,
        verbose=False,
        background_value=0.0,
        background_distance=None,
        average_with_blocks=False,
        interpolation_target: str = 'value',
        avoid_visited_threshold_enabled: bool = False,
        avoid_visited_threshold: int = 100,
        ants_sampling_percentage: float = 100.0,
        pheromone_decay_rate: int = 1,
        explore_bias: float = 2.5,
        trail_bias: float = 1.5,
        same_mark_bias: float = 2.0,
        directional_alignment_bias: float = 0.0,
        directional_cross_class_cap: float = 0.3,
        return_bias: float = 4.0,
        revisit_penalty: float = 0.2,
        background_return_enabled: bool = True,
    ):
        super().__init__(
            range_size=range_size,
            max_pheromone=max_pheromone,
            ants_per_sample=ants_per_sample,
            verbose=verbose,
            background_value=background_value,
            background_distance=background_distance,
            average_with_blocks=average_with_blocks,
            interpolation_target=interpolation_target,
            prioritize_blank_unvisited_neighbors=True,
            avoid_visited_threshold_enabled=avoid_visited_threshold_enabled,
            avoid_visited_threshold=avoid_visited_threshold,
            ants_sampling_percentage=ants_sampling_percentage,
            pheromone_decay_rate=pheromone_decay_rate,
        )
        self.explore_bias = float(explore_bias)
        self.trail_bias = float(trail_bias)
        self.same_mark_bias = float(same_mark_bias)
        self.directional_alignment_bias = float(directional_alignment_bias)
        self.directional_cross_class_cap = min(max(float(directional_cross_class_cap), 0.0), 1.0)
        self.return_bias = float(return_bias)
        self.revisit_penalty = float(revisit_penalty)
        self.background_return_enabled = bool(background_return_enabled)

    def get_algorithm_name(self) -> str:
        return "Ant Colony II"

    def get_metadata(self) -> Dict:
        metadata = super().get_metadata()
        metadata.update({
            'explore_bias': self.explore_bias,
            'trail_bias': self.trail_bias,
            'same_mark_bias': self.same_mark_bias,
            'directional_alignment_bias': self.directional_alignment_bias,
            'directional_cross_class_cap': self.directional_cross_class_cap,
            'return_bias': self.return_bias,
            'revisit_penalty': self.revisit_penalty,
            'background_return_enabled': self.background_return_enabled,
        })
        return metadata

    def _get_directional_alignment_multiplier(self, ant, current_block: Block, npos: Tuple[int, int, int]) -> float:
        bias = max(self.directional_alignment_bias, 0.0)
        previous_pos = ant.previous_pos
        if bias <= 0.0 or previous_pos is None:
            return 1.0

        forward_vector = tuple(curr - prev for prev, curr in zip(previous_pos, ant.current_pos))
        candidate_vector = tuple(candidate - curr for curr, candidate in zip(ant.current_pos, npos))
        forward_norm_sq = sum(component * component for component in forward_vector)
        candidate_norm_sq = sum(component * component for component in candidate_vector)
        if forward_norm_sq <= 0 or candidate_norm_sq <= 0:
            return 1.0

        dot_product = sum(a * b for a, b in zip(forward_vector, candidate_vector))
        if dot_product <= 0:
            return 1.0

        forward_norm = forward_norm_sq ** 0.5
        candidate_norm = candidate_norm_sq ** 0.5
        alignment = dot_product / (forward_norm * candidate_norm)
        alignment = min(max(alignment, 0.0), 1.0)
        if alignment <= 0.0:
            return 1.0

        target_mark_class = ant.mark_class
        current_mark_class = getattr(current_block, 'mark_class', None)
        if current_mark_class is not None:
            target_mark_class = current_mark_class

        class_similarity = 1.0
        candidate_block = self.blocks.get(npos)
        if candidate_block is not None:
            mark_delta = abs(float(candidate_block.mark_class) - float(target_mark_class))
            class_similarity = 1.0 / (1.0 + mark_delta / max(float(self.range_size), 1e-9))
            if mark_delta > 0.0:
                class_similarity = min(class_similarity, self.directional_cross_class_cap)

        return 1.0 + (bias * alignment * class_similarity)

    def create_ants(self):
        super().create_ants()
        for ant in self.ants:
            setattr(ant, 'returning', False)

    def calculate_block_value(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int], ant) -> float:
        value = float(super().calculate_block_value(pos, dims, ant))
        block = self.blocks.get(pos)
        if block is None or block.is_sample:
            return value
        if self.interpolation_target != 'value':
            return value
        if self.background_distance is None:
            return value

        distance = max(float(block.distance_to_sample), 0.0)
        background_distance = float(self.background_distance)
        if background_distance <= 0.0:
            return float(self.background_value)
        if distance >= background_distance:
            return float(self.background_value)

        # Ant Colony II should progressively relax toward background instead of
        # only snapping to it at the cutoff distance.
        blend = min(distance / background_distance, 1.0)
        blend *= blend
        blended_value = ((1.0 - blend) * value) + (blend * float(self.background_value))
        return float(round(blended_value, 2))

    def _is_returning(self, ant, current_block: Block) -> bool:
        returning = bool(getattr(ant, 'returning', False))
        if self.background_return_enabled and self.background_distance is not None:
            if current_block.distance_to_sample >= self.background_distance:
                returning = True
            elif returning and current_block.distance_to_sample <= 1:
                returning = False
        setattr(ant, 'returning', returning)
        return returning

    def _get_candidate(self, ant, current_block: Block, npos: Tuple[int, int, int], is_returning: bool) -> Optional[Dict]:
        if npos not in self.allowed_positions:
            return None

        domain = ant.domain or "default"
        if hasattr(self, 'domain_mapping'):
            domain = self.domain_mapping.get(npos, "Undomained")
            if domain != ant.domain:
                return None
            if self.avoid_visited_threshold_enabled and domain in self.domains_frozen:
                return None

        block = self.blocks.get(npos)
        score = 1.0
        previous_pos = ant.previous_pos

        if block is None:
            if is_returning:
                return None
            score *= max(self.explore_bias, 0.01)
            score *= 1.0 / (1.0 + current_block.distance_to_sample + 1)
            if previous_pos is not None and npos == previous_pos:
                score *= 0.5
            score *= self._get_directional_alignment_multiplier(ant, current_block, npos)
            return {'pos': npos, 'score': max(score, 1e-6), 'domain': domain}

        if block.domain != ant.domain:
            return None
        if self.avoid_visited_threshold_enabled and block.domain in self.domains_frozen:
            return None

        if previous_pos is not None and npos == previous_pos:
            score *= 1.25 if is_returning else 0.45

        pheromone_ratio = float(block.pheromone) / float(max(self.max_pheromone, 1))
        score *= 1.0 + self.trail_bias * max(0.0, pheromone_ratio)

        if block.mark_class == ant.mark_class:
            score *= max(self.same_mark_bias, 0.01)
        else:
            mark_delta = abs(block.mark_class - ant.mark_class)
            score *= 1.0 / (1.0 + mark_delta / max(float(self.range_size), 1e-9))

        score *= self._get_directional_alignment_multiplier(ant, current_block, npos)

        if block.nearest_sample_value == current_block.nearest_sample_value:
            score *= 1.15

        if not block.is_sample:
            score *= 1.0 / (1.0 + max(0, block.visit_count - 1) * max(self.revisit_penalty, 0.0))

        if is_returning:
            distance_gain = current_block.distance_to_sample - block.distance_to_sample
            if distance_gain > 0:
                score *= 1.0 + self.return_bias * distance_gain
            else:
                score *= 0.15
            if not block.visited and not block.is_sample:
                score *= 0.2
        else:
            if not block.visited:
                score *= max(self.explore_bias, 0.01)
            elif not block.is_sample:
                score *= 1.05
            if self.background_distance is not None and current_block.distance_to_sample >= max(self.background_distance - 1, 0):
                if block.distance_to_sample > current_block.distance_to_sample:
                    score *= 0.5

        return {'pos': npos, 'score': max(score, 1e-6), 'domain': domain}

    def _select_next_position(self, ant, current_block: Block, dims: Tuple[int, int, int]) -> Optional[Tuple[int, int, int]]:
        neighbors = self.get_neighbors(ant.current_pos, dims)
        random.shuffle(neighbors)
        is_returning = self._is_returning(ant, current_block)

        candidates = []
        for npos in neighbors:
            candidate = self._get_candidate(ant, current_block, npos, is_returning)
            if candidate is not None:
                candidates.append(candidate)

        if not candidates and ant.previous_pos in neighbors and ant.previous_pos in self.blocks:
            previous_block = self.blocks[ant.previous_pos]
            if previous_block.domain == ant.domain:
                return ant.previous_pos
        if not candidates:
            return None

        positions = [candidate['pos'] for candidate in candidates]
        weights = [candidate['score'] for candidate in candidates]
        return random.choices(positions, weights=weights, k=1)[0]

    def move_ants(self, dims: Tuple[int, int, int], quiet=False):
        if self.interpolation_target == 'domain':
            return super().move_ants(dims, quiet=quiet)
        if not quiet:
            print("Moving ants (Ant Colony II)...")

        changes_made = False
        for block in self.blocks.values():
            block.ant_count = 0

        if self.avoid_visited_threshold_enabled and self.domains_frozen:
            self.ants = [ant for ant in self.ants if ant.domain not in self.domains_frozen]

        for ant in self.ants:
            if ant.current_pos in self.blocks:
                self.blocks[ant.current_pos].ant_count += 1

        ant_iterator = self.ants if quiet else tqdm(self.ants)
        for ant in ant_iterator:
            if self.avoid_visited_threshold_enabled and ant.domain in self.domains_frozen:
                continue
            if ant.current_pos not in self.blocks:
                continue

            current_block = self.blocks[ant.current_pos]
            next_pos = self._select_next_position(ant, current_block, dims)
            if next_pos is None:
                continue

            if next_pos not in self.blocks:
                domain = ant.domain or "default"
                if hasattr(self, 'domain_mapping'):
                    domain = self.domain_mapping.get(next_pos, "Undomained")
                self.blocks[next_pos] = Block(
                    value=None,
                    is_sample=False,
                    mark_class=current_block.mark_class,
                    block_id=self.next_block_id,
                    pheromone=max(0, min(self.max_pheromone, current_block.pheromone + 1)),
                    visited=False,
                    visit_count=0,
                    ant_count=0,
                    distance_to_sample=current_block.distance_to_sample + 1,
                    nearest_sample_value=current_block.nearest_sample_value,
                    domain=domain,
                )
                self.domain_created_counts[domain] = self.domain_created_counts.get(domain, 0) + 1
                next_block = self.blocks[next_pos]
                if self.avoid_visited_threshold_enabled and next_block.visit_count >= self.avoid_visited_threshold and not next_block.heavy_reached:
                    next_block.heavy_reached = True
                    self.domain_heavy_counts[domain] = self.domain_heavy_counts.get(domain, 0) + 1
                    self._maybe_freeze_domain(domain)
                self.next_block_id += 1
                next_block.value = self.calculate_block_value(next_pos, dims, ant)
                next_block.mark_class = self.get_mark_class(next_block.value)
                next_block.visited = True
                next_block.visit_count = 1
                next_block.ant_count = 1
                changes_made = True
            else:
                next_block = self.blocks[next_pos]
                if not next_block.is_sample:
                    reinforced = max(current_block.pheromone, next_block.pheromone) + 1
                    next_block.pheromone = min(self.max_pheromone, reinforced)
                    if not current_block.is_sample:
                        current_block.pheromone = min(self.max_pheromone, max(current_block.pheromone, next_block.pheromone - 1))

                    if current_block.distance_to_sample + 1 < next_block.distance_to_sample:
                        next_block.distance_to_sample = current_block.distance_to_sample + 1
                        next_block.nearest_sample_value = current_block.nearest_sample_value
                    elif next_block.distance_to_sample + 1 < current_block.distance_to_sample and not current_block.is_sample:
                        current_block.distance_to_sample = next_block.distance_to_sample + 1
                        current_block.nearest_sample_value = next_block.nearest_sample_value

                    should_update_value = True
                    if self.avoid_visited_threshold_enabled and next_block.visit_count >= self.avoid_visited_threshold:
                        should_update_value = False
                    if should_update_value:
                        next_block.value = self.calculate_block_value(next_pos, dims, ant)
                        next_block.mark_class = self.get_mark_class(next_block.value)
                    next_block.visited = True
                    next_block.visit_count += 1
                    if self.avoid_visited_threshold_enabled and not next_block.heavy_reached and next_block.visit_count >= self.avoid_visited_threshold:
                        next_block.heavy_reached = True
                        self.domain_heavy_counts[next_block.domain] = self.domain_heavy_counts.get(next_block.domain, 0) + 1
                        self._maybe_freeze_domain(next_block.domain)
                    next_block.ant_count += 1
                    changes_made = True

            if ant.current_pos in self.blocks:
                self.blocks[ant.current_pos].ant_count = max(0, self.blocks[ant.current_pos].ant_count - 1)
            ant.previous_pos = ant.current_pos
            ant.current_pos = next_pos
            ant.steps += 1
            ant.steps_history.append(next_pos)
            self._is_returning(ant, self.blocks[next_pos])

        return changes_made