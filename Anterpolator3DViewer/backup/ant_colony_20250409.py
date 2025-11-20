import numpy as np
from collections import defaultdict
from functools import lru_cache  # Fix import location
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set
import random
from tqdm import tqdm

@dataclass
class Block:
    value: float
    is_sample: bool
    mark_class: int
    block_id: int = 0  # Add unique ID
    pheromone: int = 0
    visited: bool = False  # Add visited flag
    visit_count: int = 0  # Add visit counter
    ant_count: int = 0  # Make sure ant_count is initialized
    distance_to_sample: int = 999999  # Initialize with a large number
    nearest_sample_value: float = 0.0  # Add nearest sample value tracking
    domain: str = ""  # NEW: Domain of the block

@dataclass
class Ant:
    origin_block: Tuple[int, int, int]
    current_pos: Tuple[int, int, int]
    mark_class: int
    steps: int = 0
    origin_sample_id: int = 0  # Add origin sample ID tracking
    domain: str = ""  # NEW: Track ant's domain

class AntColonyInterpolator:
    def __init__(self, range_size=10, max_pheromone=150, ants_per_sample=3, verbose=False):
        self.range_size = range_size
        self.max_pheromone = max_pheromone
        self.ants_per_sample = ants_per_sample
        self.blocks: Dict[Tuple[int, int, int], Block] = {}
        self.ants: List[Ant] = []
        self.verbose = verbose
        self.next_block_id = 1  # Track block IDs
        self._neighbor_cache = {}  # Cache for neighbor calculations
        self.allowed_positions = set()  # NEW: Allowed grid positions from blocks_file
        
    def get_mark_class(self, value: float) -> float:
        decimals = len(str(self.range_size).split('.')[-1]) if '.' in str(self.range_size) else 0
        mark_class = value - (value % self.range_size) + self.range_size / 2
        return round(mark_class, decimals)
    
    def initialize_blocks(self, sample_blocks: Dict[Tuple[int, int, int], float],
                          dims: Tuple[int, int, int],
                          min_bounds, block_size,
                          use_domain_mapping: bool = False):
        print(f"Total volume: {dims[0]*dims[1]*dims[2]} blocks")
        print(f"Initializing {len(sample_blocks)} sample blocks...")
        self.dims = dims  # Store grid dimensions
        if use_domain_mapping:
            if hasattr(self, 'allowed_grid_override'):
                self.allowed_positions = self.allowed_grid_override
            else:
                self.allowed_positions = set(sample_blocks.keys())
        else:
            self.allowed_positions = {(i, j, k) for i in range(dims[0])
                                                  for j in range(dims[1])
                                                  for k in range(dims[2])}
        self.blocks = {}
        for pos, value in sample_blocks.items():
            mark_class = self.get_mark_class(value)
            # If using domain mapping and a domain mapping is attached, assign the domain from it
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
            self.next_block_id += 1
            if self.verbose:
                print(f"Sample at {pos}: value={value:.2f}, mark_class={mark_class}, domain={domain}")
    
    def get_or_create_block(self, pos: Tuple[int, int, int]) -> Block:
        if pos not in self.blocks:
            # Initialize new block with zero values
            self.blocks[pos] = Block(0.0, False, 0, 0, False)
        return self.blocks[pos]
    
    @lru_cache(maxsize=1024)
    def get_neighbors(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int]) -> List[Tuple[int, int, int]]:
        """Cached version of neighbor calculation"""
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
                    ant.domain = block.domain  # NEW: Set ant's domain from sample block
                    self.ants.append(ant)
    
    def calculate_block_value(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int]) -> float:
        """Calculate block value based on weighted average of neighbors"""
        neighbors = self.get_neighbors(pos, dims)
        values = []
        weights = []
        current_block = self.blocks.get(pos)
        current_mark_class = current_block.mark_class if current_block else None
        for npos in neighbors:
            if npos in self.blocks:
                nblock = self.blocks[npos]
                values.append(nblock.value)

                # pheromone factor
                use_pheromone = False
                weight = 1.0 + (nblock.pheromone / self.max_pheromone) if use_pheromone else 1.0
                
                # mark of class factor
                if current_mark_class is not None:
                    mark_class_diff = abs(nblock.value - nblock.mark_class)
                    similarity_factor = 1.0 / (1.0 + mark_class_diff/self.range_size)
                    weight *= similarity_factor
                
                # distance factor
                distance_factor = 1.0 / (1.0 + nblock.distance_to_sample)
                weight *= distance_factor
                
                # nearest sample value factor
                if not nblock.is_sample:
                    sample_diff = abs(nblock.value - nblock.nearest_sample_value)
                    sample_similarity = 1.0 / (1.0 + sample_diff/self.range_size)
                    weight *= sample_similarity
                else:
                    weight *= 3.0
                
                # nearest sample proximity factor
                nearest_sample_diff = abs(nblock.nearest_sample_value - current_block.nearest_sample_value)
                proximity_factor = 1.0 / (1.0 + nearest_sample_diff/self.range_size)
                weight *= proximity_factor
                
                # NEW: value closeness to nearest sample factor
                value_to_sample_diff = abs(nblock.value - nblock.nearest_sample_value)
                value_closeness_factor = 1.0 / (1.0 + value_to_sample_diff/self.range_size)
                weight *= value_closeness_factor
                
                weights.append(weight)
        if values:
            weights = np.array(weights)
            values = np.array(values)
            return round(np.average(values, weights=weights), 2)  # Round to 2 decimals
        return 0.0

    def move_ants(self, dims: Tuple[int, int, int], quiet=False):
        if not quiet:
            print("Moving ants...")
        changes_made = False
        for block in self.blocks.values():
            block.ant_count = 0
        for ant in self.ants:
            if ant.current_pos in self.blocks:
                self.blocks[ant.current_pos].ant_count += 1
        ant_iterator = self.ants if quiet else tqdm(self.ants)
        for ant in ant_iterator:
            if ant.current_pos not in self.blocks:
                continue
            current_block = self.blocks[ant.current_pos]
            neighbors = self.get_neighbors(ant.current_pos, dims)
            
            # Shuffle neighbors to remove bias
            random.shuffle(neighbors)

            # Priority 1: Choose the first found unvisited block (including not yet created blocks)
            for npos in neighbors:
                if npos not in self.blocks and npos in self.allowed_positions:
                    # Create the block as unvisited
                    domain = "default"
                    if hasattr(self, 'domain_mapping'):
                        domain = self.domain_mapping.get(npos, "Undomained")
                        if domain != ant.domain:
                            continue
                    self.blocks[npos] = Block(
                        value=0.0,
                        is_sample=False,
                        mark_class=current_block.mark_class,
                        block_id=self.next_block_id,
                        pheromone=max(0, current_block.pheromone - 1),
                        visited=False,  # Mark as unvisited
                        visit_count=0,
                        ant_count=0,
                        distance_to_sample=current_block.distance_to_sample + 1,
                        nearest_sample_value=current_block.nearest_sample_value,
                        domain=domain
                    )
                    self.next_block_id += 1
                    next_pos = npos
                    break
                elif npos in self.blocks:
                    nblock = self.blocks[npos]
                    if not nblock.visited and nblock.domain == ant.domain:
                        next_pos = npos
                        break
            else:
                # Priority 2: Factor for different nearest sample and greater distance to sample
                different_sample_greater_distance = [
                    (npos, self.blocks[npos]) for npos in neighbors
                    if npos in self.blocks and
                    self.blocks[npos].nearest_sample_value != current_block.nearest_sample_value and
                    self.blocks[npos].distance_to_sample > current_block.distance_to_sample and
                    self.blocks[npos].domain == ant.domain
                ]
                if different_sample_greater_distance:
                    next_pos = random.choice([npos for npos, _ in different_sample_greater_distance])
                else:
                    # Priority 3: Among neighbors with the same mark_class but different nearest sample, choose the one with the shortest distance to a sample
                    same_mark_class_different_sample = [
                        (npos, self.blocks[npos]) for npos in neighbors
                        if npos in self.blocks and
                        self.blocks[npos].mark_class == ant.mark_class and
                        self.blocks[npos].nearest_sample_value != current_block.nearest_sample_value and
                        self.blocks[npos].domain == ant.domain
                    ]
                    if same_mark_class_different_sample:
                        next_pos = min(same_mark_class_different_sample, key=lambda x: x[1].distance_to_sample)[0]
                    else:
                        # Priority 4: Among neighbors, choose the one with the value closer to the ant's origin sample
                        closest_to_origin_sample = [
                            (npos, self.blocks[npos]) for npos in neighbors
                            if npos in self.blocks and self.blocks[npos].domain == ant.domain
                        ]
                        if closest_to_origin_sample:
                            next_pos = min(
                                closest_to_origin_sample,
                                key=lambda x: abs(x[1].value - current_block.nearest_sample_value)
                            )[0]
                        else:
                            # Priority 5: Random valid block
                            valid_neighbors = [npos for npos in neighbors if npos in self.blocks and self.blocks[npos].domain == ant.domain]
                            if valid_neighbors:
                                next_pos = random.choice(valid_neighbors)
                            else:
                                continue  # No valid neighbors, skip this ant's movement

            # Move the ant to the chosen position
            if next_pos not in self.blocks:
                if next_pos in self.allowed_positions:
                    domain = "default"
                    if hasattr(self, 'domain_mapping'):
                        domain = self.domain_mapping.get(next_pos, "Undomained")
                        if domain != ant.domain:
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
                    self.next_block_id += 1
                    self.blocks[next_pos].value = self.calculate_block_value(next_pos, dims)
                    self.blocks[next_pos].mark_class = self.get_mark_class(self.blocks[next_pos].value)  # Update mark_class
                    changes_made = True
                else:
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
                    next_block.value = self.calculate_block_value(next_pos, dims)
                    next_block.mark_class = self.get_mark_class(next_block.value)  # Update mark_class
                    next_block.visited = True
                    next_block.visit_count += 1
                    next_block.ant_count += 1
                    changes_made = True
            self.blocks[ant.current_pos].ant_count = max(0, self.blocks[ant.current_pos].ant_count - 1)
            ant.current_pos = next_pos
            ant.steps += 1
        return changes_made
    
    def update_block_value(self, pos: Tuple[int, int, int], dims: Tuple[int, int, int]):
        block = self.blocks[pos]
        if block.is_sample:  # Never update sample blocks
            return
            
        neighbors = self.get_neighbors(pos, dims)
        weights = []
        values = []
        
        # Collect neighbors info
        sample_neighbors = []
        visited_neighbors = []
        
        for npos in neighbors:
            nblock = self.blocks.get(npos)
            if nblock is None:
                continue
                
            if nblock.is_sample:
                weights.append(self.max_pheromone * 3.0)
                values.append(nblock.value)
                sample_neighbors.append(npos)
            elif nblock.visited:
                weights.append(1.0)
                values.append(nblock.value)
                visited_neighbors.append(npos)
        
        if values:  # Update only if we have valid values
            try:
                weights = np.array(weights, dtype=float)
                if weights.sum() > 0:
                    block.value = np.average(values, weights=weights)
                    block.mark_class = self.get_mark_class(block.value)  # Update mark_class
                    
                    # Enhanced debug output with filter status
                    if self.verbose:
                        print(f"Block update at {pos}:")
                        print(f"  Type: {'Sample' if block.is_sample else 'Regular'}")
                        print(f"  Value: {block.value:.2f}")
                        print(f"  Mark class: {block.mark_class}")
                        print(f"  Visit count: {block.visit_count}")
                        print(f"  Sample neighbors: {len(sample_neighbors)} {sample_neighbors}")
                        print(f"  Visited neighbors: {len(visited_neighbors)}")
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
            if not block.is_sample:  # Never modify sample blocks
                block.pheromone = max(0, block.pheromone - 1)
                if block.pheromone <= 0:
                    blocks_to_remove.append(pos)
        
        # Only remove non-sample blocks
        for pos in blocks_to_remove:
            if pos in self.blocks and not self.blocks[pos].is_sample:
                del self.blocks[pos]
                
        # Update ant positions if their block was removed
        for ant in self.ants:
            if ant.current_pos not in self.blocks:
                ant.current_pos = ant.origin_block
                ant.steps = 0
    
    def get_blocks_data(self) -> Dict[Tuple[int, int, int], float]:
        # Include all blocks that have been visited plus ALL sample blocks
        return {pos: block.value 
                for pos, block in self.blocks.items() 
                if block.is_sample or block.visited}
