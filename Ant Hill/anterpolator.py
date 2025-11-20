import csv
import pandas as pd
import pickle
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from scipy.spatial import cKDTree
import random
#from IPython.display import display
import os
#from google.colab import files
from tqdm import tqdm
#from tqdm.notebook import tqdm
import gc
import copy
from math import copysign
from itertools import product
from decimal import Decimal, ROUND_HALF_UP, getcontext

# Define the offsets for the 26 neighboring cells
OFFSETS = [(i, j, k) for i in [-1, 0, 1] for j in [-1, 0, 1] for k in [-1, 0, 1] if (i, j, k) != (0, 0, 0)]

def round_half_up(n): return int(n + 0.5 * np.sign(n))

def sign(a): return copysign(1,a)

def detect_delimiter(filename, num_chars=1024, default =","):
    """
    Detect the delimiter used in a CSV file.

    Args:
    - filename (str): Path to the CSV file.
    - num_chars (int): Number of characters to read from the file for detection.

    Returns:
    - str: Detected delimiter.
    """
    with open(filename, 'r', encoding='utf-8', errors='ignore') as file:
        sample = file.read(num_chars)
        sniffer = csv.Sniffer()
        try:
            delimiter = sniffer.sniff(sample).delimiter
        except csv.Error:
            # If the sniffer fails (e.g., due to lack of clear delimiting),
            # you could set a default delimiter or re-raise the exception.
            delimiter = default  # default to comma if detection fails
    return delimiter


def get_mark_of_class(value, class_size,precision='0.01'):
    if value is None:
        return None
    # Set the precision high enough to handle your calculations
    getcontext().prec = 28
    # Convert values to Decimal using strings to avoid floating-point issues
    value_dec = Decimal(str(value))
    class_size_dec = Decimal(str(class_size))
    # Perform the calculation using Decimal arithmetic
    division = value_dec / class_size_dec
    rounded = division.to_integral_value(rounding=ROUND_HALF_UP)
    result = rounded * class_size_dec
    # Round the result to 2 decimal places
    result = result.quantize(Decimal(precision), rounding=ROUND_HALF_UP)
    return float(result)

class Grid:
    def __init__(self, grid_dimension, lower_bounds, cell_size, samples_df,
                 default_domain='Unassigned',class_size=None,max_pheromone = None):
        self.grid_dimension = grid_dimension
        self.cell_size = cell_size
        self.lower_bounds = lower_bounds
        self.cells = np.empty(grid_dimension, dtype=object)
        self.total_cells = self.grid_dimension[0] * self.grid_dimension[1] * self.grid_dimension[2]

        self.domained = False
        self.default_domain = default_domain
        self.samples_per_domain = {}
        self.class_size=class_size

        empty_cell_template = GridCell(grid=self,position=None, samples=pd.DataFrame(), size=self.cell_size, class_size=self.class_size, ant_capacity=100, domain=self.default_domain)
        for i in tqdm(range(grid_dimension[0]), desc='Initializing empty grid...'):
            for j in range(grid_dimension[1]):
                for k in range(grid_dimension[2]):
                    centroid = tuple(self.lower_bounds[d] + (idx + 0.5) * self.cell_size[d] for d, idx in enumerate([i, j, k]))
                    new_cell = copy.deepcopy(empty_cell_template)
                    new_cell.position = centroid
                    new_cell.indices = (i, j, k)
                    self.cells[i, j, k] = new_cell
        self.samples_df = samples_df

        self.out_of_grid_cells = set()

        self.max_pheromone = max_pheromone #max(grid_dimension[i] for i in range(3))

        # Build a KDTree from the sample points
        self.kdtree = cKDTree(self.samples_df[['x', 'y', 'z']].values)

    def initialize_regular_grid(self):
        assigned_samples = set()
        self.domained = True  # Indicate that the grid uses domains
        self.samples_per_domain[self.default_domain] = pd.DataFrame(columns=self.samples_df.columns)

        for index, row in tqdm(self.samples_df.iterrows(), desc='Assigning samples to cells...'):
            # Skip this sample if it has already been assigned to a cell
            if index in assigned_samples:
                continue
            # Get the cell indices from the 'grid_indices' column
            i, j, k = row['grid_indices']
            centroid = tuple(self.lower_bounds[d] + (idx + 0.5) * self.cell_size[d] for d, idx in enumerate([i, j, k]))
            # Query the KDTree for samples within the cell:
            indices = self.kdtree.query_ball_point(centroid, r=np.sqrt(3)/2*self.cell_size[0])
            # Filter out samples that have already been assigned to cells
            indices = [idx for idx in indices if idx not in assigned_samples]
            samples = self.samples_df.iloc[indices].copy()
            samples['Domain'] = self.default_domain  # Assign the default domain to samples

            # Collect samples under the default domain
            self.samples_per_domain[self.default_domain] = pd.concat([self.samples_per_domain[self.default_domain], samples], ignore_index=True)

            # Assign samples to cell
            cell = GridCell(grid=self, position=centroid, samples=samples, size=self.cell_size, max_pheromone=self.max_pheromone,
                            domain=self.default_domain, class_size=self.class_size)
            self.cells[i, j, k] = cell
            # Add the indices of the assigned samples to the set
            assigned_samples.update(indices)
        print(f"Grid initialization complete. {len(assigned_samples)} samples assigned")

    def initialize_irregular_grid(self, cells_df):
        assigned_samples = set()
        self.domained = True
        # Create a dictionary to hold samples per domain
        #self.samples_per_domain = {}
        for index, row in tqdm(cells_df.iterrows(), total=cells_df.shape[0], desc=f'Processing cells'):
            centroid = (row['x'], row['y'], row['z'])
            try:
                domain = row['Domain']
                self.domained = True
            except KeyError:
                print(f"Error: no domain for cell {index}")
                domain = None
                pass
            # Calculate the cell indices for the given position
            i, j, k = tuple(min(int((pos - self.lower_bounds[d]) / self.cell_size[d]), self.grid_dimension[d] - 1) for d, pos in enumerate(centroid))
            # Query the KDTree for samples within the cell
            indices = self.kdtree.query_ball_point(centroid, r=np.sqrt(3)/2*self.cell_size[0])
            # Filter out samples that have already been assigned to cells
            indices = [index for index in indices if index not in assigned_samples]
            samples = self.samples_df.iloc[indices]

            # Assign domain to samples
            samples = samples.copy()
            samples['Domain'] = domain
            # Collect samples per domain
            if domain not in self.samples_per_domain:
                self.samples_per_domain[domain] = samples
            else:
                self.samples_per_domain[domain] = pd.concat([self.samples_per_domain[domain], samples])

            # Assign samples to cell
            self.cells[i, j, k] = GridCell(centroid, samples, size=self.cell_size, max_pheromone=self.max_pheromone, domain=domain)
            #self.cells[i, j, k].domain = domain
            assigned_samples.update(indices)
        print(f"Grid initialization complete. {len(assigned_samples)} samples assigned")

    def get_cell(self, position):
        #Return a cell for a given position
        indices = self.get_indices(position)

        # Return the cell at the calculated indices
        return self.cells[indices]

    def get_indices(self,position):
        # Return the cell indices for a given position
        return tuple(min(int((pos - self.lower_bounds[d]) / self.cell_size[d]), self.grid_dimension[d] - 1) for d, pos in enumerate(position))

    def compute_domain_averages(self, printdata=True):
        self.domain_averages = {}
        if printdata: print("\nDomain Averages:")
        for domain, samples in tqdm(self.samples_per_domain.items(), desc='Computing domain averages'):
            if not samples.empty:
                average_value = samples['Value'].mean()
                self.domain_averages[domain] = average_value
                if printdata: print(f"Domain {domain}: Average Value = {round(average_value,2)}")

    def get_average_of_neighbors(self, cell, attribute_name="value"):
        x, y, z = self.get_indices(cell.position)
        total_value = 0
        count = 0

        for offset in OFFSETS:
            nx, ny, nz = x + offset[0], y + offset[1], z + offset[2]

            # Ensure that the neighboring cell is within the grid bounds
            if 0 <= nx < self.grid_dimension[0] and 0 <= ny < self.grid_dimension[1] and 0 <= nz < self.grid_dimension[2]:
                neighbor = self.cells[nx, ny, nz]

                # Use getattr to retrieve the attribute specified by attribute_name
                neighbor_value = getattr(neighbor, attribute_name, None)

                # Check if the neighboring cell has the attribute
                if neighbor is not None and neighbor_value is not None:
                    total_value += neighbor_value
                    count += 1

        if count == 0:
            return None  # No neighbors with values
        else:
            average_value = total_value / count  # Average of the values of the neighboring cells
            return round(average_value, 2)  # Return the average rounded to two decimal places


    def fill_unvisited_cells(self, fillvalue = None):

        attribs = ["value","tracked_value"]

        # Loop through each cell in the grid.
        for x, y, z in tqdm(product(range(self.grid_dimension[0]),
                                    range(self.grid_dimension[1]),
                                    range(self.grid_dimension[2])),
                            total=self.total_cells, desc="Filling cells"):

                    cell = self.cells[x, y, z]

                    for attrib in attribs:

                        # If the cell has no value (i.e., it hasn't been visited)
                        if getattr(cell,attrib,None) is None:
                            if fillvalue is None:
                                #Fills unvisited cells with the average value of their neighbors.
                                setattr(cell,attrib,self.get_average_of_neighbors(cell,attrib))
                            else:
                                #Fills unvisited cells with fill value (for instance, a background value)
                                setattr(cell,attrib,fillvalue)

class GridCell:
    def __init__(self, grid=None, position=None, samples=None, size = (10,10,10),
                 ant_capacity=100, class_size=None, max_pheromone=None,domain=None):
        #self.grid = grid
        self.position = position
        self.indices = grid.get_indices(self.position) if self.position is not None else None
        self.size = size
        self.samples = samples #Why was this line commented?
        self.best_prediction = None #lesser difference between prediction and cell value
        self.best_coefficients = None
        self.best_direction=None
        self.best_xp=None

        self.domain = domain
        self.has_samples = not samples.empty

        self.value = round(samples['Value'].mean(),2) if self.has_samples else None

        self.num_ants=0
        self.ant_capacity=ant_capacity
        self.visits = 0

        #self.tracked_value = get_mark_of_class(self.value, class_size) if not self.value is None else None
        self.tracked_samples = {}  # Dictionary to hold pheromone levels per class
        if self.has_samples:
            tracked_sample = get_mark_of_class(self.value, class_size)
            sample_signature = self.indices #This should never change
            pheromone = max_pheromone
            sample_distance = 0 #Distance to nearest sample
            closest_sample = self.value #Closest sample
            mode = "Nest Trace" #Every cell with sample is initially a nest cell
            self.tracked_samples[tracked_sample]=Trail(signature=sample_signature,pheromone=pheromone,
                                                       closest_sample=closest_sample,sample_distance=sample_distance,mode=mode)

        self.iteration= None #Iteration when the cell was last visited

    def update_best_prediction(self, prediction, ant):

        coefficients = ant.coefficients

        # Calculate the absolute difference between the prediction and the actual value
        difference = abs(prediction - self.value)

        # If the prediction is closer to the actual value than the best prediction so far, update the best prediction and coefficients
        if self.best_prediction is None or difference < self.best_prediction:
            self.best_prediction = difference
            self.best_coefficients = coefficients
            self.best_direction = ant.direction
            ant.xp+=1
            self.best_xp=ant.xp
        else:
            ant.xp= max(ant.xp-1,0) #ant looses experience

    def get_infimum_supremum_pheromone_neighbor(self, tracked_value):
        """
        Finds neighboring cells with:
        - The maximum pheromone value less than or equal to the current cell's pheromone value (infimum).
        - The minimum pheromone value greater than or equal to the current cell's pheromone value (supremum).
        Returns a tuple: (infimum_neighbor_cell, supremum_neighbor_cell)
        """
        x, y, z = self.indices
        current_trail = self.tracked_samples.get(tracked_value, None)
        if current_trail is None:
            return (None, None)  # No pheromone for this tracked_value in current cell

        current_pheromone = current_trail.pheromone
        max_lower_pheromone = None
        max_lower_cell = None
        min_upper_pheromone = None
        min_upper_cell = None

        for offset in OFFSETS:
            nx, ny, nz = x + offset[0], y + offset[1], z + offset[2]

            # Check if neighbor indices are within grid bounds Warning: grid is being used despite being outside the class scope
            if (0 <= nx < grid.grid_dimension[0] and
                0 <= ny < grid.grid_dimension[1] and
                0 <= nz < grid.grid_dimension[2]):

                neighbor = grid.cells[nx, ny, nz]
                if neighbor is not None:
                    neighbor_trail = neighbor.tracked_samples.get(tracked_value, None)
                    if neighbor_trail is not None:
                        neighbor_pheromone = neighbor_trail.pheromone
                        # Check for infimum
                        if neighbor_pheromone <= current_pheromone:
                            if (max_lower_pheromone is None or neighbor_pheromone > max_lower_pheromone):
                                max_lower_pheromone = neighbor_pheromone
                                max_lower_cell = neighbor
                        # Check for supremum
                        if neighbor_pheromone >= current_pheromone:
                            if (min_upper_pheromone is None or neighbor_pheromone < min_upper_pheromone):
                                min_upper_pheromone = neighbor_pheromone
                                min_upper_cell = neighbor

        return (max_lower_cell, min_upper_cell)

class Trail:
    def __init__(self, signature=None, pheromone=None, closest_sample=None, sample_distance=None,mode=None):
        self.signature = signature
        self.pheromone = pheromone
        self.closest_sample = closest_sample
        self.sample_distance = sample_distance
        self.mode = mode

class Ant:
    def __init__(self, grid=None, position=None,
                 cell_size=None, lower_bounds=None, upper_bounds=None,
                 domain=None,
                 coefficients=None, tolerance=None, class_size=None,
                 index=None, life=0,
                 self_inclusion=0.1, decay_rate=1, max_pheromone=None,
                 average_distance=10, background_distance=20,
                 none_is_value=False, background = 0,weight_method='nearest_sample'):
        self.index=index
        self.life = life
        self.position = position
        self.last_position = position
        self.coefficients = coefficients
        self.domain = domain
        self.tolerance = tolerance
        self.average_distance = average_distance
        self.background_distance = background_distance
        self.none_is_value = none_is_value
        self.background = background
        self.class_size = class_size
        self.decay_rate = decay_rate
        self.self_inclusion = self_inclusion #probability of a cell being included in the prediction of it's own value. 1: always, 0: never
        self.weight_method = weight_method  # Store the weighting method
        self.cell_size = cell_size
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.direction = (0,0,0) #relative direction the ant is moving
        self.xp=0 #experience points

        current_cell = grid.get_cell(self.position)
        if current_cell.has_samples:
            self.start_sample = current_cell.value
            self.start_position = self.position
            self.sample_signature=grid.get_indices(self.start_position)
            self.tracked_sample=get_mark_of_class(current_cell.value, self.class_size) #This should never change
            self.closest_sample=current_cell.value
            self.sample_distance=0
        else:
            print(f"Error: Ant {self.index} born on a cell with no samples.")
         #if current_cell.value is not None else None

        self.max_pheromone=max_pheromone

       #Ant behaviors:
       #Search: search for a sample with the same tracked value as the ant home
       #Return: follow back the track after a same tracked value sample was found

        self.behavior = "No Food"

        current_cell.num_ants +=1
        current_cell.visits +=1

    def move(self, grid):
        # Create a list to store the positions of empty cells and cells that are out of the ant's prediction tolerance
        remaining_offsets= set(OFFSETS)
        second_choice_cells = set()
        neighbor_out_of_grid_cells = set()
        out_of_domain_cells = set()
        occupied_cells = set()
        follow_leader= False

        current_cell = grid.get_cell(self.position)
        last_cell=grid.get_cell(self.last_position)

        #Cells with Supremum and Infimum pheromome value for ant tracked_value are obtained
        infimum_cell, supremum_cell = current_cell.get_infimum_supremum_pheromone_neighbor(self.tracked_sample)
        infimum_offset = tuple(infimum_cell.indices[i] - current_cell.indices[i] for i in range(3)) if infimum_cell is not None else None
        supremum_offset = tuple(supremum_cell.indices[i] - current_cell.indices[i] for i in range(3)) if supremum_cell is not None else None
        inf_sup = infimum_offset in remaining_offsets or supremum_offset in remaining_offsets

        chkcurrent_ontrack = self.tracked_sample in current_cell.tracked_samples and current_cell.tracked_samples[self.tracked_sample].pheromone > 0
        chkleader_exist = current_cell.best_direction is not None and current_cell.best_xp > self.xp

        #while True:
        while remaining_offsets:

            # Choose a random offset from the remaining ones and remove it from the set
            if follow_leader and chkleader_exist:
                #the ant changes direction towards more experienced ant
                cLeader=current_cell.best_xp/(current_cell.best_xp+self.xp)
                cFollower = 1 - cLeader #self.xp/(current_cell.best_xp+self.xp)
                #The most difference of experience between leader and follower, the closest the new direction to that of the leader
                offset = tuple(round_half_up(cFollower*self.direction[i]+cLeader*current_cell.best_direction[i]) for i in range(3))
                if offset == (0,0,0): offset = current_cell.best_direction
                follow_leader = False #to avoid endless while
                if any(abs(coordinate) > 1 for coordinate in offset) or all(coordinate == 0 for coordinate in offset):
                    print(f"New direction error: {round(cFollower,2)} * {self.direction} + {round(cLeader,2)} * {current_cell.best_direction} = {offset}")
            else:
                offset = random.choice(list(remaining_offsets))

            try:
                remaining_offsets.remove(offset)
                #print(f"Offset {offset} removed")
            except KeyError:
                        print(f"Invalid offset error: {offset} not in {remaining_offsets}")

            # Calculate the position of neighbor cell
            #if offset in {infimum_offset,supremum_offset}:print(f"offset: {offset}")
            neighbor_position = tuple(self.position[i] + offset[i]*self.cell_size[i] for i in range(3))
            if neighbor_position == self.position:
                print(f"Offset Error: neighbor cell  =  current Cell: {current_cell.position}. Offset: {offset}")

            # Validate the position of the neighboring cell
            chkneighbor_outofgrid = (neighbor_position in grid.out_of_grid_cells)
            chkneighbor_inbounds = all(self.lower_bounds[i] <= neighbor_position[i] < self.upper_bounds[i] for i in range(3))

            # If it's a valid position, the cell is obtained
            if chkneighbor_inbounds:
                neighbor_cell = grid.get_cell(neighbor_position) # Get the cell at the new position

                chkneighbor_indomain = (not grid.domained) or (neighbor_cell.domain == self.domain)
                chkneighbor_ontrack = self.tracked_sample in neighbor_cell.tracked_samples and neighbor_cell.tracked_samples[self.tracked_sample].pheromone > 0

                #neighbor_prediction = self.predict(grid,neighbor_cell)

                #update distance to nearest sample
                if chkneighbor_indomain:
                    if chkcurrent_ontrack:
                        self.update_sample_distance(current_cell, neighbor_cell)

                #Some conditions are calculated:
                #Conditions for prioritary movement based on infimum and maximum pheromone values:
                chk_move = False
                if (chkneighbor_ontrack and
                    self.tracked_sample in current_cell.tracked_samples):
                    if (self.behavior == "No Food" and
                        (current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        neighbor_cell.tracked_samples[self.tracked_sample].mode == "Food Trace")
                        ):
                        print(f"i{current_cell.iteration}:Searching for new sample, found food trace")
                        chk_move = True
                    elif (self.behavior == "No Food" and
                        (current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        neighbor_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        current_cell.tracked_samples[self.tracked_sample].signature != neighbor_cell.tracked_samples[self.tracked_sample].signature)
                        ):
                        print(f"i{current_cell.iteration}:Searching for new sample, found foreign nest trace")
                        chk_move = True
                    elif (self.behavior == "No Food" and
                        (current_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                        neighbor_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        self.sample_signature != neighbor_cell.tracked_samples[self.tracked_sample].signature and
                        neighbor_cell == supremum_cell)
                        ):
                        print(f"i{current_cell.iteration}:Following food trace, found foreign nest trace")
                        chk_move = True
                    elif (self.behavior == "No Food" and
                          (current_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                          neighbor_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                          neighbor_cell == supremum_cell)
                        ):
                          print(f"i{current_cell.iteration}:Following path to new sample, found food trace with stronger pheromone")
                          chk_move = True
                    elif (self.behavior == "Food" and
                          (current_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                          neighbor_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                          neighbor_cell.tracked_samples[self.tracked_sample].signature == self.sample_signature)
                          ):
                          print(f"i{current_cell.iteration}:Returning to home sample, found nest trace")
                          chk_move = True
                    elif (self.behavior == "Food" and
                          (current_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                          neighbor_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                          neighbor_cell == infimum_cell)
                        ):
                          print(f"i{current_cell.iteration}:Returning to home sample, following back food trace")
                          chk_move = True
                    elif (self.behavior == "Food" and
                          (current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        neighbor_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        neighbor_cell == supremum_cell)
                      ):
                        print(f"i{current_cell.iteration}:Returning to home sample, following back nest trace")
                        chk_move = True
                    elif (self.behavior == "No Food" and False and
                        (current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                        neighbor_cell.tracked_samples[self.tracked_sample].mode is None)
                        ):
                        print(f"i{current_cell.iteration}:Searching for new sample, found unvisited cell")
                        chk_move = True


                #<COMPLETAR OTRAS ALTERNATIVAS>

                #Cell status
                chkneighbor_nonecell = neighbor_cell is None
                chkneighbor_fullcell = neighbor_cell.num_ants >= neighbor_cell.ant_capacity
                chkneighbor_nonevalue = not chkneighbor_nonecell and neighbor_cell.value is None

                #Value
                chkneighbor_closevalue = not chkneighbor_nonevalue and abs(current_cell.value - neighbor_cell.value) < (self.class_size * self.tolerance)
                chkneighbor_closervalue = not chkneighbor_nonevalue and abs(current_cell.value - neighbor_cell.value) < abs(current_cell.value - last_cell.value)
                chkneighbor_closetolast = not chkneighbor_nonevalue and abs(last_cell.value - neighbor_cell.value) < (self.class_size * self.tolerance)

                #Sample value
                chkneighbor_hassamples = neighbor_cell.has_samples
                chksample2sample = current_cell.has_samples and neighbor_cell.has_samples
                chkneighbor_closesample = chkneighbor_indomain and chkneighbor_ontrack and (neighbor_cell.tracked_samples[self.tracked_sample].sample_distance <= self.background_distance)
                chkneighbor_closertosample = chkneighbor_indomain and chkneighbor_ontrack and chkcurrent_ontrack and (current_cell.tracked_samples[self.tracked_sample].sample_distance > neighbor_cell.tracked_samples[self.tracked_sample].sample_distance)
                chkneighbor_furthertosample = chkneighbor_indomain and chkneighbor_ontrack and chkcurrent_ontrack and  (current_cell.tracked_samples[self.tracked_sample].sample_distance <= neighbor_cell.tracked_samples[self.tracked_sample].sample_distance)
                chkneighbor_worsesample = chkneighbor_indomain and chkneighbor_ontrack and chkcurrent_ontrack and (current_cell.tracked_samples[self.tracked_sample].closest_sample > neighbor_cell.tracked_samples[self.tracked_sample].closest_sample)
                chkneighbor_hassamplesignature = chkneighbor_indomain and chkneighbor_ontrack and not neighbor_cell.tracked_samples[self.tracked_sample].closest_sample is None

                #Tracked value
                neighbor_trackedvalue = get_mark_of_class(neighbor_cell.value, self.class_size) if not(chkneighbor_nonevalue) else None
                neighbor_trackedsample = get_mark_of_class(neighbor_cell.tracked_samples[self.tracked_sample].closest_sample, self.class_size) if chkneighbor_hassamplesignature else None
                chkneighbor_closetracked = not chkneighbor_nonevalue and self.tracked_sample is not None and abs(self.tracked_sample - neighbor_cell.value) < (self.class_size * self.tolerance)
                chkneighbor_closertracked = not chkneighbor_nonevalue and self.tracked_sample is not None and abs(self.tracked_sample - neighbor_cell.value) < abs(self.tracked_sample - current_cell.value)
                chklast_ontrack = self.tracked_sample in last_cell.tracked_samples and last_cell.tracked_samples[self.tracked_sample].pheromone > 0


                #Pheromone
                chkneighbor_nonepheromone = not chkneighbor_nonecell and (self.tracked_sample not in neighbor_cell.tracked_samples or neighbor_cell.tracked_samples[self.tracked_sample] is None)
                chkneighbor_strongerpheromone = chkneighbor_ontrack and (not chkcurrent_ontrack or neighbor_cell.tracked_samples[self.tracked_sample].pheromone > current_cell.tracked_samples[self.tracked_sample].pheromone)
                chkneighbor_weakerpheromone = chkneighbor_ontrack and (chkcurrent_ontrack and neighbor_cell.tracked_samples[self.tracked_sample].pheromone <= current_cell.tracked_samples[self.tracked_sample].pheromone)
                chkneighbor_strongerpath = chkneighbor_closetracked and (chkneighbor_nonepheromone or chkneighbor_strongerpheromone)
                chkneighbor_searchpath = self.behavior == "Search" and chkneighbor_strongerpheromone
                chkneighbor_returnpath = self.behavior == "Return" and chkneighbor_weakerpheromone

                chkneighbor_sametrack = chkneighbor_indomain and (neighbor_trackedsample == self.tracked_sample)
                chkneighbor_distinctsample = chkneighbor_ontrack and (neighbor_cell.tracked_samples[self.tracked_sample].signature != self.sample_signature)
                chkneighbor_newpath = chkneighbor_ontrack and chkneighbor_distinctsample
                chkneighbbor_samepath = chkneighbor_ontrack and not chkneighbor_distinctsample
                chkneighbor_samesample = chkneighbor_ontrack and (neighbor_cell.tracked_samples[self.tracked_sample].signature == self.sample_signature)

                #Misc
                chkneighbor_lessvisited = current_cell.visits > neighbor_cell.visits
                chkneighbor_farthertosample = chkneighbor_indomain and chkneighbor_ontrack and chkcurrent_ontrack and (current_cell.tracked_samples[self.tracked_sample].sample_distance > neighbor_cell.tracked_samples[self.tracked_sample].sample_distance)
                chkneighbor_sameaslast = (neighbor_position == self.last_position)
                chkneighbor_samedirection = (offset == self.direction)
                chkneighbor_sametrend = chkneighbor_samedirection and chkneighbor_closetolast
                chkneighbor_stray = (chkneighbor_indomain and chkneighbor_ontrack and not chkneighbor_nonevalue and
                                     get_mark_of_class(neighbor_cell.value, self.class_size) != get_mark_of_class(neighbor_cell.tracked_samples[self.tracked_sample].closest_sample,self.class_size) and
                                     get_mark_of_class(neighbor_cell.tracked_samples[self.tracked_sample].closest_sample, self.class_size) == self.tracked_sample
                                    )


                chkneighbor_pathfound = not chkcurrent_ontrack and chkneighbor_ontrack
                chkneighbor_closeweak = (chkneighbor_sametrack and chkneighbor_distinctsample
                                         and chkneighbor_closertosample and chkneighbor_weakerpheromone)
                chkneighbor_farstrong = (chkneighbor_sametrack and chkneighbor_distinctsample
                                         and chkneighbor_furthertosample and chkneighbor_strongerpheromone)
                chkneighbor_betterpath = (chk_move
                                         )

                #Conditions for a neighbor cell to be considered of first priority for moving
                chkneighbor_firstchoice = (not chkneighbor_sameaslast and
                                           (not chkneighbor_hassamples or True) and
                                              (
                                              (chkneighbor_nonevalue and False) or
                                              (chkneighbor_newpath and False) or
                                              (chkneighbor_betterpath and True) or
                                              (chkneighbor_stray and False) or
                                              (chkneighbor_sametrend and False) or
                                              (chkneighbor_lessvisited and False) or
                                              (chkneighbor_worsesample and False) or
                                              (chkneighbor_farthertosample and False) or
                                              (chkneighbor_closevalue and False) or
                                              (chkneighbor_closervalue and False) or
                                              (chkneighbor_closetolast and False) or
                                              (chkneighbor_closetracked and False) or
                                              (chkneighbor_closertracked and False) or
                                              (chkneighbor_closesample and False)
                                              )
                                            )
            else:
                if not chkneighbor_outofgrid:
                    grid.out_of_grid_cells.add(tuple(neighbor_position))
                    chkneighbor_outofgrid = True
                neighbor_out_of_grid_cells.add(tuple(neighbor_position))

            # Check if the new position is a known out of grid cell:
            if chkneighbor_outofgrid or chkneighbor_nonecell:
                grid.out_of_grid_cells.add(tuple(neighbor_position))
                neighbor_out_of_grid_cells.add(tuple(neighbor_position))

            # Check if the new position is inside the grid and domain
            elif chkneighbor_inbounds:

                # If the neighbor cell's domain doesn't match the ant domain
                if not chkneighbor_indomain:
                    out_of_domain_cells.add(tuple(neighbor_position))

                # if the neighbor cell is full
                elif chkneighbor_fullcell:
                    occupied_cells.add(tuple(neighbor_position))

                # If the neighbor cell is eligible
                elif chkneighbor_firstchoice:
                    self.update_position(current_cell,neighbor_cell, offset)
                    return

                # If the neighbor cell's value is not close to the ant's prediction or the current cell, add it to the list of second choice cells
                else:
                    second_choice_cells.add(tuple(neighbor_position))
            else:
                grid.out_of_grid_cells.add(tuple(neighbor_position))
                neighbor_out_of_grid_cells.add(tuple(neighbor_position))

            # If all neighboring cells have been checked, move to a random available cell
            if len(second_choice_cells) + len(neighbor_out_of_grid_cells) + len(out_of_domain_cells) +len(occupied_cells) == len(OFFSETS):
                if second_choice_cells:  # Check if the list is not empty
                    print(f"i{current_cell.iteration}: Ant moved to second choice cell")
                    neighbor_position = random.choice(list(second_choice_cells))
                    neighbor_cell = grid.get_cell(neighbor_position)
                    new_direction = tuple(round((neighbor_cell.position[i] - current_cell.position[i])/cell_size[i],0) for i in range(3))
                    self.update_position(current_cell,neighbor_cell, new_direction)
                else:
                    print(f"No second choice cells available, ant won't move: {self.position}. Current Cell: {current_cell.position}. Domain: {self.domain}")
                return

    def update_position(self, current_cell, new_cell, new_direction=(0,0,0)):

        lock_cell = False

        if current_cell.position != self.position:
            print(f"Position Error. Ant: {self.position}. Current Cell: {current_cell.position}")
        self.last_position = self.position
        self.position = new_cell.position #neighbor_position
        self.direction = new_direction #offset

        #Update number of ants after moving to a new cell
        current_cell.num_ants-=1
        new_cell.num_ants+=1

        #lock cell after leaving
        if lock_cell:
            if current_cell.num_ants == 0:
                current_cell.has_samples = True

        #increase visits to new cell
        if new_cell != current_cell:
            new_cell.visits+=1

    def predict(self, grid, cell,self_inclusion=0, use_coefficient=False):

        prediction = 0.0
        total_weight = 0.0

        if self.weight_method == 'pheromone_strenght':
            for tracked_sample in cell.tracked_samples.keys():
                #print(f"tv: {tracked_sample}. ph: {cell.tracked_samples[tracked_sample]}")
                weight = cell.tracked_samples[tracked_sample]**self.decay_rate
                prediction += tracked_sample * weight
                total_weight += weight
            # Normalize the prediction by the total weight
            if total_weight != 0:
                prediction /= total_weight
            else:
                #print(f"zero-sum for total weight of coefficients. Cell Value: {cell.value}. Samples: {cell.has_samples}. Tracked Value: {self.tracked_value}. Pheromone: {cell.tracked_samples[self.tracked_value]}")
                prediction = self.background
            return round(prediction,2)


        # Get the indices of the current cell
        i, j, k = grid.get_indices(self.position)


        last_cell=grid.get_cell(self.last_position)

        # Initialize weighted parameters
        weighted = True  # If True, coefficients are weighted

        #Initialize the prediction to the current cell's value times the last coefficient
        #Including the current cell in the prediction of its own value makes the interpolator converge to 0 too quickly away from samples
        #Not including the current cell in the prediction makes the interpolator to return values too high away from samples
        #self_inclusion parameter aims to calibrate this two behaviors into something in the middle.
        #self_inclusion: probability of a cell being included in the prediction of it's own value
        if self.self_inclusion > random.random():
            try:
                if (cell.value is None and self.none_is_value):
                    if (last_cell.sample_distance >= self.background_distance):
                        cell.value = self.background
                    elif (last_cell.sample_distance >= self.average_distance):
                        cell.value = grid.domain_averages[self.domain]
                if cell.value is not None:
                    weight = self.get_weight(last_cell, cell, grid) if weighted else 1 #and not cell.has_samples else 1

                    # cell_coefficient = self.coefficients[-1] if use_coefficient else 1
                    # weighted_coefficient = weight * cell_coefficient
                    # prediction += weighted_coefficient * cell.value
                    prediction += weight * cell.value
                    total_weight += weight #weighted_coefficient
            except AttributeError:
                print(f"\nPrediction Error")
                print(f"Current cell indices: {(i,j,k)}. Current cell position: {cell if(cell is None) else cell.position}")
                print(f"Expected ant cell indices: {grid.get_indices(self.position)}. Current ant position: {self.position}")
                exit()

        # Loop over the offsets to find the max_dTracked across all neighbors

        nNoneNeighbors = 0
        # Loop over the offsets to get the neighboring cells
        for idx, offset in enumerate(OFFSETS):
            # Compute the indices of the neighboring cell
            neighbor_i = i + offset[0]
            neighbor_j = j + offset[1]
            neighbor_k = k + offset[2]

            # Check if the indices are within the grid
            if (0 <= neighbor_i < grid.grid_dimension[0] and
                0 <= neighbor_j < grid.grid_dimension[1] and
                0 <= neighbor_k < grid.grid_dimension[2]):

                # Get the neighboring cell
                neighbor_cell = grid.cells[neighbor_i, neighbor_j, neighbor_k]

                # Check if the neighboring cell is None
                if neighbor_cell is None or neighbor_cell.value is None:
                    nNoneNeighbors += 1

                if neighbor_cell is not None:

                    # Check is the cell is within the domain
                    if (not grid.domained) or (self.domain == neighbor_cell.domain):

                        # #update distance to nearest sample in the same domain (if any) (also in move method)
                        # self.update_sample_distance(cell, neighbor_cell)

                        # If the neghbor cell has a value, add the weighted value of the cell to the prediction and add the coefficient to the list of valid coefficients
                        if (neighbor_cell.value is None and self.none_is_value):
                            if (last_cell.tracked_samples[self.tracked_sample].sample_distance >= self.background_distance):
                                neighbor_cell.value = self.background
                            elif (last_cell.tracked_samples[self.tracked_sample].sample_distance >= self.average_distance):
                                neighbor_cell.value = grid.domain_averages[self.domain]
                        if neighbor_cell.value is not None:
                            weight = self.get_weight(last_cell, neighbor_cell, grid) if weighted else 1 #and not neighbor_cell.has_samples else 1

                            # cell_coefficient = self.coefficients[idx] if use_coefficient else 1
                            # weighted_coefficient = weight * cell_coefficient
                            # prediction += weighted_coefficient * neighbor_cell.value
                            prediction += weight * neighbor_cell.value

                            # Accumulate the weighted coefficient for normalization
                            total_weight += weight #weighted_coefficient

            else:
                out_position = tuple(grid.lower_bounds[d] + (idx + 0.5) * grid.cell_size[d] for d, idx in enumerate([neighbor_i, neighbor_j,  neighbor_k]))
                grid.out_of_grid_cells.add(out_position)

        # Normalize the prediction by the total weight
        if total_weight != 0:
            prediction /= total_weight
        else:
            print(f"zero-sum for total weight of coefficients. Samples: {cell.has_samples}. None neighbors: {nNoneNeighbors}")
            prediction = 0#self.background

        # Return the prediction
        return round(prediction,2)

    def get_weight(self, last_cell, cell, grid):
        if self.weight_method == 'nearest_sample':
            fValue = get_mark_of_class(abs(last_cell.closest_sample - cell.value), self.class_size)
            weight = self.class_size / (self.class_size + fValue ** self.decay_rate)
        elif self.weight_method == 'domain_average':
            # Get the average value for the ant's domain
            domain_avg = grid.domain_averages.get(self.domain, self.background)
            fValue = get_mark_of_class(abs(domain_avg - cell.value), self.class_size)
            weight = self.class_size / (self.class_size + fValue ** self.decay_rate)
        elif self.weight_method == 'greater_value':
            weight = (cell.value**self.decay_rate)
        elif self.weight_method == 'greater_sample':
            weight = (cell.tracked_samples[self.tracked_sample].closest_sample**self.decay_rate) if self.tracked_sample in cell.tracked_samples else 0
        elif self.weight_method == 'value_consistence':
            fValue = get_mark_of_class(abs(cell.value - cell.tracked_samples[self.tracked_sample].closest_sample),self.class_size)
            #weight = self.class_size / (fValue ** self.decay_rate) if fValue==0 else 1
            weight = self.class_size / (self.class_size + fValue ** self.decay_rate)
        elif self.weight_method == 'pheromone_strenght':
            weight = (cell.tracked_samples[cell.tracked_value]**self.decay_rate) if self.tracked_sample in cell.tracked_samples else 0
        else:
            # Default weight
            weight = 1
        return weight

    def update_coefficients(self, best_coefficients):
        # If the cell has no best coefficients yet, do nothing
        if best_coefficients is None:
            return

        # Calculate the new coefficients as the average of the ant's current coefficients and the best coefficients of the cell
        self.coefficients = [(self.coefficients[i] + best_coefficients[i]) / 2 for i in range(27)]

    def get_cells_distance(self,cell1,cell2):
        # Calculate the distance between the current cell and the neighbor cell
        dx = (cell1.position[0] - cell2.position[0])/cell1.size[0]
        dy = (cell1.position[1] - cell2.position[1])/cell1.size[1]
        dz = (cell1.position[2] - cell2.position[2])/cell1.size[2]

        return np.sqrt(dx**2 + dy**2 + dz**2)

    def update_sample_distance(self, cell, neighbor_cell):

        distance = self.get_cells_distance(cell,neighbor_cell)

        # Update neighbor_cell's sample_distance
        if (self.tracked_sample not in neighbor_cell.tracked_samples or
            neighbor_cell.tracked_samples[self.tracked_sample].sample_distance is None or
            neighbor_cell.tracked_samples[self.tracked_sample].sample_distance > (cell.tracked_samples[self.tracked_sample].sample_distance + distance)
           ):
            if cell.tracked_samples[self.tracked_sample].sample_distance is None:
                print(f"Error: Current cell has None distance for closest sample {self.tracked_sample}. Samples: {cell.has_samples}")
            #signature = self.sample_signature
            signature = cell.tracked_samples[self.tracked_sample].signature
            pheromone = max(cell.tracked_samples[self.tracked_sample].pheromone-1,0) #Is this correct?
            sample_distance = cell.tracked_samples[self.tracked_sample].sample_distance + distance
            closest_sample = cell.tracked_samples[self.tracked_sample].closest_sample
            neighbor_cell.tracked_samples[self.tracked_sample]=Trail(signature=signature,pheromone=pheromone,closest_sample=closest_sample,sample_distance=sample_distance)

        # Update cell's sample_distance
        if (cell.tracked_samples[self.tracked_sample].sample_distance is None or
            cell.tracked_samples[self.tracked_sample].sample_distance > (neighbor_cell.tracked_samples[self.tracked_sample].sample_distance + distance)):
            #signature = self.sample_signature
            signature = neighbor_cell.tracked_samples[self.tracked_sample].signature
            pheromone = max(neighbor_cell.tracked_samples[self.tracked_sample].pheromone-1,0) #Is this correct?
            sample_distance = neighbor_cell.tracked_samples[self.tracked_sample].sample_distance + distance
            closest_sample = neighbor_cell.tracked_samples[self.tracked_sample].closest_sample
            mode = "Nest Trace"
            cell.tracked_samples[self.tracked_sample]=Trail(signature=signature,pheromone=pheromone,closest_sample=closest_sample,
                                                            sample_distance=sample_distance, mode = mode)

        self.closest_sample = cell.tracked_samples[self.tracked_sample].closest_sample
        self.sample_distance = cell.tracked_samples[self.tracked_sample].sample_distance


    def update_trail(self,grid,verbose=False):

        #Current cell and last cell visited are retrieved, aong with its their tracked value and sample mark for messaging purposes.
        current_cell = grid.get_cell(self.position)
        last_cell = grid.get_cell(self.last_position)

        #Cell is visited for the first time or trace faded out:
        if self.tracked_sample not in current_cell.tracked_samples:
            print("Cell visited for the first time")
            signature = self.sample_signature
            if (not current_cell.has_samples
                or self.tracked_sample != get_mark_of_class(current_cell.value,self.class_size)
               ):
                pheromone = 0
            closest_sample = self.closest_sample
            sample_distance = self.sample_distance + self.get_cells_distance(last_cell, current_cell)
            if self.behavior == "No Food":
                mode = "Nest Trace"
            elif self.behavior == "Food":
                mode = "Food Trace"
            current_cell.tracked_samples[self.tracked_sample]=Trail(signature=signature,pheromone=pheromone,closest_sample=closest_sample,sample_distance=sample_distance,mode=mode)


        if self.tracked_sample not in last_cell.tracked_samples or last_cell.tracked_samples[self.tracked_sample].mode == None:
            print("Error: last cell has no trace mode")
            #last_cell.tracked_samples[self.tracked_sample].mode = "Nest Trace"
        if current_cell.tracked_samples[self.tracked_sample].mode == None:
            if self.behavior == "No Food":
                current_cell.tracked_samples[self.tracked_sample].mode = "Nest Trace"
        if self.behavior == "No Food":
            if (last_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                current_cell.tracked_samples[self.tracked_sample].mode == "Food Trace"):
                if last_cell.tracked_samples[self.tracked_sample].pheromone <= current_cell.tracked_samples[self.tracked_sample].pheromone:
                    last_cell.tracked_samples[self.tracked_sample].mode = "Food Trace"
                    last_cell.tracked_samples[self.tracked_sample].pheromone = max(0,current_cell.tracked_samples[self.tracked_sample].pheromone-1)
                    if not last_cell.has_samples:
                        last_cell.tracked_samples[self.tracked_sample].signature = current_cell.tracked_samples[self.tracked_sample].signature
            elif (last_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                  current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace"):
                  if (not last_cell.has_samples and
                      current_cell.tracked_samples[self.tracked_sample].signature != last_cell.tracked_samples[self.tracked_sample].signature):
                      last_cell.tracked_samples[self.tracked_sample].signature = current_cell.tracked_samples[self.tracked_sample].signature
                  if last_cell.tracked_samples[self.tracked_sample].pheromone < current_cell.tracked_samples[self.tracked_sample].pheromone:
                      last_cell.tracked_samples[self.tracked_sample].pheromone = max(0,current_cell.tracked_samples[self.tracked_sample].pheromone-1)
                  elif last_cell.tracked_samples[self.tracked_sample].pheromone > current_cell.tracked_samples[self.tracked_sample].pheromone:
                      current_cell.tracked_samples[self.tracked_sample].pheromone = max(0,last_cell.tracked_samples[self.tracked_sample].pheromone-1)
            if (current_cell.has_samples and
                current_cell.tracked_samples[self.tracked_sample].signature != self.sample_signature): #Ant found new sample
                self.behavior = "Food"
                current_cell.tracked_samples[self.tracked_sample].mode = "Food Trace"
        elif self.behavior == "Food":
            if (last_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace"):
                if last_cell.tracked_samples[self.tracked_sample].pheromone >= current_cell.tracked_samples[self.tracked_sample].pheromone:
                    current_cell.tracked_samples[self.tracked_sample].mode = "Food Trace"
                    current_cell.tracked_samples[self.tracked_sample].pheromone = max(0,last_cell.tracked_samples[self.tracked_sample].pheromone-1)
                    if not current_cell.has_samples:
                        current_cell.tracked_samples[self.tracked_sample].signature = last_cell.tracked_samples[self.tracked_sample].signature
            elif (last_cell.tracked_samples[self.tracked_sample].mode == "Food Trace" and
                  current_cell.tracked_samples[self.tracked_sample].mode == "Food Trace"):
                if (current_cell.tracked_samples[self.tracked_sample].pheromone < last_cell.tracked_samples[self.tracked_sample].pheromone):
                    current_cell.tracked_samples[self.tracked_sample].pheromone = max(0,last_cell.tracked_samples[self.tracked_sample].pheromone-1)
                    if not current_cell.has_samples:
                        current_cell.tracked_samples[self.tracked_sample].signature = last_cell.tracked_samples[self.tracked_sample].signature
            elif (last_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace" and
                  current_cell.tracked_samples[self.tracked_sample].mode == "Nest Trace"):
                  if last_cell.tracked_samples[self.tracked_sample].pheromone < current_cell.tracked_samples[self.tracked_sample].pheromone:
                      last_cell.tracked_samples[self.tracked_sample].pheromone = max(0,current_cell.tracked_samples[self.tracked_sample].pheromone-1)
            if (current_cell.has_samples and
                current_cell.tracked_samples[self.tracked_sample].signature == self.sample_signature): #Ant returned to home sample
                self.behavior = "No Food"
                #current_cell.tracked_samples[self.tracked_sample].mode = "Nest Trace"


class Interpolator:
    def __init__(self, ants, grid, background=0, fade_rate = 1,max_life=None):
        self.ants = ants
        self.max_life=max_life
        self.grid = grid
        self.background=background
        self.fade_rate = fade_rate

    def run(self, num_iterations,pinclusion = 0):

        iteration = 1
        for _ in tqdm(range(num_iterations),desc='Interpolating'):

            # Shuffle the ants
            random.shuffle(self.ants)

            for ant in self.ants:

                StartCell=self.grid.get_cell(ant.position)
                StartCell.iteration = iteration - 1

                #Ant moves to a new cell
                ant.move(self.grid)
                EndCell=self.grid.get_cell(ant.position)

                if EndCell is None:
                    print("\nError: Ant moved to None cell")
                    exit()

                #Ant predicts the value of its current cell
                cell = self.get_cell(ant.position)
                if EndCell!= cell:
                    print("\nError: Cell Discrepancy.")

                prediction = ant.predict(self.grid, cell, pinclusion)  # Pass the Grid object here

                #Ant updates the cell based on it's prediction
                if cell.has_samples:  # Check if the cell has samples
                    #ant life is replenished when she finds a new sample
                    ant.life = self.max_life
                    cell.update_best_prediction(prediction, ant)
                    ant.update_coefficients(cell.best_coefficients)
                else:
                    cell.value = max(prediction, self.background)

                if (cell.best_xp is None or (ant.xp > cell.best_xp)):
                    cell.best_xp = ant.xp
                    cell.best_direction = ant.direction

                #Trail is updated
                ant.update_trail(self.grid)

                #Ant gets weaker and dies
                ant.life-=1
                if ant.life == 0:
                    cell.num_ants -= 1
                    ant_dies=True
                    ant_respawns = not ant_dies
                    if ant_dies:
                        #Remove ant
                        self.ants.remove(ant)
                    elif ant_respawns:
                        #Respawn ant
                        ant.life = self.max_life
                        ant.position = ant.start_position
                        ant.last_position = ant.position
                        respawn_cell=self.get_cell(ant.position)
                        if respawn_cell.position != cell.position:
                            respawn_cell.num_ants+=1
                            respawn_cell.visits +=1

                cell.iteration = iteration

            self.evaporate_pheromone(fade_rate = self.fade_rate)
            iteration += 1

    def get_cell(self, position):
      # Calculate the cell indices for the given position
      indices = self.grid.get_indices(position)
      #indices = tuple(min(int((pos - lower_bounds[d]) / cell_size[d]), self.grid.grid_dimension[d] - 1) for d, pos in enumerate(position))

      # Return the cell at the calculated indices
      return self.grid.cells[indices]

    def evaporate_pheromone(self, fade_rate=1):
        for cell in self.grid.cells.flat:
            if cell and cell.tracked_samples and not cell.has_samples:
                # Create a new dictionary to hold updated pheromones
                updated_tracked_samples = {}
                for ts, trail in cell.tracked_samples.items():
                    # Decrease the pheromone level
                    trail.pheromone = max(trail.pheromone - fade_rate, 0)
                    # Keep the Trail object if pheromone > 0
                    if trail.pheromone >= 0: #Removed cells may be problematic
                        updated_tracked_samples[ts] = trail
                # Update the cell's pheromones
                cell.tracked_samples = updated_tracked_samples

#Input parameters
samplesfile_columns = ['x', 'y', 'z', 'Value'] #Column definition for samples file
#SamplesFile= "/content/Samples_CNN_ABA_BED.csv" #Samples filename
#SamplesFile= "/content/Samples_CNN 03.csv"
#SamplesFile= "/content/Leda_assay.csv"
#SamplesFile= "/content/Samples2D.csv"
#SamplesFile="/content/Sondeos_muestras.csv"
#SamplesFile="/content/Alati_Leyes.csv"
SamplesFile="/content/ANT-Samples.csv"
SamplesFile="/content/Cubo-Samples.csv"
SamplesFile="Samples_CNN_ABA_BED.csv"
samplesfile_delimiter = detect_delimiter(SamplesFile) #Inferred column delimiter
base_name, extension = os.path.splitext(SamplesFile)

cellsfile_columns = ['x', 'y', 'z','Domain'] #Column definition for domained cells
#CellsFile = "/content/Leda_DomainedCells.csv" #Domained cells filename. None: undomained
#CellsFile = "/content/Leda_5x5x5_DomainedCells.csv"
#CellsFile = "/content/Sondeos_muestras_5x5x5_DomainedCells.csv"
#CellsFile = "/content/Alati 20x20x20_Domained.csv"
CellsFile="/content/ANT-Domains.csv"
CellsFile = None
if CellsFile is not None:
    try:
        cellsfile_delimiter = detect_delimiter(CellsFile)
        base_name=base_name + "_Domained"
    except FileNotFoundError:
        print(f"Error: File {CellsFile} not found. Using undomained interpolation.")
        CellsFile = None
        #exit()

InitGridFile = f"{base_name}_InitGrid.pkl" #f"{base_name}_InitGrid.csv"

#Grid parameters
cell_size = (20,20,20) #block size
#cell_size = (1,1,1)
grid_rebuild = True #False: uses de grid generated in previous code execution
default_domain = 'Unassigned' if CellsFile is None else None

# Interpolation parameters
ants_per_sample = 20#nAnts = 3000  #number of ants
max_life=200
max_pheromone = 150
fade_rate = 1
#pheromone_factor = 10
iterations = 200 #number of times ants and cells visited are updated
weight_method = 'greater_sample'#'nearest_sample'#'pheromone_strenght'#'domain_average'#'greater_value'#'value_consistence'#
class_size = 1 #size of mark of class
tolerance = 1 #factor of class size
background=0 #background value
average_distance = 7 #Distance to samples where for blank cells is assumed domain's average value
background_distance = 10 #Distance to samples where for blank cells is assumed background value (out of the anomaly)
self_inclusion = 0 #probability of cell value considered in prediction. O: only neighbors are used for prediction
decay_rate = 2 #Influence of neighbor cell value on weighted prediction
none_is_value=True #True: background value is assumed for blank cells.
fillbackground=False
fillaverage = False #True: fill unvisited cells with average value of neighbors

#Output parameters
plotscatter=False #A plot is generated after interpolation ends
download = True #True: file with interpolated cells is downloaded automatically
export_pheromones = True

#Load samples
print(f"Loading samples from file {SamplesFile}. Delimiter: '{samplesfile_delimiter}'...")
samples_df = pd.read_csv(SamplesFile,delimiter=samplesfile_delimiter)
print(f'Samples loaded: {len(samples_df)}')
samples_df.to_csv("SamplesRaw.csv", index=False)

# Check if the required columns are present in the DataFrame
if not all(column in samples_df.columns for column in samplesfile_columns):
    # If not, try to rename the columns based on their positions
    samples_df = samples_df.iloc[:, :4].copy()
    samples_df.columns = samplesfile_columns
else:
    # If the required columns are already present, select them
    samples_df = samples_df[samplesfile_columns].copy()

# Print the first few rows to verify
print(samples_df[samplesfile_columns].head())

#Rows with invalid samples are removed
samples_df['Value'] = pd.to_numeric(samples_df['Value'], errors='coerce') #Converts column Value to numeric
samples_df = samples_df.dropna() #remove rows with NaN values
print(f'Valid (numeric) samples: {len(samples_df)}')
samples_df.to_csv("SamplesCheck.csv", index=False)

# Calculate the grid bounds
inputcells_df=None
if(CellsFile is None):

    lower_bounds = samples_df[['x', 'y', 'z']].min().values - 0.5 * np.array(cell_size)
    upper_bounds = samples_df[['x', 'y', 'z']].max().values + 0.5 * np.array(cell_size)

else:
    print(f"Loading cells definition from file {CellsFile}. Delimiter: '{cellsfile_delimiter}'...")
    inputcells_df = pd.read_csv(CellsFile,delimiter=cellsfile_delimiter)
    print(f"Cells loaded: {len(inputcells_df)}")
    # Convert non-numeric values in the 'x', 'y', and 'z' columns to NaN
    inputcells_df['x'] = pd.to_numeric(inputcells_df['x'], errors='coerce')
    inputcells_df['y'] = pd.to_numeric(inputcells_df['y'], errors='coerce')
    inputcells_df['z'] = pd.to_numeric(inputcells_df['z'], errors='coerce')

    # Remove rows with NaN values in the 'x', 'y', or 'z' columns
    inputcells_df = inputcells_df.dropna(subset=['x', 'y', 'z'])
    print(f"Cells with valid (numeric) coordinates: {len(inputcells_df)}")

    # Cell size is calculated
    lensubset=100
    subset_df = inputcells_df.head(lensubset)
    while subset_df['x'].nunique()<=1 or subset_df['y'].nunique()<=1 or subset_df['z'].nunique()<=1:
        lensubset+=lensubset
        subset_df = inputcells_df.head(lensubset)
    subset_df = subset_df.sort_values(['x', 'y', 'z'])

    differences = subset_df[['x', 'y', 'z']].diff().abs() # Calculate the difference between consecutive rows
    cell_size_series = differences[differences > 0].min() # The cell size in each dimension is the minimum non-zero difference
    cell_size = tuple(int(round(size)) for size in cell_size_series)  # Convert the Series to a tuple of integers

    #print(f"Cell Size: {cell_size}")

    lower_bounds = inputcells_df[['x', 'y', 'z']].min().values - 0.5 * np.array(cell_size)
    upper_bounds = inputcells_df[['x', 'y', 'z']].max().values + 0.5 * np.array(cell_size)

    #print(f"min Values: {inputcells_df[['x', 'y', 'z']].min().values}")
    #print(f"max Values: {inputcells_df[['x', 'y', 'z']].max().values}")

    # Print the first few rows to verify
    print(inputcells_df[cellsfile_columns].head())

print("Lower bounds: " + str(lower_bounds))
print("Upper bounds: " + str(upper_bounds))

# Calculate the grid size
grid_size = upper_bounds - lower_bounds

print(grid_size)

# Calculate the grid dimension
grid_dimension = [int(grid_size[i] / cell_size[i]) for i in range(3)]
total_cells = grid_dimension[0]*grid_dimension[1]*grid_dimension[2]

#Initialize the grid
print("Grid size: " + str(grid_size))
print("Grid dimension: " + str(grid_dimension))
print("Cells size: " + str(cell_size))
print(f"Total cells: {total_cells}")

print("Assigning cell indices to samples...")
samples_df['grid_indices'] = list(zip(
    (np.minimum((samples_df['x'] - lower_bounds[0]) / cell_size[0], grid_dimension[0] - 1)).astype(int),
    (np.minimum((samples_df['y'] - lower_bounds[1]) / cell_size[1], grid_dimension[1] - 1)).astype(int),
    (np.minimum((samples_df['z'] - lower_bounds[2]) / cell_size[2], grid_dimension[2] - 1)).astype(int)
))
print("Indices assigned")

#Grid initialization
if (grid_rebuild):
    print("Initializing grid...")
    grid = Grid(grid_dimension=grid_dimension, lower_bounds=lower_bounds, cell_size=cell_size, samples_df=samples_df,
                class_size=class_size, max_pheromone = max_pheromone)
    if(CellsFile is None):
        grid.initialize_regular_grid()
    else:
        grid.initialize_irregular_grid(inputcells_df)
    # Compute domain averages
    grid.compute_domain_averages()

    # Serialize the grid object to a file
    with open(InitGridFile, 'wb') as f:
        pickle.dump(grid, f)
    print(f"Initial grid saved to {InitGridFile}")
else:
    print(f"Loading initial grid from {InitGridFile}...")
    with open(InitGridFile, 'rb') as f:
        grid = pickle.load(f)
    print("Grid successfully loaded from serialized file.")

del inputcells_df
gc.collect()

# Get the list of cells that contain samples
sample_cells = [cell for cell in grid.cells.flatten() if cell is not None and cell.value is not None]
print(f'Cells with samples: {len(sample_cells)}')

# Choose cells for the ants
nAnts = int(len(sample_cells)*ants_per_sample)
if nAnts >= len(sample_cells):
    ant_cells = [cell for cell in sample_cells for _ in range(ants_per_sample)]
else:
    ant_cells = np.random.choice(sample_cells, size=nAnts, replace=True)

#max_pheromone = max([int((upper_bounds[i]-lower_bounds[i])/cell_size[i]) for i in range(3)])*pheromone_factor
print(f"Max pheromone: {max_pheromone}")

# Initialize the ants at the centroids of the chosen cells
ants = [
        Ant(grid = grid, position = cell.position, lower_bounds = lower_bounds, upper_bounds = upper_bounds, cell_size = cell_size,
            domain=cell.domain if cell.domain is not None else default_domain,
            index = i, life=max_life, max_pheromone = max_pheromone,
            coefficients = np.random.dirichlet(np.ones(27),size=1)[0],
            class_size = class_size, tolerance=tolerance, decay_rate = decay_rate,
            average_distance=average_distance, background_distance=background_distance,
            self_inclusion=self_inclusion, none_is_value=none_is_value, background=background,
            weight_method = weight_method
            ) for i,cell in enumerate(ant_cells)
       ]
print(f'Ants: {len(ants)}')
for ant in ants:
    current_cell = grid.get_cell(ant.position)
    ant.update_position(current_cell, current_cell)
    ant.update_trail(grid)

# Create a DataFrame from the list of ant coordinates
ants_df = pd.DataFrame([ant.position for ant in ants], columns=['x', 'y', 'z'])

# Save the DataFrame as a CSV file
ants_df.to_csv('ants_Initial_coordinates.csv', index=False)

print(f"Class size: {class_size} - Tolerance: {tolerance}")
print(f"Background value: {background}")

# Create the interpolator and run it
print("Starting interpolation...")
interpolator = Interpolator(ants = ants, max_life = max_life, grid = grid, background = background, fade_rate=fade_rate) # Interpolator(ants, grid) uses a background value of 0 as default
interpolator.run(num_iterations=iterations, pinclusion=self_inclusion)
print("Interpolation complete.")

# Calculate cells visited directly from the grid object
cells_visited = np.sum(np.vectorize(lambda cell: cell.value is not None)(grid.cells))
print(f"Cells visited: {cells_visited} ({(cells_visited / total_cells) * 100:.2f}%)")

# Fill unvisited cells (optional)
if fillaverage:
    print("Filling unvisited cells with average value of neighbors...")
    grid.fill_unvisited_cells()
elif fillbackground:
    print(f"Filling unvisited cells with background value of {background}...")
    grid.fill_unvisited_cells(background)

# Create the DataFrame for the grid cells
print("Creating csv file...")
cells_df = pd.DataFrame([(cell.position[0], cell.position[1], cell.position[2],
                          cell.value, cell.domain, cell.has_samples,
                          cell.num_ants, cell.best_xp,
                          cell.visits, cell.iteration)
                         for cell in grid.cells.flatten() if cell is not None and cell.value is not None],
                        columns=['x', 'y', 'z',
                                 'Value', 'Domain', 'Has_Samples',
                                 'Num_Ants', 'Best_XP',
                                 'Visits','Iteration'])

cells_df['size'] = cells_df['Has_Samples'].apply(lambda has_samples: 20 if has_samples else 8)

# Save the DataFrame to a CSV file
if weight_method == "domain_average":
    wm = "da"
elif weight_method ==  "value_consistence":
    wm = "vc"
elif weight_method ==  "greater_value":
    wm = "gv"
elif weight_method ==  "nearest_sample":
    wm = "ns"
elif weight_method ==  "greater_sample":
    wm = "gs"
elif weight_method ==  "pheromone_strenght":
    wm = "ps"
else:
    wm = weight_method
cells_file = (
    f"{base_name}_{cell_size[0]}x{cell_size[1]}x{cell_size[2]}_a{nAnts}_cs{class_size}_"
    f"b{background}_ad{average_distance}_bd{background_distance}_{wm}_dr{decay_rate}_i{iterations}.csv"
)
cells_df.to_csv(cells_file, index=False)
print("Cells saved in " + cells_file)

if export_pheromones:
    print("Exporting pheromone data...")
    # Collect all possible tracked_values (marks of class)
    tracked_samples = set()
    for cell in tqdm(grid.cells.flatten(), desc="Collecting tracked samples", total=grid.total_cells):
        if cell is not None and cell.value is not None:
            tracked_samples.update(cell.tracked_samples.keys())

    # Convert the set to a sorted list
    sorted_tracked_samples = sorted(tracked_samples,reverse=True)

    print("Tracked samples (marks of class):", sorted_tracked_samples)

    # For each tracked_value, create a CSV file with pheromone intensities > 0
    pheromone_files = []
    print("Processing tracked values...")
    for tracked_sample in tqdm(sorted_tracked_samples, desc="Processing tracked samples"):
        # Initialize an empty list to collect data
        data = []
        # Loop over the cells with a progress bar
        for cell in tqdm(grid.cells.flatten(), desc=f"Processing cells for tracked sample {tracked_sample}", total=grid.total_cells, leave=False):
            if cell is not None and cell.value is not None and tracked_sample in cell.tracked_samples:
                sample_signature = cell.tracked_samples[tracked_sample].signature
                pheromone_intensity = cell.tracked_samples[tracked_sample].pheromone
                closest_sample=cell.tracked_samples[tracked_sample].closest_sample
                sample_distance=cell.tracked_samples[tracked_sample].sample_distance
                cell_mode = cell.tracked_samples[tracked_sample].mode
                if pheromone_intensity > 0:
                    data.append((
                        cell.position[0], cell.position[1], cell.position[2],
                        cell.iteration,
                        sample_signature, pheromone_intensity,closest_sample,sample_distance, cell_mode  # Only pheromone intensity
                    ))

        # Create the DataFrame from the collected data
        cells_pheromone_df = pd.DataFrame(data, columns=[
            'x', 'y', 'z',
            'Iteration',
            'Sample_Signature', 'Pheromone','Closest_Sample','Sample_Distance','Mode'
        ])

        # Optionally, skip creating empty DataFrames
        if cells_pheromone_df.empty:
            print(f"No cells with pheromone >0 for tracked value {tracked_sample}. Skipping file creation.")
            continue

        # Save the DataFrame to a CSV file
        pheromone_file = (
            f"{base_name}_pheromone_{tracked_sample}.csv"
        )
        cells_pheromone_df.to_csv(pheromone_file, index=False)
        print(f"Pheromone data for tracked value {tracked_sample} saved in {pheromone_file}")

        pheromone_files.append(pheromone_file)


# After saving the cells_file, print the required data
print("8<----------------------------------------------------")
print(f"Samples: {len(samples_df)}")
print(f"Cell Size: {cell_size[0]}x{cell_size[1]}x{cell_size[2]}")
print(f"Total cells: {total_cells}")
print(f"Cells with samples: {len(sample_cells)}")
print(f"Starting ants: {nAnts} - Surviving ants: {len(ants)}")
print(f"Max pheromone: {max_pheromone}")
print(f"Class size: {class_size} - Tolerance: {tolerance}")
print(f"Domain average distance: {average_distance} - Background distance: {background_distance}")
print(f"Background value: {background}")
print(f"Decay rate: {decay_rate}")
print(f"Probability of including cell value during prediction: {self_inclusion}")
print(f"Value for None cells: {none_is_value}")
print(f"Weight method: {weight_method}")
print(f"Fill unvisited cells with neighbors average: {fillaverage}")
print(f"Iterations: {iterations}")
print(f"Cells visited: {cells_visited} ({(cells_visited / total_cells) * 100:.2f}%)")
print("8<----------------------------------------------------")

if plotscatter:

    # Create a 3D scatter plot of the interpolated values
    fig1 = px.scatter_3d(cells_df, x='x', y='y', z='z', color='Value', size='size',
                         color_continuous_scale=['green','yellow','orange','red'],
                         render_mode='webgl')
    fig1.update_traces(marker=dict(line=dict(width=2, color='Black')))
    fig1.update_layout(scene=dict(
        xaxis=dict(range=[lower_bounds[0], upper_bounds[0]], tickmode='linear', tick0=lower_bounds[0], dtick=cell_size[0],backgroundcolor="black"),
        yaxis=dict(range=[lower_bounds[1], upper_bounds[1]], tickmode='linear', tick0=lower_bounds[1], dtick=cell_size[1],backgroundcolor="black"),
        zaxis=dict(range=[lower_bounds[2], upper_bounds[2]], tickmode='linear', tick0=lower_bounds[2], dtick=cell_size[2],backgroundcolor="black"),
        camera=dict(projection=dict(type='orthographic'))
    ))
    fig1.show()

""" if download:
  files.download(cells_file)
  if export_pheromones:
    for pheromone_file in pheromone_files:
        files.download(pheromone_file) """