import pandas as pd
import pyvista as pv
import numpy as np
import csv
from tqdm import tqdm
import sys
import os
from datetime import datetime
sys.path.append("C:/Projects/Anterpolator")
from ant_colony import AntColonyInterpolator

def format_point_info(point, value):
    return f"Position: ({point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f})\nValue: {value:.2f}"

def middle_click_callback(plotter, event_id):
    if hasattr(plotter, 'picked_point'):
        point = plotter.picked_point
        if point is not None:
            # Update camera focus
            plotter.camera_position = [
                plotter.camera.position,
                point,
                plotter.camera.up
            ]
            plotter.renderer.GetActiveCamera().SetFocalPoint(point)
            plotter.render()

def create_blocks(points, values, block_size=10, verbose=False, range_size=10, max_pheromone=150,
                  ants_per_sample=3, blocks_file=None):
    if blocks_file is not None:
        print("Loading predefined cells from blocks_file...")
        df_blocks = pd.read_csv(blocks_file)
        # Use ALL rows from df_blocks for boundaries
        centroids_all = df_blocks[['x','y','z']].values
        all_min_bounds = np.min(centroids_all, axis=0)
        all_max_bounds = np.max(centroids_all, axis=0)
        # Calculate unified block dimensions using all df_blocks:
        dims_all = []
        for col in ['x','y','z']:
            uniq = np.sort(df_blocks[col].unique())
            diffs = np.diff(uniq)
            positive_diffs = diffs[diffs > 0]
            dim = positive_diffs.min() if len(positive_diffs) > 0 else 10
            dims_all.append(dim)
        unified_dims = np.array(dims_all)
        dims_grid = np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int)
        print("Calculated unified block dimensions:", unified_dims)
        print("Calculated grid dimensions:", dims_grid)
        # Build mapping from grid index to Domain using ALL rows from df_blocks
        domain_mapping = {}
        for idx, row in df_blocks.iterrows():
            centroid = np.array([row['x'], row['y'], row['z']])
            grid_idx = tuple(np.floor((centroid - all_min_bounds) / unified_dims + 1e-6).astype(int))
            domain = row['Domain'] if 'Domain' in row and pd.notnull(row['Domain']) and str(row['Domain']).strip() != "" else "Undomained"
            domain_mapping[grid_idx] = domain
        allowed_grid = set(domain_mapping.keys())
        # Group sample points using the same all_min_bounds
        sample_blocks_dict = {}
        block_values = {}
        for point, value in tqdm(zip(points, values), total=len(values), desc="Assigning points to blocks"):
            block_idx = tuple(np.floor((point - all_min_bounds) / unified_dims + 1e-6).astype(int))
            if block_idx in allowed_grid:
                if block_idx not in sample_blocks_dict:
                    sample_blocks_dict[block_idx] = []
                    block_values[block_idx] = []
                sample_blocks_dict[block_idx].append(point)
                block_values[block_idx].append(value)
        block_data = []
        for idx in tqdm(sample_blocks_dict.keys(), desc="Creating blocks"):
            corner = all_min_bounds + np.array(idx) * unified_dims
            cell = pv.Box(bounds=(
                corner[0], corner[0] + unified_dims[0],
                corner[1], corner[1] + unified_dims[1],
                corner[2], corner[2] + unified_dims[2]
            ))
            avg_value = np.mean(block_values[idx])
            cell.cell_data['Value'] = np.full(cell.n_cells, avg_value)
            cell.cell_data['Is_Sample'] = np.full(cell.n_cells, True)
            cell.cell_data['Block_ID'] = np.full(cell.n_cells, 0)  # to be set later
            domain = domain_mapping.get(idx, "Undomained")
            cell.cell_data['Domain'] = np.full(cell.n_cells, domain)
            block_data.append(cell)
        block_info = {
            'min_bounds': all_min_bounds,
            'dims': np.ceil((all_max_bounds - all_min_bounds) / unified_dims).astype(int),
            'block_size': unified_dims.tolist(),
            'allowed_grid': list(allowed_grid)
        }
        multiblock = pv.MultiBlock(block_data)
        multiblock.block_info = block_info
        multiblock.sample_blocks = {idx: np.mean(vals) for idx, vals in block_values.items()}
        interpolator = AntColonyInterpolator(range_size=range_size, max_pheromone=max_pheromone,
                                               ants_per_sample=ants_per_sample, verbose=verbose)
        # Attach allowed_grid and domain_mapping to the interpolator
        interpolator.allowed_grid_override = allowed_grid
        interpolator.domain_mapping = domain_mapping  # Pass domain mapping to interpolator
        interpolator.initialize_blocks(multiblock.sample_blocks, tuple(block_info['dims']),
                                       all_min_bounds, unified_dims.tolist(), use_domain_mapping=True)
        interpolator.create_ants()
        multiblock.ant_colony = interpolator
        return multiblock
    else:
        print("Assigning points to blocks...")
        min_bounds = np.min(points, axis=0)
        max_bounds = np.max(points, axis=0)
        dims = np.ceil((max_bounds - min_bounds) / np.array(block_size)).astype(int)
        blocks = {}
        block_values = {}
        block_info = {
            'min_bounds': min_bounds,
            'dims': dims,
            'block_size': block_size
        }
        for point, value in tqdm(zip(points, values), total=len(values)):
            block_idx = tuple(((point - min_bounds) // np.array(block_size)).astype(int))
            if block_idx not in blocks:
                blocks[block_idx] = []
                block_values[block_idx] = []
            blocks[block_idx].append(point)
            block_values[block_idx].append(value)
        block_data = []
        next_block_id = 1
        for idx in tqdm(blocks.keys(), desc="Creating blocks"):
            corner = min_bounds + np.array(idx) * np.array(block_size)
            cell = pv.Box(bounds=(
                corner[0], corner[0] + block_size[0],
                corner[1], corner[1] + block_size[1],
                corner[2], corner[2] + block_size[2]
            ))
            avg_value = np.mean(block_values[idx])
            cell.cell_data['Value'] = np.full(cell.n_cells, avg_value)
            cell.cell_data['Is_Sample'] = np.full(cell.n_cells, True)
            cell.cell_data['Block_ID'] = np.full(cell.n_cells, next_block_id)
            next_block_id += 1
            block_data.append(cell)
        sample_blocks = {idx: np.mean(vals) for idx, vals in block_values.items()}
        multiblock = pv.MultiBlock(block_data)
        multiblock.block_info = block_info
        multiblock.sample_blocks = sample_blocks
        interpolator = AntColonyInterpolator(range_size=range_size, max_pheromone=max_pheromone,
                                               ants_per_sample=ants_per_sample, verbose=verbose)
        interpolator.initialize_blocks(sample_blocks, tuple(block_info['dims']),
                                       min_bounds, block_size, use_domain_mapping=False)
        interpolator.create_ants()
        multiblock.ant_colony = interpolator
        return multiblock

def toggle_blocks(plotter):
    if hasattr(plotter, '_blocks_actor'):
        is_visible = not plotter._blocks_actor.GetVisibility()
        plotter._blocks_actor.SetVisibility(is_visible)
        plotter.render()

def update_interpolation(plotter):
    print("Checking ant colony data...")
    if not hasattr(plotter, 'blocks_data'):
        print("No blocks data found")
        return
        
    blocks = plotter.blocks_data
    if not hasattr(blocks, 'ant_colony'):
        print("No ant colony found in blocks")
        return
        
    print("Found ant colony, updating interpolation...")
    interpolator = blocks.ant_colony
    dims = tuple(blocks.block_info['dims'])
    min_bounds = blocks.block_info['min_bounds']
    block_size = blocks.block_info['block_size']
    
    # Move ants and track if changes were made
    print("Moving ants...")
    changes_made = interpolator.move_ants(dims)
    interpolator.decay_pheromone()
    
    if changes_made:
        colony_blocks = interpolator.get_blocks_data()
        
        # Track blocks by both position and ID
        block_mapping = {}
        id_mapping = {}
        
        # First pass: map existing blocks and track IDs - no value filtering here
        for i, block in enumerate(blocks):
            corner = block.bounds[::2]
            pos = tuple(np.floor((corner - min_bounds) / np.array(block_size)).astype(int))
            block_id = block.cell_data['Block_ID'][0]
            
            if pos in blocks.sample_blocks:
                block_mapping[pos] = (i, block)
                id_mapping[block_id] = pos
            else:
                # Keep track of all non-sample blocks regardless of value
                block_mapping[pos] = (i, block)
        
        # Process blocks
        new_blocks = []
        modified = False
        
        # Update blocks - create all blocks, filter only for display
        for pos, value in colony_blocks.items():
            try:
                if pos in blocks.sample_blocks:
                    continue
                
                if pos not in block_mapping:
                    # Create new interpolated block regardless of value
                    corner = min_bounds + np.array(pos) * np.array(block_size)
                    half_size = np.array(block_size) / 2
                    center = corner + half_size
                    
                    new_block = pv.Box(bounds=(
                        center[0] - half_size[0]/2, center[0] + half_size[0]/2,
                        center[1] - half_size[1]/2, center[1] + half_size[1]/2,
                        center[2] - half_size[2]/2, center[2] + half_size[2]/2
                    ))
                    new_block.cell_data['Value'] = np.full(new_block.n_cells, value)
                    new_block.cell_data['Is_Sample'] = np.full(new_block.n_cells, False)
                    new_block.cell_data['Block_ID'] = np.full(new_block.n_cells, interpolator.next_block_id)
                    interpolator.next_block_id += 1
                    new_blocks.append(new_block)
                elif pos in block_mapping:
                    _, block = block_mapping[pos]
                    if abs(block.cell_data['Value'][0] - value) > 0.0001:
                        block.cell_data['Value'][:] = value
                        modified = True
            except Exception as e:
                print(f"Error processing position {pos}: {str(e)}")
                continue

        # Add all new blocks to the multiblock
        if new_blocks:
            print(f"Adding {len(new_blocks)} new blocks...")
            for block in new_blocks:
                blocks.append(block)
        
        # Filter blocks for display only
        visible_blocks = pv.MultiBlock()
        for block in blocks:
            if block.cell_data['Is_Sample'][0]:  # Always show sample blocks
                visible_blocks.append(block)
            elif block.cell_data['Value'][0] >= plotter.value_filter:
                visible_blocks.append(block)
        
        # Create fresh mesh actor with filtered blocks
        plotter.remove_actor(plotter._blocks_actor)
        plotter._blocks_actor = plotter.add_mesh(
            visible_blocks,
            style='surface',
            scalars='Value',
            opacity=0.5,
            show_edges=True,
            cmap='rainbow',
            clim=[0, max(colony_blocks.values())]
        )
        plotter.render()
        print("Visualization updated")
    else:
        print("No blocks to update")

def export_blocks_to_csv(blocks, filepath):
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    print(f"Exporting blocks to {filepath}...")
    
    # Prepare data for export
    data = []
    min_bounds = blocks.block_info['min_bounds']
    block_size = blocks.block_info['block_size']
    
    for pos, block in tqdm(blocks.ant_colony.blocks.items(), desc="Processing blocks"):
        # Calculate block centroid
        corner = min_bounds + np.array(pos) * np.array(block_size)
        centroid = corner + np.array(block_size)/2
        
        # Get block data including domain
        data.append({
            'x': centroid[0],
            'y': centroid[1],
            'z': centroid[2],
            'Value': block.value,
            'Visits': block.visit_count,
            'Pheromone': block.pheromone,
            'Is_Sample': block.is_sample,
            'Mark_Class': block.mark_class,
            'Block_ID': block.block_id,
            'Distance_To_Sample': block.distance_to_sample,
            'Nearest_Sample_Value': block.nearest_sample_value,  # Already included
            'Ant_Count': block.ant_count,
            'Domain': block.domain   # NEW: Export domain
        })
    
    # Export to CSV
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Exported {len(data)} blocks to {filepath}")

def silent_interpolation(plotter, iterations, interpolation_file):
    blocks = plotter.blocks_data
    interpolator = blocks.ant_colony
    dims = tuple(blocks.block_info['dims'])
    
    for i in tqdm(range(iterations), desc="Interpolating"):
        interpolator.move_ants(dims, quiet=True)
        interpolator.decay_pheromone(quiet=True)
    
    # Export results
    export_blocks_to_csv(blocks, interpolation_file)

def load_and_visualize_samples(samples_file, block_size=10, value_filter=60, verbose=False, iterations=100, range_size=10, max_pheromone=150, ants_per_sample=3, blocks_file=None):
    try:
        # Read CSV
        df = pd.read_csv(samples_file)
        points = df[['x', 'y', 'z']].values
        values = df['Value'].values

        # Setup visualization first
        plotter = pv.Plotter()
        plotter.set_background('white')
        
        # Create blocks and store in plotter using new parameters
        blocks = create_blocks(points, values, block_size, verbose, range_size, max_pheromone, ants_per_sample, blocks_file)
        
        # Store settings in plotter
        plotter.blocks_data = blocks
        plotter.verbose = verbose
        plotter.value_filter = value_filter  # Store value_filter
        
        # Add point cloud with scalar bar settings
        cloud = pv.PolyData(points)
        cloud['Values'] = values
        sargs = dict(
            vertical=True,
            position_x=0.93,
            position_y=0.2,
            height=0.6,
            label_font_size=10,
            title='Value'
        )
        plotter.add_mesh(cloud, 
                        render_points_as_spheres=True,
                        point_size=5,
                        scalars='Values',
                        cmap='rainbow',
                        scalar_bar_args=sargs)

        # Add blocks without scalar bar
        plotter._blocks_actor = plotter.add_mesh(
            blocks,
            style='surface',
            scalars='Value',
            opacity=0.5,
            show_edges=True,
            cmap='rainbow',
            show_scalar_bar=False  # Hide second scalar bar
        )
        
        # Add controls
        plotter.add_key_event('b', lambda: toggle_blocks(plotter))
        plotter.add_key_event('i', lambda: update_interpolation(plotter))
        plotter.add_key_event('I', lambda: silent_interpolation(plotter, iterations, interpolation_file))  # Shift+I
        
        # Store iterations setting
        plotter.iterations = iterations
        
        # Show visualization
        plotter.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    
    #File paths and settings
    samples_file = "Data\\ANT-Samples.csv" # Path to the CSV file with sample data
    blocks_file = "Data\\ANT-Domains.csv" # Defintion of cells and domains
    output_folder = "output"
    interpolation_file = os.path.join(output_folder, f"{os.path.splitext(os.path.basename(samples_file))[0]}_interpolation.csv")
    
    # Interpolation settings
    block_size = (10,10,10)  # Block size in each dimension
    range_size = .2         # Size used to define mark of class 
    max_pheromone = 150     # Maximum pheromone value
    ants_per_sample = 3     # Number of ants per sample
    iterations = 500        # Number of iterations for silent interpolation

    # Display settings
    value_filter = 0
    verbose = False
    
    load_and_visualize_samples(
        samples_file, 
        block_size=block_size,
        value_filter=value_filter,
        verbose=verbose,
        iterations=iterations,
        range_size=range_size,
        max_pheromone=max_pheromone,
        ants_per_sample=ants_per_sample,
        blocks_file=blocks_file
    )