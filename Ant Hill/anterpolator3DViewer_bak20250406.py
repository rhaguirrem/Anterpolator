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

def create_blocks(points, values, block_size=10, verbose=False):
    # Calculate bounding box
    min_bounds = np.min(points, axis=0)
    max_bounds = np.max(points, axis=0)
    
    # Calculate number of blocks in each dimension
    dims = np.ceil((max_bounds - min_bounds) / block_size).astype(int)
    
    # Initialize block values
    blocks = {}
    block_values = {}
    
    # Store extra info for ant colony
    block_info = {
        'min_bounds': min_bounds,
        'dims': dims,
        'block_size': block_size
    }
    
    # Assign points to blocks with progress bar
    print("Assigning points to blocks...")
    for point, value in tqdm(zip(points, values), total=len(values)):
        block_idx = tuple(((point - min_bounds) // block_size).astype(int))
        if block_idx not in blocks:
            blocks[block_idx] = []
            block_values[block_idx] = []
        blocks[block_idx].append(point)
        block_values[block_idx].append(value)
    
    # Create block geometries
    block_data = []
    print("\nCreating blocks...")
    for idx in tqdm(blocks.keys(), desc="Creating blocks"):
        # Create block geometry
        corner = min_bounds + np.array(idx) * block_size
        block = pv.Box(bounds=(
            corner[0], corner[0] + block_size,
            corner[1], corner[1] + block_size,
            corner[2], corner[2] + block_size
        ))
        
        # Set block value
        avg_value = np.mean(block_values[idx])
        block['Value'] = np.full(block.n_points, avg_value)
        block_data.append(block)
    
    # Store sample information
    sample_blocks = {}
    for idx in blocks.keys():
        avg_value = np.mean(block_values[idx])
        sample_blocks[idx] = avg_value
    
    # Create multiblock with extra attributes
    multiblock = pv.MultiBlock(block_data)
    multiblock.block_info = block_info
    multiblock.sample_blocks = sample_blocks
    
    # Initialize ant colony without value filter
    interpolator = AntColonyInterpolator(range_size=10, max_pheromone=150, ants_per_sample=3, verbose=verbose)
    interpolator.initialize_blocks(sample_blocks, tuple(dims))
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
        
        # Track existing block positions
        block_mapping = {}
        for i, block in enumerate(blocks):
            corner = block.bounds[::2]
            pos = tuple(np.floor((corner - min_bounds) / block_size).astype(int))
            block_mapping[pos] = (i, block)
        
        print(f"Processing blocks (value filter: {value_filter})...")
        new_blocks = []
        modified = False
        
        # Process blocks
        for pos, value in tqdm(colony_blocks.items(), desc="Updating blocks"):
            if pos in blocks.sample_blocks:
                continue
            
            if pos not in block_mapping:
                # Create new block if it meets value threshold
                if value >= value_filter:
                    corner = min_bounds + np.array(pos) * block_size
                    new_block = pv.Box(bounds=(
                        corner[0], corner[0] + block_size,
                        corner[1], corner[1] + block_size,
                        corner[2], corner[2] + block_size
                    ))
                    new_block['Value'] = np.full(new_block.n_points, value)
                    new_blocks.append(new_block)
            else:
                # Update existing block value if needed
                _, block = block_mapping[pos]
                if abs(block.point_data['Value'][0] - value) > 0.0001:
                    block.point_data['Value'][:] = value
                    modified = True

        # Handle updates efficiently
        if new_blocks or modified:
            if new_blocks:
                print(f"Adding {len(new_blocks)} new blocks...")
                # Add new blocks to multiblock
                for block in new_blocks:
                    blocks.append(block)
            
            # Create fresh mesh actor with updated blocks
            plotter.remove_actor(plotter._blocks_actor)
            plotter._blocks_actor = plotter.add_mesh(
                blocks,
                style='surface',
                scalars='Value',
                opacity=0.5,
                show_edges=True,
                cmap='rainbow',
                clim=[0, max(colony_blocks.values())]  # Set color range
            )
            plotter.render()
            print("Visualization updated")
            
    else:
        print("No blocks to update")

def export_blocks_to_csv(blocks, filepath=None):
    # Generate default filepath if none provided
    if filepath is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        filepath = os.path.join("C:/Projects/Anterpolator", f"interpolation_{timestamp}.csv")
    
    print(f"Exporting blocks to {filepath}...")
    
    # Prepare data for export
    data = []
    min_bounds = blocks.block_info['min_bounds']
    block_size = blocks.block_info['block_size']
    
    for pos, block in tqdm(blocks.ant_colony.blocks.items(), desc="Processing blocks"):
        # Calculate block centroid
        corner = min_bounds + np.array(pos) * block_size
        centroid = corner + block_size/2
        
        # Get block data
        data.append({
            'x': centroid[0],
            'y': centroid[1],
            'z': centroid[2],
            'Value': block.value,
            'Visits': block.visit_count,
            'Pheromone': block.pheromone,
            'Is_Sample': block.is_sample,
            'Mark_Class': block.mark_class
        })
    
    # Export to CSV
    df = pd.DataFrame(data)
    df.to_csv(filepath, index=False)
    print(f"Exported {len(data)} blocks to {filepath}")

def silent_interpolation(plotter, iterations):
    print(f"Running {iterations} silent interpolation iterations...")
    blocks = plotter.blocks_data
    interpolator = blocks.ant_colony
    dims = tuple(blocks.block_info['dims'])
    
    for i in tqdm(range(iterations), desc="Interpolating"):
        interpolator.move_ants(dims)
        interpolator.decay_pheromone()
    
    # Export results
    export_path = export_blocks_to_csv(blocks)
    print("Silent interpolation complete")

def load_and_visualize_samples(csv_file, block_size=10, value_filter=60, verbose=False, iterations=100):
    try:
        # Read CSV
        df = pd.read_csv(csv_file)
        points = df[['x', 'y', 'z']].values
        values = df['Value'].values

        # Setup visualization first
        plotter = pv.Plotter()
        plotter.set_background('white')
        
        # Create blocks and store in plotter
        blocks = create_blocks(points, values, block_size, verbose)
        plotter.blocks_data = blocks  # Store after plotter creation
        plotter.verbose = verbose
        
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
        plotter.add_key_event('I', lambda: silent_interpolation(plotter, iterations))  # Shift+I
        
        # Store iterations setting
        plotter.iterations = iterations
        
        # Show visualization
        plotter.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    csv_file = "C:\Projects\Anterpolator\Samples_CNN_ABA_BED.csv"
    block_size = 20
    value_filter = 0
    verbose = False
    iterations = 10  # Number of iterations for silent interpolation
    
    load_and_visualize_samples(
        csv_file, 
        block_size=block_size,
        value_filter=value_filter,
        verbose=verbose,
        iterations=iterations
    )