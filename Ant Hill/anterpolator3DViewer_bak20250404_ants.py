import pandas as pd
import pyvista as pv
import numpy as np
import csv
from tqdm import tqdm
import sys
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

def create_blocks(points, values, block_size=10):
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
    
    # Initialize ant colony
    interpolator = AntColonyInterpolator(range_size=10, max_pheromone=150, ants_per_sample=3)
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
    block_size = blocks.block_info['block_size']
    min_bounds = blocks.block_info['min_bounds']
    
    # Move ants and update pheromones
    print("Moving ants...")
    interpolator.move_ants(dims)
    interpolator.decay_pheromone()
    
    # Get all blocks including newly visited ones
    all_blocks = interpolator.get_blocks_data()
    changes = 0
    
    # Create new blocks for newly visited positions
    for coords, value in all_blocks.items():
        if coords not in blocks.sample_blocks:
            block_exists = False
            for existing_block in blocks:
                corner = existing_block.bounds[::2]  # Get lower corner
                block_coords = tuple(np.floor((corner - min_bounds) / block_size).astype(int))
                if block_coords == coords:
                    # Update existing block
                    if abs(existing_block.point_data['Value'][0] - value) > 0.0001:
                        existing_block.point_data['Value'][:] = value
                        changes += 1
                    block_exists = True
                    break
            
            if not block_exists:
                # Create new block
                corner = min_bounds + np.array(coords) * block_size
                new_block = pv.Box(bounds=(
                    corner[0], corner[0] + block_size,
                    corner[1], corner[1] + block_size,
                    corner[2], corner[2] + block_size
                ))
                new_block['Value'] = np.full(new_block.n_points, value)
                blocks.append(new_block)
                changes += 1
    
    if changes > 0:
        print(f"Updated/created {changes} blocks, refreshing display")
        plotter._blocks_actor.mapper.dataset.Modified()
        plotter._blocks_actor.mapper.Update()
        plotter.render()
    else:
        print("No blocks needed updating")

def load_and_visualize_samples(csv_file, block_size=10):
    try:
        # Read CSV
        df = pd.read_csv(csv_file)
        points = df[['x', 'y', 'z']].values
        values = df['Value'].values

        # Setup visualization first
        plotter = pv.Plotter()
        plotter.set_background('white')

        # Create blocks and store in plotter
        blocks = create_blocks(points, values, block_size)
        plotter.blocks_data = blocks  # Store after plotter creation
        
        # Add point cloud
        cloud = pv.PolyData(points)
        cloud['Values'] = values
        plotter.add_mesh(cloud, 
                        render_points_as_spheres=True,
                        point_size=5,
                        scalars='Values',
                        cmap='rainbow')

        # Add blocks
        plotter._blocks_actor = plotter.add_mesh(
            blocks,
            style='surface',
            scalars='Value',
            opacity=0.5,
            show_edges=True,
            cmap='rainbow'
        )
        
        # Add controls
        plotter.add_key_event('b', lambda: toggle_blocks(plotter))
        plotter.add_key_event('i', lambda: update_interpolation(plotter))
        
        # Show visualization
        plotter.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    csv_file = "C:\Projects\Anterpolator\Samples_CNN_ABA_BED.csv"
    block_size = 50  # Adjust block size as needed
    # Load and visualize samples
    load_and_visualize_samples(csv_file, block_size=block_size)