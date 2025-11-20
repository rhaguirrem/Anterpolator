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
        # Get all blocks including newly created ones
        colony_blocks = interpolator.get_blocks_data()
        
        # Track existing block positions
        existing_positions = set()
        new_blocks = []
        
        # Disable rendering during updates
        plotter.ren_win.SetOffScreenRendering(1)
        
        # Collect all blocks to add before updating mesh
        for pos, value in tqdm(colony_blocks.items(), desc="Processing blocks"):
            if pos in blocks.sample_blocks:
                continue
                
            if pos not in existing_positions and (pos in blocks.sample_blocks or value >= value_filter):
                corner = min_bounds + np.array(pos) * block_size
                new_block = pv.Box(bounds=(
                    corner[0], corner[0] + block_size,
                    corner[1], corner[1] + block_size,
                    corner[2], corner[2] + block_size
                ))
                new_block['Value'] = np.full(new_block.n_points, value)
                new_blocks.append(new_block)
                
                if plotter.verbose:
                    print(f"Creating block at {pos} with value {value:.2f}")
        
        if new_blocks:
            print(f"\nAdding {len(new_blocks)} new blocks to visualization...")
            
            # Create single multiblock from all new blocks
            new_multiblock = pv.MultiBlock(new_blocks)
            blocks.combine(new_multiblock)
            
            # Update display with combined mesh
            plotter.remove_actor(plotter._blocks_actor)
            plotter._blocks_actor = plotter.add_mesh(
                blocks,
                style='surface',
                scalars='Value',
                opacity=0.5,
                show_edges=True,
                cmap='rainbow',
                immediate=True,  # Force immediate update
                culling='back'   # Enable back-face culling
            )
            
            # Re-enable rendering and update
            plotter.ren_win.SetOffScreenRendering(0)
            plotter.reset_camera_clipping_range()
            plotter.render()
            
            print("Visualization updated with new blocks")
    else:
        print("No new blocks created or updated")

def load_and_visualize_samples(csv_file, block_size=10, value_filter=60, verbose=False):  # Add parameter
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
        plotter.verbose = verbose  # Store verbose setting in plotter
        
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
    block_size = 10
    value_filter = 40  # Already added
    verbose = False  # Added for verbosity control
    # Update function call with new parameter
    load_and_visualize_samples(csv_file, block_size=block_size, value_filter=value_filter, verbose=verbose)