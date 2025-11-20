import pandas as pd
import pyvista as pv
import numpy as np
import csv
from tqdm import tqdm

def create_blocks(points, values, block_size=10):
    # Calculate bounding box
    min_bounds = np.min(points, axis=0)
    max_bounds = np.max(points, axis=0)
    
    # Calculate number of blocks in each dimension
    dims = np.ceil((max_bounds - min_bounds) / block_size).astype(int)
    
    # Initialize block values
    blocks = {}
    block_values = {}
    
    # Assign points to blocks with progress bar
    print("Assigning points to blocks...")
    for point, value in tqdm(zip(points, values), total=len(values)):
        block_idx = tuple(((point - min_bounds) // block_size).astype(int))
            plotter.renderer.GetActiveCamera().SetFocalPoint(point)
            
            # Get point data
            dataset = plotter.picked_data
            if dataset is not None:
                point_id = plotter.picked_point_id
                if point_id >= 0:
                    value = dataset.point_data['Values'][point_id]
                    # Update point info display
                    plotter.point_info_text.set_text(text=format_point_info(point, value), position='lower_left')
                    
            plotter.renderer.ResetCamera()
            plotter.render()

def toggle_grid(plotter):
    if hasattr(plotter, '_grid_actor'):
        plotter._grid_actor.SetVisibility(not plotter._grid_actor.GetVisibility())
        plotter.render()

def move_camera(plotter, direction):
    camera = plotter.camera
    speed = 0.5  # Adjust speed as needed
    
    # Calculate forward vector from position to focal point
    pos = np.array(camera.position)
    focal = np.array(camera.focal_point)
    forward = focal - pos
    forward = forward / np.linalg.norm(forward)
    
    # Calculate right vector
    right = np.cross(forward, camera.up)
    right = right / np.linalg.norm(right)
    
    # Calculate movement based on direction
    if direction == 'w':
        movement = forward * speed
    elif direction == 's':
        movement = -forward * speed
    elif direction == 'a':
        movement = -right * speed
    elif direction == 'd':
        movement = right * speed
    
    # Update camera position and focal point
    new_pos = pos + movement
    new_focal = focal + movement
    camera.position = new_pos
    camera.focal_point = new_focal
    plotter.render()

def create_blocks(points, values, block_size=10):
    # Calculate bounding box
    min_bounds = np.min(points, axis=0)
    max_bounds = np.max(points, axis=0)
    
    # Calculate number of blocks in each dimension
    dims = np.ceil((max_bounds - min_bounds) / block_size).astype(int)
    
    # Store min_bounds for coordinate mapping
    block_info = {
        'min_bounds': min_bounds,
        'dims': dims,
        'block_size': block_size
    }
    
    # Initialize block values
    blocks = {}
    block_values = {}
    
    # Assign points to blocks with progress bar
    print("Assigning points to blocks...")
    for point, value in tqdm(zip(points, values), total=len(values)):
        block_idx = tuple(((point - min_bounds) // block_size).astype(int))
        if block_idx not in blocks:
            blocks[block_idx] = []
            block_values[block_idx] = []
        blocks[block_idx].append(point)
        block_values[block_idx].append(value)
    
    # Create block geometries and prepare for ant colony
    block_data = []
    sample_blocks = {}
    block_map = {}
    
    # Create 3D array to store block references
    dims = tuple(dims)
    block_array = np.empty(dims, dtype=object)
    
    # Create a template block at (min_bounds - 1)
    template_position = min_bounds - 1
    template_block = pv.Box(bounds=(
        template_position[0], template_position[0] + 1,
        template_position[1], template_position[1] + 1,
        template_position[2], template_position[2] + 1
    ))
    template_block['Value'] = np.full(template_block.n_points, 0.0)
    
    print("\nTemplate block info:")
    print(f"Position: {template_position}")
    print(f"Bounds: {template_block.bounds}")
    print(f"Center: {template_block.center}")
    print(f"Size: 1 (will be scaled to {block_size})\n")
    
    # Create only sample blocks
    print("Creating sample blocks...")
    for idx in tqdm(blocks.keys(), desc="Creating sample blocks"):
        # Scale and translate template for this position
        block = template_block.copy()
        block.scale(block_size)  # Scale unit block to desired size
        corner = min_bounds + np.array(idx) * block_size
        block.translate(corner - template_position)  # Adjust translation based on template position
        
        # Set block value
        avg_value = np.mean(block_values[idx])
        block['Value'] = np.full(block.n_points, avg_value)
        
        # Store block
        block_data.append(block)
        block_map[idx] = block
        block_array[idx] = block
        sample_blocks[idx] = avg_value
        print(f"Sample at {idx}: value={avg_value:.2f}, mark_class={int(avg_value)}")

    # Delete template
    del template_block
    
    print("\nInitializing ant colony...")
    print(f"Number of sample blocks: {len(sample_blocks)}")
    print(f"Grid dimensions: {dims}")
    
    # Initialize ant colony interpolator with debug info
    interpolator = AntColonyInterpolator(range_size=10, max_pheromone=150, ants_per_sample=3)
    print("Created interpolator")
    
    interpolator.initialize_blocks(sample_blocks, tuple(dims))
    print("Initialized blocks")
    
    print("Creating ants...")
    interpolator.create_ants()  # Create initial ants
    print(f"Created {len(interpolator.ants)} ants")
    
    # Create multiblock with debug info
    print("\nCreating multiblock dataset...")
    multiblock = pv.MultiBlock(block_data) if block_data else None
    if multiblock:
        print("Attaching data to multiblock...")
        multiblock.ant_colony = interpolator
        multiblock.block_info = block_info
        multiblock.block_map = block_map
        multiblock.block_array = block_array
        multiblock.sample_blocks = set(sample_blocks.keys())
        
        # Add method to create blocks on demand
        def create_block(coords):
            if coords in block_map:
                return block_map[coords]
            
            # Create new block from template geometry
            x, y, z = coords
            block = pv.Box(bounds=(-1, 0, -1, 0, -1, 0))
            block.scale(block_size)
            corner = min_bounds + np.array((x, y, z)) * block_size
            block.translate(corner)
            block['Value'] = np.full(block.n_points, 0.0)
            
            # Store block
            block_map[coords] = block
            block_array[x, y, z] = block
            block_data.append(block)
            return block
            
        multiblock.create_block = create_block
        print("Multiblock creation complete")
    
    return multiblock

def toggle_blocks(plotter):
    if hasattr(plotter, 'block_actors'):
        print("\nToggling blocks visibility...")
        is_visible = not plotter.block_actors[0].GetVisibility()
        
        # Toggle all block actors
        for actor in plotter.block_actors:
            actor.SetVisibility(is_visible)
            actor.GetProperty().SetOpacity(0.5 if is_visible else 0.0)
            actor.Modified()
        
        plotter.renderer.ResetCamera()
        plotter.render()
        
        print(f"Visibility: {is_visible}")
        print(f"Number of actors toggled: {len(plotter.block_actors)}")

def update_blocks(plotter):
    print("Update blocks called")  # Debug info
    if not hasattr(plotter, '_blocks_actor'):
        print("No blocks actor found")
        return
        
    if not hasattr(plotter, 'blocks_data'):
        print("No blocks data found")
        return
        
    blocks = plotter.blocks_data
    
    if not hasattr(blocks, 'ant_colony'):
        print("No ant colony found")
        return
    
    interpolator = blocks.ant_colony
    info = blocks.block_info
    
    # Show iteration info and progress
    total_ants = len(interpolator.ants)
    active_ants = sum(1 for ant in interpolator.ants if ant.steps < 100)
    progress = (total_ants - active_ants) / total_ants * 100
    
    # Update progress text
    if hasattr(plotter, 'progress_text'):
        plotter.progress_text.set_text(
            text=f"Progress: {progress:.1f}% ({active_ants} active ants)",
            position='upper_left'
        )
    
    if active_ants == 0:
        print("No active ants remaining - interpolation complete")
        if hasattr(plotter, 'progress_text'):
            plotter.progress_text.set_text(
                text="Interpolation Complete!",
                position='upper_left'
            )
        return
    
    # Debug ant positions before movement
    print("\nAnt positions before movement:")
    for i, ant in enumerate(interpolator.ants):
        if ant.steps < 100:  # Only show active ants
            print(f"Ant {i}: pos={ant.current_pos}, steps={ant.steps}")
    
    # Move ants and update pheromones
    interpolator.move_ants(info['dims'])
    interpolator.decay_pheromone()
    
    # Debug ant positions after movement
    print("\nAnt positions after movement:")
    for i, ant in enumerate(interpolator.ants):
        if ant.steps < 100:  # Only show active ants
            print(f"Ant {i}: pos={ant.current_pos}, steps={ant.steps}")
    
    # Update visualization
    print(f"Updating block values... Progress: {progress:.1f}%")
    blocks = plotter.blocks_data
    block_array = blocks.block_array
    sample_blocks = blocks.sample_blocks
    
    # Track changes
    changes_made = 0
    errors = 0
    
    # Debug info
    print(f"Dimensions: {block_array.shape}")
    print(f"Sample blocks: {len(sample_blocks)}")
    print(f"Total blocks in interpolator: {len(interpolator.blocks)}")
    
    # Print some sample block values
    print("\nSample block values:")
    for coords, block_data in list(interpolator.blocks.items())[:5]:
        print(f"Block {coords}: value={block_data.value}")
    
    # Update all blocks
    for coords, block_data in interpolator.blocks.items():
        x, y, z = coords
        try:
            if 0 <= x < block_array.shape[0] and \
               0 <= y < block_array.shape[1] and \
               0 <= z < block_array.shape[2]:
                
                # Get or create block on demand
                current_block = block_array[x, y, z]
                if current_block is None and coords not in sample_blocks:
                    current_block = blocks.create_block(coords)
                
                if current_block is not None:
                    old_value = current_block.point_data['Value'][0]
                    new_value = block_data.value
                    
                    # Update if value changed
                    if coords not in sample_blocks:
                        if abs(new_value - old_value) > 0.0001:
                            current_block.point_data['Value'][:] = new_value
                            changes_made += 1
                            if changes_made < 10:  # Print first 10 updates
                                print(f"Block update: {coords} from {old_value:.2f} to {new_value:.2f}")
        except Exception as e:
            print(f"Error updating block {coords}: {str(e)}")
            errors += 1
    
    print(f"Updated {changes_made} blocks with {errors} errors")
    if changes_made > 0:
        plotter._blocks_actor.mapper.dataset.Modified()
        plotter._blocks_actor.mapper.Update()
    plotter.render()

def load_and_visualize_samples(csv_file, block_size=10):
    # Read CSV file
    try:
        # Try to detect delimiter, default to comma
        with open(csv_file, 'r') as file:
            first_line = file.readline()
            sniffer = csv.Sniffer()
            dialect = sniffer.sniff(first_line)
            delimiter = dialect.delimiter if dialect.delimiter else ','

        df = pd.read_csv(csv_file, sep=delimiter)
        required_columns = ['x', 'y', 'z', 'Value']
        if not all(col in df.columns for col in required_columns):
            raise ValueError("CSV must contain columns: x, y, z, Value")

        # Create point cloud from coordinates
        points = df[['x', 'y', 'z']].values
        values = df['Value'].values

        # Create PyVista point cloud
        cloud = pv.PolyData(points)
        cloud['Values'] = values

        # Create discretized blocks
        print("\nStarting block creation...")
        blocks = create_blocks(points, values, block_size)
        print("Block creation complete")
        
        print("\nInitializing visualization...")
        # Create plotter
        plotter = pv.Plotter()
        plotter.set_background('white')
        print("Created plotter")
        
        # Add text elements with debug info
        print("Adding UI elements...")
        plotter.progress_text = plotter.add_text(
            text="Progress: 0% (waiting to start)",
            position='upper_left',
            font_size=12,
            font='arial'
        )
        # Add point info display
        plotter.point_info_text = plotter.add_text(
            text="Click a point to see info",
            position='lower_left',
            font_size=12,
            font='arial'
        )
        
        # Add grid toggle callback
        plotter.add_key_event('g', lambda: toggle_grid(plotter))
        
        # Add WASD navigation callbacks
        plotter.add_key_event('w', lambda: move_camera(plotter, 'w'))
        plotter.add_key_event('a', lambda: move_camera(plotter, 'a'))
        plotter.add_key_event('s', lambda: move_camera(plotter, 's'))
        plotter.add_key_event('d', lambda: move_camera(plotter, 'd'))
        
        # Initialize grid and store its actor
        plotter._grid_actor = plotter.show_bounds(grid='front', location='outer', color='darkblue')
        plotter._grid_actor.SetVisibility(False)  # Hide grid by default
        
        # Set text colors for all elements
        plotter.add_title("Point Cloud Visualization", color='darkblue', font='arial', font_size=8)
        
        # Setup picker and callbacks
        plotter.enable_point_picking(callback=middle_click_callback, 
                                   show_message=False,
                                   use_picker=True,
                                   left_clicking=False,
                                   show_point=True,
                                   tolerance=0.1)  # Increased tolerance for easier picking
        
        # Create picker for middle button
        plotter.enable_3_lights()  # Improve point visibility
        
        print("\nAdding visualization data...")
        # Add point cloud to plotter with scalar coloring
        plotter.add_mesh(cloud, 
                        render_points_as_spheres=True,
                        point_size=5,
                        scalars='Values',
                        cmap='rainbow',
                        show_scalar_bar=False)  # Hide individual scalar bar
        print("Added point cloud")

        # Add blocks to the scene (initially hidden)
        if blocks is not None:
            print(f"\nProcessing {len(blocks)} blocks for visualization...")
            plotter.blocks_data = blocks
            
            # Initialize list to store all block actors
            plotter.block_actors = []
            
            # Add each block individually for better control
            for i in range(len(blocks)):
                block = blocks[i]
                actor = plotter.add_mesh(
                    block,
                    style='surface',
                    scalars='Value',
                    cmap='rainbow',
                    show_scalar_bar=(i == 0),  # Only show scalar bar for first block
                    reset_camera=False,
                    pickable=True,
                    opacity=0.5,
                    show_edges=True
                )
                
                # Set VTK properties
                actor.GetProperty().SetEdgeColor(0.0, 0.0, 0.0)
                actor.GetProperty().SetAmbient(0.1)
                actor.GetProperty().SetDiffuse(0.9)
                actor.GetProperty().SetSpecular(0.1)
                
                # Store actor reference
                plotter.block_actors.append(actor)
            
            print("Added blocks to scene")
            print(f"Number of block actors: {len(plotter.block_actors)}")
            
            # Reset camera after adding all blocks
            plotter.renderer.ResetCamera()
            
            # Add callbacks with debug info
            def trigger_update():
                print("Update triggered by key press")
                update_blocks(plotter)
            
            def trigger_toggle():
                print("Block toggle triggered by key press")
                toggle_blocks(plotter)
            
            plotter.add_key_event('b', trigger_toggle)
            plotter.add_key_event('i', trigger_update)
        
        # Add axes with specific text color
        plotter.add_axes(color='darkblue')
        
        # Show a single scalar bar for all elements in dark blue
        plotter.add_scalar_bar(title='Values', 
                             fmt='%.1f',
                             color='darkblue',
                             vertical=True,
                             position_x=0.93,
                             position_y=0.2,
                             height=0.6,
                             width=0.08,
                             title_font_size=12,
                             label_font_size=12,
                             font_family='arial')

        # Enable camera controls
        plotter.camera.zoom(1.5)
        
        print("\nFinishing visualization setup...")
        print("Starting interactive viewer...")
        plotter.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    # Example usage
    csv_file = "C:\Projects\Anterpolator\Samples_CNN_ABA_BED.csv"  # Replace with your CSV file path
    block_size=100
    load_and_visualize_samples(csv_file, block_size=block_size)