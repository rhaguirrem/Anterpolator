import pandas as pd
import pyvista as pv
import numpy as np
import csv
from tqdm import tqdm

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
    
    # Create multiblock
    return pv.MultiBlock(block_data)

def toggle_blocks(plotter):
    if hasattr(plotter, '_blocks_actor'):
        is_visible = not plotter._blocks_actor.GetVisibility()
        plotter._blocks_actor.SetVisibility(is_visible)
        plotter.render()

def load_and_visualize_samples(csv_file, block_size=10):
    try:
        # Read CSV
        df = pd.read_csv(csv_file)
        points = df[['x', 'y', 'z']].values
        values = df['Value'].values

        # Create blocks
        blocks = create_blocks(points, values, block_size)
        
        # Setup visualization
        plotter = pv.Plotter()
        plotter.set_background('white')
        
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
        
        # Show visualization
        plotter.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    csv_file = "C:\Projects\Anterpolator\Samples_CNN_ABA_BED.csv"
    block_size = 20
    load_and_visualize_samples(csv_file, block_size=block_size)