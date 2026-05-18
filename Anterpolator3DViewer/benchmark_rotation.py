import sys
import os
import numpy as np
import io
from math import degrees, atan2
from anterpolator3DViewer import detect_grid_rotation

def run_benchmark():
    grid_sizes = {'small': (40, 40, 6), 'large': (100, 100, 25)}
    angles = [0, 2, 4, 8, 15]
    sigmas = [0.0, 0.2, 0.4, 0.6, 1.0]
    seeds = range(20)
    spacing = 10.0

    print("grid,angle,sigma,detected_rate,median_abs_err,detected_count")

    results = []

    for size_name, (nx, ny, nz) in grid_sizes.items():
        for angle_deg in angles:
            for sigma in sigmas:
                detected_flags = []
                angle_errors = []
                
                for seed in seeds:
                    np.random.seed(seed)
                    
                    # Create grid
                    x = np.arange(nx) * spacing
                    y = np.arange(ny) * spacing
                    z = np.arange(nz) * spacing
                    xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')
                    points = np.stack([xv, yv, zv], axis=-1).reshape(-1, 3)
                    
                    # Rotate around center
                    center = np.mean(points, axis=0)
                    rad = np.radians(angle_deg)
                    c, s = np.cos(rad), np.sin(rad)
                    rmat = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
                    points = (points - center) @ rmat.T + center
                    
                    # Add noise
                    noise = np.random.normal(0, sigma, (points.shape[0], 2))
                    points[:, :2] += noise
                    
                    # Call detector (suppress stdout)
                    old_stdout = sys.stdout
                    sys.stdout = io.StringIO()
                    try:
                        rot_mat, rot_center, is_rotated = detect_grid_rotation(points, block_size_hint=(10, 10, 10))
                    finally:
                        sys.stdout = old_stdout
                        
                    detected_flags.append(is_rotated)
                    if is_rotated:
                        # Angle estimation from R: [cos -sin; sin cos] or [sin cos; -cos sin]? 
                        # detect_grid_rotation uses axes. We expect the matrix to represent the transformation.
                        # The request says: theta_deg = degrees(arctan2(R[1,1], R[1,0])) - 90
                        # R[1,0] is sin, R[1,1] is cos. atan2(cos, sin) = 90 - atan2(sin, cos).
                        # So atan2(R[1,1], R[1,0]) - 90 is effectively -theta if theta is rotation.
                        # Let's verify the wrap.
                        est_theta = degrees(atan2(rot_mat[1,1], rot_mat[1,0])) - 90
                        est_theta = (est_theta + 180) % 360 - 180
                        angle_errors.append(abs(est_theta - angle_deg))
                
                det_rate = np.mean(detected_flags)
                med_err = np.median(angle_errors) if angle_errors else float('nan')
                det_count = sum(detected_flags)
                print(f"{size_name},{angle_deg},{sigma},{det_rate:.2f},{med_err:.4f},{det_count}")
                results.append({
                    'size': size_name, 'angle': angle_deg, 'sigma': sigma, 
                    'det_rate': det_rate, 'med_err': med_err, 'det_count': det_count
                })

    # Summary prints
    print("\nSummary:")
    for sigma in sigmas:
        fpr = np.mean([r['det_rate'] for r in results if r['angle'] == 0 and r['sigma'] == sigma])
        print(f"FPR at angle=0, sigma={sigma}: {fpr:.2f}")
    
    for sigma in sigmas:
        tpr = np.mean([r['det_rate'] for r in results if r['angle'] >= 8 and r['sigma'] == sigma])
        print(f"TPR at angle>=8, sigma={sigma}: {tpr:.2f}")

    print("\nNotes on ambiguous range (angle=2, 4):")
    for a in [2, 4]:
        rates = [r['det_rate'] for r in results if r['angle'] == a]
        print(f"Angle {a} detection rates: {min(rates):.2f} to {max(rates):.2f}")

if __name__ == '__main__':
    run_benchmark()
