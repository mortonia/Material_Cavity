import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

def visualize_displacements(before_xyz, after_xyz, displacement_threshold=0.1):
    atoms_before = read(before_xyz)
    atoms_after = read(after_xyz)

    pos_before = atoms_before.get_positions()
    pos_after = atoms_after.get_positions()

    displacements = np.linalg.norm(pos_after - pos_before, axis=1)

    # Plot setup
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter points colored by displacement magnitude
    sc = ax.scatter(pos_after[:, 0], pos_after[:, 1], pos_after[:, 2],
                    c=displacements, cmap='inferno', s=50, alpha=0.8)

    # Add colorbar
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label('Displacement (Å)', fontsize=12)

    # Highlight atoms displaced beyond threshold with black edge
    displaced_mask = displacements > displacement_threshold
    ax.scatter(pos_after[displaced_mask, 0], pos_after[displaced_mask, 1], pos_after[displaced_mask, 2],
               facecolors='none', edgecolors='black', s=100, linewidths=1.5, label='Displaced > threshold')

    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Atomic Displacements After LJ Relaxation')

    ax.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python visualize_displacements.py before_lj.xyz after_lj.xyz [displacement_threshold]")
        sys.exit(1)

    before_file = sys.argv[1]
    after_file = sys.argv[2]
    threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1

    visualize_displacements(before_file, after_file, threshold)

