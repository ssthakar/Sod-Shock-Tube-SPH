import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import glob
import os


# Create figure and subplots
fig, axs = plt.subplots(5, 1, figsize=(10, 12))
fig.suptitle('SPH Simulation Results', fontsize=14)

# Labels for each subplot
ylabels = ['Position', 'Internal Energy', 'Density', 'Pressure', 'Velocity']

# Initialize empty lists to store lines
lines = []
for ax in axs:
    line = ax.plot([], [], 'b.', markersize=2)[0]
    lines.append(line)
    ax.set_xlabel('Position')
    ax.grid(True)

# Set y-axis labels
for ax, label in zip(axs, ylabels):
    ax.set_ylabel(label)

def init():
    """Initialize animation"""
    for line in lines:
        line.set_data([], [])
    return lines

def animate(frame):
    """Animation function"""
    # Read data file
    filename = f'./output/data_{frame}.dat'
    try:
        data = np.loadtxt(filename)
    except:
        print(f"Could not read file: {filename}")
        return lines
    
    # Position array (x-axis for all plots)
    positions = data[:, 0]
    
    # Update each subplot
    for i, (line, ax) in enumerate(zip(lines, axs)):
        line.set_data(positions, data[:, i])
        
        # Set x limits based on position domain
        ax.set_xlim(-0.6, 0.6)
        
        # Dynamically adjust y limits with some padding
        ymin = min(data[:, i])
        ymax = max(data[:, i])
        padding = (ymax - ymin) * 0.1
        ax.set_ylim(ymin - padding, ymax + padding)
    
    # Add frame number to title
    fig.suptitle(f'SPH Simulation Results (Frame {frame})', fontsize=14)
    return lines

# Get number of frames from available data files
data_files = glob.glob('./output/data_*.dat')
num_frames = len(data_files)

if num_frames == 0:
    print("No data files found in ./serial/ directory!")
    exit()

# Create animation
anim = FuncAnimation(fig, animate, frames=num_frames,
                    init_func=init, blit=True,
                    interval=50)  # 50ms between frames

# Adjust layout to prevent overlap
plt.tight_layout()

# Save the animation
print("Saving animation...")
writer = PillowWriter(fps=20)
anim.save('sph_simulation.gif', writer=writer)
print("Animation saved as 'sph_simulation.gif'")

# Also display the animation
plt.show()
