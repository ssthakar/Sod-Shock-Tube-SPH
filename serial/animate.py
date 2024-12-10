import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import glob
import os

# Set up monochrome style
plt.style.use('grayscale')
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = 'white'

# Create figure and subplots
fig, axs = plt.subplots(4, 1, figsize=(10, 10))
fig.suptitle('SPH Simulation Results', fontsize=14, color='black')

# Labels for each subplot
ylabels = ['Internal Energy', 'Density', 'Pressure', 'Velocity']

# Initialize empty lists to store lines
lines = []
for ax in axs:
    # Black dots with white edge for better visibility
    line = ax.plot([], [], 'k.', markersize=2, markeredgecolor='black')[0]
    lines.append(line)
    ax.set_xlabel('Position', color='black')
    ax.grid(True, color='gray', linestyle='--', alpha=0.3)
    # Set tick colors to black
    ax.tick_params(colors='black')
    # Set spine colors to black
    for spine in ax.spines.values():
        spine.set_color('black')

# Set y-axis labels
for ax, label in zip(axs, ylabels):
    ax.set_ylabel(label, color='black')

def init():
    """Initialize animation"""
    for line in lines:
        line.set_data([], [])
    return lines

def animate(frame):
    """Animation function"""
    filename = f'./output/data_{frame}.dat'
    try:
        data = np.loadtxt(filename)
    except:
        print(f"Could not read file: {filename}")
        return lines
    
    positions = data[:, 0]
    
    for i, (line, ax) in enumerate(zip(lines, axs)):
        line.set_data(positions, data[:, i+1])
        
        ax.set_xlim(-0.6, 0.6)
        
        ymin = min(data[:, i+1])
        ymax = max(data[:, i+1])
        padding = (ymax - ymin) * 0.1
        ax.set_ylim(ymin - padding, ymax + padding)
    
    fig.suptitle(f'SPH Simulation Results (Frame {frame})', fontsize=14, color='black')
    return lines

# Get number of frames
data_files = glob.glob('./output/data_*.dat')
num_frames = len(data_files)

if num_frames == 0:
    print("No data files found in ./serial/ directory!")
    exit()

# Create animation
anim = FuncAnimation(fig, animate, frames=num_frames,
                    init_func=init, blit=True,
                    interval=50)

plt.tight_layout()

# Save the animation
print("Saving animation...")
writer = PillowWriter(fps=20)
anim.save('sph_simulation.gif', writer=writer)
print("Animation saved as 'sph_simulation.gif'")

plt.show()
