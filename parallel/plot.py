import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import re

# Function to extract number from filename
def get_frame_number(filename):
    return int(re.findall(r'data_(\d+)\.dat', filename)[0])

# Get all data files and sort them by frame number
data_files = sorted(glob.glob('data_*.dat'), key=get_frame_number)

# Read first file to get number of particles
first_frame = np.loadtxt(data_files[0])
num_particles = len(first_frame)

# Initialize figure and axis with black background
plt.style.use('dark_background')
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 14), 
                                        gridspec_kw={'height_ratios': [1, 1, 1, 1],
                                                    'hspace': 0.3})
fig.patch.set_facecolor('black')

# Set up each subplot
for ax in [ax1, ax2, ax3, ax4]:
    ax.set_facecolor('black')
    ax.grid(True, color='gray', alpha=0.3)
    ax.set_xlim(-0.6, 0.6)  # Updated x-axis limits
    ax.tick_params(colors='white', which='both')  # Make ticks white
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')

# Set y-limits based on the physics
ax1.set_ylim(0, 3)      # Energy
ax2.set_ylim(0, 1.2)    # Density
ax3.set_ylim(0, 1.2)    # Pressure
ax4.set_ylim(-0.5, 0.5) # Velocity

# Labels with larger font size
fontsize = 12
ax1.set_ylabel('Energy', color='white', fontsize=fontsize)
ax2.set_ylabel('Density', color='white', fontsize=fontsize)
ax3.set_ylabel('Pressure', color='white', fontsize=fontsize)
ax4.set_ylabel('Velocity', color='white', fontsize=fontsize)
ax4.set_xlabel('Position', color='white', fontsize=fontsize)

# Initialize plots with empty data and larger markers
line1, = ax1.plot([], [], 'w.', markersize=2, alpha=0.8)
line2, = ax2.plot([], [], 'c.', markersize=2, alpha=0.8)  # Cyan for density
line3, = ax3.plot([], [], 'r.', markersize=2, alpha=0.8)  # Red for pressure
line4, = ax4.plot([], [], 'y.', markersize=2, alpha=0.8)  # Yellow for velocity

# Set title with white text and larger font
title = fig.suptitle('Time: 0.000', color='white', fontsize=14, y=0.95)

def init():
    """Initialize animation"""
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return line1, line2, line3, line4, title

def animate(frame):
    """Animation function"""
    # Read data file
    data = np.loadtxt(data_files[frame])
    
    # Extract columns
    x = data[:, 0] - 0.6  # Shift x values to center the plot
    e = data[:, 1]      # Energy
    rho = data[:, 2]    # Density
    P = data[:, 3]      # Pressure
    v = data[:, 4]      # Velocity
    
    # Update data
    line1.set_data(x, e)
    line2.set_data(x, rho)
    line3.set_data(x, P)
    line4.set_data(x, v)
    
    # Update title with current time
    time = frame * 0.005  # dt = 0.005 from original code
    title.set_text(f'Time: {time:.3f}')
    
    return line1, line2, line3, line4, title

# Create animation
anim = animation.FuncAnimation(fig, animate, init_func=init,
                             frames=len(data_files), 
                             interval=50, blit=True)

# Save animation with high quality
anim.save('sph_simulation.mp4', fps=20, 
          extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'],
          savefig_kwargs={'facecolor': 'black'})

plt.close()
