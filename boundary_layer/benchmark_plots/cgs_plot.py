import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function names
functions = ['iter', 'coeff', 'solve_2d', 'conv', 'dphidy']

cpu_tottime_grid1 = [5.635, 0.310, 0.161, 0.092, 0.109]
gpu_tottime_grid1 = [6.595, 0.288, 0.172, 0.099, 0.100]
cpu_tottime_grid5 = [37.571, 28.532, 9.762, 7.486, 7.493]
gpu_tottime_grid5 = [2.764, 0.436, 0.441, 0.146, 0.087]
cpu_tottime_grid10 = [1306.906, 118.773, 38.007, 30.872, 29.660]
gpu_tottime_grid10 = [69.800, 1.502, 94.464, 0.861, 0.133]

# Set a maximum value for the z-axis to handle outliers
z_max = 100

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


# Scatter plot for CPU values for grid 10
ax.scatter([1] * len(functions), [i + 2 for i in range(len(functions))], cpu_tottime_grid1, marker='^', facecolors='none', edgecolors='g', label='CPU_1')

# Scatter plot for GPU values for grid 10
ax.scatter([2] * len(functions), [i + 2 for i in range(len(functions))], gpu_tottime_grid1, marker='^', facecolors='c', edgecolors='c', label='GPU_1')

# Scatter plot for CPU values for grid 5
ax.scatter([3] * len(functions), [i + 1 for i in range(len(functions))], cpu_tottime_grid5, marker='o', facecolors='none', edgecolors='r', label='CPU_5')

# Scatter plot for GPU values for grid 5
ax.scatter([4] * len(functions), [i + 1 for i in range(len(functions))], gpu_tottime_grid5, marker='o', facecolors='b', edgecolors='b', label='GPU_5')

# Scatter plot for CPU values for grid 10
ax.scatter([5] * len(functions), [i + 1.5 for i in range(len(functions))], cpu_tottime_grid10, marker='s', facecolors='none', edgecolors='g', label='CPU_10')

# Scatter plot for GPU values for grid 10
ax.scatter([6] * len(functions), [i + 1.5 for i in range(len(functions))], gpu_tottime_grid10, marker='s', facecolors='c', edgecolors='c', label='GPU_10')

# Set labels and title with smaller font size
ax.set_xlabel('Grid Size', fontsize=10)
ax.set_ylabel('Function', fontsize=10)
ax.set_zlabel('Total Time (s)', fontsize=10)
ax.set_title('Exec time comparison of functions between CPU and GPU (log)', fontsize=12)

# Set ticks and labels with smaller font size
ax.set_xticks([1, 2, 3, 4, 5, 6])
ax.set_xticklabels(['CPU1', 'GPU1', 'CPU5', 'GPU5', 'CPU10', 'GPU10'], fontsize=8)
ax.set_yticks([1, 2, 3, 4, 5])
ax.set_yticklabels(functions, fontsize=8)

# Add legend with smaller font size
ax.legend(fontsize=8)

# Set z-axis limits to handle outliers
ax.set_zlim(0, z_max)

plt.show()
