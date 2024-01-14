import matplotlib.pyplot as plt
import numpy as np

cells = np.array([10, 20, 30, 40])
cpu_times = np.array([9.038, 101.077, 262.092, 560.075])
gpu_times = np.array([3.094, 4.659, 11.954, 22.240])

bar_width = 0.35
bar_positions_cpu = cells - bar_width / 2
bar_positions_gpu = cells + bar_width / 2

plt.bar(bar_positions_cpu, cpu_times, width=bar_width, label='CPU', log=True)
plt.bar(bar_positions_gpu, gpu_times, width=bar_width, label='GPU', log=True)

plt.xlabel('Number of Cells')
plt.ylabel('Time (seconds, log scale)')
plt.title('Comparison of CPU and GPU Performance (log scale)')
plt.legend()

plt.savefig('lorena.png', bbox_inches='tight')

plt.show()
