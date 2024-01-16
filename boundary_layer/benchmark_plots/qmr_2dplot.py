import matplotlib.pyplot as plt
import numpy as np

functions_selected = ['coeff', 'solve_2d', 'conv']


cpu_cumtime_size1 = [0.376, 10.778, 0.197]
cpu_cumtime_size5 = [32.576, 191.601, 13.353]
cpu_cumtime_size10 = [136.959, 984.549, 56.651]

bar_width = 0.2 
index = np.arange(len(functions_selected))  

fig, ax = plt.subplots()

ax.bar(index - bar_width, cpu_cumtime_size1, width=bar_width, label='Size 1', alpha=0.7)
ax.bar(index, cpu_cumtime_size5, width=bar_width, label='Size 5', alpha=0.7)
ax.bar(index + bar_width, cpu_cumtime_size10, width=bar_width, label='Size 10', alpha=0.7)

ax.set_xlabel('Functions')
ax.set_ylabel('Cumulative Time (s)')
ax.set_title('CPU Cumtime for Selected Functions and Different Grid Sizes (qmr)')

ax.set_xticks(index)
ax.set_xticklabels(functions_selected)

ax.set_yscale('log')

ax.legend()

plt.savefig('qmr_2dplot.pdf', bbox_inches='tight')

plt.show()
