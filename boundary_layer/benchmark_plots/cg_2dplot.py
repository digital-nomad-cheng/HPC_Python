import matplotlib.pyplot as plt

# coeff = Coefficients
# solve_2d = Solving Navier-Strokes (?)
# tocoo = to COOrdinate format
functions = ['coeff', 'solve_2d', 'tocoo']


cpu_cumtime_grid1 = [0.539, 33.452, 0.717]
cpu_cumtime_grid5 = [32.938, 184.976, 17.195]
cpu_cumtime_grid10 = [137.658, 1303.395, 73.657]

fig, ax = plt.subplots()

bar_width = 0.2
bar_positions_grid1 = range(len(functions))
bar_positions_grid5 = [pos + bar_width for pos in bar_positions_grid1]
bar_positions_grid10 = [pos + 2 * bar_width for pos in bar_positions_grid1]

ax.bar(bar_positions_grid1, cpu_cumtime_grid1, width=bar_width, label='Grid 1')
ax.bar(bar_positions_grid5, cpu_cumtime_grid5, width=bar_width, label='Grid 5')
ax.bar(bar_positions_grid10, cpu_cumtime_grid10, width=bar_width, label='Grid 10')

ax.set_xlabel('Functions')
ax.set_ylabel('Cumulative Time (s)')
ax.set_title('CPU Cumtime for Different Grid Sizes')

ax.set_xticks([pos + bar_width for pos in bar_positions_grid1])
ax.set_xticklabels(functions)
ax.set_yscale('log')


ax.legend()

plt.show()
