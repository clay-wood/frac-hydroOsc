from matplotlib.lines import Line2D

def legend_stuff(color1, color2):
    legend_elements1N = [Line2D([0], [0], marker='o', color = 'w', label='T3R3', markerfacecolor=color1, mec = 'w', mew = 0, ms=10),
                          Line2D([0], [0], marker='o', color = 'w', label='T5R5', markerfacecolor=color2, mec = 'w', mew = 0, ms=10)]

    legend_elements1P = [Line2D([0], [0], marker='o', color = 'w', label='T3R3', markerfacecolor=color1, mec = 'w', mew = 0, ms=10),
                          Line2D([0], [0], marker='o', color = 'w', label='T5R5', markerfacecolor=color2, mec = 'w', mew = 0, ms=10)]

    legend_elements2N = [Line2D([0], [0], marker='d', color = 'w', label='0.1 Hz', markerfacecolor='w', mec = 'k', mew = 1, ms=10), 
                          Line2D([0], [0], marker='o', color = 'w', label='1.0 Hz', markerfacecolor='w', mec = 'k', mew = 1, ms=10),
                          Line2D([0], [0], marker='*', color = 'w', label='10 Hz', markerfacecolor='w', mec = 'k', mew = 1, ms=12),
                          Line2D([0], [0], marker='s', color = 'w', label='40 Hz', markerfacecolor='w', mec = 'k', mew = 1, ms=9)]

    legend_elements2P = [Line2D([0], [0], marker='o', color = 'w', label='Set 1', markerfacecolor='w', mec = 'k', mew = 1, ms=10), 
                          Line2D([0], [0], marker='s', color = 'w', label='Set 2', markerfacecolor='w', mec = 'k', mew = 1, ms=10)]

    legend_elements2Q = [Line2D([0], [0], marker='o', color = 'w', label='Set 1', markerfacecolor='w', mec = 'k', mew = 1, ms=10), 
                          Line2D([0], [0], marker='s', color = 'w', label='Set 2', markerfacecolor='w', mec = 'k', mew = 1, ms=9),
                          Line2D([0], [0], marker='d', color = 'w', label='Set 3', markerfacecolor='w', mec = 'k', mew = 1, ms=10),
                          Line2D([0], [0], marker='*', color = 'w', label='Set 4', markerfacecolor='w', mec = 'k', mew = 1, ms=12)]

    legend_elements2R = [Line2D([0], [0], marker='o', color = 'w', label='Set 1', markerfacecolor='w', mec = 'k', mew = 1, ms=10), 
                          Line2D([0], [0], marker='s', color = 'w', label='Set 2', markerfacecolor='w', mec = 'k', mew = 1, ms=9),
                          Line2D([0], [0], marker='d', color = 'w', label='Set 3', markerfacecolor='w', mec = 'k', mew = 1, ms=10)]
    
    return legend_elements1N, legend_elements1P, legend_elements2N, legend_elements2P, legend_elements2Q, legend_elements2R