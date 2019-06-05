#!/usr/bin/python

__author__ = 'Ralf Hauenschild, Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '1.0'

################################################
# Library imports
################################################
import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import matplotlib
matplotlib.use('PS')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as grid_spec
import numpy as np
import pylab as py
from matplotlib.ticker import ScalarFormatter

table_paths = []
table_paths.append(sys.argv[1])

visual_range = 5
threshold = 0.1
target_parameter = 'mm_patterns'
aspect_ratio = '1500_800'
epsilon = 0.1
show_base = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
show_arrest = True
show_sequence = True
show_legend = True
topx = 0

width = 1.0 * int((aspect_ratio.split('_'))[0])
height = 1.0 * int((aspect_ratio.split('_'))[1])

x_height = 1.0 * height / width
x_width = 1.0

start = int(sys.argv[3])
end = int(sys.argv[4])

difference_analysis = False

middle_index = visual_range

line_list = []

colors_dict = {'K': 'k', 'R': 'k', 'Y': 'k', 'S': 'k', 'A': 'g', 'G': 'orange', 'C': 'b', 'T': 'r', 't': 'r', 'N': 'k',
               '': 'k', '_': 'k', '.': 'k', 'X': 'k'}

fig = py.figure(0, figsize=(2 * (width / 300.0), 2 * (height / 300.0)), dpi=800)
py.rc('font', size=12)

r = 0.92 - min(max(87016 * width ** (-2.083), 0.0), 0.17)
l = 0.05 + min(max(87016 * width ** (-2.083), 0.0), 0.17)
t = 1.0
b = 0.1 + min(max(18397 * height ** (-2.071), 0.0), 0.11)

py.subplots_adjust(hspace=0.1, bottom=b, top=t, left=l, right=r)

# ax = fig.add_subplot(1,1,1) # first two digits mean 2*2 = 4 subplot cells, last digit is index. Index may not be > #cells
gs = grid_spec.GridSpec(2, 1, height_ratios=[2.8, 1 + 0.015 * height])
legend_ax = plt.subplot(gs[0])
legend_ax.axis('off')
ax = plt.subplot(gs[1])
ax.locator_params(axis='x', tight=True, nbins=int((width / 100.0)))
ax.locator_params(axis='y', tight=True, nbins=6)
for axis in [ax.yaxis]:
    axis.set_major_formatter(ScalarFormatter())

ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 4))
ax.ticklabel_format(axis='x', style='sci', scilimits=(-2, 4))

prop_dicts = {}

skip_lines = start - visual_range - 5

# Parse positional data table
for t in range(0, len(table_paths)):

    prop_dicts[t] = {'pos': [], 'ref_base': [], 'cov': [], 'a_mism': [], 'g_mism': [], 't_mism': [], 'c_mism': [],
                     'mm_ratio': [], 'arrest_rate': []}

    infile = open(table_paths[t], 'r')

    line = infile.readline()
    line = infile.readline()

    while len(line) > 0 and int((line.split('\t'))[0]) <= skip_lines:
        line = infile.readline()
    pos = 0

    while (pos <= visual_range * 2) and len(line) > 0:
        pos += 1

        line_list.append((line[:-1]).split('\t'))
        if len(line_list) > visual_range * 2 + 1:
            del line_list[0]

        splitlist = (line[:-1]).split('\t')
        if pos <= middle_index + 1:
            # table format:
            # position    refbase    coverage    mm_ratio    a_mism    g_mism    t_mism    c_mism    arrest_rate
            # 1            G            16           0.0        0        16        0        0        0.36

            ((prop_dicts[t])['pos']).append(int(splitlist[0]))
            ((prop_dicts[t])['ref_base']).append(splitlist[1])
            ((prop_dicts[t])['cov']).append(int(splitlist[2]))
            ((prop_dicts[t])['mm_ratio']).append(float(splitlist[3]))
            ((prop_dicts[t])['a_mism']).append(int(splitlist[4]))
            ((prop_dicts[t])['g_mism']).append(int(splitlist[5]))
            ((prop_dicts[t])['t_mism']).append(int(splitlist[6]))
            ((prop_dicts[t])['c_mism']).append(int(splitlist[7]))
            ((prop_dicts[t])['arrest_rate']).append(float(splitlist[8]))

        line = infile.readline()

    while len(line) > 0 and int(splitlist[0]) <= end:

        splitlist = (line[:-1]).split("\t")

        # table format:
        # position    refbase    coverage    mm_ratio    a_mism    g_mism    t_mism    c_mism    arrest_rate
        # 1            G            16           0.0        0        16        0        0        0.36

        del line_list[0]
        line_list.append(splitlist)

        ((prop_dicts[t])['pos']).append(int((line_list[middle_index])[0]))
        ((prop_dicts[t])['ref_base']).append((line_list[middle_index])[1])
        ((prop_dicts[t])['cov']).append(int((line_list[middle_index])[2]))
        ((prop_dicts[t])['mm_ratio']).append(float((line_list[middle_index])[3]))
        ((prop_dicts[t])['a_mism']).append(int((line_list[middle_index])[4]))
        ((prop_dicts[t])['g_mism']).append(int((line_list[middle_index])[5]))
        ((prop_dicts[t])['t_mism']).append(int((line_list[middle_index])[6]))
        ((prop_dicts[t])['c_mism']).append(int((line_list[middle_index])[7]))
        ((prop_dicts[t])['arrest_rate']).append(float((line_list[middle_index])[8]))

        arrestrate = float((line_list[middle_index])[8])
        arrestratemedian = np.median(
            [float(line2[8]) for line2 in line_list[:middle_index]] + [float(line2[8]) for line2 in
                                                                       line_list[middle_index + 1:]])

        line = infile.readline()

    while middle_index < len(line_list) - 1:
        del line_list[0]

        ((prop_dicts[t])['pos']).append(int((line_list[middle_index])[0]))
        ((prop_dicts[t])['ref_base']).append((line_list[middle_index])[1])
        ((prop_dicts[t])['cov']).append(int((line_list[middle_index])[2]))
        ((prop_dicts[t])['mm_ratio']).append(float((line_list[middle_index])[3]))
        ((prop_dicts[t])['a_mism']).append(int((line_list[middle_index])[4]))
        ((prop_dicts[t])['g_mism']).append(int((line_list[middle_index])[5]))
        ((prop_dicts[t])['t_mism']).append(int((line_list[middle_index])[6]))
        ((prop_dicts[t])['c_mism']).append(int((line_list[middle_index])[7]))
        ((prop_dicts[t])['arrest_rate']).append(float((line_list[middle_index])[8]))
    infile.close()

prop_dict = prop_dicts[0]
pos_list = [pos for pos in (prop_dicts[0])['pos']]

markers = []
marker_values = []


maximum = 0

saturation = 10.0
norm1 = matplotlib.colors.Normalize(vmin=1.0, vmax=saturation)

bottomlevel = 0.7

color_dict1 = {'red': ((0.0, bottomlevel, bottomlevel),
                       (1.0, 1.0, 1.0)),
          'green': ((0.0, bottomlevel, bottomlevel),
                    (1.0, 0.0, 0.0)),
          'blue': ((0.0, bottomlevel, bottomlevel),
                   (1.0, 0.0, 0.0))
               }

cmap1 = matplotlib.colors.LinearSegmentedColormap('graurot', color_dict1)

m1 = cm.ScalarMappable(norm=norm1, cmap=cmap1)

color_dict2 = {'red': ((0.0, bottomlevel, bottomlevel),
                       (1.0, 0.0, 0.0)),
          'green': ((0.0, bottomlevel, bottomlevel),
                    (1.0, 1.0, 1.0)),
          'blue': ((0.0, bottomlevel, bottomlevel),
                   (1.0, 0.0, 0.0))
               }

cmap2 = matplotlib.colors.LinearSegmentedColormap('graugruen', color_dict2)

norm2 = matplotlib.colors.Normalize(vmin=1.0, vmax=saturation)
m2 = cm.ScalarMappable(norm=norm2, cmap=cmap2)

cmap3 = cm.coolwarm

norm3 = matplotlib.colors.LogNorm(vmin=0.1, vmax=saturation)
m3 = cm.ScalarMappable(norm=norm3, cmap=cm.coolwarm)

if end - start <= 1000:
    if target_parameter == 'mm_patterns':
        mism_comps = ['a_mism', 'g_mism', 't_mism', 'c_mism']

        for i in range(pos_list.index(start) - 1, min(pos_list.index(end) + 1, len(pos_list))):
            srhbar_1 = ax.bar(pos_list[i], (prop_dict['cov'])[i], width=1.0, bottom=0.0, color='lightgrey',
                              edgecolor='none')

        if topx == 0:  # report all
            for i in range(pos_list.index(start) - 1,
                           min(pos_list.index(end) + 1, len(pos_list))):

                if (prop_dict['mm_ratio'])[i] >= threshold and (prop_dict['ref_base'])[i] in show_base:
                    markers.append(pos_list[i])
                    marker_values.append((prop_dict['mm_ratio'])[i])

                    last_bottom = 0
                    n_height = (prop_dict['a_mism'])[i]
                    a_height = n_height
                    bar = ax.barh(last_bottom, 1, left=pos_list[i] - 0.5, align='edge', height=n_height, color='g',
                                  edgecolor='none')
                    last_bottom += n_height

                    n_height = (prop_dict['g_mism'])[i]
                    bar = ax.barh(last_bottom, 1, left=pos_list[i] - 0.5, align='edge', height=n_height, color='orange',
                                  edgecolor='none')
                    last_bottom += n_height

                    n_height = (prop_dict['t_mism'])[i]
                    bar = ax.barh(last_bottom, 1, left=pos_list[i] - 0.5, align='edge', height=n_height, color='r',
                                  edgecolor='none')
                    last_bottom += n_height

                    n_height = (prop_dict['c_mism'])[i]
                    bar = ax.barh(last_bottom, 1, left=pos_list[i] - 0.5, align='edge', height=n_height, color='b',
                                  edgecolor='none')
                    last_bottom += n_height

            py.ylabel('coverage')
            plt.ylim(ymin=max(0.0, threshold))
	    plt.ticklabel_format(useOffset=False, style='plain')
            ax2 = ax.twinx()

            ax2.tick_params(axis='y', which='major')

            if show_arrest:
                arrest_rate, = ax2.plot(pos_list, prop_dict['arrest_rate'], 'white', lw=3)
                arrest_rate, = ax2.plot(pos_list, prop_dict['arrest_rate'], 'red')
            ax2.set_ylim(0.0, 1.0)

            mma = ax2.scatter([1000000], [10000000], marker='s', color='green', edgecolor='none', s=50, zorder=10,
                              label='A', alpha=1.0)
            mmg = ax2.scatter([1000000], [10000000], marker='s', color='orange', edgecolor='none', s=50, zorder=10,
                              label='G', alpha=1.0)
            mmt = ax2.scatter([1000000], [10000000], marker='s', color='red', edgecolor='none', s=50, zorder=10,
                              label='T', alpha=1.0)
            mmc = ax2.scatter([1000000], [10000000], marker='s', color='blue', edgecolor='none', s=50, zorder=10,
                              label='C', alpha=1.0)
            ngsreads = ax2.scatter([1000000], [10000000], marker='s', color='lightgrey', edgecolor='lightgrey', s=50,
                                   zorder=10, label='C', alpha=1.0)
            ar = ax2.scatter([1000000], [10000000], marker='_', color='red', edgecolor='red', linewidth=1.5, s=50,
                             zorder=10, label='Arrest rate', alpha=1.0)
            legendevents = ax2.scatter([1000000], [10000000], marker='v', color='yellow', edgecolor='black', s=50,
                                       zorder=10, label='Above-threshold events')

            ax2.set_ylim(0.0, 1.0)

            py.ylabel('mism. rate,\narrest rate')

            plt.ylim(ymin=0, ymax=1.0)

            maximum = max(
                (prop_dict['cov'])[pos_list.index(start):pos_list.index(end) + 1 - 1])
            print(maximum)
            if topx == 0:
                events = ax.scatter(markers, [maximum for marker in markers], marker='v', color='yellow',
                                    edgecolor='black', s=100, zorder=10, label='Above-threshold events')
                mismatches = ax2.scatter(markers, marker_values, marker='x', color='black', edgecolor='black', s=100,
                                         zorder=10, label='Mismatch ratio')

            if show_legend:
                legend_ax.legend((mma, mmg, mmt, mmc, legendevents, mismatches, ngsreads, ar), (
                'A', 'G', 'T', 'C', 'Mism.' + r'$\geq$' + " " + str(threshold), 'Mismatch rate', 'Reads',
                'Arrest rate'), loc='upper center', bbox_to_anchor=(0.5, 0.95), bbox_transform=plt.gcf().transFigure,
                                 scatterpoints=1, ncol=4, frameon=False, labelspacing=0.2, columnspacing=0.8,
                                 handlelength=0.4, prop={'size': 11})

    else:
        pass

    scope = min(pos_list.index(end), len(pos_list)) - pos_list.index(start) + 1
    if show_sequence:
        if (0.1 + (1.0 - min(1.0, scope / 120.0))) * (0.1 + min(1.0, (width / 1000.0))) >= 0.5:
            for i in range(0, scope):
                ax.annotate((prop_dict['ref_base'])[pos_list.index(start) + i],
                            xy=((float(i) / (scope - 1)), 0.03), xycoords=('axes fraction', 'figure fraction'),
                            color=colors_dict[(prop_dict['ref_base'])[pos_list.index(start) + i]], ha='center')
        else:
            py.xlabel('position')
    else:
        py.xlabel('position')
else:  # low detail
    x = [pos_list[i] for i in range(0, len(pos_list)) if pos_list[i] % (100) == 1]

    maximum = max(
        (prop_dict['cov'])[pos_list.index(start):pos_list.index(end) + 1 - 1])

    y = np.row_stack(([(prop_dict['cov'])[i] for i in range(0, len(pos_list)) if pos_list[i] % (100) == 1],))

    ax.stackplot(x, y, edgecolor='none')



    py.ylabel('coverage')
    py.title('Coverage profile')


ax.xaxis.grid('minor', linestyle='dotted', linewidth=0.5)
ax.xaxis.grid('major', linestyle='dotted', linewidth=0.5)
ax.yaxis.grid('major', linestyle='dotted', linewidth=0.5)
ax.yaxis.grid('minor', linestyle='dotted', linewidth=0.5)

py.xlim(min(start, pos_list[-1]), min(end, pos_list[-1]))
py.xlabel('position')
py.savefig(sys.argv[2], dpi=800, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches='tight', pad_inches=0.1)

py.close()

'''
/home/akhelm/Lukas/Ergebnisse/CovAn_Test/positional_tables/test/my_Table_test.txt
/home/akhelm/Lukas/Ergebnisse/CovAn_Test/positional_tables/test/myplot
50
70
'''

# The End

