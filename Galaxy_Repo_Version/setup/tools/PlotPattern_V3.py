#!/usr/bin/python

__author__ = 'Ralf Hauenschild, Lukas Schmidt'
__maintainer__ = 'Lukas Schmidt'
__email__ = 'sclukas@students.uni-mainz.de'
__version__ = '2.0'

################################################
# Library imports
################################################
import sys
sys.path.append('/usr/local/lib/python2.7/dist-packages')
import matplotlib
matplotlib.use('PS')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as grid_spec
import numpy as np
import pylab as py
from matplotlib.ticker import ScalarFormatter


########################################################################################################################

def build_positional_table(profile):
    """
    This function writes the information for each position into a dictionary.
    :param profile: Profile file
    :return: Dictionary with all info
    """
    prop_dict = {'pos': [], 'ref_base': [], 'cov': [], 'mismatch_rate': [], 'a_mism': [], 'g_mism': [], 't_mism': [],
                 'c_mism': [], 'arrest_rate': []}

    ref = sys.argv[3]
    print(ref.replace('__tt__', '|'))
    for line in profile:
        line1 = line.strip().split()
        if line1[0] == ref.replace('__tt__', '|') and start <= int(line1[1]) <= end:
            prop_dict['pos'].append(int(line1[1]))
            prop_dict['ref_base'].append(line1[2])
            prop_dict['cov'].append(int(line1[3]))
            prop_dict['mismatch_rate'].append(float(line1[5]))
            prop_dict['a_mism'].append(int(line1[6]) + int(line1[11]))
            prop_dict['g_mism'].append(int(line1[7]) + int(line1[12]))
            prop_dict['t_mism'].append(int(line1[8]) + int(line1[13]))
            prop_dict['c_mism'].append(int(line1[9]) + int(line1[14]))
            prop_dict['arrest_rate'].append(float(line1[-1]))

    return prop_dict

########################################################################################################################


# input: Profile-file
start = int(sys.argv[4])
end = int(sys.argv[5])
with open(sys.argv[1], 'r') as fin:
    prop_dict = build_positional_table(fin)

font_size = int(sys.argv[8])
line_width = float(sys.argv[9])
marker_size = float(sys.argv[10])


visual_range = 5
threshold = 0.1
target_parameter = 'mm_patterns'
width = 1.0 * int(sys.argv[6])
height = 1.0 * int(sys.argv[7])

#aspect_ratio = aspect_x + '_' + '390_390'  # 1500_800
epsilon = 0.1
show_base = {'A': 1, 'C': 1, 'G': 1, 'T': 1}
show_arrest = True
show_sequence = True
show_legend = True if sys.argv[11] == 'yes' else False
padding = float(sys.argv[12])
topx = 0

line_list = []

colors_dict = {'K': 'k', 'R': 'k', 'Y': 'k', 'S': 'k', 'A': 'g', 'G': 'orange', 'C': 'b', 'T': 'r', 't': 'r', 'N': 'k',
               '': 'k', '_': 'k', '.': 'k', 'X': 'k'}

fig = py.figure(0, figsize=(2 * (width / 300.0), 2 * (height / 300.0)), dpi=800)
py.rc('font', size=font_size)  # 12


r = 0.92 - min(max(87016 * width ** (-2.083), 0.0), 0.17)
l = 0.05 + min(max(87016 * width ** (-2.083), 0.0), 0.17)
t = 0.90  # 1.0
b = 0.1 + min(max(18397 * height ** (-2.071), 0.0), 0.11)

py.subplots_adjust(hspace=0.1, bottom=b, top=t, left=l, right=r)

# ax = fig.add_subplot(2,2,1) # first two digits mean 2*2 = 4 subplot cells, last digit is index.
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


pos_list = [pos for pos in prop_dict['pos']]


markers = []
marker_values = []

if end - start <= 1000:
    if target_parameter == 'mm_patterns':
        mism_comps = ['a_mism', 'g_mism', 't_mism', 'c_mism']

        # Plot coverage bars
        for i in range(pos_list.index(start) - 1, min(pos_list.index(end) + 1, len(pos_list))):
            srhbar_1 = ax.bar(pos_list[i], (prop_dict['cov'])[i], width=1.0, bottom=0.0, color='lightgrey',
                              edgecolor='none')

        if topx == 0:  # report all
            for i in range(pos_list.index(start) - 1,
                           min(pos_list.index(end) + 1, len(pos_list))):

                # Plot the mismatch
                if (prop_dict['mismatch_rate'])[i] >= threshold and (prop_dict['ref_base'])[i] in show_base:
                    markers.append(pos_list[i])
                    marker_values.append((prop_dict['mismatch_rate'])[i])

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
                arrest_rate, = ax2.plot(pos_list, prop_dict['arrest_rate'], 'white', lw=line_width)  # lw=3
                arrest_rate, = ax2.plot(pos_list, prop_dict['arrest_rate'], 'red', lw=line_width * 0.5)
            ax2.set_ylim(0.0, 1.0)

            mma = ax2.scatter([1000000], [10000000], marker='s', color='green', edgecolor='none', s=marker_size,
                              zorder=10, label='A', alpha=1.0)  # s=50
            mmg = ax2.scatter([1000000], [10000000], marker='s', color='orange', edgecolor='none', s=marker_size,
                              zorder=10, label='G', alpha=1.0)
            mmt = ax2.scatter([1000000], [10000000], marker='s', color='red', edgecolor='none', s=marker_size,
                              zorder=10, label='T', alpha=1.0)
            mmc = ax2.scatter([1000000], [10000000], marker='s', color='blue', edgecolor='none', s=marker_size,
                              zorder=10, label='C', alpha=1.0)
            ngsreads = ax2.scatter([1000000], [10000000], marker='s', color='lightgrey', edgecolor='lightgrey',
                                    s=marker_size, zorder=10, label='C', alpha=1.0)
            ar = ax2.scatter([1000000], [10000000], marker='_', color='red', edgecolor='red', linewidth=line_width,
                             s=marker_size, zorder=10, label='Arrest rate', alpha=1.0)
            legendevents = ax2.scatter([1000000], [10000000], marker='v', color='yellow', edgecolor='black',
                                       s=marker_size, zorder=10, label='Above-threshold events')

            ax2.set_ylim(0.0, 1.0)

            py.ylabel('mism. rate,\narrest rate')

            plt.ylim(ymin=0, ymax=1.0)

            maximum = max(
                (prop_dict['cov'])[pos_list.index(start):pos_list.index(end) + 1 - 1])
            if topx == 0:
                events = ax.scatter(markers, [maximum for marker in markers], marker='v', color='yellow',
                                    edgecolor='black', s=marker_size * 0.75, zorder=10,
                                    label='Above-threshold events')  # s=30
                mismatches = ax2.scatter(markers, marker_values, marker='x', color='black', edgecolor='black',
                                         s=marker_size * 0.75, zorder=10, label='Mismatch ratio')

            if show_legend:
                legend_ax.legend((mma, mmg, mmt, mmc, legendevents, mismatches, ngsreads, ar), (
                'A', 'G', 'T', 'C', 'Mism.' + r'$\geq$' + " " + str(threshold), 'Mismatch rate', 'Reads',
                'Arrest rate'), loc='upper center', bbox_to_anchor=(0.5, 0.86), bbox_transform=plt.gcf().transFigure,
                                 scatterpoints=1, ncol=4, frameon=False, labelspacing=0.2, columnspacing=0.8,
                                 handlelength=0.4, prop={'size': font_size - 1})  # 11
    else:
        pass

    # Add sequence to graphic
    scope = min(pos_list.index(end), len(pos_list)) - pos_list.index(start) + 1
    
    '''
    if width > 450:
        padding = 0.03
    else:
	padding = 0.06
    '''
    if show_sequence:
        if (0.1 + (1.0 - min(1.0, scope / 120.0))) * (0.1 + min(1.0, (width / 1000.0))) >= 0.05:  # 0.5
            for i in range(0, scope):
                ax.annotate((prop_dict['ref_base'])[pos_list.index(start) + i],
                            xy=((float(i) / (scope - 1)), padding), xycoords=('axes fraction', 'figure fraction'),
                            color=colors_dict[(prop_dict['ref_base'])[pos_list.index(start) + i]], ha='center')
                #  xy=((float(i) / (scope - 1)), 0.03  # for smaller graphs 0.06
        else:
            py.xlabel('position')
    else:
        py.xlabel('position')
else:  # low detail
    x = [pos_list[i] for i in range(0, len(pos_list)) if pos_list[i] % 100 == 1]

    maximum = max(
        (prop_dict['cov'])[pos_list.index(start):pos_list.index(end) + 1 - 1])

    y = np.row_stack(([(prop_dict['cov'])[i] for i in range(0, len(pos_list)) if pos_list[i] % 100 == 1],))

    ax.stackplot(x, y, edgecolor='none')

    py.ylabel('coverage')
    py.title('Coverage profile')


# Add grid
ax.xaxis.grid('minor', linestyle='dotted', linewidth=0.5)
ax.xaxis.grid('major', linestyle='dotted', linewidth=0.5)
ax.yaxis.grid('major', linestyle='dotted', linewidth=0.5)
ax.yaxis.grid('minor', linestyle='dotted', linewidth=0.5)


py.xlim(min(start, pos_list[-1]), min(end, pos_list[-1]))
py.xlabel('position')
py.savefig(sys.argv[2], dpi=1200, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png',
           transparent=False, bbox_inches='tight', pad_inches=0.1)

py.close()

'''
/home/akhelm/Lukas/Ergebnisse/all_tRNAs/new_Workflow/Yeast/Raw/MH1512.profile
/home/akhelm/Downloads/myplot
tdbR00000369__tt__Saccharomyces_cerevisiae__tt__4932__tt__Arg__tt__ACG
50
60
390
180
5
2
30
no
0.06
'''

# The End

