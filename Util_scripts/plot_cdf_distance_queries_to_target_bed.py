#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import numpy as np
from pybedtools import BedTool
import pybedtools
from scipy.stats import mannwhitneyu
from statsmodels.distributions.empirical_distribution import ECDF
from reddylab_utils.reddylab_plotting_utils import discrete_cmap

font = {'size' : 8}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import palettable.colorbrewer.sequential

import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
description="""

plot_cdf_distance_queries_to_target_bed.py

This script plots the cumulative distribution function (CDF) of the distances
of all sites in each of a list of query bed files to a target bed file.

""")

##################################################
# required args:

parser.add_argument("--queries", nargs='+', type=str,
                    help="""required, file paths to beds of query sites of interest, e.g.:
/path/to/site_group1.bed /path/to/site_group2.bed etc.
""", dest="queries", required=True)

parser.add_argument("--names", nargs='+', type=str,
                    help="""required, list of names corresponding to query bed files, e.g.:
"Group 1" "Group 2" etc.

""", required=True)

parser.add_argument("--target", type=str,
                    help="""required, file path to bed of target sites:
""", required=True)

parser.add_argument("-o", '--outplot', dest='outplot', type=str,
                    help='required, outplot', required=True)

##################################################
# optional args:

parser.add_argument("--min_log_dist", type=float, default=0.,
help="min_log_dist", dest="min_log_dist")
parser.add_argument("--max_log_dist", type=float, default=6.,
help="max_log_dist", dest="max_log_dist")

parser.add_argument("--min_y", type=float, default=0.,
help="min_y", dest="min_y")
parser.add_argument("--max_y", type=float, default=1.,
help="max_y", dest="max_y")

parser.add_argument("--line_width", type=int, default=2,
help="line_width", dest="line_width")

parser.add_argument("--target_name", type=str,
default="target site", help="target_name", dest="target_name")

parser.add_argument("--colors", nargs="+")
parser.add_argument("--cmap", default="rainbow", 
                    help="""optional, name of matplotlib colormap, see:
                    http://matplotlib.org/examples/color/colormaps_reference.html
                    (default: %(default)s)""")

parser.add_argument("--fig_size_x", type=float,
default=3, help="fig_size_x", dest="fig_size_x")
parser.add_argument("--fig_size_y", type=float,
default=2, help="fig_size_y", dest="fig_size_y")

parser.add_argument("--ignore_overlapping",
                    help="""if you would like to ignore sites that are overlapping between
query and target, then indicate with flag. (Basically, rather than the closest site if
a given site has zero distance, the next closest site will be found with non-zero distance.)

""", dest="ignore_overlapping", action='store_true', default=False)

###################################################

args = parser.parse_args()

###################################################

tableau20 = [(31/255., 119/255., 180/255.), 
             (174/255., 199/255., 232/255.),
             (255/255., 127/255., 14/255.),
             (255/255., 187/255., 120/255.),    
             (44/255., 160/255., 44/255.),
             (152/255., 223/255., 138/255.),
             (214/255., 39/255., 40/255.),
             (255/255., 152/255., 150/255.),    
             (148/255., 103/255., 189/255.),
             (197/255., 176/255., 213/255.),
             (140/255., 86/255., 75/255.),
             (196/255., 156/255., 148/255.),    
             (227/255., 119/255., 194/255.),
             (247/255., 182/255., 210/255.),
             (127/255., 127/255., 127/255.),
             (199/255., 199/255., 199/255.),    
             (188/255., 189/255., 34/255.),
             (219/255., 219/255., 141/255.),
             (23/255., 190/255., 207/255.),
             (158/255., 218/255., 229/255.)]    
tableau10 = tableau20[::2]

if not args.colors:
    if args.cmap == "tableau":
        if len(args.queries) <= 10:
            colors = tableau10[:len(args.queries)+1]
        else:
            colors = tableau20
    else:
        colors = discrete_cmap(len(args.queries), args.cmap)
else:
    colors = args.colors

##########
# Get distance distribution and ECDF
##########

def extract_dist_from_BedTool(bed):
    dists = np.array([int(dist.split("\t")[-1]) for dist in str(bed).split('\n') if dist != ''])
    return(dists)

dx = 0.01
min_log_dist = args.min_log_dist
max_log_dist = args.max_log_dist

def find_ecdf_and_supporting(query, target, io):
    
    dists = extract_dist_from_BedTool(query.closest(target, io=io, d=True))
    dists[dists == -1] = 999999999
    
    median = np.median(dists)    
    ecdf = ECDF(np.log10(dists + 1))        
    
    results_dict = {'dists':dists, \
                    'median':median, \
                    'ecdf':ecdf}
    
    return(results_dict)

results_dict = {}
for name, query in zip(args.names, args.queries):
    results_dict[name] = find_ecdf_and_supporting(BedTool(query).sort(), BedTool(args.target).sort(), args.ignore_overlapping)

x = np.arange(int(np.floor(min_log_dist)), int(np.ceil(max_log_dist)), dx)

##########
# Plot
##########

def create_ecdf_plot(ax, ecdf, label, color):
    
    ax.plot(x, ecdf(x), color=color, lw=args.line_width, label=label)
    ax.axhline(0.5, color='gray', ls='--')
    ax.set_yticks([args.min_y,(args.min_y+args.max_y)/2,args.max_y])
    ax.set_xticks(np.arange(min_log_dist, max_log_dist+1))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.legend(loc='upper left', fontsize=8, frameon=False)

fig, ax = plt.subplots(figsize=(args.fig_size_x,args.fig_size_y))

for name, color in zip(args.names, colors):
    create_ecdf_plot(ax, results_dict[name]['ecdf'], label=name, color=color)

lgd = ax.legend(loc="upper left", frameon=False, fontsize=8)
# plt.gca().add_artist(lgd)
ax.set_ylabel('Cumulative probability', fontsize=8)
ax.set_xlabel('Log10 dist. to nearest %s'%(args.target_name), fontsize=8)
ax.set_ylim((args.min_y+-0.01,args.max_y+0.01))
ax.set_xlim((args.min_log_dist, args.max_log_dist*1.01))
ax.yaxis.set_tick_params(width=args.line_width)
ax.spines['left'].set_linewidth(args.line_width)
ax.xaxis.set_tick_params(width=args.line_width)
ax.spines['bottom'].set_linewidth(args.line_width)
plt.tight_layout()
plt.savefig(args.outplot)

##########
# Print median distance and Mann-Whitney U-test pairwise comparisons
##########

# for name in args.names:
#     print "Median distance %s to nearest target site = %s Kb" % (name, float(results_dict[name]['median'])/1000)

# for i in range(len(args.names)):
#     name1 = args.names[i]    
#     for j in range(i):
#         name2 = args.names[j]
#         _,p = mannwhitneyu(results_dict[name1]['dists'], results_dict[name2]['dists'])
#         print 'Mann-whitney U-test, %s vs. %s, distance to nearest target site, p = %.3e' % (name1, name2, p)
