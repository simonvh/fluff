#!/usr/bin/env python
# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

### Standard imports ###
from optparse import OptionParser,OptionGroup
import sys
import os

### External imports ###
from pylab import plot, show, ylim, yticks,savefig
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,NullLocator
from matplotlib.font_manager import fontManager, FontProperties
from numpy import array,mean,amax,max,std,median,arange
from scipy.interpolate import splprep, splev
from scipy.stats import scoreatpercentile
import pybedtools

### My imports ###
from fluff.fluffio import *
from fluff.plot import *
from fluff.util import *
from fluff.color import DEFAULT_COLORS,parse_colors
from fluff.config import *

######## EDIT CONSTANTS TO CHANGE BEHAVIOUR OF THE SCRIPT #############
# Sizes of the plots (in inches)
PLOTWIDTH = 0.6
PLOTHEIGHT = 0.6
PAD = 0.05
PADLEFT = 1.5                        # For labels
PADTOP = 0.4                         # For labels
PADBOTTOM = 0.05
PADRIGHT = 0.05
FONTSIZE = 8
BINS = 100                                # Number of bins for profile
#########################################################################
font = FontProperties(size=FONTSIZE / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])

VERSION = str(FL_VERSION)

usage = "Usage: %prog -i <bedfile> -d <file1>[,<file2>,...] -o <out> [options]"
version = "%prog " + str(VERSION)
parser = OptionParser(version=version, usage=usage)
group1 = OptionGroup(parser, 'Optional')
parser.add_option("-i", dest="clust_file", help="BED file with cluster in 4th column", metavar="FILE")
parser.add_option("-d", dest="datafiles", help="data files (reads in BAM or BED format)", metavar="FILE(S)")
parser.add_option("-o", dest="outfile", help="output file (type determined by extension)", metavar="FILE")
group1.add_option("-S", dest="summary", help="create summary graphs", default=False, action="store_true")
group1.add_option("-c", dest="colors", help="color(s) (name, colorbrewer profile or hex code)", metavar="NAME(S)", default=DEFAULT_COLORS)
group1.add_option("-b", dest="bins", help="number of bins", metavar="INT", default=BINS, type=int)
group1.add_option("-s", dest="scalegroups", help="scale groups", metavar="GROUPS")
group1.add_option("-p", dest="percs", help="range of percentiles (default 50,90)", metavar="INT,INT", default="50,90")
group1.add_option("-F", dest="fragmentsize", help="fragment length (default: read length)",type="int",  default=None)
group1.add_option("-r", dest="rpkm", help="use RPKM instead of read counts", metavar="", action="store_true", default=False)
group1.add_option("-D", dest="rmdup", help="keep duplicate reads (removed by default)", metavar="", default=True, action="store_false")
group1.add_option("-R", dest="rmrepeats", help="keep repeats (removed by default, bwa only) ", metavar="", action="store_false", default=True)


parser.add_option_group(group1)
(options, args) = parser.parse_args()

for opt in [options.clust_file, options.datafiles, options.outfile]:
    if not opt:
        parser.print_help()
        sys.exit()
clust_file = options.clust_file
for x in options.datafiles.split(","):
  if '.bam' in x and not os.path.isfile("{0}.bai".format(x)):
    print "Data file '{0}' does not have an index file".format(x)
    print "Creating an index file for {0}".format(x)
    pysam.index(x)
    print "Done!"
datafiles = [x.strip() for x in options.datafiles.split(",")]
fragmentsize = options.fragmentsize
tracks = [os.path.basename(x) for x in datafiles]
titles = [os.path.splitext(x)[0] for x in tracks]
colors = parse_colors(options.colors)
scalegroups = process_groups(options.scalegroups)
percs = [int(x) for x in options.percs.split(",")]
rmdup = options.rmdup
rpkm = options.rpkm
rmrepeats = options.rmrepeats
bins = options.bins
summary = options.summary

# Calculate the profile data
data = load_cluster_data(clust_file, datafiles, bins, rpkm, rmdup, rmrepeats, fragmentsize=fragmentsize)
# Get cluster information
cluster_data = load_bed_clusters(clust_file)
clusters = cluster_data.keys()
#Init x-axis
t = arange(bins)
rows = len(tracks)
cols = len(clusters)
if summary:
    rows += 1
    cols += 1
# Get a figure with a lot of subplots
fig, axes = create_grid_figure(rows, cols, plotwidth=PLOTWIDTH, plotheight=PLOTHEIGHT, padleft=PADLEFT, padtop=PADTOP, pad=PAD, padright=PADRIGHT, padbottom=PADBOTTOM) 
track_max = []
for track_num, track in enumerate(tracks):
    percentiles = [scoreatpercentile([data[track][x] for x in cluster_data[cluster]], 90) for cluster in clusters]
    track_max.append(max(percentiles))
for track_num, track in enumerate(tracks):
    for i,cluster in enumerate(clusters):
        # Retrieve axes
        ax = axes[track_num][i]
        # Get the data
        vals = array([data[track][x] for x in cluster_data[cluster]])
        # Make the plot
        coverage_plot(ax, t, vals, colors[track_num % len(colors)], percs)
        # Get scale max
        maxscale = track_max[track_num]
        if scalegroups and len(scalegroups) > 0:
            for group in scalegroups:
                if (track_num + 1) in group:
                    maxscale = max([track_max[j - 1] for j in group])
                    break
        # Set scale    
        ax.set_ylim(0, maxscale)
        ax.set_xlim(0, bins - 1)    
        # Cluster titles
        if track_num == 0:
            ax.set_title("%s\nn=%s" % (cluster, len(cluster_data[cluster])), font_properties=font)
        # Track title and scale
        if i == 0:
            pos = axes[track_num][0].get_position().get_points()
            text_y = (pos[1][1] + pos[0][1]) / 2
            text_x = pos[0][0] - (PAD / fig.get_figwidth())
            plt.figtext(text_x, text_y, titles[track_num], clip_on=False, horizontalalignment="right", verticalalignment="center", font_properties=font)
            plt.figtext(text_x,  pos[1][1], "%.4g" % maxscale, clip_on=False, horizontalalignment="right", verticalalignment="top", font_properties=font)
            plt.figtext(text_x,  pos[0][1], 0, clip_on=False, horizontalalignment="right", verticalalignment="bottom", font_properties=font)
if summary:
    for i,track in enumerate(tracks):
        ax = axes[i][cols - 1]
        l = len(clusters)
        min_alpha = 0.3
        max_alpha = 0.9
        if l > 1:
            step = (max_alpha - min_alpha) / (l - 1)
            alphas = arange(min_alpha, max_alpha + step, step)
        else:
            alphas = [max_alpha]
        for j,cluster in enumerate(clusters):
            vals = array([data[track][x] for x in cluster_data[cluster]])
            m = median(vals, axis=0)
            ax.plot(arange(len(m)), m, color=colors[i % len(colors)], alpha=alphas[j])
        ax.set_ylim(0, track_max[i])
    for i,cluster in enumerate(clusters):
        ax = axes[rows - 1][i]
        max_max = 0
        for j,track in enumerate(tracks):
            vals = array([data[track][x] for x in cluster_data[cluster]])
            m = median(vals, axis=0)
            ax.plot(arange(len(m)), m, color=colors[j % len(colors)], alpha=0.8)
            if track_max[j] > max_max:
                max_max = track_max[j]
        ax.set_ylim(0, max_max)
    ax = axes[rows - 1][cols - 1]
    ax.set_frame_on(False)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
print "Saving figure"
savefig(options.outfile, dpi=600)
