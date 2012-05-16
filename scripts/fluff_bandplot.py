#!/usr/bin/env python
# Copyright (c) 2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

### Standard imports ###
from optparse import OptionParser
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

### Other imports ###
from solexatools.track import SimpleTrack
from solexatools.peak_stats import bin_formatter, bam_binned_peak_stats, peak_stats
from fluff.plot import *

######## EDIT CONSTANTS TO CHANGE BEHAVIOUR OF THE SCRIPT #############
# Sizes of the plots (in inches)
PLOTWIDTH = 0.6
PLOTHEIGHT = 0.6
PAD = 0.05
PADLEFT = 1.5						# For labels
PADTOP = 0.4 						# For labels
PADBOTTOM = 0.05
PADRIGHT = 0.05
FONTSIZE = 8
BINS = 21								# Number of bins for profile
RPKM = False						# If False, use absolute count
REMOVE_DUP = True				# Remove duplicate reads, if present
DEFAULT_COLORS = "red,blue,green,orange,brown,purple,yellow" 
#########################################################################
font = FontProperties(size=FONTSIZE / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])

VERSION = "1.0"

def load_data(clust_file, datafiles):
	data = {}
	for datafile in datafiles:
		result = []
		clust_track = SimpleTrack(clust_file)
		if datafile.endswith("bam"):
			if not os.path.exists(datafile + ".bai"):
				print "Please provide indexed bam file(s) %s" % (datafile + ".bai")
				sys.exit()
			result = bam_binned_peak_stats(clust_track, datafile, BINS, RPKM, REMOVE_DUP)
		else:
			result = peak_stats(clust_track, SimpleTrack(datafile), bin_formatter, {"bins": BINS})
		result =  [row.split("\t") for row in result]
		data[os.path.basename(datafile)] = dict([["%s:%s-%s" % (vals[0], vals[1], vals[2]), [float(x) for x in vals[3:]]] for vals in result])
	return data

def load_clusters(clust_file):
	cluster_data = {}
	clust_track = SimpleTrack(clust_file)
	f = clust_track.get_next_feature()
	while f:
		cluster_data.setdefault(f[3], []).append("%s:%s-%s" % (f[0], f[1], f[2]))
		f = clust_track.get_next_feature()
	return cluster_data

parser = OptionParser(version="%prog " + str(VERSION))
parser.add_option("-c", "--clusterfile", dest="clust_file", help="File containing clusters", metavar="FILE")
parser.add_option("-d", "--datafiles", dest="datafiles", help="Data files (reads in BED/BAM format)", metavar="FILE(S)")
parser.add_option("-l", "--colors", dest="colors", help="Colors", metavar="NAME(S)", default=DEFAULT_COLORS)
parser.add_option("-o", "--outfile", dest="outfile", help="Output file (type determined by extension)", metavar="FILE")

(options, args) = parser.parse_args()

for opt in [options.clust_file, options.datafiles, options.outfile]:
	if not opt:
		parser.print_help()
		sys.exit()

clust_file = options.clust_file
datafiles = [x.strip() for x in options.datafiles.split(",")]
tracks = [os.path.basename(x) for x in datafiles]
colors = [x.strip() for x in options.colors.split(",")]

# Calculate the profile data
data = load_data(clust_file, datafiles)
# Get cluster information
cluster_data = load_clusters(clust_file)
clusters = [int(x) for x in cluster_data.keys()]

#Init x-axis
t = arange(BINS)

# Get a figure with a lot of subplots
fig, axes = create_grid_figure(len(tracks), len(clusters), plotwidth=PLOTWIDTH, plotheight=PLOTHEIGHT, padleft=PADLEFT, padtop=PADTOP, pad=PAD, padright=PADRIGHT, padbottom=PADBOTTOM) 

for track_num, track in enumerate(tracks):
	percentiles = [scoreatpercentile([data[track][x] for x in cluster_data[cluster]], 90) for cluster in clusters]
	track_max = max(percentiles)
	for i,cluster in enumerate(clusters):
		# Retrieve axes
		ax = axes[track_num][i]
		
		# Get the data
		vals = array([data[track][x] for x in cluster_data[cluster]])
		
		# Make the plot
		coverage_plot(ax, t, vals, colors[track_num % len(colors)])
		
		# Set scale	
		ax.set_ylim(0, track_max)
		ax.set_xlim(0, BINS - 1)	
		
		# Cluster titles
		if track_num == 0:
			ax.set_title("%s\nn=%s" % (cluster, len(cluster_data[cluster])), font_properties=font)
		
		# Track title and scale
		if i == 0:
			pos = axes[track_num][0].get_position().get_points()
			text_y = (pos[1][1] + pos[0][1]) / 2
			text_x = pos[0][0] - (PAD / fig.get_figwidth())
			plt.figtext(text_x, text_y, track, clip_on=False, horizontalalignment="right", verticalalignment="center", font_properties=font)
			plt.figtext(text_x,  pos[1][1], track_max, clip_on=False, horizontalalignment="right", verticalalignment="top", font_properties=font)
			plt.figtext(text_x,  pos[0][1], 0, clip_on=False, horizontalalignment="right", verticalalignment="bottom", font_properties=font)

print "Saving figure"
savefig(options.outfile)
