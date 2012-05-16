#!/usr/bin/env python
# Copyright (c) 2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

### Standard imports ###
from optparse import OptionParser
from tempfile import NamedTemporaryFile
import sys
import os

from numpy import array,hstack
#import numpy
from scipy.stats import scoreatpercentile
from scipy.stats.mstats import rankdata
#from scipy.cluster.vq import kmeans2
import Pycluster
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import colorConverter, LinearSegmentedColormap
from matplotlib.ticker import NullFormatter,NullLocator
import pp

### Other imports ###
from solexatools.track import SimpleTrack
from solexatools.peak_stats import bin_formatter, bam_binned_peak_stats, peak_stats

#from kmeans import kmeanssample, Lqmetric

BINS = 100
RPKM = False
REMOVE_DUP = True
VERSION = "1.0"
DEFAULT_COLORS = "red,blue,green,orange,brown,purple,yellow"
DEFAULT_COLORS = "#e41a1c,#377eb8,#4daf4a,#984ea3,#ff7f00,#ffff33,#a65628"
METRIC = "e"		# Euclidian, PyCluster

def load_data(featurefile, datafile, bins=100, up=5000, down=5000, remove_dup=True, rpkm=False):
	tmp = tempfile.NamedTemporaryFile(delete=False)
	t = solexatools.track.SimpleTrack(featurefile)
	f = t.get_next_feature()
	regions = []
	while f:
		strand = f[4]
		if not strand:
			strand = "+"
		middle = (f[2] + f[1]) / 2
		start,end = middle, middle
		if strand == "+":
			start -= up
			end += down
		else:
			start -= down
			end += up
		if start >= 0:
			regions.append([f[0], start, end, strand])
			tmp.write("%s\t%s\t%s\t0\t0\t%s\n" % (f[0], start, end, strand))
		f = t.get_next_feature()
	tmp.flush()
	
	data = {}
	result = []
	if datafile.endswith("bam"):
		if not os.path.exists(datafile + ".bai"):
			print "Please provide indexed bam file(s) %s" % (datafile + ".bai")
			sys.exit()
		result = solexatools.peak_stats.bam_binned_peak_stats(solexatools.track.SimpleTrack(tmp.name), datafile, bins, rpkm, remove_dup)
	else:
		result = solexatools.peak_stats.peak_stats(solexatools.track.SimpleTrack(tmp.name), solexatools.track.SimpleTrack(datafile), solexatools.peak_stats.bin_formatter, {"bins": bins}, )
	return os.path.basename(datafile), regions, numpy.array([[float(x) for x in row.split("\t")[3:]] for row in result])

def normalize_data(data, percentile=75):
	norm_data = {}
	for track,ar in data.items():
		x =  ar / scoreatpercentile(ar.flatten(), percentile)
		#x[x <= 0.5] = 0
		x[x >= 1.0] = 1
		norm_data[track] = x
	return norm_data

def create_colormap(col1, col2):
	c1 = colorConverter.to_rgb(col1)
	c2 = colorConverter.to_rgb(col2)

	cdict = {
		'red': ((0.,c1[0], c1[0]),(1.,c2[0], c2[0])),
		'green': ((0.,c1[1], c1[1]),(1.,c2[1], c2[1])),
		'blue': ((0.,c1[2], c1[2]),(1.,c2[2], c2[2]))
	}
	return LinearSegmentedColormap('custom', cdict, 256)


parser = OptionParser(version="%prog " + str(VERSION))
parser.add_option("-f", "--featurefile", dest="featurefile", help="File containing features", metavar="FILE")
parser.add_option("-d", "--datafiles", dest="datafiles", help="Data files (reads in BED/BAM format)", metavar="FILE(S)")
parser.add_option("-k", "--numclusters", dest="numclusters", help="Number of clusters", metavar="INT", type="int")
parser.add_option("-l", "--colors", dest="colors", help="Colors", metavar="NAME(S)", default=DEFAULT_COLORS)
parser.add_option("-o", "--outfile", dest="outfile", help="Output file name", metavar="FILE")

(options, args) = parser.parse_args()

for opt in [options.featurefile, options.datafiles, options.outfile, options.numclusters]:
	if not opt:
		parser.print_help()
		sys.exit()

featurefile = options.featurefile
datafiles = [x.strip() for x in options.datafiles.split(",")]
tracks = [os.path.basename(x) for x in datafiles]
colors = [x.strip() for x in options.colors.split(",")]
outfile = options.outfile
extend_up = 5000
extend_down = 5000

# Calculate the profile data
# Load data in parallel
print "Loading data"
job_server = pp.Server()
jobs = []
for datafile in datafiles:
	jobs.append(job_server.submit(load_data, (featurefile, datafile, BINS, extend_up, extend_down, REMOVE_DUP, RPKM), (), ("tempfile","sys","os","solexatools.peak_stats","solexatools.track","numpy")))

data = {}
regions = []
for job in jobs:
	track,regions,profile = job()
	data[track] = profile

# Normalize
norm_data = normalize_data(data, 75)
clus = hstack([norm_data[t] for t in tracks])

print "Clustering"
## K-means clustering
# PyCluster
labels, error, nfound = Pycluster.kcluster(clus, options.numclusters, dist=METRIC)
# Other cluster implementation
#	centres, labels, dist = kmeanssample(clus, options.numclusters, len(clus) / 10,  metric=cl, maxiter=200, verbose=1, delta=0.00001)
## Hierarchical clusterling
#tree = Pycluster.treecluster(clus, method="m", dist=METRIC)
#labels = tree.cut(options.numclusters)

f = open("%s_clusters.bed" % outfile, "w")
for (chrom,start,end,strand), cluster in zip(regions, labels):
	f.write("%s\t%s\t%s\t%s\t0\t%s\n" % (chrom, start, end, cluster, strand))
f.close()

ind = labels.argsort()
fig = plt.figure(figsize=(10,5))

for i, track in enumerate(tracks):
	c = create_colormap('white', colors[i])
	ax = fig.add_subplot(1,len(tracks),i + 1)
	ax.pcolormesh(data[track][ind], cmap=c, vmin=0, vmax=15)
	ax.set_title("bla")
	for x in [ax.xaxis, ax.yaxis]:
		x.set_major_formatter(NullFormatter())
		x.set_major_locator(NullLocator())
	for loc,spine in ax.spines.iteritems():
		spine.set_color('none')
fig.subplots_adjust(wspace=0, hspace=0)

print "Saving image"
plt.savefig("%s.png" % (outfile))
