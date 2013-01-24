#!/usr/bin/env python
# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

### Standard imports ###
from optparse import OptionParser
from tempfile import NamedTemporaryFile
import sys
import os

### External imports ###
from numpy import array,hstack,arange,median,mean,zeros
from scipy.stats import scoreatpercentile
from scipy.stats.mstats import rankdata
import Pycluster
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import colorConverter, LinearSegmentedColormap
from matplotlib.ticker import NullFormatter,NullLocator
from matplotlib.font_manager import fontManager, FontProperties
import pp
from math import sqrt,log

### My imports ###
from fluff.util import *
from fluff.fluffio import *

#from kmeans import kmeanssample, Lqmetric

BINS = 100
RPKM = False
RMDUP = True
RMREPEATS = True
VERSION = "1.0"
COLOR_MAP = {"red":"#e41a1c","blue":"#377eb8","green":"#4daf4a","orange":"#ff7f00","brown":"#a65628","purple":"#984ea3","yellow":"#ffff33"}
DEFAULT_COLORS = "#e41a1c,#377eb8,#4daf4a,#984ea3,#ff7f00,#ffff33,#a65628"
METRIC = "e"		# Euclidian, PyCluster
FONTSIZE = 8
DEFAULT_SCALE = 15 
DEFAULT_EXTEND = 5000
DEFAULT_PERCENTILE = 99
DEFAULT_CLUSTERING = "kmeans"

def load_data(featurefile, datafile, bins=100, up=5000, down=5000, rmdup=True, rpkm=False, rmrepeats=True):
	tmp = tempfile.NamedTemporaryFile(delete=False)
	regions = []
	order = {}
	count = 0
	for line in open(featurefile):
		if line.startswith("#") or line[:5] == "track":
			continue
		vals = line.strip().split("\t")
		strand = "+"
		if len(vals) >= 6:
			strand = vals[5]
		middle = (int(vals[2]) + int(vals[1])) / 2
		start,end = middle, middle
		if strand == "+":
			start -= up
			end += down
		else:
			start -= down
			end += up
		if start >= 0:
			regions.append([vals[0], start, end, strand])
			order["%s:%s-%s" % (vals[0], start, end)] = count
			count += 1
			tmp.write("%s\t%s\t%s\t0\t0\t%s\n" % (vals[0], start, end, strand))
	tmp.flush()
	
	result = fluff.fluffio.get_binned_stats(tmp.name, datafile, bins, rpkm=rpkm, rmdup=rmdup, rmrepeats=rmrepeats)
	
	# Retrieve original order
	#r_regions = ["{}:{}-{}".format(*row.split("\t")[:3]) for row in result]
	#r_order = numpy.array([order[region] for region in r_regions]).argsort()[::-1]
	r_data = numpy.array([[float(x) for x in row.split("\t")[3:]] for row in result])

	return os.path.basename(datafile), regions, r_data#[r_order]

def normalize_data(data, percentile=75):
	norm_data = {}
	for track,ar in data.items():
		s = scoreatpercentile(ar.flatten(), percentile)
		if s == 0:
			sys.stderr.write("Error normalizing track %s as score at percentile %s is 0, normalizing to maximum value instead\n" % (track, percentile))
			x =  ar / max(ar.flatten())
		else:
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
parser.add_option("-c", "--clustering", dest="clustering", help="kmeans, hierarchical or none", default=DEFAULT_CLUSTERING)
parser.add_option("-k", "--numclusters", dest="numclusters", help="Number of clusters", metavar="INT", type="int")
parser.add_option("-l", "--colors", dest="colors", help="Colors", metavar="NAME(S)", default=DEFAULT_COLORS)
parser.add_option("-o", "--outfile", dest="outfile", help="Output file name", metavar="FILE")
parser.add_option("-s", "--scale", dest="scale", help="Scale", metavar="INT", type="int", default=DEFAULT_SCALE)
parser.add_option("-e", "--extend", dest="extend", help="Extend", metavar="INT", type="int", default=DEFAULT_EXTEND)

(options, args) = parser.parse_args()

for opt in [options.featurefile, options.datafiles, options.outfile, options.numclusters]:
	if not opt:
		parser.print_help()
		sys.exit()

featurefile = options.featurefile
datafiles = [x.strip() for x in options.datafiles.split(",")]
tracks = [os.path.basename(x) for x in datafiles]
colors = [x.strip() for x in options.colors.split(",")]
for i,color in enumerate(colors):
	if COLOR_MAP.has_key(color):
		colors[i] = COLOR_MAP[color]

outfile = options.outfile
extend_up = options.extend
extend_down = options.extend
cluster_type = options.clustering[0].lower()

if not cluster_type in ["k", "h", "n"]:
	sys.stderr.write("Unknown clustering type!\n")
	sys.exit(1)

## Get scale for each track
scale = [1.0 for track in datafiles]

# Calculate the profile data
# Load data in parallel
print "Loading data"
job_server = pp.Server(ncpus=4)
jobs = []
for datafile in datafiles:
	jobs.append(job_server.submit(load_data, (featurefile, datafile, BINS, extend_up, extend_down, RMDUP, RPKM, RMREPEATS),  (), ("tempfile","sys","os","fluff.fluffio","numpy")))

data = {}
regions = []
for job in jobs:
	track,regions,profile = job()
	#print "##### %s " % track
	#for row in profile:
	#	print row
	data[track] = profile

# Normalize
norm_data = normalize_data(data, DEFAULT_PERCENTILE)
clus = hstack([norm_data[t] for t in tracks])

if cluster_type == "k":
	print "K-means clustering"
	## K-means clustering
	# PyCluster
	labels, error, nfound = Pycluster.kcluster(clus, options.numclusters, dist=METRIC)
	ind = labels.argsort()
	# Other cluster implementation
	#	centres, labels, dist = kmeanssample(clus, options.numclusters, len(clus) / 10,  metric=cl, maxiter=200, verbose=1, delta=0.00001)
elif cluster_type == "h":
	print "Hierarchical clustering"
	tree = Pycluster.treecluster(clus, method="m", dist=METRIC)
	labels = tree.cut(options.numclusters)
	ind = sort_tree(tree, arange(len(regions)))
else:
	ind = arange(len(regions))
	#print ind
	labels = zeros(len(regions))

font = FontProperties(size=FONTSIZE / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])

f = open("%s_clusters.bed" % outfile, "w")
for (chrom,start,end,strand), cluster in zip(array(regions)[ind], array(labels)[ind]):
	f.write("%s\t%s\t%s\t%s\t0\t%s\n" % (chrom, start, end, cluster, strand))
f.close()

fig = plt.figure(figsize=(10,5))

axes = []
for i, track in enumerate(tracks):
	c = create_colormap('white', colors[i % len(colors)])
	ax = fig.add_subplot(1,len(tracks),i + 1)
	ax.set_title(track.replace(".bam",""),  fontproperties=font)
	axes.append(ax)
	ax.pcolormesh(data[track][ind], cmap=c, vmin=0, vmax=options.scale * scale[i])
	#print "%s\t%s\t%s\t%s" % (track, scale[i] * options.scale, mean(data[track][ind][:,0:20]), median(data[track][ind]))
	for x in [ax.xaxis, ax.yaxis]:
		x.set_major_formatter(NullFormatter())
		x.set_major_locator(NullLocator())
	for loc,spine in ax.spines.iteritems():
		spine.set_color('none')
fig.subplots_adjust(wspace=0, hspace=0)

#for i, track in enumerate(tracks):
	#axes[i].set_title(track.replace(".bam",""), verticalalignment={0:"top",1:"bottom"}[i % 2], zorder=100000)

ext = outfile.split(".")[-1]
if not ext in ["png", "svg", "ps", "eps", "pdf"]:
	outfile += ".png"
print "Saving image"
if outfile.endswith("png"):
	plt.savefig(outfile, dpi=600)
else:
	plt.savefig(outfile)
