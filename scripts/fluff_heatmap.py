#!/usr/bin/env python
# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

### Standard imports ###
from optparse import OptionParser,OptionGroup
from tempfile import NamedTemporaryFile
import sys
import os

### External imports ###
from numpy import array,hstack,arange,median,mean,zeros
from scipy.stats.mstats import rankdata
import Pycluster
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter,NullLocator
from matplotlib.font_manager import fontManager, FontProperties
import pp
from math import sqrt,log

### My imports ###
from fluff.util import *
from fluff.fluffio import *
from fluff.color import create_colormap, COLOR_MAP, DEFAULT_COLORS, parse_colors

VERSION = "1.2"

DEFAULT_BINSIZE = 100
METRIC = "e"        # Euclidian, PyCluster
FONTSIZE = 8
DEFAULT_SCALE = "90%" 
DEFAULT_EXTEND = 5000
DEFAULT_PERCENTILE = 99
DEFAULT_CLUSTERING = "none"
DEFAULT_BG = "white"

usage = "Usage: %prog -f <bedfile> -d <file1>[,<file2>,...] -o <out> [options]"
version = "%prog " + str(VERSION)
parser = OptionParser(version=version, usage=usage)
group1 = OptionGroup(parser, 'Optional')

parser.add_option("-f", dest="featurefile", help="BED file containing features", metavar="FILE")
parser.add_option("-d", dest="datafiles", help="data files (reads in BAM or BED format)", metavar="FILE(S)")
parser.add_option("-o", dest="outfile", help="output file (type determined by extension)", metavar="FILE")

# Optional arguments
group1.add_option("-C", dest="clustering", help="kmeans, hierarchical or none", default=DEFAULT_CLUSTERING, metavar="METHOD")
group1.add_option("-k", dest="numclusters", help="number of clusters", metavar="INT", type="int", default=1)
group1.add_option("-m", dest="merge_mirrored", help="merge mirrored clusters (only with kmeans)", metavar="", default=False, action="store_true")
group1.add_option("-c", dest="colors", help="color(s) (name, colorbrewer profile or hex code)", metavar="NAME(S)", default=DEFAULT_COLORS)
group1.add_option("-B", dest="bgcolors", help="background color(s) (name, colorbrewer profile or hex code)", metavar="NAME(S)", default=DEFAULT_BG)
group1.add_option("-e", dest="extend", help="extend (in bp)", metavar="INT", type="int", default=DEFAULT_EXTEND)
group1.add_option("-b", dest="binsize", help="bin size (default %s)" % DEFAULT_BINSIZE, metavar="INT", type="int", default=DEFAULT_BINSIZE)
group1.add_option("-s", dest="scale", help="scale (absolute or percentage)", metavar="", type="string", default=DEFAULT_SCALE)
group1.add_option("-F", dest="fragmentsize", help="Fragment length (default: read length)",type="int",  default=None)
group1.add_option("-r", dest="rpkm", help="use RPKM instead of read counts", metavar="", action="store_true", default=False)
group1.add_option("-D", dest="rmdup", help="keep duplicate reads (removed by default)", metavar="", default=True, action="store_false")
group1.add_option("-R", dest="rmrepeats", help="keep repeats (removed by default, bwa only) ", metavar="", action="store_false", default=True)


parser.add_option_group(group1)
(options, args) = parser.parse_args()

for opt in [options.featurefile, options.datafiles, options.outfile]:
    if not opt:
        parser.print_help()
        sys.exit()

featurefile = options.featurefile
datafiles = [x.strip() for x in options.datafiles.split(",")]
tracks = [os.path.basename(x) for x in datafiles]
titles = [os.path.splitext(x)[0] for x in tracks]
colors = parse_colors(options.colors)
bgcolors = parse_colors(options.bgcolors)

outfile = options.outfile
extend_up = options.extend
extend_down = options.extend
fragmentsize = options.fragmentsize
cluster_type = options.clustering[0].lower()
merge_mirrored = options.merge_mirrored
bins = (extend_up + extend_down) / options.binsize
rmdup = options.rmdup
rpkm = options.rpkm
rmrepeats = options.rmrepeats

if not cluster_type in ["k", "h", "n"]:
    sys.stderr.write("Unknown clustering type!\n")
    sys.exit(1)

if cluster_type == "k" and not options.numclusters >= 2:
    sys.stderr.write("Please provide number of clusters!\n")
    sys.exit(1)

## Get scale for each track
tscale = [1.0 for track in datafiles]

# Calculate the profile data
# Load data in parallel
print "Loading data"
job_server = pp.Server(ncpus=4)
jobs = []
for datafile in datafiles:
    jobs.append(job_server.submit(load_heatmap_data, (featurefile, datafile, bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize),  (), ("tempfile","sys","os","fluff.fluffio","numpy")))

data = {}
regions = []
for job in jobs:
    track,regions,profile = job()
    #print "##### %s " % track
    #for row in profile:
    #    print row
    data[track] = profile

scale = get_absolute_scale(options.scale, [data[track] for track in tracks])

# Normalize
norm_data = normalize_data(data, DEFAULT_PERCENTILE)
clus = hstack([norm_data[t] for t in tracks])

if cluster_type == "k":
    print "K-means clustering"
    ## K-means clustering
    # PyCluster
    labels, error, nfound = Pycluster.kcluster(clus, options.numclusters, dist=METRIC)
    
    if merge_mirrored:
        (i,j) = mirror_clusters(data, labels)
        while j:
            for track in data.keys():
                data[track][labels == j] = [row[::-1] for row in data[track][labels == j]]
            for k in range(len(regions)):
                if labels[k] == j:
                    (chrom,start,end,strand) = regions[k]
                    if strand == "+":
                        strand = "-"
                    else:
                        strand = "+"
                    regions[k] = (chrom, start, end, strand)
            n = len(set(labels))
            labels[labels == j] = i
            for k in range(j + 1, n):
                labels[labels == k] = k - 1
            (i,j) = mirror_clusters(data, labels)
            
    ind = labels.argsort()
    # Other cluster implementation
    #    centres, labels, dist = kmeanssample(clus, options.numclusters, len(clus) / 10,  metric=cl, maxiter=200, verbose=1, delta=0.00001)
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
for (chrom,start,end,strand), cluster in zip(array(regions, dtype="object")[ind], array(labels)[ind]):
    f.write("%s\t%s\t%s\t%s\t0\t%s\n" % (chrom, start, end, cluster, strand))
f.close()

fig = plt.figure(figsize=(3,1 * len(tracks)))

axes = []
for i, track in enumerate(tracks):
    c = create_colormap(bgcolors[i % len(bgcolors)], colors[i % len(colors)])
    ax = fig.add_subplot(1,len(tracks),i + 1)
    ax.set_title(titles[i],  fontproperties=font)
    axes.append(ax)
    ax.pcolormesh(data[track][ind], cmap=c, vmin=0, vmax=scale * tscale[i])
    print "%s\t%s\t%s\t%s" % (track, tscale[i] * scale, mean(data[track][ind][:,0:20]), median(data[track][ind]))
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
<<<<<<< HEAD
	plt.savefig(outfile, dpi=300, bbox_inches='tight')
=======
    plt.savefig(outfile, dpi=600)
>>>>>>> 3af74a46210a0d466a5042eb785a759d4e6fa5ac
else:
    plt.savefig(outfile)
