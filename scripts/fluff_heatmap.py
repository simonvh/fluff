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
import multiprocessing

### External imports ###
from numpy import array,hstack,arange,median,mean,zeros
from scipy.stats.mstats import rankdata
import Pycluster
import matplotlib.cm as cm
from math import sqrt,log

### My imports ###
from fluff.util import *
from fluff.fluffio import *
from fluff.color import DEFAULT_COLORS, parse_colors
from fluff.plot import heatmap_plot

VERSION = "1.3"

DEFAULT_BINSIZE = 100
DEFAULT_METRIC = "e"        # Euclidian, PyCluster
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

parser.add_option("-f", 
                  dest="featurefile", 
                  help="BED file containing features", 
                  metavar="FILE")
parser.add_option("-d", 
                  dest="datafiles", 
                  help="data files (reads in BAM or BED format)", 
                  metavar="FILE(S)")
parser.add_option("-o", 
                  dest="outfile", 
                  help="output file (type determined by extension)", 
                  metavar="FILE")

# Optional arguments
group1.add_option("-p", 
                  dest="pick", 
                  help="pick specific data files to use for clustering", 
                  default=None,
                  type="string")
group1.add_option("-C", 
                  dest="clustering", 
                  help="kmeans, hierarchical or none", 
                  default=DEFAULT_CLUSTERING, 
                  metavar="METHOD")
group1.add_option("-k", 
                  dest="numclusters", 
                  help="number of clusters", 
                  metavar="INT", 
                  type="int", 
                  default=1)
group1.add_option("-m", 
                  dest="merge_mirrored", 
                  help="merge mirrored clusters (only with kmeans)", 
                  metavar="", 
                  default=False, 
                  action="store_true")
group1.add_option("-c", 
                  dest="colors", 
                  help="color(s) (name, colorbrewer profile or hex code)", 
                  metavar="NAME(S)", 
                  default=DEFAULT_COLORS)
group1.add_option("-B", 
                  dest="bgcolors", 
                  help="background color(s) (name, colorbrewer profile or hex code)", 
                  metavar="NAME(S)", 
                  default=DEFAULT_BG)
group1.add_option("-e", 
                  dest="extend", 
                  help="extend (in bp)", 
                  metavar="INT", 
                  type="int", 
                  default=DEFAULT_EXTEND)
group1.add_option("-b", 
                  dest="binsize", 
                  help="bin size (default {0})".format(DEFAULT_BINSIZE), 
                  metavar="INT", 
                  type="int", 
                  default=DEFAULT_BINSIZE)
group1.add_option("-s", 
                  dest="scale", 
                  help="scale (absolute or percentage)", 
                  metavar="", 
                  type="string", 
                  default=DEFAULT_SCALE)
group1.add_option("-F", 
                  dest="fragmentsize", 
                  help="Fragment length (default: read length)",
                  type="int",
                  default=None)
group1.add_option("-r", 
                  dest="rpkm", 
                  help="use RPKM instead of read counts", 
                  metavar="", 
                  action="store_true", 
                  default=False)
group1.add_option("-D", 
                  dest="rmdup", 
                  help="keep duplicate reads (removed by default)", 
                  metavar="", 
                  default=True, 
                  action="store_false")
group1.add_option("-R", 
                  dest="rmrepeats", 
                  help="keep repeats (removed by default, bwa only) ", 
                  metavar="", 
                  action="store_false", 
                  default=True)
group1.add_option("-P", 
                  dest="cpus", 
                  help="number of CPUs (default: 4)", 
                  metavar="INT", 
                  type="int", 
                  default=4)
group1.add_option("-M", 
                  dest="distancefunction", 
                  help="Euclidean or Pearson (default: Euclidean)", 
                  default=DEFAULT_METRIC,
                  metavar="METHOD")
group1.add_option("-g", 
                  dest="graphdynamics", 
                  help="Cluster as 1 bin, diplay as original number of bins", 
                  metavar="", 
                  action="store_true", 
                  default=False)
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
ncpus = options.cpus
distancefunction = options.distancefunction[0].lower()
dynam = options.graphdynamics

if (ncpus>multiprocessing.cpu_count()):
  print "Warning: You can use only up to {0} processors!".format(multiprocessing.cpu_count())
  sys.exit(1)
  
if (len(tracks)>4):
  print "Warning: Running fluff with too many files might make you system use enormous amount of memory!"

if (options.pick != None):
  pick = [i - 1 for i in split_ranges(options.pick)]
else:
  pick = range(len(datafiles))

if not cluster_type in ["k", "h", "n"]:
    sys.stderr.write("Unknown clustering type!\n")
    sys.exit(1)

if cluster_type == "k" and not options.numclusters >= 2:
    sys.stderr.write("Please provide number of clusters!\n")
    sys.exit(1)

if not distancefunction in ["e", "p"]: 
  sys.stderr.write("Unknown distance function!\n")
  sys.exit(1)
else:
  if distancefunction == "e":
    METRIC = DEFAULT_METRIC
    print "Euclidean distance method"
  else:
    METRIC = "c"
    print "Pearson distance method"
## Get scale for each track
tscale = [1.0 for track in datafiles]

def load_data(featurefile, amount_bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize):
  # Calculate the profile data
  data = {}
  regions = []
  
  print "Loading data"
  try:
    # Load data in parallel
    import pp

    job_server = pp.Server(ncpus)
    jobs = []
    for datafile in datafiles:
        jobs.append(job_server.submit(load_heatmap_data, (featurefile, datafile, amount_bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize),  (), ("tempfile","sys","os","fluff.fluffio","numpy")))

    for job in jobs:
        track,regions,profile = job()
        data[track] = profile
  except:
    sys.stderr.write("Parallel Python (pp) not installed, can't load data in parallel\n")
    for datafile in datafiles:
        track,regions,profile = load_heatmap_data(featurefile, datafile, amount_bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize)
        data[track] = profile

  return data, regions

if dynam:
  amount_bins = 1
else:
  amount_bins = bins
  
data, regions = load_data(featurefile, amount_bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize)

# Normalize
norm_data = normalize_data(data, DEFAULT_PERCENTILE)
clus = hstack([norm_data[t] for i,t in enumerate(tracks) if (i in pick or not pick)])

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
                    (chrom,start,end,gene,strand) = regions[k]
                    if strand == "+":
                        strand = "-"
                    else:
                        strand = "+"
                    regions[k] = (chrom, start, end, gene, strand)
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
    labels = zeros(len(regions))

f = open("{0}_clusters.bed".format(outfile), "w")
for (chrom,start,end,gene,strand), cluster in zip(array(regions, dtype="object")[ind], array(labels)[ind]):
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, start, end, gene, cluster, strand))
f.close()

if not cluster_type == "k":
    labels = None

#Save read count matrix
input_file = open('{0}_readcount.txt'.format(outfile), 'w')
for k, v in data.items():
  for row in v:
    for x in row:
      input_file.write('{0}\t'.format(str(x)))
    input_file.write('\n')

if dynam:
  data, regions = load_data(featurefile, bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize)

scale = get_absolute_scale(options.scale, [data[track] for track in tracks])
heatmap_plot(data, ind, outfile, tracks, titles, colors, bgcolors, scale, tscale, labels)