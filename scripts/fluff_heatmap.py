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
import os.path
import multiprocessing

### External imports ###
from numpy import array,hstack,arange,median,mean,zeros
from scipy.stats.mstats import rankdata
import Pycluster
from math import sqrt,log

### My imports ###
from fluff.util import *
from fluff.fluffio import *
from fluff.color import DEFAULT_COLORS, parse_colors
from fluff.plot import heatmap_plot
from fluff.config import *

VERSION = str(FL_VERSION)

DEFAULT_BINSIZE = 100
DEFAULT_METRIC = "e"        # Euclidian, PyCluster
FONTSIZE = 8
DEFAULT_SCALE = "90%" 
DEFAULT_EXTEND = 5000
DEFAULT_PERCENTILE = 99
DEFAULT_CLUSTERING = "none"
DEFAULT_BG = "white"
DEFAULT_DYN_EXTEND = 5000

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
                  help="merge mirrored clusters (only with kmeans and without -g option)", 
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
                  help="extend (in bp. Default: {0})".format(DEFAULT_EXTEND), 
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
                  help="Identify dynamics by extending features 1kb up/down stream(just for clustering), cluster as 1 bin, diplay as original number of bins and with the default extend values",  
                  metavar="", 
                  action="store_true", 
                  default=False)
parser.add_option_group(group1)
(options, args) = parser.parse_args()

#Check if files exist
for opt in [options.featurefile, options.datafiles, options.outfile]:
    if not opt:
        parser.print_help()
        sys.exit()
if not os.path.isfile(options.featurefile):
  print "ERROR: The BED files which contains the features file does not exist!"
  sys.exit(1)
for x in options.datafiles.split(","):
  if not os.path.isfile(x):
    print "ERROR: Data file '{0}' does not exist".format(x)
    sys.exit(1)
for x in options.datafiles.split(","):
  if '.bam' in x and not os.path.isfile("{0}.bai".format(x)):
    print "ERROR: Data file '{0}' does not have an index file!".format(x)
    sys.exit(1)
#Options Parser
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

#Check for given number of processors
if (ncpus>multiprocessing.cpu_count()):
  print "ERROR: You can use only up to {0} processors!".format(multiprocessing.cpu_count())
  sys.exit(1)
#Check for mutually exclusive parameters
if merge_mirrored and dynam:
  print "ERROR: -m and -g option CANNOT be used together"
  sys.exit(1)
#Warning about too much files
if (len(tracks)>4):
  print "Warning: Running fluff with too many files might make you system use enormous amount of memory!"
#Method of clustering
if (options.pick != None):
  pick = [i - 1 for i in split_ranges(options.pick)]
else:
  pick = range(len(datafiles))
if not cluster_type in ["k", "h", "n"]:
    sys.stderr.write("Unknown clustering type!\n")
    sys.exit(1)
#Number of clusters
if cluster_type == "k" and not options.numclusters >= 2:
    sys.stderr.write("Please provide number of clusters!\n")
    sys.exit(1)
#Distance function
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

#Function to load heatmap data
def load_data(featurefile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam, guard=[]):
  # Calculate the profile data
  data = {}
  regions = []
  if guard or not dynam:
    print "Loading data"
  try:
    # Load data in parallel
    import pp
    job_server = pp.Server(ncpus)
    jobs = []
    for datafile in datafiles:
        jobs.append(job_server.submit(load_heatmap_data, (featurefile, datafile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam, guard),  (), ("tempfile","sys","os","fluff.fluffio","numpy")))
    for job in jobs:
        track,regions,profile,guard = job()
        data[track] = profile
  except:
    sys.stderr.write("Parallel Python (pp) not installed, can't load data in parallel\n")
    for datafile in datafiles:
        track,regions,profile,guard = load_heatmap_data(featurefile, datafile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam, guard)
        data[track] = profile
  return data, regions, guard

#-g : Option to try and get dynamics
#Extend features 1kb up/down stream
#Cluster them in one bin
guard=[]
if dynam:
  amount_bins = 1
  extend_dyn_up = 1000
  extend_dyn_down = 1000
  data, regions, guard = load_data(featurefile, bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam, guard)
else:
  amount_bins = bins
  extend_dyn_up = extend_up
  extend_dyn_down = extend_down

#Load data for clustering
data, regions, guard = load_data(featurefile, amount_bins, extend_dyn_up, extend_dyn_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam, guard)
# Normalize
norm_data = normalize_data(data, DEFAULT_PERCENTILE)
clus = hstack([norm_data[t] for i,t in enumerate(tracks) if (not pick or i in pick)])

#Clustering
if cluster_type == "k":
    print "K-means clustering"
    ## K-means clustering
    # PyCluster
    labels, error, nfound = Pycluster.kcluster(clus, options.numclusters, dist=METRIC)
    if not dynam and merge_mirrored:
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
  if not gene:
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chrom, start, end, cluster+1, strand))
  else: 
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrom, start, end, gene, cluster+1, strand))
f.close()

if not cluster_type == "k":
    labels = None

#Save read count matrix
input_file = open('{0}_readcount.txt'.format(outfile), 'w')
for i, track in enumerate(tracks):
  for k, v in data.items():
    if track == k:
      order = []
      for row in v:
        for x in row:
          order.append(x)
      for j in ind:
        input_file.write('{0}\t'.format(str(order[j])))
      input_file.write('\n')

#Load data for visualization if -g option was used
if dynam:
  data, regions, guard = load_data(featurefile, bins, extend_up, extend_down, rmdup, rpkm, rmrepeats, fragmentsize, dynam, guard)

scale = get_absolute_scale(options.scale, [data[track] for track in tracks])
heatmap_plot(data, ind, outfile, tracks, titles, colors, bgcolors, scale, tscale, labels)
