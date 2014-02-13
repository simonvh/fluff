#!/usr/bin/env python
# Copyright (c) 2012-2014 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

### Standard imports ###
from optparse import OptionParser,OptionGroup
import sys
import os

### My imports ###
from fluff.fluffio import load_cluster_data,load_bed_clusters
from fluff.plot import bandplot
from fluff.color import DEFAULT_COLORS,parse_colors
from fluff.config import *
from fluff.util import process_groups

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
group1.add_option("-b", dest="bins", help="number of bins", metavar="INT", default=BANDPLOT_BINS, type=int)
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

# Data
clust_file = options.clust_file
datafiles = [x.strip() for x in options.datafiles.split(",")]
fragmentsize = options.fragmentsize
tracks = [os.path.basename(x) for x in datafiles]
rmdup = options.rmdup
rpkm = options.rpkm
rmrepeats = options.rmrepeats
bins = options.bins

# Plot
titles = [os.path.splitext(x)[0] for x in tracks]
colors = parse_colors(options.colors)
scalegroups = process_groups(options.scalegroups)
percs = [int(x) for x in options.percs.split(",")]
summary = options.summary

# Calculate the profile data
data = load_cluster_data(clust_file, datafiles, bins, rpkm, rmdup, rmrepeats, fragmentsize=fragmentsize)
# Get cluster information
cluster_data = load_bed_clusters(clust_file)

# Make the plot
bandplot(options.outfile, tracks, cluster_data, data, bins=bins, summary=summary, colors=colors, scalegroups=scalegroups, percs=percs, titles=titles) 
