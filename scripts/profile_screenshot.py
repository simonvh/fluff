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
from numpy import array

### My imports ###
from fluff.plot import profile_screenshot

VERSION = 1.0
DEFAULT_COLORS = ["#e41a1c","#4daf4a","#377eb8"]

def process_groups(groups):
	if not groups:
		return None
	pg = []
	for group in groups.split(","):
		ids = [int(x) for x in group.split(":")]
		if len(ids) == 2:
			pg.append(range(ids[0], ids[1] + 1))
		else:
			pg.append(ids)
	return pg

def split_interval(interval):
	chrom,coords = interval.split(":")
	start,end = [int(x) for x in coords.replace(",","").split("-")]
	return chrom,start,end

parser = OptionParser(version="%prog " + str(VERSION))
parser.add_option("-i", "--intervals", dest="intervals", help="Intervals (chrom:start-end)")
parser.add_option("-d", "--datafiles", dest="datafiles", help="Data files (reads in BAM format)", metavar="FILE(S)")
parser.add_option("-o", "--outfile", dest="outfile", help="Output file name (type determined by extension)", metavar="FILE")
parser.add_option("-l", "--colors", dest="colors", help="Colors", metavar="NAME(S)", default=DEFAULT_COLORS)
parser.add_option("-a", "--annotation", dest="annotation", help="Annotation in BED12 format", metavar="FILE")
parser.add_option("-t", "--trackgroups", dest="trackgroups", help="Track groups", metavar="GROUPS")
parser.add_option("-s", "--scalegroups", dest="scalegroups", help="Scale groups", metavar="GROUPS")

(options, args) = parser.parse_args()

for opt in [options.intervals, options.datafiles, options.outfile]:
	if not opt:
		parser.print_help()
		sys.exit()

intervals = [x.strip() for x in options.intervals.split(",")]
datafiles = [x.strip() for x in options.datafiles.split(",")]
annotation = options.annotation
outfile = options.outfile
colors = [x.strip() for x in options.colors.split(",")]

trackgroups = process_groups(options.trackgroups)	
if not trackgroups:
	trackgroups = [[x] for x in range(1, len(datafiles) + 1)]

scalegroups = process_groups(options.scalegroups)	

if trackgroups and scalegroups:
	if len(trackgroups) != sum([len(x) for x in scalegroups]):
		sys.stderr.write("Track groups and scales do not match!\n")
		sys.exit()

# Group the tracks according to track_groups
tracks = []
for group in trackgroups:
	tracks.append([datafiles[i - 1] for i in group])

# Intervals
intervals = [split_interval(x) for x in intervals]

# Create the image
profile_screenshot(outfile, intervals, tracks, annotation=annotation, scalegroups=scalegroups, colors=colors, bgmode="white")
