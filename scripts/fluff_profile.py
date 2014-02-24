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
from numpy import array

### My imports ###
from fluff.plot import profile_screenshot
from fluff.util import *
from fluff.color import DEFAULT_COLORS, parse_colors
from fluff.config import *

VERSION = str(FL_VERSION)

BACKGROUNDS = ["white", "stripes", "color"]
FRAGMENTLENGTH = 200

usage = "Usage: %prog -i <loc1>[,<loc2>,...] -d <file1>[,<file2>,...] -o <out> [options]"
version = "%prog " + str(VERSION)
parser = OptionParser(version=version, usage=usage)
group1 = OptionGroup(parser, 'Optional')
parser.add_option("-i", dest="intervals", help="one or more genomic intervals (chrom:start-end)", metavar="INTERVAL(S)")
parser.add_option("-d", dest="datafiles", help="data files (reads in BAM or BED format)", metavar="FILE(S)")
parser.add_option("-o", dest="outfile", help="output file name (type determined by extension)", metavar="FILE")
group1.add_option("-a", dest="annotation", help="annotation in BED12 format", metavar="FILE")
group1.add_option("-c", dest="colors", help="color(s) (name, colorbrewer profile or hex code)", metavar="NAME(S)", default=DEFAULT_COLORS)
group1.add_option("-t", dest="trackgroups", help="track groups", metavar="GROUPS")
group1.add_option("-s", dest="scalegroups", help="scale groups", metavar="GROUPS")
group1.add_option("-S", dest="scale", help="scale: 'auto' (default), 'off' or int for each track", metavar="SCALE", default="auto")
group1.add_option("-b", dest="background", help="background color: white | color | stripes", default="white")
group1.add_option("-f", dest="fragmentsize", help="fragment length (default: %s)" % FRAGMENTLENGTH,type="int",  default=FRAGMENTLENGTH)

parser.add_option_group(group1)
(options, args) = parser.parse_args()

for opt in [options.intervals, options.datafiles, options.outfile]:
    if not opt:
        parser.print_help()
        sys.exit()

if not options.background in BACKGROUNDS:
    print "Please specify a correct background!"
    sys.exit(1)


intervals = [x.strip() for x in options.intervals.split(",")]
datafiles = [x.strip() for x in options.datafiles.split(",")]
annotation = options.annotation
outfile = options.outfile
colors = parse_colors(options.colors)

trackgroups = process_groups(options.trackgroups)    
if not trackgroups:
    trackgroups = [[x] for x in range(1, len(datafiles) + 1)]

scalegroups = process_groups(options.scalegroups)    
scale = options.scale
if scale == "auto":
    scale = True
elif scale == "off":
    scale = False
elif scale:
    try:
        scale = [int(x) for x in scale.split(",")]
    except:
        print "Error in scale argument"
        sys.exit(1)

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
profile_screenshot(outfile, intervals, tracks, annotation=annotation, scalegroups=scalegroups, colors=colors, bgmode=options.background, fragmentsize=options.fragmentsize, scale=scale)
