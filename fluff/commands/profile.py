__author__ = 'george'

import os

### My imports ###
from fluff.color import parse_colors
from fluff.plot import profile_screenshot
from fluff.util import *


def profile(args):
    intervals = [x.strip() for x in args.intervals.split(",")]
    datafiles = [x.strip() for x in args.datafiles]
    annotation = args.annotation
    outfile = args.outfile
    colors = parse_colors(args.colors)

    for x in args.datafiles:
        if '.bam' in x and not os.path.isfile("{0}.bai".format(x)):
            print "Data file '{0}' does not have an index file. Creating an index file for {0}.".format(x)
            pysam.index(x)

    trackgroups = process_groups(args.trackgroups)
    if not trackgroups:
        trackgroups = [[x] for x in range(1, len(datafiles) + 1)]

    scalegroups = process_groups(args.scalegroups)
    scale = args.scale
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
    profile_screenshot(outfile, intervals, tracks,
                       annotation=annotation,
                       scalegroups=scalegroups,
                       fontsize=args.textfontsize,
                       colors=colors,
                       bgmode=args.background,
                       fragmentsize=args.fragmentsize,
                       scale=scale,
                       rmdup=args.rmdup,
                       rmrepeats=args.rmrepeats,
                       reverse=args.reverse
                       )
