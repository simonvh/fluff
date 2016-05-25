__author__ = 'george'

import os
import sys
import pysam

### My imports ###
from fluff.color import parse_colors
from fluff.plot import profile_screenshot
from fluff.util import process_groups

def profile(args):
    interval = args.interval
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
    if scale:
        try:
            scale = [int(x) for x in scale.split(",")]
        except Exception:
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

    # Create the image
    profile_screenshot(outfile, interval, tracks,
                       annotation=annotation,
                       scalegroups=scalegroups,
                       fontsize=args.textfontsize,
                       colors=colors,
                       bgmode=args.background,
                       fragmentsize=args.fragmentsize,
                       scale=scale,
                       show_scale=args.show_scale,
                       rmdup=args.rmdup,
                       rmrepeats=args.rmrepeats,
                       reverse=args.reverse,
                       adjscale=args.adjscale
                       )
