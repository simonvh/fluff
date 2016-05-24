__author__ = 'george'

import os
import sys
import pysam

### External imports ###
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from scipy.stats import scoreatpercentile

### My imports ###
from fluff.color import parse_colors
from fluff.fluffio import load_read_counts, load_bed_clusters, load_cluster_data
from fluff.util import process_groups
from fluff.plot import create_grid_figure,coverage_plot
import fluff.config as cfg

def bandplot(args):
    if (0 > args.scalar) or (args.scalar > 100):
      print "ERROR: -P value has to be between 0 and 100"
      sys.exit(1)
    else:
      scalar = args.scalar

    if not args.datafiles and not args.readCount:
        print 'You should provide data file(s) or the read counts file.'
        sys.exit()
    if args.datafiles and args.readCount:
        print 'You should choose only ONE option. Either data file(s) or the read counts file.'
        sys.exit()


    clust_file = args.clust_file

    if args.datafiles:
        for x in args.datafiles:
          if '.bam' in x and not os.path.isfile("{0}.bai".format(x)):
              print "Data file '{0}' does not have an index file. Creating an index file for {0}.".format(x)
              pysam.index(x)
        datafiles = [x.strip() for x in args.datafiles]


    fragmentsize = args.fragmentsize
    colors = parse_colors(args.colors)
    scalegroups = process_groups(args.scalegroups)
    percs = [int(x) for x in args.percs.split(",")]
    rmdup = args.rmdup
    rpkm = args.rpkm
    rmrepeats = args.rmrepeats
    bins = args.bins
    summary = args.summary
    fontsize = args.textfontsize
    font = FontProperties(size=fontsize / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])
    # Calculate the profile data
    if args.datafiles:
        data = load_cluster_data(clust_file, datafiles, bins, rpkm, rmdup, rmrepeats, fragmentsize=fragmentsize)
        tracks = [os.path.basename(x) for x in datafiles]
        titles = [os.path.splitext(x)[0] for x in tracks]
    else:
        titles, data = load_read_counts(args.readCount)
        tracks = titles
        for x in data:
            for i in data[x]:
                bins = len(data[x][i])
                break
            break

    # Get cluster information
    cluster_data = load_bed_clusters(clust_file)
    clusters = cluster_data.keys()
    #Init x-axis
    t = np.arange(bins)
    rows = len(tracks)
    cols = len(clusters)
    if summary:
        rows += 1
        cols += 1
    # Get a figure with a lot of subplots
    fig, axes = create_grid_figure(rows, cols, 
            plotwidth=cfg.PLOTWIDTH, 
            plotheight=cfg.PLOTHEIGHT, 
            padleft=cfg.PADLEFT, 
            padtop=cfg.PADTOP, 
            pad=cfg.PAD, 
            padright=cfg.PADRIGHT, 
            padbottom=cfg.PADBOTTOM)
    track_max = []
    for track_num, track in enumerate(tracks):
        percentiles = [scoreatpercentile([data[track][x] for x in cluster_data[cluster]], scalar) for cluster in clusters]
        track_max.append(max(percentiles))
    for track_num, track in enumerate(tracks):
        for i,cluster in enumerate(clusters):
            # Retrieve axes
            ax = axes[track_num][i]
            # Get the data
            vals = np.array([data[track][x] for x in cluster_data[cluster]])
            # Make the plot
            coverage_plot(ax, t, vals, colors[track_num % len(colors)], percs)
            # Get scale max
            maxscale = track_max[track_num]
            if scalegroups and len(scalegroups) > 0:
                for group in scalegroups:
                    if (track_num + 1) in group:
                        maxscale = max([track_max[j - 1] for j in group])
                        break
            # Set scale
            ax.set_ylim(0, maxscale)
            ax.set_xlim(0, bins - 1)
            # Cluster titles
            if track_num == 0:
                ax.set_title("%s\nn=%s" % (cluster, len(cluster_data[cluster])), font_properties=font)
            # Track title and scale
            if i == 0:
                pos = axes[track_num][0].get_position().get_points()
                text_y = (pos[1][1] + pos[0][1]) / 2
                text_x = pos[0][0] - (cfg.PAD / fig.get_figwidth())
                plt.figtext(text_x, text_y, titles[track_num], clip_on=False, horizontalalignment="right", verticalalignment="center", font_properties=font)
                plt.figtext(text_x,  pos[1][1], "%.4g" % maxscale, clip_on=False, horizontalalignment="right", verticalalignment="top", font_properties=font)
                plt.figtext(text_x,  pos[0][1], 0, clip_on=False, horizontalalignment="right", verticalalignment="bottom", font_properties=font)
    if summary:
        for i,track in enumerate(tracks):
            ax = axes[i][cols - 1]
            l = len(clusters)
            min_alpha = 0.3
            max_alpha = 0.9
            if l > 1:
                step = (max_alpha - min_alpha) / (l - 1)
                alphas = np.arange(min_alpha, max_alpha + step, step)
            else:
                alphas = [max_alpha]
            for j,cluster in enumerate(clusters):
                vals = np.array([data[track][x] for x in cluster_data[cluster]])
                m = np.median(vals, axis=0)
                ax.plot(np.arange(len(m)), m, color=colors[i % len(colors)], alpha=alphas[j])
            ax.set_ylim(0, track_max[i])
        for i,cluster in enumerate(clusters):
            ax = axes[rows - 1][i]
            max_max = 0
            for j,track in enumerate(tracks):
                vals = np.array([data[track][x] for x in cluster_data[cluster]])
                m = np.median(vals, axis=0)
                ax.plot(np.arange(len(m)), m, color=colors[j % len(colors)], alpha=0.8)
                if track_max[j] > max_max:
                    max_max = track_max[j]
            ax.set_ylim(0, max_max)
        ax = axes[rows - 1][cols - 1]
        ax.set_frame_on(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
    print "Saving figure"
    plt.savefig(args.outfile, dpi=600)
