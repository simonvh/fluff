import os
import re
import sys

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.offsetbox import HPacker, TextArea, AnnotationBbox
from matplotlib.patches import FancyArrowPatch, ArrowStyle, Polygon
from matplotlib.ticker import NullFormatter, NullLocator
from scipy.stats import scoreatpercentile

from fluff.color import create_colormap
from fluff.config import FONTSIZE
from fluff.fluffio import load_annotation
from fluff.track import Track

DEFAULT_COLORS = ["#e41a1c", "#4daf4a", "#377eb8"]
GENE_ARROW = "->"
GENE_ARROW = ArrowStyle._Curve(beginarrow=False, endarrow=True, head_length=.4, head_width=.4)

def colortext(x, y, texts, colors, **kwargs):
    pos = {
            "right": 1,
            "center": 0.5,
            "left": 0,
            "top": 0,
            "bottom": 1
            }
    
    ax = kwargs.get("ax")
    verticalalignment = pos[kwargs.get("verticalalignment", "center")]
    horizontalalignment = pos[kwargs.get("horizontalalignment", "center")]
    annotation_clip = kwargs.get("clip_on", False)
    fontproperties = kwargs.get("fontproperties", None)
    textprops = {"fontproperties":fontproperties}
    transform = kwargs.get("transform", None)

    areas = []
    for t,c in zip(texts, colors):
        textprops["color"] = c    
        text = TextArea(t, textprops=textprops)
        areas.append(text)
        
    txt = HPacker(children=areas,
                    align="baseline",
                    pad=0, sep=0)
    
    bbox =  AnnotationBbox(txt, xy=(x, y),
                            xycoords='data',
                            annotation_clip=annotation_clip,
                            frameon=False,
                            boxcoords=("axes fraction"),
                            box_alignment=(
                                horizontalalignment, 
                                verticalalignment), # alignment center, center
                            #bboxprops={"bbox_transmuter":transform},
                            )

    ax.add_artist(bbox)

def hide_axes(ax):
    for x in [ax.xaxis, ax.yaxis]:
        x.set_major_formatter(NullFormatter())
        x.set_major_locator(NullLocator())
    for _, spine in ax.spines.iteritems():
        spine.set_color('none')

def heatmap_plot(data, ind, outfile, tracks, titles, colors, bgcolors, scale, tscale, labels, fontsize):
    font = FontProperties(size=fontsize / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])
    label_ratio = 4.0
    # space between heatmaps
    btw_space = 0.02
    plot_width = 1.75 * len(tracks) + btw_space * len(tracks)
    plot_height = 6
    width_ratios = [label_ratio] * len(tracks)
    numplots = len(tracks)
    if labels is not None and len(labels) == len(ind):
        plot_width += 1 / label_ratio
        numplots += 1
        width_ratios += [1]

    # Create figure
    fig = plt.figure(figsize=(plot_width, plot_height))
    # Create subplot layout
    gs = gridspec.GridSpec(1, numplots, width_ratios=width_ratios)

    axes = []
    for i, track in enumerate(tracks):
        c = create_colormap(bgcolors[i % len(bgcolors)], colors[i % len(colors)])
        ax = plt.subplot(gs[i])
        ax.set_title(titles[i], fontproperties=font, y=1)
        axes.append(ax)
        ax.pcolormesh(data[track][ind], cmap=c, vmin=0, vmax=scale * tscale[i])
        hide_axes(ax)
        ylim = ax.get_ylim()

    if labels is not None and len(labels) == len(ind):
        ax = plt.subplot(gs[len(tracks)])
        ax.axis('off')
        min_y, max_y = ylim
        s = 0

        plt.axhline(y=0,
                    color="grey",
                    linewidth=0.5,
                    alpha=0.5
                    )
        labels = np.array(labels)
        # Smaller cluster on the top ([::-1])
        for i in range(max(labels) + 1)[::-1]:
            prev = s
            s += sum(labels == i)
            plt.axhline(y=s + 1 - 1,
                        color="grey",
                        linewidth=0.5,
                        alpha=0.5
                        )
            plt.text(0.5, (prev + s) / 2,
                     str(i + 1),
                     verticalalignment="center",
                     horizontalalignment="center",
                     fontproperties=font)
        
        ax.set_ylim(ylim)

    fig.subplots_adjust(wspace=btw_space, hspace=0)
    ext = outfile.split(".")[-1]
    if ext not in ["png", "svg", "ps", "eps", "pdf"]:
        outfile += ".png"
    sys.stderr.write("Saving figure\n")
    if outfile.endswith("png"):
        plt.savefig(outfile, dpi=600, bbox_inches='tight')
    else:
        plt.savefig(outfile)

def coverage_plot(ax, x, data, color="red", percs=None):
    """
    ax = matplotlib axes instance
    x = x-axis coordinates
    data = profile data
    color = color in any way matplotlib accepts
    """

    # Might change this into an argument for the function
    if percs is None:
        percs = [50, 90]
    percs = [(100 - float(p)) / 2 for p in percs[::-1]]
    alphas = [0.1, 0.4]

    # Convert to numpy array
    vals = np.array(data)

    # Get the median
    m = np.median(vals, axis=0)

    # Draw the minimum percentiles
    lines = [np.array([scoreatpercentile(vals[:, i], perc) for i in range(len(vals[0]))]) for perc in percs] + [m]
    for (line_min, line_max), alpha in zip([(lines[i], lines[i + 1]) for i in range(len(percs))], alphas):
        ax.fill_between(x, line_min, line_max, facecolor=color, alpha=alpha, edgecolor='face')

        # Draw the maximum percentiles
    lines = [m] + [np.array([scoreatpercentile(vals[:, i], 100 - perc) for i in range(len(vals[0]))]) for perc in
                   percs[::-1]]
    for (line_min, line_max), alpha in zip([(lines[i], lines[i + 1]) for i in range(len(percs))], alphas[::-1]):
        ax.fill_between(x, line_min, line_max, facecolor=color, alpha=alpha, edgecolor='face')

        # Draw the median
        ax.plot(x, m, color="black", alpha=0.95, linewidth=0.8)
        # ax.plot(x, mean(vals, axis = 0), color = "purple", alpha = 0.95, linewidth = 0.8)


def create_grid_figure(nrows, ncolumns, plotwidth=2.0, plotheight=2.0, pad=0.1, padleft=0.1, padright=0.1, padtop=0.1,
                       padbottom=0.1, clean=True):
    wsize = padleft + (ncolumns * plotwidth) + (pad * (ncolumns - 1)) + padright
    hsize = padtop + (nrows * plotheight) + (pad * (nrows - 1)) + padbottom
    fig = plt.figure(figsize=(wsize, hsize))
    wpadfraction = pad / wsize
    hpadfraction = pad / hsize
    wplotsize = plotwidth / wsize
    hplotsize = plotheight / hsize
    axes = {}
    # Create all the subplots
    for row in range(nrows):
        axes[row] = {}
        for col in range(ncolumns):
            axes[row][col] = plt.subplot(nrows, ncolumns, row * ncolumns + col + 1)

            # No labels, ticks, etc.
            if clean:
                for ax in [axes[row][col].xaxis, axes[row][col].yaxis]:
                    ax.set_major_formatter(NullFormatter())
                    ax.set_major_locator(NullLocator())

    # Resize all the subplots
    for row in range(nrows):
        for col in range(ncolumns):
            x0 = (padleft / wsize) + (wplotsize + wpadfraction) * col
            x1 = wplotsize
            y0 = (padbottom / hsize) + (nrows - row - 1) * (hplotsize + hpadfraction)
            y1 = hplotsize
            coords = [x0, y0, x1, y1]
            axes[row][col].set_position(coords)

            for s in axes[row][col].spines.values():
                s.set_linewidth(0.8)

    return fig, axes

def profile_screenshot(fname, interval, tracks, fontsize=None, colors=None, scalegroups=None, scale=None, show_scale=True, annotation=None, bgmode="color", fragmentsize=200,  dpi=600, rmdup=False, rmrepeats=False, reverse=False, adjscale=False):
    """
    Plot a genome browser like profile
    
    Parameters
    ----------
    fname: string
        output file name
    
    interval: string
        interval to plot in "chrom:start-end" format
    
    tracks: list
        list of filenames
    """
    if scalegroups is None:
        scalegroups = []

    if not fontsize:
        fontsize = FONTSIZE

    if not colors:
        colors = DEFAULT_COLORS
  
    
    # Plot size and padding definition
    plotwidth = 6
    plotheight = 0.3
    pad = {
            "left": 1.5,
            "right": 0.05,
            "top": 0.05,
            "bottom": 0.05,
            "row": 0,
            "column": 3,
            }
    
    # adjust width for track names if they are to long
    # kind of a quick hack
    max_len = 0
    for group in tracks:
        names = [os.path.splitext(os.path.basename(t))[0].strip() for t in group]
        l = sum([len(name) for name in names])
        if l > max_len:
            max_len = l
    if max_len > 27:
        pad["left"] = 3
    
    # Genomic scale
    scale_height = 0.1
    # Annotation track height
    annotation_height = 0.01
    
    chrom, start, end = re.split(r'[-:]', interval)
    start, end = int(start), int(end)
    
    if annotation:
        ann = load_annotation([chrom,start,end], annotation)
        if ann:
            annotation_height = 0.2 * len(ann.keys())
        else:
            annotation = False
    
    nrows = len(tracks)

    wsize = pad["left"] + plotwidth + pad["right"]
    hsize = pad["top"] + (nrows * plotheight) + (pad["row"] * (nrows - 1)) + pad["bottom"]
    hsize += scale_height + pad["row"] + annotation_height + pad["row"]

    # initialize figure 
    fig = plt.figure(figsize=(wsize, hsize))
 
    # initialize profile figure
    pfig = ProfileFigure(fig=fig, fontsize=fontsize, pad=pad)

    # add the genomic scale
    pfig.add_panel(ScalePanel())
   
    if type(scale) != type([]):
        scale = [scale]

    # add the signal tracks
    c = 0
    for group in tracks:
        for i,track in enumerate(group):
            panel = pfig.add_panel(
                    BamProfilePanel(track,
                        color = colors[c % len(colors)], 
                        bgmode = bgmode,
                        name = os.path.splitext(os.path.split(track)[-1])[0],
                        fragmentsize = fragmentsize,
                        rmrepeats = rmrepeats,
                        rmdup = rmdup,
                        adjscale = adjscale,
                        show_scale = show_scale,
                        ),
                    overlay= i != 0
                    )
            panel.ymax = scale[c % len(scale)]
            c += 1
    
    # add the annotation panel
    if annotation:
        pfig.add_panel(AnnotationPanel(annotation))
    
    pfig.plot([chrom, start, end], scalegroups=scalegroups, reverse=reverse)
    plt.savefig(fname, dpi=dpi)

class ProfileFigure(object):
    def __init__(self, fig=None, gs=None, fontsize=FONTSIZE, pad=None):
        self._panels = []
        if not fig:
            fig = plt.figure()
        self.fig = fig

        self.pad = {}
        if pad:
            self.pad.update(pad)

        relpad = {}
        for k in ["left", "right"]:
            relpad[k] = float(self.pad.get(k,0)) / fig.get_figwidth()        
        for k in ["top", "bottom"]:
            relpad[k] = float(self.pad.get(k,0)) / fig.get_figheight()        

        if gs:
            self.gs = gs
        else:
            gs = gridspec.GridSpec(1, 1)
            gs.update(
                    left=relpad["left"], 
                    right=1 - relpad["right"], 
                    top=1 - relpad["top"], 
                    bottom=relpad["bottom"], 
                    wspace=0, 
                    hspace=0
                    )
            self.gs = gs[0]

        self.font = FontProperties(size=fontsize / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])

    def _plot_panel_names(self, ax, panels):
        names = [p.name for p in panels]
        colors = ["black"]

        if len(names) > 1:
            tmp_names = []
            colors = []
            for name,color in zip(names, [p.color for p in panels]):
                tmp_names.append("= ")
                tmp_names.append(name + ", ")
                colors += [color,"black"]
            names = tmp_names
            names[-1] = names[-1].strip(", ")

        colortext(-0.01, 0.5,
                names,
                colors,
                ax=ax,
                horizontalalignment='right',
                verticalalignment="center",
#                    transform=ax.transAxes,
                clip_on=False,
                fontproperties=self.font)

    def plot(self, interval, scalegroups=None, reverse=False, **kwargs):
        if scalegroups is None:
            scalegroups = []
        
        for panels in self._panels:
            for panel in panels:
                panel._load_data(interval)

        gs0 = gridspec.GridSpecFromSubplotSpec(
                len(self._panels),
                1,
                subplot_spec=self.gs,
                height_ratios=[max([p.height for p in panels]) for panels in self._panels]
        )
        
        for panels in self._panels:
            if isinstance(panels[-1], BamProfilePanel):
                ymax = max([p.ymax for p in panels])
                for panel in panels:
                    panel.ymax = ymax
        
        if scalegroups and len(scalegroups) > 0:
            for group in scalegroups:
                ymax = max([self._panels[g][-1].ymax for g in group])
                for g in group:
                    for panel in self._panels[g]:
                        panel.ymax = ymax
        
        # These are quick hacks to to get the track groups to work
        for panels in self._panels:
            if len(panels) > 1:
                # Set the alpha for overlapping tracks
                for panel in panels:
                    panel.alpha = 0.5
                
        for i, panels in enumerate(self._panels):
            ax = plt.Subplot(self.fig, gs0[i])
            plt.subplots_adjust(bottom=0, top=1, left=0, right=1, hspace=0)
            
            # add track labels
            self._plot_panel_names(ax, panels)
            
            for panel in panels:
                panel._plot(ax, interval, fig=self.fig, reverse=reverse, odd=i % 2, font=self.font, **kwargs)
            
            self.fig.add_subplot(ax)

    def add_panel(self, panel, overlay=False):
        if overlay and len(self._panels) > 0:
            self._panels[-1].append(panel)
        else:
            self._panels.append([panel])
        return panel

class ProfilePanel(object):
    name = ""
    
    def hide_axes(self, axes):
        for ax in [axes.xaxis, axes.yaxis]:
            ax.set_major_formatter(NullFormatter())
            ax.set_minor_formatter(NullFormatter())
            ax.set_major_locator(NullLocator())
            ax.set_minor_locator(NullLocator())

        for s in axes.spines.values():
            s.set_color('none')

class BamProfilePanel(ProfilePanel):
    def __init__(self, bamfile, height=1, color=None, bgmode=None, alpha=None, fragmentsize=200, rmdup=True,
                 rmrepeats=True, **kwargs):
        self.height = height
        self.track = Track.load(bamfile, fragmentsize=fragmentsize, rmdup=rmdup, rmrepeats=rmrepeats)

        self.ymax = None
        self.bgmode = bgmode

        self.scalepm = kwargs.get("adjscale", False)
        self.show_scale = kwargs.get("show_scale", True)

        if color:
            self.color = color
        else:
            self.color = "#a7004b"

        if alpha:
            self.alpha = alpha
        else:
            self.alpha = 1

        self.fragmentsize = fragmentsize
        self.rmdup = rmdup
        self.rmrepeats = rmrepeats

        self.name = kwargs.get('name')

    def _load_data(self, interval):
        self.profile = self.track.get_profile(interval, 
                scalepm=self.scalepm)

        if not self.ymax:
            self.ymax = np.nanmax(self.profile) * 1.10

    def _plot(self, ax, interval, reverse=False, fig=None, odd=False, font=None, **kwargs):

        # Background of profile
        if self.bgmode == "stripes":
            bgcol = {0: "white", 1: (0.95, 0.95, 0.95)}[int(odd)]
            ax.set_facecolor(bgcol)
        elif self.bgmode == "color":
            ax.set_facecolor(self.color)
            ax.patch.set_alpha(0.07)

        # get interval
        chrom, start, end = interval
        profile = self.profile
        if reverse:
            profile = profile[::-1]

        # plot data
        ax.fill_between(
                range(start, end),
                np.zeros(len(profile)),
                profile,
                edgecolor='face',
                facecolor=self.color,
                linewidth=0.5,
                alpha=self.alpha)
        
        # set the y-limit
        ax.set_ylim(0, self.ymax)

        # add y-limit label
        if self.show_scale:
            ax.text(0.005, 0.90,
                int(ax.get_ylim()[-1] + 0.5),
                horizontalalignment='left',
                verticalalignment="top",
                transform=ax.transAxes,
                clip_on=False,
                fontproperties=font)
        
        ax.set_xlim(start, end)

        self.hide_axes(ax)

class AnnotationPanel(ProfilePanel):
    def __init__(self, annofile, height=0.3, vis="stack", color="black"):
        self.annofile = annofile
        self.height = height
        self.vis = vis
        self.color = color

    def _load_data(self, interval):
        self.gene_track = load_annotation(interval, self.annofile, vis=self.vis)
        self.max_tracks = len(self.gene_track.keys())
        self.height *= self.max_tracks

    def _plot(self, ax, interval, reverse=False, fig=None, odd=False, font=None, **kwargs):
        
        chrom, start, end = interval
        ax.set_ylim(- 1 * self.max_tracks, 0)
        for track_id, genes in self.gene_track.items():
            for gene in genes:
                h_gene = -1 * track_id - 0.5

                genestart = gene[1]
                geneend = gene[2]
                genename = gene[3]
                
                if len(gene) >= 6:
                    genestrand = gene[5]
                else:
                    genestrand = "+"

                # BED12 format
                if len(gene) == 12:

                    exonstarts = [int(x) for x in gene[11].split(",") if x]
                    exonsizes = [int(x) for x in gene[10].split(",") if x]
                else:
                    exonstarts = [0]
                    exonsizes = [geneend - genestart]

                x1 = (genestart - start)
                x2 = (geneend - start)
                if reverse:
                    x1 = end - genestart
                    x2 = end - geneend

                gstart = x1 / float(end - start)
                gend = x2 / float(end - start)

                # Horizontal line for complete gene
                ax.axhline(h_gene,
                           gstart,
                           gend,
                           color=self.color,
                           solid_capstyle="butt",
                           )
                # Exons 
                for exonstart, exonsize in zip(exonstarts, exonsizes):
                    estart = (genestart + exonstart - start)
                    eend = (genestart + exonstart + exonsize - start)
                    if reverse:
                        estart = end - (genestart + exonstart)
                        eend = end - (genestart + exonstart + exonsize)

                    ax.axhspan(
                            h_gene - 0.35,
                            h_gene + 0.35,
                            estart / float(end - start),
                            eend / float(end - start),
                            linewidth=0.1,
                            color=self.color)

                # Only draw arrows for BED12 entries
                if len(gene) == 12:
                    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                    figwidth, figheight = bbox.width, bbox.height

                    # Scale with absolute width of figure
                    step = 0.08 / figwidth

                    if reverse:
                        step = -step
                    for i in np.arange(gstart + step, gend - step, step):
                        if genestrand == "-":
                            astart = (i + step, h_gene)
                            aend = (i, h_gene)
                        else:
                            astart = (i, h_gene)
                            aend = (i + step, h_gene)

                        arr = FancyArrowPatch(
                                astart,
                                aend,
                                arrowstyle=GENE_ARROW,
                                mutation_scale=(figheight * fig.dpi) / 2 / self.max_tracks * 1.5,
                                linewidth=0.5,
                                color=self.color,
                        )
                        ax.add_patch(arr)
                if gstart > 0:
                    ax.text(gstart - 0.01, h_gene, genename, 
                            horizontalalignment="right",
                            verticalalignment="center",
                            fontproperties=font)
        
        self.hide_axes(ax)


class ScalePanel(ProfilePanel):
    def __init__(self, height=0.3, color=None, alpha=None):
        self.height = height

        if color:
            self.color = color
        else:
            self.color = "black"

        if alpha:
            self.alpha = alpha
        else:
            self.alpha = 1

    def _load_data(self, interval):
        pass

    def _plot(self, ax, interval, reverse=False, fig=None, odd=False, font=None, **kwargs):

        chrom, start, end = interval

        # Formatting
        for s in ax.spines.values():
            s.set_color('none')

        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_locator(NullLocator())
        ax.set_xlim(start, end)
        # ax.set_ylim(0,1)
        # Set font

        # Plot the numbers
        ticks = [s for s in ax.xaxis.get_ticklocs()[:-1] if s > start and s < end]
        xcoords = [(s - start) / (end - start) + 0.01 for s in ticks]

        if reverse:
            ticks = ticks[::-1]

        for s, x in zip(ticks[:-1], xcoords[:-1]):
            ax.text(
                    x,
                    0.5,
                    str(int(s)),
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes,
                    fontproperties=font,
                    color=self.color)
        ax.text(
                0,
                0.5,
                chrom,
                horizontalalignment='left',
                verticalalignment='center',
                transform=ax.transAxes,
                fontproperties=font,
                color=self.color)


class ConservationPanel(ProfilePanel):
    def __init__(self, track, target, height=1):
        self.track = track
        self.height = height
        self.data = []
        self.target = target

    def _load_data(self, ival1):
        for line in open(self.track):
            vals = line.strip().split("\t")
            for i in [1, 2, 4, 5]:
                vals[i] = int(vals[i])
            self.data.append(vals)

    def _plot(self, ax, interval, reverse=False, fig=None, odd=False, font=None, **kwargs):
        reverse_other = reverse
        reverse_self = kwargs.get("reverse_self", False)

        chrom, start, end = interval
        c2, s2, e2 = self.target
        span1 = float(end - start)
        span2 = float(e2 - s2)
        for [chrom1, start1, end1, chrom2, start2, end2] in self.data:
            if reverse_self:
                if reverse_other:
                    coords = [
                        [1 - (end1 - start) / span1, 1],
                        [1 - (end2 - s2) / span2, 0],
                        [1 - (start2 - s2) / span2, 0],
                        [1 - (start1 - start) / span1, 1]
                    ]
                else:
                    coords = [
                        [1 - (end1 - start) / span1, 1],
                        [(start2 - s2) / span2, 0],
                        [(end2 - s2) / span2, 0],
                        [1 - (start1 - start) / span1, 1]
                    ]
            else:
                if reverse_other:
                    coords = [
                        [(start1 - start) / span1, 1],
                        [1 - (end2 - s2) / span2, 0],
                        [1 - (start2 - s2) / span2, 0],
                        [(end1 - start) / span1, 1]
                    ]
                else:
                    coords = [
                        [(start1 - start) / span1, 1],
                        [(start2 - s2) / span2, 0],
                        [(end2 - s2) / span2, 0],
                        [(end1 - start) / span1, 1]
                    ]

            poly = Polygon(coords,
                           facecolor="black",
                           edgecolor='none',
                           alpha=0.2,
                           )
            ax.add_patch(poly)
        self.hide_axes(ax)
