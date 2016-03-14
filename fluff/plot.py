import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch, ArrowStyle, Polygon
from matplotlib.ticker import NullFormatter, NullLocator
from numpy import *
from scipy.stats import scoreatpercentile

from fluff.color import create_colormap
from fluff.config import FONTSIZE
from fluffio import *

DEFAULT_COLORS = ["#e41a1c", "#4daf4a", "#377eb8"]
PROFILE_MIN_Y = 75
GENE_ARROW = "->"
GENE_ARROW = ArrowStyle._Curve(beginarrow=False, endarrow=True, head_length=.4, head_width=.4)


def hide_axes(ax):
    for x in [ax.xaxis, ax.yaxis]:
        x.set_major_formatter(NullFormatter())
        x.set_major_locator(NullLocator())
    for loc, spine in ax.spines.iteritems():
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
        min_y, max_y = ylim
        s = 0
        plt.axhline(y=0,
                    color="grey",
                    linewidth=0.5,
                    alpha=0.5
                    )
        labels = array(labels)
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
        hide_axes(ax)
        ax.set_ylim(ylim)

    fig.subplots_adjust(wspace=btw_space, hspace=0)
    ext = outfile.split(".")[-1]
    if not ext in ["png", "svg", "ps", "eps", "pdf"]:
        outfile += ".png"
    sys.stderr.write("Saving figure\n")
    if outfile.endswith("png"):
        plt.savefig(outfile, dpi=600, bbox_inches='tight')
    else:
        plt.savefig(outfile)


def coverage_plot(ax, x, data, color="red", percs=[50, 90]):
    """
    ax = matplotlib axes instance
    x = x-axis coordinates
    data = profile data
    color = color in any way matplotlib accepts
    """

    # Might change this into an argument for the function
    percs = [(100 - float(p)) / 2 for p in percs[::-1]]
    alphas = [0.1, 0.4]

    # Convert to numpy array
    vals = array(data)

    # Get the median
    m = median(vals, axis=0)

    # Draw the minimum percentiles
    lines = [array([scoreatpercentile(vals[:, i], perc) for i in range(len(vals[0]))]) for perc in percs] + [m]
    for (line_min, line_max), alpha in zip([(lines[i], lines[i + 1]) for i in range(len(percs))], alphas):
        ax.fill_between(x, line_min, line_max, facecolor=color, alpha=alpha, edgecolor='face')

        # Draw the maximum percentiles
    lines = [m] + [array([scoreatpercentile(vals[:, i], 100 - perc) for i in range(len(vals[0]))]) for perc in
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


def profile_screenshot(fname, intervals, tracks, fontsize, colors=None, scalegroups=[], annotation=None, bgmode="color",
                       fragmentsize=200, scale=False, dpi=600, rmdup=False, rmrepeats=False, reverse=False):
    # Colors
    if not colors:
        colors = DEFAULT_COLORS
    font = FontProperties(size=fontsize / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])

    # Sizes
    #Adjust width based on titles length
    if max([len(os.path.splitext(os.path.basename(i[0]))[0].strip()) for i in tracks]) > 10:
        padleft = 1
    else:
        padleft = 0.1
    plotwidth = 6
    plotheight = 0.3
    padh = 0
    padw = 3
    #padleft = 0.1
    padright = 0.1
    padtop = 0.1
    padbottom = 0.1
    clean = True

    # Genomic scale
    scale_height = 0.1
    # Annotation track height

    annotation_height = 0.01
    gene_tracks = []
    if annotation:
        for interval in intervals:
            ann = load_annotation(interval, annotation)
            gene_tracks.append(ann)
        if gene_tracks[0]:
            max_tracks = max([len(x.keys()) for x in gene_tracks])
            annotation_height = 0.2 * max_tracks
        else:
            annotation = False
    #

    ncolumns = len(intervals)
    nrows = len(tracks)

    wsize = padleft + (ncolumns * plotwidth) + (padw * (ncolumns - 1)) + padright
    # Profile plots
    hsize = padtop + (nrows * plotheight) + (padh * (nrows - 1)) + padbottom
    hsize += scale_height + padh + annotation_height + padh

    fig = plt.figure(figsize=(wsize, hsize))
    padtopfraction = padtop / hsize
    padbottomfraction = padbottom / hsize
    wpadfraction = padw / wsize
    hpadfraction = padh / hsize
    wplotsize = plotwidth / wsize
    hplotsize = plotheight / hsize
    hannotation = annotation_height / hsize
    hscale = scale_height / hsize

    rellabelsize = 10.0

    for int_num, interval in enumerate(intervals):
        all_axes = []

        # Scale ax
        ax = plt.axes([padleft / wsize + (wplotsize + wpadfraction) * int_num + wplotsize / rellabelsize,
            1 - hscale - padtopfraction,
            wplotsize - wplotsize / rellabelsize,
            hscale], axisbg="white")
        all_axes.append(ax)

        # Profile axes
        for i in range(nrows):
            bgcol = "white"
            if bgmode == "stripes":
                bgcol = {0: "white", 1: (0.95, 0.95, 0.95)}[i % 2]
            elif bgmode == "color":
                bgcol = colors[i % len(colors)]

            xleft = padleft / wsize + (wplotsize + wpadfraction) * int_num
            ax = plt.axes([xleft,
                           1 - hscale - padtopfraction - (i + 1) * (hplotsize + hpadfraction),
                           wplotsize / rellabelsize,
                           hplotsize],
                          axisbg=bgcol)
            hide_axes(ax)
            ax.text(0.95, 0.5, os.path.splitext(os.path.basename(tracks[i][0]))[0].strip(), fontproperties=font, horizontalalignment="right", verticalalignment="center")
            if bgmode == "color":
                ax.patch.set_alpha(0.07)

            ax = plt.axes([padleft / wsize + (wplotsize + wpadfraction) * int_num + wplotsize / rellabelsize,
                           1 - hscale - padtopfraction - (i + 1) * (hplotsize + hpadfraction), wplotsize - wplotsize / rellabelsize, hplotsize],
                          axisbg=bgcol)
            if bgmode == "color":
                ax.patch.set_alpha(0.07)
            all_axes.append(ax)

        # Annotation axes
        ax = plt.axes(
                [padleft / wsize + (wplotsize + wpadfraction) * int_num + wplotsize / rellabelsize,
                    0 + padbottomfraction,
                    wplotsize - wplotsize / rellabelsize,
                    hannotation],
                axisbg="white")
        all_axes.append(ax)

        # No labels, ticks, etc.
        for axes in all_axes[1:]:
            for ax in [axes.xaxis, axes.yaxis]:
                ax.set_major_formatter(NullFormatter())
                ax.set_major_locator(NullLocator())

        for axes in all_axes:
            for s in axes.spines.values():
                s.set_color('none')

        chrom, start, end = interval

        # Format the genomic scale
        ax = all_axes[0]
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_locator(NullLocator())
        ax.set_xlim(start, end)

        for s in [s for s in ax.xaxis.get_ticklocs()[:-1] if s > start and s < end][:-1]:
            plt.text((s - start) / (end - start) + 0.01, 0.5, str(int(s)), horizontalalignment='left',
                     verticalalignment='center', transform=ax.transAxes, fontproperties=font)

        # Load the actual data
        profiles = load_profile(interval, tracks, fragmentsize=fragmentsize, rmdup=rmdup, rmrepeats=rmrepeats,
                                reverse=reverse)

        # Plot the profiles
        color_index = 0
        track_maxes = []
        for i, profile_group in enumerate(profiles):
            ax = all_axes[i + 1]
            if type(profile_group) == type([]):
                maxes = []
                for profile in profile_group:
                    if len(profile_group) > 1:
                        alpha = 0.4
                    else:
                        alpha = 1
                    ax.fill_between(range(start, end), zeros(len(profile)), profile, edgecolor='face',
                                    facecolor=colors[color_index % len(colors)], linewidth=0, alpha=alpha)
                    color_index += 1
                    maxes.append(numpy.nanmax(profile) * 1.1)
                track_maxes.append(max(maxes))
                ax.set_ylim(0, max(maxes))
            else:
                # Single track
                profile = profile_group
                ax.fill_between(range(start, end), zeros(len(profile)), profile, edgecolor='face',
                                facecolor=colors[color_index % len(colors)], linewidth=0)
                color_index += 1
                track_maxes.append(max(profile) * 1.1)
                ax.set_ylim(0, max(profile) * 1.1)

        for i, profile_group in enumerate(profiles):
            # Get maximum for this track based on scalegroups
            ylim_max = track_maxes[i]
            if scalegroups and len(scalegroups) > 0:
                for group in scalegroups:
                    if (i + 1) in group:
                        ylim_max = max([track_maxes[g - 1] for g in group]) * 1.1
                        break

            if ylim_max < PROFILE_MIN_Y:
                ylim_max = PROFILE_MIN_Y

            if type(scale) == type([]):
                ylim_max = scale[i]
            # Set maximum
            all_axes[i + 1].set_ylim(0, ylim_max)

            # Label scale
            if scale:
                all_axes[i + 1].text(0.005, 0.90, int(ylim_max + 0.5), horizontalalignment='left',
                                     verticalalignment="top", transform=all_axes[i + 1].transAxes, clip_on=False,
                                     fontproperties=font)

        # Plot the gene annotation
        if annotation:
            ax = all_axes[-1]
            ax.set_ylim(- 4 * max_tracks, 0)
            for track_id, genes in gene_tracks[int_num].items():
                for gene in genes:
                    h_gene = -4 * track_id - 2

                    genestart = gene[1]
                    geneend = gene[2]
                    exonstarts = [int(x) for x in gene[11].split(",") if x]
                    exonsizes = [int(x) for x in gene[10].split(",") if x]
                    genestrand = gene[5]
                    genename = gene[3]

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
                               color="black")

                    ax.text(gstart, h_gene - 3, genename, fontsize=fontsize)

                    # Exons 
                    for exonstart, exonsize in zip(exonstarts, exonsizes):
                        estart = (genestart + exonstart - start)
                        eend = (genestart + exonstart + exonsize - start)
                        if reverse:
                            estart = end - (genestart + exonstart)
                            eend = end - (genestart + exonstart + exonsize)

                        ax.axhspan(
                                h_gene - 0.5,
                                h_gene + 0.5,
                                estart / float(end - start),
                                eend / float(end - start),
                                color="black")

                    step = plotwidth / 600.0
                    if reverse:
                        step = -step
                    for i in arange(gstart + step, gend - step, step):
                        if genestrand == "-":
                            arr = FancyArrowPatch(
                                    (i + step, h_gene),
                                    (i, h_gene),
                                    arrowstyle=GENE_ARROW,
                                    mutation_scale=14,
                                    linewidth=0.5,
                            )

                        else:
                            arr = FancyArrowPatch(
                                    (i, h_gene),
                                    (i + step, h_gene),
                                    arrowstyle=GENE_ARROW,
                                    mutation_scale=14,
                                    linewidth=0.5,
                            )
                        ax.add_patch(arr)
    print 'Saving figure'




    plt.savefig(fname, dpi=dpi)
    plt.close()


class ProfileFigure():
    def __init__(self, fig=None, gs=None):
        self._panels = []
        if not fig:
            fig = plt.figure()
        self.fig = fig

        if gs:
            self.gs = gs
        else:
            gs = gridspec.GridSpec(1, 1)
            gs.update(left=0, right=1, wspace=0, hspace=0)
            self.gs = gs[0]

        self.font = FontProperties(size=FONTSIZE / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])

    def plot(self, interval, scalegroups=[], reverse=False, **kwargs):
        for panel in self._panels:
            panel._load_data(interval)

        gs0 = gridspec.GridSpecFromSubplotSpec(
                len(self._panels),
                1,
                subplot_spec=self.gs,
                height_ratios=[p.height for p in self._panels]
        )

        if scalegroups and len(scalegroups) > 0:
            for group in scalegroups:
                ymax = max([self._panels[g].ymax for g in group])
                for g in group:
                    self._panels[g].ymax = ymax

        # for p in self._panels:
        #    if hasattr(p, "ymax") and p.ymax < PROFILE_MIN_Y:
        #        p.ymax = PROFILE_MIN_Y

        for i, panel in enumerate(self._panels):
            ax = plt.Subplot(self.fig, gs0[i])
            plt.subplots_adjust(bottom=0, top=1, left=0, right=1, hspace=0)
            panel._plot(ax, interval, fig=self.fig, reverse=reverse, odd=i % 2, font=self.font, **kwargs)
            self.fig.add_subplot(ax)

    def add_panel(self, panel):
        self._panels.append(panel)
        return panel


class ProfilePanel():
    def hide_axes(self, axes):
        for ax in [axes.xaxis, axes.yaxis]:
            ax.set_major_formatter(NullFormatter())
            ax.set_major_locator(NullLocator())

        for s in axes.spines.values():
            s.set_color('none')


class BamProfilePanel(ProfilePanel):
    def __init__(self, bamfile, height=1, color=None, bgmode=None, alpha=None, fragmentsize=200, rmdup=True,
                 rmrepeats=True, **kwargs):
        self.height = height
        self.track = TrackWrapper(bamfile)

        self.ymax = None
        self.bgmode = bgmode

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

        self.profile = self.track.get_profile(
                interval,
                self.fragmentsize,
                self.rmdup,
                self.rmrepeats
        )

        if not self.ymax:
            self.ymax = numpy.nanmax(self.profile) * 1.10

    def _plot(self, ax, interval, reverse=False, fig=None, odd=False, font=None, **kwargs):

        # Background of profile
        if self.bgmode == "stripes":
            bgcol = {0: "white", 1: (0.95, 0.95, 0.95)}[int(odd)]
            ax.set_axis_bgcolor(bgcol)
        elif self.bgmode == "color":
            ax.set_axis_bgcolor(self.color)
            ax.patch.set_alpha(0.07)

        chrom, start, end = interval
        profile = self.profile
        if reverse:
            profile = profile[::-1]

        # r = np.split(range(start, end), np.where(profile == 0)[0])
        # p = np.split(profile, np.where(profile == 0)[0])
        ax.fill_between(
                range(start, end),
                zeros(len(profile)),
                profile,
                edgecolor='face',
                facecolor=self.color,
                linewidth=0.5,
                alpha=self.alpha)

        # for i,pos in zip(profile, range(start, end)):
        #    if i ! = 0:
        #        ax.axvline(
        #                pos, 
        #                ymin = 0, 
        #                ymax = i,
        #                linewidth = 0.5,
        #                color = self.color,
        #                alpha = self.alpha,
        #                )


        ax.set_ylim(0, self.ymax)

        ax.text(0.005, 0.90,
                int(ax.get_ylim()[-1] + 0.5),
                horizontalalignment='left',
                verticalalignment="top",
                transform=ax.transAxes,
                clip_on=False,
                fontproperties=font)

        if self.name:
            ax.text(0.005, 0.5,
                    self.name,
                    horizontalalignment='left',
                    verticalalignment="center",
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
                    for i in arange(gstart + step, gend - step, step):
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


if __name__ == "__main__":
    #    import matplotlib.pyplot as plt
    #
    #    data = [
    #        [10, 7, 5, 7, 10],
    #        [9, 7, 5, 7, 9],
    #        [8, 6, 5, 6, 8],
    #        [7, 6, 5, 6, 7],
    #        [6, 5, 5, 5, 6],
    #        [5, 5, 5, 5, 5],
    #        [4, 4, 5, 4, 4],
    #        [3, 4, 5, 4, 3],
    #        [2, 3, 5, 3, 2],
    #        [1, 3, 5, 3, 1],
    #    ]
    #    x = arange(5)
    #    plt.figure()
    #
    #    ax = plt.subplot(111)
    #
    #    coverage_plot(ax, x, data)
    #    plt.show()

    tracks = [
        "/home/simon/prj/xenopus/xenTro3b/bam/EZH2_stage9.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/Jarid2_stage9.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage9.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage12.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage16.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage30.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K4me1_stage9.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K4me3_stage9.bam",
        "/home/simon/prj/xenopus/xenTro3b/bam/H3K4me3_stage12.bam",
    ]
    intervals = [
        ["scaffold_3", 112047153, 112091309],
        #        ["scaffold_1", 20081849, 20104636],
        #        ["scaffold_1", 126564001, 126606826],
        #        ["scaffold_9", 20687459, 20911523],
    ]

    colors = ["#a65628", "#ff7f00", "#e41a1c", "#e41a1c", "#e41a1c", "#e41a1c", "#377eb8", "#4daf4a", "#4daf4a"]
    annotation = "/home/simon/prj/xenopus/xenTro3b/annotation/Xentr7_2_Stable_name.bed"
    scalegroups = [[1], [2], [3, 4, 5, 6], [7], [8, 9]]
    profile_screenshot("test.png", intervals, tracks, colors=colors, annotation=annotation, scalegroups=scalegroups,
                       bgmode="color")
