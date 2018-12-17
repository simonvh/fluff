# __author__ = 'george'
import argparse, sys

from fluff.commands.heatmap import heatmap
from fluff.commands.bandplot import bandplot
from fluff.commands.profile import profile
from fluff.color import DEFAULT_COLORS
from fluff.config import *
from fluff.fluffio import *

def parse_cmds():

    description = """
    fluff v{0}
    """.format(FL_VERSION)

    epilog = """
    commands:
        heatmap      Produce a heatmap
        bandplot     Show the profiles as bandplots
        profile      Genome Browser screenshot

    type `fluff <command> -h` for more details
    """

    usage = "%(prog)s [-h] [options]"

    parser = argparse.ArgumentParser(
            usage=usage,
            description=description,
            epilog=epilog,
            formatter_class=argparse.RawDescriptionHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='subcommand_name')
    # heatmap subparser
    p = subparsers.add_parser('heatmap', add_help=False)
    # heatmap Required arguments
    req_grp = p.add_argument_group(title='Required arguments')
    req_grp.add_argument("-f",
                         required=True,
                         dest="featurefile",
                         help="BED file containing features",
                         metavar="FILE",
                         default=None)
    req_grp.add_argument("-d",
                         required=True,
                         dest="datafiles",
                         help="data files (reads in BAM or BED format)",
                         metavar="FILE",
                         nargs='*')
    req_grp.add_argument("-o",
                         required=True,
                         dest="outfile",
                         help="output file (type determined by extension)",
                         metavar="name",
                         default=None)
    # heatmap Optional arguments
    clu_grp = p.add_argument_group(title='Clustering')
    clu_grp.add_argument("-C",
                         dest="clustering",
                         help="kmeans, hierarchical or none",
                         default=DEFAULT_CLUSTERING,
                         metavar="METHOD")
    clu_grp.add_argument("-k",
                         dest="numclusters",
                         help="number of clusters",
                         metavar="INT",
                         type=int,
                         default=1)
    clu_grp.add_argument("-g",
                         dest="graphdynamics",
                         help="identify dynamics",
                         action="store_true",
                         default=False)
    clu_grp.add_argument("-M",
                         dest="distancefunction",
                         help="cluster metric: Euclidean or Pearson (default: Euclidean)",
                         default=DEFAULT_METRIC,
                         metavar="METHOD")
    clu_grp.add_argument("-p",
                         dest="pick",
                         help="pick specific data files to use for clustering",
                         default=None,
                         type=str)
    clu_grp.add_argument("-S",
                         dest="seed",
                         help="random seed (int) to use for K-means clustering",
                         default=None,
                         type=int)
    dap_grp = p.add_argument_group(title='Data processing')
    dap_grp.add_argument("-r",
                         dest="rpkm",
                         help="normalize using RPKM instead of read counts",
                         action="store_true",
                         default=False)
    dap_grp.add_argument("-e",
                         dest="extend",
                         help="extend (in bp, default: {0})".format(DEFAULT_EXTEND),
                         metavar="INT",
                         type=int,
                         default=DEFAULT_EXTEND)
    dap_grp.add_argument("-B",
                         dest="binsize",
                         help="bin size (default: {0})".format(DEFAULT_BINSIZE),
                         metavar="INT",
                         type=int,
                         default=DEFAULT_BINSIZE)
    dap_grp.add_argument("-F",
                         dest="fragmentsize",
                         help="fragment length (default: read length)",
                         type=int,
                         default=None)
    dap_grp.add_argument("-D",
                         dest="rmdup",
                         help="keep duplicate reads (removed by default)",
                         default=True,
                         action="store_false")
    dap_grp.add_argument("-R",
                         dest="rmrepeats",
                         help="keep reads with mapq 0 (removed by default) ",
                         action="store_false",
                         default=True)
    dap_grp.add_argument("-m",
                         dest="merge_mirrored",
                         help="merge mirrored clusters (only with kmeans and without -g option)",
                         default=False,
                         action="store_true")
    dap_grp.add_argument("-s",
                         dest="scale",
                         help="scale (absolute or percentage)",
                         type=str,
                         default=DEFAULT_SCALE)

    vis_grp = p.add_argument_group(title='Visualization')
    vis_grp.add_argument("-c",
                         dest="colors",
                         help="color(s) (name, colorbrewer profile or hex code), separate each color name by comma",
                         metavar="NAME(S)",
                         default=DEFAULT_COLORS)
    vis_grp.add_argument("-b",
                         dest="bgcolors",
                         help="background color(s) (name, colorbrewer profile or hex code)",
                         metavar="NAME(S)",
                         default=DEFAULT_BG)
    vis_grp.add_argument("-T",
                         dest="textfontsize",
                         help="text font size(default: {0})".format(FONTSIZE),
                         type=int,
                         default=FONTSIZE)
    vis_grp.add_argument("--no-colorbar",
                         dest="colorbar",
                         help="don't show colorbars of each heatmap",
                         action="store_false",
                         default=True)
                         
    opt_grp = p.add_argument_group(title='Other')
    opt_grp.add_argument("-P",
                         dest="cpus",
                         help="number of CPUs (default: 4)",
                         metavar="INT",
                         type=int,
                         default=4)
    opt_grp.add_argument("-h", "--help",
                         dest="help",
                         help="show this help message and exit",
                         action="help")

    p.set_defaults(func=heatmap)

    # bandplot subparser
    p = subparsers.add_parser('bandplot', add_help=False)
    # bandplot Required arguments
    req_grp = p.add_argument_group(title='Required arguments')
    req_grp.add_argument("-f",
                         required=True,
                         dest="clust_file",
                         help="BED file with cluster in 5th column",
                         metavar="FILE",
                         default=None)
    req_grp.add_argument("-d",
                         dest="datafiles",
                         help="data files (reads in BAM or BED format)",
                         metavar="FILE",
                         nargs='*',
                         default=None)

    req_grp.add_argument("-counts",
                         dest="readCount",
                         help="read counts table (instead of data files)",
                         metavar="FILE",
                         default=None)

    req_grp.add_argument("-o",
                         required=True,
                         dest="outfile",
                         help="output file (type determined by extension)",
                         metavar="name",
                         default=None)

    # bandplot Optional arguments
    dap_grp = p.add_argument_group(title='Data processing')
    dap_grp.add_argument("-r", dest="rpkm",
                         help="normalize using RPKM instead of read counts",
                         action="store_true",
                         default=False)
    dap_grp.add_argument("-S",
                         dest="summary",
                         help="create summary graphs",
                         default=False,
                         action="store_true")
    dap_grp.add_argument("-b", dest="bins",
                         help="number of bins",
                         metavar="INT",
                         default=BINS,
                         type=int)
    dap_grp.add_argument("-F", dest="fragmentsize",
                         help="fragment length (default: read length)",
                         type=int,
                         default=None)
    dap_grp.add_argument("-D", dest="rmdup",
                         help="keep duplicate reads (removed by default)",
                         default=True,
                         action="store_false")
    dap_grp.add_argument("-R", dest="rmrepeats",
                         help="keep repeats with mapq 0 (removed by default) ",
                         action="store_false",
                         default=True)
    dap_grp.add_argument("-s", dest="scalegroups",
                         help="scale groups",
                         metavar="GROUPS")
    dap_grp.add_argument("-p", dest="percs",
                         help="range of percentiles (default 50,90)",
                         metavar="INT,INT", default="50,90")
    dap_grp.add_argument("-P", dest="scalar",
                         help="percentile at which to extract score (should be in range [0,100], default: 99)",
                         metavar="INT",
                         default=SCALAR,
                         type=float)
    vis_grp = p.add_argument_group(title='Visualization')
    vis_grp.add_argument("-c", dest="colors",
                         help="color(s) (name, colorbrewer profile or hex code), separate each color name by comma",
                         metavar="NAME(S)",
                         default=DEFAULT_COLORS)
    vis_grp.add_argument("-T",
                         dest="textfontsize",
                         help="text font size(default: {0})".format(FONTSIZE),
                         type=int,
                         default=FONTSIZE)
    opt_grp = p.add_argument_group(title='Other')
    opt_grp.add_argument("-h", "--help",
                         dest="help",
                         help="show this help message and exit",
                         action="help")

    p.set_defaults(func=bandplot)

    # profile subparser
    p = subparsers.add_parser('profile', add_help=False)
    # profile Required arguments
    req_grp = p.add_argument_group(title='Required arguments')
    req_grp.add_argument("-i",
                         dest="interval",
                         help="genomic interval (chrom:start-end)",
                         metavar="INTERVAL")
    req_grp.add_argument("-d",
                         required=True,
                         dest="datafiles",
                         help="data files (reads in BAM or BED format)",
                         metavar="FILE",
                         nargs='*')
    req_grp.add_argument("-o",
                         required=True,
                         dest="outfile",
                         help="output file (type determined by extension)",
                         metavar="name",
                         default=None)
    # profile Optional arguments
    dap_grp = p.add_argument_group(title='Data processing')
    dap_grp.add_argument("-n",
                         dest="adjscale",
                         help="normalize to per million mapped reads",
                         default=False,
                         action="store_true")
    dap_grp.add_argument("-a", dest="annotation",
                         help="annotation in BED12 format",
                         metavar="FILE")
    dap_grp.add_argument("-t",
                         dest="trackgroups",
                         help="track groups",
                         metavar="GROUPS")
    dap_grp.add_argument("-s",
                         dest="scalegroups",
                         help="scale groups",
                         metavar="GROUPS")
    dap_grp.add_argument("-S", dest="scale",
                         help="scale, one value or comma-separated values for each track",
                         metavar="INT",
                         default=None
                         )
    dap_grp.add_argument("-f", dest="fragmentsize",
                         help="fragment length (default: %s)" % FRAGMENTLENGTH,
                         type=int,
                         default=FRAGMENTLENGTH)
    dap_grp.add_argument("-D",
                         dest="rmdup",
                         help="keep duplicate reads (removed by default)",
                         default=True,
                         action="store_false")
    dap_grp.add_argument("-R",
                         dest="rmrepeats",
                         help="keep repeats with mapq 0 (removed by default) ",
                         action="store_false",
                         default=True)
    dap_grp.add_argument("-r",
                         dest="reverse",
                         help="reverse ",
                         action="store_true",
                         default=False)
    vis_grp = p.add_argument_group(title='Visualization')
    vis_grp.add_argument("-c", dest="colors",
                         help="color(s) (name, colorbrewer profile or hex code), separate each color name by comma",
                         metavar="NAME(S)",
                         default=DEFAULT_COLORS)
    vis_grp.add_argument("-b", dest="background",
                         help="background color: stripes | white | color",
                         default="stripes")
    vis_grp.add_argument("-T",
                         dest="textfontsize",
                         help="text font size(default: {0})".format(FONTSIZE),
                         type=int,
                         default=FONTSIZE)
    vis_grp.add_argument("-H",
                         dest="show_scale",
                         help="hide track scale label",
                         action="store_false",
                         default=True)
    opt_grp = p.add_argument_group(title='Optional arguments')
    opt_grp.add_argument("-h", "--help",
                         dest="help",
                         help="show this help message and exit",
                         action="help")
    p.set_defaults(func=profile)

    return parser

def main():
    """entry point"""
    parser = parse_cmds()
    args = parser.parse_args()

    subcommand = args.subcommand_name
    if subcommand == 'heatmap':
        heatmap(args)
    elif subcommand == 'bandplot':
        bandplot(args)
    elif subcommand == 'profile':
        profile(args)
    else:
        parser.print_help()

    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me! ;-) Bye!\n")
        sys.exit(0)
