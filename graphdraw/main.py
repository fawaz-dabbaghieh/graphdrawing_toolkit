"""
The scripts I need to put together
1- Taking a BAM file of the graph alignments and producing keeping one alignment per node with coloring for chromosomes (Done)
2- Integrating that big script that takes a table of contig chrom start end and produce all kind of interesting plots (plot_contig_alignment_info) (Done)
3- Color graph with 2 samples/colors (color_graph_2_samples) (Done)
4- Output junction points between two different chromosomes (chromosome_junction_points)
5- Generating WC segment, make Python call Rscript plotting_wc_segment
6- Subgraph extraction (extract_subgraph_frmo_ref_coord) (I think some changes were made, should check the server version) (Done)
7- A subcommand that takes variation table (like the somatic BNDs or Insertions), and takes some flagging distance between these and outputs all of the graphs areas around such SV also as PNG (extract_sections.sh using extract_subgraph_coords)
7- Generating that SVG plot which shows nodes against the reference with edges
    7.1- need another script that filters out a BAM file into a TSV that chooses one best alignment for a node (done)
    7.2- Add subcommand for drawing contigs (done)
    7.3- Add subcommand for drawing graph against reference (on going)
    7.4- make it work for both BAM and PAF, and both with or without chaining
    7.5- Add ticks to reference with some small numbers
    7.6- Highlight overlaps in contigs (Maybe not, because the overlaps can be seen)

8- Needs a setup.py and documentation

- Check a graph's overlap (low priority)
"""

import sys
import os
import pysam
import argparse
import logging
import pdb
from matplotlib.backends.backend_pdf import PdfPages
from graphdraw.color_graph_with_chromosomes import filter_and_output
from graphdraw.extract_subgraph_from_ref_coord import extract_subgraph
from graphdraw.draw_contig_graph import draw_chained_contigs, draw_unchained_contigs
from graphdraw.plot_depth import output_windowed_table, plot_chromosome
import graphdraw.draw_subgraph_to_ref_svg


################# helper functions
def sys_exit():
    print("Error Happened! Please check the log file")
    sys.exit()


def check_in_file(in_file):
    if not os.path.exists(in_file):
        logging.error(f"File {in_file} does not exist, please check the input")
        sys_exit()


def check_out_file(out_file):
    if os.path.exists(out_file):
        logging.error(f"Output file {out_file} already exist, please give another output, will not overwrite")
        sys_exit()
    else:
        try:
            outfile = open(out_file, "w")
            outfile.close()
            os.remove(out_file)
        except MemoryError:
            logging.error(f"There was a memory error writing {out_file}")
            sys_exit()
        except PermissionError:
            logging.error(f"Permission denied to write {out_file}")
            sys_exit()


def exit_error(message):
    logging.error(message)
    sys_exit()


#################################  Setting up parser

parser = argparse.ArgumentParser(description='Graph Toolkit for cancer graphs', add_help=True)
subparsers = parser.add_subparsers(help='Available subcommands', dest="subcommands")

parser._positionals.title = 'Subcommands'
parser._optionals.title = 'Global Arguments'

parser.add_argument("--log_file", dest="log_file", type=str, default="log.log",
                    help="The name/path of the log file. Default: log.log")

parser.add_argument("--log", dest="log_level", type=str, default="DEBUG",
                    help="The logging level [DEBUG, INFO, WARNING, ERROR, CRITICAL]")

########################## Color graph with chromosomes ###############################
#         print("Please give the input sam or bam file and output csv file name")
color_graph = subparsers.add_parser('color_graph', help='Colors each node based on the chromosome it aligns to')

color_graph.add_argument("-b", "--bam", dest="in_bam", default=None, type=str,
                         help="The input alignments of the graph nodes, BAM or SAM")

color_graph.add_argument("--out_csv", dest="out_csv", default=None, type=str,
                         help="The Bandage compatible CSV for coloring the graph")

color_graph.add_argument("--out_tsv", dest="out_tsv", default=None, type=str,
                         help="Outputs a TSV that works for drawing graph against reference subcommand")

########################## Plot contig alignment info ###############################
# print("You need to give input tsv with contigID-chr-start-end and a prefix for plots produced")

plot_contigs = subparsers.add_parser('plot_contigs', help='Generates different plots for the contig alignments')

plot_contigs.add_argument("-b", "--bam", dest="in_bam", default=None, type=str,
                          help="Alignments of contigs in BAM or SAM")

plot_contigs.add_argument("--bam2", dest="in_bam2", default=None, type=str,
                          help="In case of 2 haplotypes, give the first haplotype with --bam and second "
                               "haplotype with --bam2")

plot_contigs.add_argument("--in_table", dest="in_table", default=None, type=str,
                          help="table for plotting, e.g. CNV or depth values with header containing "
                               "chromosome, start, end, and value")

plot_contigs.add_argument("--out_prefix", dest="out_prefix", default="", type=str,
                          help="Prefix for the output PDF")

########################## Plot contig alignment info ###############################
# Color graph with 2 samples/colors (color_graph_2_samples)
# print("You need to give the <gaf1> <color> <gaf2> <color> <preifx>")

two_colors = subparsers.add_parser('two_colors', help='Color graph with 2 different alignments')

two_colors.add_argument("--gaf1", dest="in_gaf1", default=None, type=str,
                        help="The input alignment gaf of the first sample")

two_colors.add_argument("--gaf2", dest="in_gaf2", default=None, type=str,
                        help="The input alignment gaf of the second sample")

two_colors.add_argument("--color1", dest="color1", default=None, type=str,
                        help="The color of the first sample, compatible with X11 colors or HEX")

two_colors.add_argument("--color2", dest="color2", default=None, type=str,
                        help="The color of the second sample, compatible with X11 colors or HEX")

two_colors.add_argument("--out_prefix", dest="out_prefix", default="", type=str,
                        help="The output Bandage compatible CSV prefix")

########################## extract subgraph ###############################
# Color graph with 2 samples/colors (color_graph_2_samples)
#     print("you need to give the <gfa> <bam> <chr:start-end> <neighborhood_size> <output.gfa>")
subgraph = subparsers.add_parser('subgraph', help='Extract subgraph using chromosome coordinates')

subgraph.add_argument("-g", "--gfa", dest="in_gfa", default=None, type=str,
                      help="The input graph that you want to extract the subgraph from")

subgraph.add_argument("-b", "--bam", dest="in_bam", default=None, type=str,
                      help="The alignment of the graph nodes against a reference")

subgraph.add_argument("--coord", dest="in_coord", default=None, type=str,
                      help="The reference coordinate to extract the subgraph, in samtools format, e.g. chr1:1000000-2000000")

subgraph.add_argument("--sg_size", dest="sg_size", default=20, type=int,
                      help="The size of the neighborhood of the subgraph, i.e. the nodes around the nodes in the coordinates given")

subgraph.add_argument("-o", "--out_gfa", dest="out_gfa", default=None, type=str,
                      help="The output subgraph in GFA format")

########################## draw_contig ###############################
# Color graph with 2 samples/colors (color_graph_2_samples)
#     print("you need to give the <gfa> <bam> <chr:start-end> <neighborhood_size> <output.gfa>")
draw_contig = subparsers.add_parser('draw_contig',
                                    help='Visualizes contig alignments against a reference using BAM, SAM, or PAF information')

draw_contig.add_argument("-a", "--in_align", dest="in_a", default=None, type=str,
                         help="Input alignment file, SAM, BAM, or PAF")

draw_contig.add_argument("-c", "--contig", dest="ctg_name", default=None, type=str,
                         help="The contig name present in the alignment file that you would like to visualize")

draw_contig.add_argument("--out_prefix", dest="svg_out", default=None, type=str,
                         help="Output svg name, default: [contig_name]")

draw_contig.add_argument("--max_dist", dest="max_dist", default=10000, type=int,
                         help="Maximum distance allowed between alignment to chain them together, default: 10,000 bp")

draw_contig.add_argument("--max_overlap", dest="max_overlap", default=30, type=int,
                         help="Maximum allowed overlap percentage between two alignments to chain, default: 20 (0-100)")

draw_contig.add_argument("--svg_width", dest="svg_width", default=600, type=int,
                         help="SVG width, default: 600")

draw_contig.add_argument("--svg_height", dest="svg_height", default=300, type=int,
                         help="SVG height, default: 300")

draw_contig.add_argument("--svg_legend", dest="svg_legend", default=10, type=int,
                         help="SVG legend, default: 10")

draw_contig.add_argument("--no_chaining", dest="no_chaining", action="store_true",
                         help="Doesn't use the chaining method and draws all the alignments")

########################## subgraph ###############################
# Color graph with 2 samples/colors (color_graph_2_samples)
#     print("you need to give the <gfa> <bam> <chr:start-end> <neighborhood_size> <output.gfa>")
draw_subgraph = subparsers.add_parser('draw_subgraph',
                                      help='visualizes a subgraph node alignments against the reference')

draw_subgraph.add_argument("-a", "--in_align", dest="in_a", default=None, type=str,
                           help="Input alignment file, SAM, BAM")

draw_subgraph.add_argument("-g", "--graph", dest="in_gfa", default=None, type=str,
                           help="The graph file in GFA format")

draw_subgraph.add_argument("--coord", dest="coord", default=None, type=str,
                           help="Coordinates in the Ref to extract the subgraph, same as samtools view format "
                                "chrX:start-end")

draw_subgraph.add_argument("--out_svg", dest="svg_out", default=None, type=str,
                           help="Output svg file name, default ")

draw_subgraph.add_argument("--svg_width", dest="svg_width", default=600, type=int,
                           help="SVG width, default: 600")

draw_subgraph.add_argument("--svg_height", dest="svg_height", default=300, type=int,
                           help="SVG height, default: 300")

draw_subgraph.add_argument("--svg_legend", dest="svg_legend", default=10, type=int,
                           help="SVG legend, default: 10")


########################## subgraph ###############################
# Color graph with 2 samples/colors (color_graph_2_samples)
#     print("you need to give the <gfa> <bam> <chr:start-end> <neighborhood_size> <output.gfa>")
depth = subparsers.add_parser('depth', help='Produce a windows depth table and depth plots')

depth.add_argument("-d", "--depth_table", dest="in_depth", default=None, type=str,
                   help="Input Samtools depth table")

depth.add_argument("--bin_size", dest="bin_size", default=1, type=float,
                   help="The bin size as a percentage of the chromosome length , default: 1")

depth.add_argument("--prefix", dest="depth_prefix", default=None, type=str,
                   help="The prefix of the output table or figure")

depth.add_argument("--output_table", dest="out_depth_table", action="store_true",
                   help="If given, then the it will output a TSV table with the windows and their values")

depth.add_argument("--output_plot", dest="out_depth_plot", action="store_true",
                   help="If given, then it will output a PDF file with depth plots for each chromosome")

depth.add_argument("--reference_name", dest="reference_name", default="t2t", type=str,
                   help="Reference used to show centromere and gaps, values: t2t, hg19, hg38. Default: T2T")

######################################################################################################
args = parser.parse_args()
# log_file = "log_" + str(time.clock_gettime(1)).split(".")[0] + ".log"
log_file = args.log_file


def main():
    if len(sys.argv) == 1:
        print("You did not provide any arguments\n"
              "Try to use -h or --help for help\n")
        sys.exit(0)

    logging.basicConfig(filename=log_file, filemode='w',
                        format='[%(asctime)s] %(message)s',
                        level=getattr(logging, args.log_level.upper()))
    logging.info(" ".join(["argument given:"] + sys.argv))

    ########################## Color graph with chromosomes ###############################
    if args.subcommands == "color_graph":

        if args.in_bam == None:
            exit_error("You need to give an input bam, -b, --bam")
        else:
            check_in_file(args.in_bam)

        if args.out_csv == None and args.out_tsv == None:
            exit_error("You need to give an output csv, -o, --out_csv, --out_tsv")
        else:
            if args.out_csv:
                check_out_file(args.out_csv)
            elif args.out_tsv:
                check_out_file(args.out_tsv)

        if args.out_csv:
            filter_and_output(args.in_bam, args.out_csv, tsv=False)
        elif args.out_tsv:
            filter_and_output(args.in_bam, args.out_tsv, tsv=True)

    ########################## Plot contig alignment info ###############################

    if args.subcommands == "plot_contigs":
        # plotting has heavier dependencies and having the import at the beginning gives some delays
        # even when plotting is not required
        from graphdraw.plot_contig_alignment_info import prepare_bam, plotting

        if args.in_bam is None:
            exit_error("You need to give an input bam, -b, --bam")
        else:
            check_in_file(args.in_bam)

        out_file = args.out_prefix + "_alignment_plots.pdf"
        if args.out_prefix is None:
            exit_error("You need to give a prefix with --out_prefix")
        else:
            check_out_file(out_file)

        # haplotype1 = prepare_bam(args.in_bam)
        print("finished loading the data from the bam, and now plotting")
        plotting(args.in_bam, out_file, args.out_prefix, in_table=args.in_table, in_bam2=args.in_bam2)
        # if args.in_bam2 is None:
        #     plotting(args.in_bam, out_file, args.out_prefix)
        # else:
        #     # haplotype2 = prepare_bam(args.in_bam2)
        #     plotting(args.in_bam, out_file, args.out_prefix, in_bam2=args.in_bam2)

    ############################ subgraph extraction ############################################
    # subgraph = subparsers.add_parser('subgraph', help='Extract subgraph using chromosome coordinates')

    # subgraph.add_argument("--sg_size", dest="sg_size", default=20, type=int,
    #     help="The size of the neighborhood of the subgraph, i.e. the nodes around the nodes in the coordinates given")

    # subgraph.add_argument("-o", "--out_gfa", dest="out_gfa", default=None, type=str,
    #     help="The output subgraph in GFA format")

    if args.subcommands == "subgraph":
        if args.in_gfa == None:
            exit_error("You need to give an input gfa, -g, --gfa")
        else:
            check_in_file(args.in_gfa)

        if args.in_bam == None:
            exit_error("You need to give an input BAM, -b, --bam")
        else:
            check_in_file(args.in_bam)

        if args.in_coord == None:
            exit_error("You need to give coordinates chrX:start-end, --coord")
        else:
            chrom, coord = args.in_coord.split(":")
            start, end = coord.split("-")
            start = int(start.replace(",", ""))
            end = int(end.replace(",", ""))

        # if args.out_gfa == None:
        #     exit_error("You need to give an output csv, -o, --out_csv")
        # else:
        #     check_out_file(args.out_csv)

        extract_subgraph(args.in_gfa, args.in_bam, chrom, start, end, args.sg_size, args.out_gfa)

    ########################## draw contig ###############################
    if args.subcommands == "draw_contig":

        if args.in_a is None:
            exit_error("You need to give an input bam, -b, --bam")
        else:
            check_in_file(args.in_a)

        if args.ctg_name is None:
            exit_error("You need to give a contig name that is present in the alignment file")

        if args.max_dist < 0:
            exit_error("Maximum distance cannot be smaller than 0")

        if args.max_overlap > 100 or args.max_overlap < 0:
            exit_error("Maximum overlap should be between 0 and 100")

        if args.svg_out is None:
            args.svg_out = args.ctg_name

        if args.no_chaining:
            draw_unchained_contigs(args.in_a, args.ctg_name, args.svg_out, args.svg_width, args.svg_height,
                                   args.svg_legend)
        else:
            draw_chained_contigs(args.in_a, args.ctg_name, args.svg_out, args.max_dist, args.max_overlap,
                                 args.svg_width, args.svg_height, args.svg_legend)
        # if args.in_a.endswith("paf"):
        #     draw_chained_contigs(args.in_a, args.ctg_name, args.svg_out, args.max_dist, args.max_overlap,
        #                          args.svg_width, args.svg_height, args.svg_legend)
        #
        # if args.no_chaining and args.in_a.endswith("bam"):
        #     draw_unchained_contigs(args.in_a, args.ctg_name, args.svg_out, args.svg_width, args.svg_height,
        #                            args.svg_legend)

    ########################## draw subgraph ###############################
    if args.subcommands == "draw_subgraph":

        if args.in_a is None:
            exit_error("You need to give an input bam or sam, -b, --bam")
        else:
            check_in_file(args.in_a)

        if args.in_gfa is None:
            exit_error("You need to give an input Graph in GFA format")
        else:
            check_in_file(args.in_gfa)

        if args.svg_out is None:
            args.svg_out = args.coord + ".svg"

        if args.coord is None:
            exit_error("You need to give coordinates in format chrx:start-end that exist in the alignment file given")

        graphdraw.draw_subgraph_to_ref_svg.draw_subgraph(args.in_a, args.in_gfa, args.coord, args.svg_out,
                                                           args.svg_width, args.svg_height, args.svg_legend)

        ########################## depth plots and table ###############################
        # depth.add_argument("-d", "--depth_table", dest="in_depth", default=None, type=str,
        #                    help="Input Samtools depth table")
        #
        # depth.add_argument("--bin_size", dest="bin_size", default=1, type=float,
        #                    help="The bin size as a percentage of the chromosome length , default: 1")
        #
        # depth.add_argument("--prefix", dest="depth_prefix", default=None, type=str,
        #                    help="The prefix of the output table or figure")
        #
        # depth.add_argument("--output_table", dest="out_depth_table", action="store_true",
        #                    help="If this option given then the a TSV with windows and depth "
        #                         "is outputted, otherwise a plot is outputted")
        #
        # depth.add_argument("--reference_name", dest="reference_name", default="t2t", type=str,
        #                    help="Reference used to show centromere and gaps")
    if args.subcommands == "depth":
        from graphdraw.constants import REF_LEN
        CHROMOSOMES = REF_LEN.keys()

        if (not args.out_depth_table) and (not args.out_depth_plot):
            print("You need to either choose to output a plot of a table")
            sys.exit(1)

        if args.in_depth is None:
            exit_error("You need to give an input samtools depth table")
        else:
            check_in_file(args.in_depth)

        if args.depth_prefix is None:  # using the output file for the prefix
            args.depth_prefix = args.in_depth.split(os.sep).split(".")[0]

        try:
            bin_perc = float(args.bin_size)
            bin_perc = bin_perc / 100
        except ValueError:
            logging.error("The percentage needs to be an integer")
            sys.exit()
        counter = 0
        if args.out_depth_table:
            out_file = args.depth_prefix + ".tsv"
            check_out_file(out_file)
            with open(out_file, "w") as outfile:  # adding header to the file
                outfile.write("\t".join(["chromosome", "start", "end", "value"]) + "\n")
        else:
            out_file = args.depth_prefix + ".pdf"
            pdf_file = PdfPages(out_file)
        counts = dict()
        with open(args.in_depth, "r") as infile:
            for l in infile:
                l = l.strip().split()
                counter += 1
                if counter % 10000000 == 0:
                    logging.info(f"Processed 10 million lines for chromosome {l[0]}")
                try:
                    # dealing with positions that don't have counts
                    while len(counts[l[0]]) < int(l[1]) - 1:
                        counts[l[0]].append(0)
                    counts[l[0]].append(int(l[2]))
                except KeyError:
                    if counts:  # finished a chromosome and now plotting
                        key = list(counts.keys())[0]
                        # chromosome with no counts
                        if len(counts[key]) == 0:
                            continue
                        # number of counts is too small
                        if int(len(counts[key]) * bin_perc) == 0:
                            continue
                        if key not in CHROMOSOMES:
                            continue
                        if args.out_depth_table:
                            output_windowed_table(counts[key], key, bin_perc, out_file)
                        else:  # plotting
                            fig = plot_chromosome(counts[key], key, args.reference_name, bin_perc, args.depth_prefix)
                            pdf_file.savefig(fig)
                        counts = dict()
                        counts[l[0]] = [0] * int(l[1])
                        counts[l[0]].append(int(l[2]))
                    else:
                        counts[l[0]] = [0] * int(l[1])
                        counts[l[0]].append(int(l[2]))
        if counts:
            key = list(counts.keys())[0]
            if key in CHROMOSOMES:
                if args.out_depth_table:
                    output_windowed_table(counts[key], key, bin_perc, out_file)
                else:
                    fig = plot_chromosome(counts[key], key, args.reference_name, bin_perc, args.depth_prefix)
                    pdf_file.savefig(fig)

        if not args.out_depth_table:
            pdf_file.close()


if __name__ == "__main__":
    main()
