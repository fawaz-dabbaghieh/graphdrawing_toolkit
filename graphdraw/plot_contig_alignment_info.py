"""
Part of this script is adapted from Bernie's plotting scripts
"""
import os
import sys
import pdb
import pysam
import logging
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import seaborn as sns
from graphdraw.utils import rgb_picker
from graphdraw.constants import *


logger = logging.getLogger(__name__)
# CHM13 T2T ref lengths


# I have implemented a "better" color picker in utils
# COLORS = ['aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black', 'blanchedalmond',
#           'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue',
#           'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 'darkgreen', 'darkgrey',
#           'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon',
#           'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 'darkviolet', 'deeppink',
#           'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'floralwhite', 'forestgreen', 'fuchsia',
#           'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow', 'grey', 'honeydew', 'hotpink',
#           'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon',
#           'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightgrey',
#           'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray', 'lightslategrey',
#           'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 'magenta', 'maroon', 'mediumaquamarine',
#           'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 'mediumspringgreen',
#           'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 'moccasin', 'navajowhite',
#           'navy', 'oldlace', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 'palegoldenrod', 'palegreen',
#           'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru', 'pink', 'plum', 'powderblue', 'purple',
#           'rebeccapurple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen',
#           'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 'slategray', 'slategrey', 'springgreen',
#           'steelblue', 'tan', 'teal', 'thistle', 'tomato', 'turquoise', 'violet', 'wheat', 'yellow', 'yellowgreen']


def read_tsv(in_tsv):
    """
    Reads in the TSV for generating plots like CNV or depth coverage
    The table needs to be tab separated as has a head with chromosome, start, end and a value
    """
    if not os.path.isfile(in_tsv):
        logger.error(f"File {in_tsv} does not exist")
        sys.exit(1)
    data = dict()
    with open(in_tsv, 'r') as infile:
        line = next(infile)
        if "chromosome" not in line:
            logger.error("The header should contain the following information chromosome, start, end, value")
            sys.exit(1)
        line = line.strip().split('\t')
        chrom_loc = line.index('chromosome')
        start_loc = line.index('start')
        end_loc = line.index('end')
        value_loc = line.index('value')
        for line in infile:
            line = line.strip().split('\t')
            if line[chrom_loc] not in data:
                data[line[chrom_loc]] = [(int(line[start_loc]), int(line[end_loc]), float(line[value_loc]))]
            else:
                data[line[chrom_loc]].append((int(line[start_loc]), int(line[end_loc]), float(line[value_loc])))
    return data


def get_chrom_counts(chrom):
    # chrom here is a dict of contig ids and alignment coordinates
    counts = 0
    for c in chrom:
        counts += len(chrom[c])
    return counts


def plot_contigs_per_chrom(data, prefix, pdf, two_haps=-1):
    chroms = list(data.keys())
    counts = [len(value) for value in data.values()]
    fig, axs = plt.subplots(2)
    fig.set_figheight(10)
    fig.set_figwidth(8)

    if two_haps != -1:
        fig.suptitle(prefix + f" Haplotype {two_haps}")
    # plt.xticks(rotation=45, ha='right')
    plt.subplots_adjust(hspace=0)
    axs[0].set(ylabel='Number of Contigs')
    axs[0].set_title(f"Number of contigs distribution in {prefix}")
    axs[0].tick_params(labelrotation=45)
    g = sns.barplot(x=chroms, y=counts, color=(0.2, 0.4, 0.6, 0.6), ax=axs[0])
    # Save figure
    # out_file = f"{outdir}/{prefix}_contigs_histogram.png"
    # plt.savefig(out_file, dpi = DPI)
    # plt.clf()

    # plotting number of alignment per chromosome
    counts = [get_chrom_counts(chrom) for chrom in data.values()]
    # fig = plt.figure(figsize=(15,5))
    # plt.xticks(rotation=45, ha='right')
    plt.subplots_adjust(hspace=0.5)
    axs[1].set(ylabel="Number of Alignments")
    axs[1].set_title(f"Number of Contig Alignments distribution in {prefix}")
    axs[1].tick_params(labelrotation=45)

    g = sns.barplot(x=chroms, y=counts, color=(0.2, 0.4, 0.6, 0.6), ax=axs[1])

    # out_file = f"{outdir}/{prefix}_histogram.png"
    pdf.savefig(fig, dpi=DPI)


def collect_intervals(data):
    """
    Collect all intervals of different contigs for all chromosomes
    """
    all_data = dict()
    for chr_id, contigs in data.items():
        intervals = set()
        for contig_id, coords in contigs.items():
            for coord in coords:
                intervals.add((contig_id, coord[0], coord[1]))
        intervals = list(intervals)
        intervals.sort(key=lambda x: x[2], reverse=False)
        all_data[chr_id] = intervals
    return all_data


def plot_contig_alignment_intervals(haplotype1, prefix, pdf, haplotype2=None, cnv=None):
    """
    plots the continuous contigs against the chromosome
    all I need is a query name, a start and an end
    and to be split into chromosome name then the hits for that chrom
    """
    if cnv is not None:
        print(f"Getting values and intervals from {cnv}")
        cnv_table = read_tsv(cnv)

    logger.info(f"Getting alignments for {haplotype1}")
    haplotype1_data = collect_intervals(haplotype1)
    if haplotype2 is not None:
        logger.info(f"Getting alignments for {haplotype2}")
        haplotype2_data = collect_intervals(haplotype2)
        if cnv is None:
            chromosomes = set(haplotype1_data.keys()).intersection(haplotype2_data.keys())
        else:
            chromosomes = set(haplotype1_data.keys()).intersection(haplotype2_data.keys()).intersection(set(cnv_table.keys()))
    else:
        chromosomes = set(haplotype1_data.keys())

    for chr_id in CENTRO_T2T_UCSC_IDEOGRAM.keys():
        # as chromosome is a set, it's ordered alphabetically, I am using instead the order in that dictionary
        # taken from constants.py and skipping any chromosome that doesn't exist in my data
        if chr_id not in chromosomes:
            continue

        intervals_haplotype1 = haplotype1_data[chr_id]
        if haplotype2 is not None:
            intervals_haplotype2 = haplotype2_data[chr_id]
        if cnv is not None:
            intervals_cnv = cnv_table[chr_id]

        # Plotting
        if haplotype2 is not None:
            if cnv is None:
                fig, axs = plt.subplots(2)
                fig.set_figheight(11)
                fig.set_figwidth(8)
            else:
                fig, axs = plt.subplots(3)
                fig.set_figheight(11)
                fig.set_figwidth(8)

        else:
            fig = plt.figure(figsize=(10, 6))
            ax = plt.subplot(1, 1, 1)

        # centromere
        if haplotype2 is not None:
            sns.lineplot(x=CENTRO_T2T_UCSC_IDEOGRAM[chr_id], y=[2, 2], linewidth=2, color="black", markers=False,
                         errorbar=None, dashes=True, ax=axs[0])

            sns.lineplot(x=CENTRO_T2T_UCSC_IDEOGRAM[chr_id], y=[2, 2], linewidth=2, color="black", markers=False,
                         errorbar=None, dashes=True, ax=axs[1])

        else:
            sns.lineplot(x=CENTRO_T2T_UCSC_IDEOGRAM[chr_id], y=[2, 2], linewidth=2, color="black", markers=False,
                         errorbar=None, dashes=True)

        y_coord = 3
        color_counter = 0
        seen_contigs = dict()
        previous_color = "#ffffff"

        # plotting lines of haplotype 1, or just one bam file given
        for interval in intervals_haplotype1:
            # should give the same color to the same contig
            # so if a contig is fragmented it will still be recognizable
            if interval[0] in seen_contigs:
                color = seen_contigs[interval[0]]
            else:
                # rolling color counter
                # seen_contigs[interval[0]] = COLORS[color_counter % (len(COLORS) - 1)]
                seen_contigs[interval[0]] = rgb_picker(previous_color, color_counter % 3)
                previous_color = seen_contigs[interval[0]]
                color = seen_contigs[interval[0]]
                color_counter += 1

            current_contig, begin, end = interval
            coords = [begin, end]
            # I can maybe add palette with a color
            if haplotype2 is not None:
                sns.lineplot(x=coords, y=[y_coord, y_coord], color=color, markers=False,
                             errorbar=None, dashes=False, ax=axs[0])
            else:
                sns.lineplot(x=coords, y=[y_coord, y_coord], color=color, markers=False,
                             errorbar=None, dashes=False)
            y_coord += 1

        # in case a second bam file given then I plot that as haplotype 2
        if haplotype2 is not None:
            y_coord = 3
            color_counter = 0
            seen_contigs = dict()
            previous_color = "#ffffff"

            for interval in intervals_haplotype2:
                # should give the same color to the same contig
                # so if a contig is fragmented it will still be recognizable
                if interval[0] in seen_contigs:
                    color = seen_contigs[interval[0]]
                else:
                    # rolling color counter
                    # seen_contigs[interval[0]] = COLORS[color_counter % (len(COLORS) - 1)]
                    seen_contigs[interval[0]] = rgb_picker(previous_color, color_counter % 3)
                    previous_color = seen_contigs[interval[0]]
                    color = seen_contigs[interval[0]]
                    color_counter += 1

                current_contig, begin, end = interval
                coords = [begin, end]
                # I can maybe add palette with a color
                sns.lineplot(x=coords, y=[y_coord, y_coord], color=color, markers=False, errorbar=None,
                             dashes=False, ax=axs[1])
                y_coord += 1

        if cnv is not None:
            for x1, x2, value in intervals_cnv:
                sns.lineplot(x=(x1, x2), y=[value, value], linewidth=2, color="black", markers=False,
                             errorbar=None, dashes=True, ax=axs[2])
        # plt.xlim(0, REF_LEN[chr_id])
        # plt.ylim(0, y_coord)

        if haplotype2 is not None:
            # common xlabel
            fig.text(0.5, 0.04, "Coord (Mb)", ha='center', fontsize=10)
            # axs[0].set_ylabel("Coordinates", fontsize=10)
            # axs[1].set_ylabel("Coordinates", fontsize=10)

            axs[0].spines['right'].set_visible(False)
            axs[0].spines['top'].set_visible(False)

            axs[1].spines['right'].set_visible(False)
            axs[1].spines['top'].set_visible(False)

            axs[0].set_title(f"{prefix} alignment in {chr_id} Haplotype 1")
            axs[1].set_title(f"{prefix} alignment in {chr_id} Haplotype 2")
            if cnv is not None:
                axs[2].spines['right'].set_visible(False)
                axs[2].spines['top'].set_visible(False)
                axs[2].set_ylabel("Copy Number", fontsize=10)
                axs[2].set_title(f"{prefix} Copy Number in {chr_id}")

            # common title
            # fig.suptitle(f"{prefix} alignment in {chr_id}")

        else:
            plt.xlabel('Coordinates', fontsize=12)
            plt.ylabel('Contigs', fontsize=12)
            plt.title(f"{prefix} alignment in {chr_id}")

            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
        # common ylable
        # fig.text(0.04, 0.5, 'Contigs', va='center', rotation='vertical')

        # xlabel for one subplot
        # axs[0].set_xlabel('Coord (Mb)', fontsize=10)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        ## Save figure
        # outFile = f"{outdir}/alignment_intervals_{chr_id}.png"
        pdf.savefig(fig, dpi=DPI)


def plot_counts_per_bin(counts, prefix, sizes_dict, out_dir):
    """
    for plotting that continuous counts per bin

    I need to define some bins and count the contigs in those bins
    one way to do it is to just look at the beginning of the alignment
    and count that as alignment in that bin

    seaborn just takes bins (I guess int or char) and counts for that bin
    1M, 2M, 3M and so on I guess
    """
    chr_list = set(counts['chr'].tolist())

    for chrom in chr_list:
        # counts_chr is a dict with 'bin' as list of bin names (location)
        # and 'counts' as a list of counts corresponding to the bins

        counts_chr = counts[counts['chr'] == chrom]
        fig = plt.figure(figsize=(10, 2))
        plt.ylabel('# contigs', fontsize=12)
        ax = fig.add_subplot(111)
        sns.lineplot(data=counts_chr, x='bin', y='counts')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # sizes_dict[chrom] gives the size of that chromosome
        plt.xlim(0, sizes_dict[chrom] / 1000000)
        ##Save figure
        outFile = f"{outdir}/alignment_bins_{chrId}.png"
        plt.savefig(outFile, dpi=DPI)


def prepare_bam(in_bam):
    if not os.path.exists(in_bam):
        logger.error(f"The given input file {in_bam} does not exist!")
        sys.exit(1)

    if in_bam.endswith(".bam"):
        alignment_file = pysam.AlignmentFile(in_bam, "rb")
    else:
        alignment_file = pysam.AlignmentFile(in_bam, "r")

    data = dict()
    for a in alignment_file:
        if a.mapq == 60:
            # create entry for the chromosome if there are none yet
            if a.reference_name not in data:
                data[a.reference_name] = dict()
            if a.qname in data[a.reference_name]:
                data[a.reference_name][a.qname].append((a.reference_start, a.reference_end))
            else:
                data[a.reference_name][a.qname] = [(a.reference_start, a.reference_end)]
        to_remove = set()
        for chrom in data.keys():
            if not data[chrom]:
                to_remove.add(chrom)
        for c in to_remove:
            del data[c]
    return data


def prepare_tsv(in_tsv):
    data = dict()
    with open(in_tsv, "r") as infile:
        header = infile.readline()
        header = header.strip().split()
        START = header.index("ref_start")
        END = header.index("ref_end")
        CHROM = header.index("chr")
        MAPQ = header.index("mapq")
        for l in infile:
            l = l.strip().split()
            mapq = int(l[MAPQ])
            if mapq != 60:  # skipping all non-unique alignments
                continue
            l[START] = int(l[START])  # start
            l[END] = int(l[END])  # end
            contig = l[0]
            if l[CHROM] not in data:
                data[l[CHROM]] = dict()
            else:
                if l[0] not in data[l[CHROM]]:
                    data[l[CHROM]][l[0]] = [(l[START], l[END])]
                else:
                    data[l[CHROM]][l[0]].append((l[START], l[END]))
        to_remove = set()
        for chrom in data.keys():
            if not data[chrom]:
                to_remove.add(chrom)
        for c in to_remove:
            del data[c]

    return data


def plotting(in_bam1, outfile, prefix, in_table=None, in_bam2=None):
    haplotype1 = prepare_bam(in_bam1)
    pdf = matplotlib.backends.backend_pdf.PdfPages(outfile)
    if not in_bam2:
        plot_contigs_per_chrom(haplotype1, prefix, pdf)
        plot_contig_alignment_intervals(haplotype1, prefix, pdf, cnv=in_table)
    else:
        haplotype2 = prepare_bam(in_bam2)
        plot_contigs_per_chrom(haplotype1, prefix, pdf, two_haps=1)
        plot_contigs_per_chrom(haplotype2, prefix, pdf, two_haps=2)
        plot_contig_alignment_intervals(haplotype1, prefix, pdf, haplotype2=haplotype2, cnv=in_table)
    pdf.close()


# if __name__ == "__main__":
#     if len(sys.argv) < 3:
#         print("You need to give input tsv with contigID-chr-start-end and a prefix for plots produced")
#         sys.exit()
#
#     # the bed file in the align folder from PAV pipeline
#     in_table = sys.argv[1]
#     prefix = sys.argv[2]
#
#     out_dir = prefix + "_alignment_plots"
#     # if not os.path.exists(out_dir):
#     #	 os.makedirs(out_dir)
#
#     haplotype1 = prepare_tsv(in_table)
#     # test_data = dict()
#     # test_data['chr17'] = data['chr17']
#     outfile = out_dir + ".pdf"
#     plotting(haplotype1, outfile, prefix)
#     # pdb.set_trace()
#     # pdf = matplotlib.backends.backend_pdf.PdfPages(f"{out_dir}.pdf")
#     # plot_contigs_per_chrom(data, prefix, pdf)
#     # plot_contig_alignment_intervals(data, prefix, pdf)
#
#     # pdf.close()
    # prepare the table into a dict divided by chromosome as keys
