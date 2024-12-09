import os
import sys
import pdb
import math
from graphdraw.constants import REF_LEN
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D


CHROMOSOMES = list(REF_LEN.keys())


def roundup(x):
    return int(math.ceil(x / 100.0)) * 100


def set_xlabels(x_coordinates, n_xlabel):
    """
    I take all the x coordinates and subset some n_xlabel of them to visualize
    """
    step = int(len(x_coordinates) / n_xlabel)
    labels = x_coordinates[0::step]
    if labels[-1] != x_coordinates[-1]:
        labels[-1] = x_coordinates[-1]
    return labels


def get_centromer_coordinates(reference, chromosome):
    if reference == "t2t":
        from graphdraw.constants import CENTRO_T2T_UCSC_IDEOGRAM as centromers
    elif reference == "hg19":
        from graphdraw.constants import CENTRO_HG19_UCSC_IDEOGRAM as centromers
    elif reference == "hg38":
        from graphdraw.constants import CENTRO_HG38_UCSC_IDEOGRAM as centromers
    else:
        print(
            f"Warning: The reference {reference} given was not recognized, centromeric locations will not be added to the plot")
        return (0, 0)

    if chromosome in centromers:
        return centromers[chromosome]
    else:
        print(f"Warning: The chromosome {chromosome} was not found in the reference list of centromers")
        return (0, 0)


def get_gaps_coordinates(reference, chromosome):
    if reference == "hg19":
        from graphdraw.constants import HG19_GAPS as gaps
    elif reference == "hg38":
        from graphdraw.constants import HG38_GAPS as gaps
    else:
        print(
            f"Warning: The reference {reference} given was not recognized, gap locations will not be added to the plot")
        return [(0, 0)]

    if chromosome in gaps:
        return gaps[chromosome]
    else:
        print(f"Warning: The chromosome {chromosome} was not found in the reference list of gaps")
        return [(0, 0)]


def plot_chromosome(chrom_count, chrom_name, reference, bin_perc, title_info):
    # calculating average for each bin
    bin_size = int(len(chrom_count) * bin_perc)
    if bin_size == 0:
        n_bins = 1
    else:
        n_bins = int(len(chrom_count) / bin_size)

    print(f"Length of chromosome {chrom_name} is {len(chrom_count)}, and with bin percentage "
          f"of {bin_perc * 100}% we'll have a bin size of {bin_size} and {n_bins} bins")

    y = []
    x = []
    start = 0
    mid_point = int(bin_size / 2)
    for i in range(n_bins):
        if not i == n_bins - 1:
            y.append(sum(chrom_count[start:start + bin_size]) / bin_size)
            x.append(start + mid_point)
        else:
            if start > len(chrom_count):
                continue
            y.append(sum(chrom_count[start:]) / len(chrom_count[start:]))
            x.append(start + mid_point)
        start = start + bin_size + 1
    fig = plt.figure(figsize=(16, 8))
    # plt.hist(chrom_count, bin_size, histtype="step")
    plt.plot(x, y, linestyle="-", color='b', label='Alignment Depth')

    # adding centromere
    x1, x2 = get_centromer_coordinates(reference, chrom_name)
    if (x1, x2) == (0, 0):
        centro_exists = False
    else:
        centro_exists = True

    if x1 < x[-1]:
        if x2 > x[-1]:
            x2 = x[-1]
        plt.plot([x1, x2], [-2, -2], linestyle='-', color='r', label='Centromere')
    # adding gaps
    gaps = get_gaps_coordinates(reference, chrom_name)
    if gaps == [(0, 0)]:
        gap_exists = False
    else:
        gap_exists = True

    for x1, x2 in gaps:
        if x1 < x[-1]:
            if x2 > x[-1]:
                x2 = x[-1]
            plt.plot([x1, x2], [-1, -1], linestyle='-', color='g', label='Gap')

    # plot information (title, legend, lables and so on)
    if centro_exists is True and gap_exists is True:
        custom_legend = [Line2D([0], [0], color='b', lw=1.5), Line2D([0], [0], color='r', lw=1.5),
                         Line2D([0], [0], color='g', lw=1.5)]
        plt.legend(custom_legend, ["Alignment Depth", "Centromere", "Gaps"])
    elif centro_exists is False and gap_exists is True:
        custom_legend = [Line2D([0], [0], color='b', lw=1.5), Line2D([0], [0], color='g', lw=1.5)]
        plt.legend(custom_legend, ["Alignment Depth", "Gaps"])
    elif centro_exists is True and gap_exists is False:
        custom_legend = [Line2D([0], [0], color='b', lw=1.5), Line2D([0], [0], color='r', lw=1.5)]
        plt.legend(custom_legend, ["Alignment Depth", "Centromere"])
    elif centro_exists is False and gap_exists is False:
        custom_legend = [Line2D([0], [0], color='b', lw=1.5)]
        plt.legend(custom_legend, ["Alignment Depth"])

    plt.xlabel(f"Chromsome {chrom_name}")
    plt.ylabel("Alignment counts")
    n_steps = 30
    if len(x) > n_steps:
        plt.xticks(set_xlabels(x, n_steps), rotation=45, ha="center", fontsize=8)

    plt.title(f"Alignment counts, {title_info}")
    fig.tight_layout()
    print(f"Info: finished saving plot for chromsome {chrom_name}")
    return fig


def output_windowed_table(counts, chrom, bin_perc, out_table):
    """
    Takes in the samtools depth table, the size of bin percentage and an output file to produce a table with
    header chromosome start end value which can be given to the contig plotting subcommand for plotting
    """

    bin_size = int(len(counts) * bin_perc)
    if bin_size == 0:
        n_bins = 1
    else:
        n_bins = int(len(counts) / bin_size)

    with open(out_table, "a") as outfile:
        start = 0
        for i in range(n_bins):
            if not i == n_bins - 1:
                end = start + bin_size
                value = sum(counts[start:start + bin_size]) / bin_size
            else:
                if start > len(counts):
                    continue
                end = len(counts) - 1
                value = sum(counts[start:]) / len(counts[start:])
            outfile.write(f"{chrom}\t{start}\t{end}\t{value}\n")
            start = start + bin_size + 1


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("You need to give the input <samtools_depth_table> <output_plot> <size_of_bin_perc> "
              "<title_info> <reference_name>")
        sys.exit()

    in_table = sys.argv[1]
    if not os.path.exists(in_table):
        print(f"the file {in_table} does not exist")
        sys.exit()

    out_plot = sys.argv[2]
    try:
        bin_perc = float(sys.argv[3])
        bin_perc = bin_perc / 100
    except ValueError:
        print("The percentage needs to be an integer")
        sys.exit()

    title_info = sys.argv[4]

    reference = sys.argv[5]
    counts = dict()

    out_file = out_plot.split(".")[0] + ".pdf"
    pdf_file = PdfPages(out_file)

    counter = 0
    with open(in_table, "r") as infile:
        for l in infile:
            l = l.strip().split()
            counter += 1
            if counter % 10000000 == 0:
                print(f"processed 10 million lines for chromosome {l[0]}")
            # basically filling 0's where there's no counts for that position
            try:
                while len(counts[l[0]]) < int(l[1]) - 1:
                    counts[l[0]].append(0)
                counts[l[0]].append(int(l[2]))

            except KeyError:
                if counts:  # already finished a chromosome, plot it and empty memory
                    key = list(counts.keys())[0]
                    # chromosome with no counts
                    if len(counts[key]) == 0:
                        continue
                    # number of counts is too small
                    if int(len(counts[key]) * bin_perc) == 0:
                        continue
                    if key not in CHROMOSOMES:
                        continue
                    pdb.set_trace()
                    fig = plot_chromosome(counts[key], key, reference, bin_perc, title_info)
                    # import pickle
                    # with open("chr1_counts.pickle", "wb") as outfile:
                    # pickle.dump(counts[key], outfile)
                    pdf_file.savefig(fig)
                    counts = dict()
                    counts[l[0]] = [0] * int(l[1])
                    counts[l[0]].append(int(l[2]))
                # pdf_file.close()
                # sys.exit()
                else:  # first chromosome
                    counts[l[0]] = [0] * int(l[1])
                    counts[l[0]].append(int(l[2]))
    if counts:
        key = list(counts.keys())[0]
        if key in CHROMOSOMES:
            fig = plot_chromosome(counts[key], key, reference, bin_perc, title_info)
            pdf_file.savefig(fig)

    pdf_file.close()
