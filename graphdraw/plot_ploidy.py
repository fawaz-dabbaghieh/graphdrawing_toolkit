import os
import sys
import pdb
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

CHROMOSOMES = ['chr1',
               'chr2',
               'chr3',
               'chr4',
               'chr5',
               'chr6',
               'chr7',
               'chr8',
               'chr9',
               'chr10',
               'chr11',
               'chr12',
               'chr13',
               'chr14',
               'chr15',
               'chr16',
               'chr17',
               'chr18',
               'chr19',
               'chr20',
               'chr21',
               'chr22',
               'chrY',
               'chrX',
               'chrM']


def set_xlabels(x_coordinates, n_xlabel):
    """
    I take all the x coordinates and subset some n_xlabel of them to visualize
    """
    step = int(len(x_coordinates) / n_xlabel)
    labels = x_coordinates[0::step]
    if labels[-1] != x_coordinates[-1]:
        labels[-1] = x_coordinates[-1]
    return labels


def plot_ploidy(ploidy, chrom_name, title_info, max_cnv):
    # ploidy is a list of tuples with (ploidy, start, end)
    # calculating average for each bin

    fig = plt.figure(figsize=(16, 8))
    #	plt.hist(chrom_count, bin_size, histtype="step")
    # plt.plot([x[1] for x in ploidy], [x[0] for x in ploidy], linestyle="-", color='b', label='Alignment Depth')
    x = []
    y = []
    for p in ploidy:
        x.append(p[1])
        x.append(p[2])
        y.append(p[0])
        y.append(p[0])
    # cnv = [x[0] for x in ploidy]
    # start = [x[1] for x in ploidy]
    # cnv.append(ploidy[-1][0])
    # start.append(ploidy[-1][2])
    plt.plot(x, y, linestyle='-', color='black', label='CNV', linewidth=2)
    # for cnv, x1, x2 in ploidy:
    # plt.plot([x1, x2], [cnv, cnv], linestyle='-', color='b', label='CNV')

    plt.plot([ploidy[0][1], ploidy[-1][2]], [2, 2], linestyle='--', color='r', linewidth=1)
    plt.xlabel(f"Chromsome {chrom_name}")

    plt.yticks(list(range(max_cnv)), ha="center", fontsize=8)

    plt.ylabel("Ploidy")
    # n_steps = 30
    # if len(x) > n_steps:
    # 	plt.xticks(set_xlabels(x, n_steps), rotation=45, ha="center", fontsize=8)

    plt.title(f"Ploidy, {title_info}")
    fig.tight_layout()
    print(f"Info: finished saving plot for chromsome {chrom_name}")
    return fig


if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("You need to give the input <ploidy_counts> <output_plot> <title_info> "
              "<chrom_col_n> <start_col_n> <end_col_n> <cnv_col_n>")
        sys.exit()

    in_table = sys.argv[1]
    if not os.path.exists(in_table):
        print(f"the file {in_table} does not exist")
        sys.exit()

    out_plot = sys.argv[2]
    title_info = sys.argv[3]
    chrom_loc = int(sys.argv[4])
    chrom_start = int(sys.argv[5])
    chrom_end = int(sys.argv[6])
    cnv_value = int(sys.argv[7])

    ploidy = dict()

    out_file = out_plot.split(".")[0] + ".pdf"
    pdf_file = PdfPages(out_file)
    max_cnv = 0
    with open(in_table, "r") as infile:
        for l in infile:
            if not l.startswith("chr"):
                continue
            l = l.strip().split()
            if max_cnv < int(l[cnv_value]):
                max_cnv = int(l[cnv_value])
            if l[0] not in ploidy:
                ploidy[l[chrom_loc]] = [
                    (int(l[cnv_value]), int(l[chrom_start]), int(l[chrom_end]))]  # ploidy, start, end
            else:
                ploidy[l[chrom_loc]].append((int(l[cnv_value]), int(l[chrom_start]), int(l[chrom_end])))
        # basically filling 0's where there's no counts for that position

    for chrom in CHROMOSOMES:
        try:
            value = ploidy[chrom]
        except KeyError:
            print(f"chromosome {chrom} does not exist in the data table")
            continue
        fig = plot_ploidy(value, chrom, title_info, max_cnv)
        pdf_file.savefig(fig)

    pdf_file.close()
