import os
import sys
import pysam
import pdb
from graphdraw.color_graph_with_chromosomes import Alignment
from graphdraw.draw_subgraph_to_ref_svg import *
from graphdraw.SVG import *
from graphdraw.PAF import *
from graphdraw.gRanges import *


def find_ctg_in_bam(in_bam, ctg_name):
    """
    return a list of alignments as a list of pysam alignment objects
    if paf is true, then it returns a populated or empty PAF object
    """
    alignment_file = None
    try:
        if in_bam.endswith(".bam"):
            alignment_file = pysam.AlignmentFile(in_bam, "rb")
        else:
            alignment_file = pysam.AlignmentFile(in_bam, "r")
    except FileNotFoundError:
        logging.error(f"The file {in_bam} does not exist")
        sys.exit(1)
    except ValueError:
        logging.error(f"The file {in_bam} is not a bam file")
        sys.exit(1)

    contig_alignments = []
    if alignment_file:
        for a in alignment_file:
            if a.qname == ctg_name:
                contig_alignments.append(Alignment(a))

    return contig_alignments


def find_ctg_in_paf(in_paf, ctg_name):
    """
    extracts the alignment of ctg_name in the paf file, generates a PAF object and return it populated or empty
    """
    alignments = []
    with open(in_paf, "r") as infile:
        for l in infile:
            if l.startswith(ctg_name):
                alignments.append(l)
    paf = PAF()
    if alignments:
        paf.load_alignments(alignments, type='paf')
        return paf
    return paf


def get_alignment_obj(in_alignment, ctg_name):
    """
    extracts the alignment of ctg_name in the paf file, generates a PAF object and return it populated or empty
    """
    paf = PAF()
    alignments = []
    if in_alignment.endswith("paf"):
        type = "paf"
        with open(in_alignment, "r") as infile:
            for l in infile:
                if l.startswith("#") or not l:
                    continue
                if l.startswith(ctg_name):
                    alignments.append(l.strip().split())

    else:
        type = "bam"
        if in_alignment.endswith(".bam"):
            alignment_file = pysam.AlignmentFile(in_alignment, "rb")
        else:
            alignment_file = pysam.AlignmentFile(in_alignment, "r")
        for a in alignment_file:
            if a.query_name == ctg_name:
                alignments.append(a)

    if alignments:
        paf.load_alignments(alignments, type)
        return paf

    return paf


def filter_chrom_paf(paf):
    """
    splits the paf object to several paf objects for each chromosome
    """
    chromosomes = dict()
    for a in paf.alignments:
        if a.tName not in chromosomes:
            new_paf = PAF()
            new_paf.alignments.append(a)
            chromosomes[a.tName] = new_paf
        else:
            chromosomes[a.tName].alignments.append(a)
    return chromosomes


def filter_chrom_bam(alignments):
    """
    split the bam/sam alignments list into different lists based on chromosome
    """
    chromosomes = dict()
    for a in alignments:
        if a.tName not in chromosomes:
            chromosomes[a.tName] = [a]
        else:
            chromosomes[a.tName].append(a)
    return chromosomes


def remove_redundant(alignments):
    """
    removes alignments with the same start and end on the ref
    keeps the one with the higher mapq
    """
    to_remove = set()
    for a in alignments:
        for aa in alignments:
            if a != aa:
                if (a.ref_start, a.ref_end) == (aa.ref_start, aa.ref_end):
                    if a.mapq < aa.mapq:
                        to_remove.add(a)
                    elif a.mapq > aa.mapq:
                        to_remove.add(a)
                    else:
                        # have equal mapping quality choose randomly
                        to_remove.add(a)
    for a in to_remove:
        alignments.remove(a)


def remove_inside(alignments):
    """
    remove the alignments that are completely inside other ones
    """
    all_intervals = []
    for a in alignments:
        all_intervals.append((a.ref_start, a.ref_end, a))
    to_remove = set()
    for s1, e1, a1 in all_intervals:
        for s2, e2, a2 in all_intervals:
            if (s1, e1) != (s2, e1):
                # if interval start1-end1 is inside interval start2-end2, remove
                if (s2 < s1 < e2) and (s2 < e1 < e2):
                    to_remove.add(a1)
    for a in to_remove:
        alignments.remove(a)


def keep_one_chrom(alignments):
    """
    take a majority vote on the chromosomes and keep the most represented one
    """
    chromosomes = dict()
    for a in alignments:
        if a.ref_name in chromosomes:
            chromosomes[a.ref_name] += 1
        else:
            chromosomes[a.ref_name] = 1
    majority_vote = ("", 0)
    for chrom, count in chromosomes.items():
        if majority_vote[1] < count:
            majority_vote = (chrom, count)
    to_remove = []
    for a in alignments:
        if a.ref_name != majority_vote[0]:
            to_remove.append(a)
    for a in to_remove:
        alignments.remove(a)

    return majority_vote[0]


def remove_outliers(alignments):
    """
    need to remove outliers, where it would totally distort the figure
    Not sure what the best way to do this, but for now, I'll check if the interval
    between one node and the next is over 40% of the reference, I remove the first node
    """
    # alignments are sorted by start positions here
    reference_length = alignments[-1].ref_end - alignments[0].ref_start

    to_remove = set()
    for i in range(1, len(alignments)):
        interval = alignments[i].ref_start - alignments[i - 1].ref_start
        if interval <= 0:
            continue
        if (interval / reference_length) > 0.4:
            to_remove.add(alignments[i - 1])
    for a in to_remove:
        alignments.remove(a)


# def alignments_to_nodes(alignments, contig_ref=False):
#     nodes = dict()
#     for idx, a in enumerate(alignments):
#         if not contig_ref:
#             nodes[f"{a.q_name}_{idx}"] = {"start": a.tBeg, "end": a.tEnd, "orientation": a.strand,
#                                           "chromosome": a.tName, "mapq": a.MAPQ}
#         else:
#             nodes[f"{a.q_name}_{idx}"] = {"start": a.qBeg, "end": a.qEnd, "orientation": a.strand,
#                                           "chromosome": a.tName, "mapq": a.MAPQ}
#     return nodes


def paf_chain_to_nodes(alignments, contig_ref=False):
    nodes = dict()
    for idx, a in enumerate(alignments):
        if not contig_ref:
            nodes[f"{a.qName}_{idx}"] = {"start": a.tBeg, "end": a.tEnd, "orientation": a.strand, "chromosome": a.tName,
                                         "mapq": a.MAPQ}
        else:
            nodes[f"{a.qName}_{idx}"] = {"start": a.qBeg, "end": a.qEnd, "orientation": a.strand, "chromosome": a.tName,
                                         "mapq": a.MAPQ}
    return nodes


def paf_chain_edges(nodes, output_svg, show_overlap=False):
    """
    Adds edges to the output svg object, the edges are just based on ordering the nodes after each other
    if show_overlap is True, then the overlap section will be darkened
    """
    ordered_nodes = sorted(list(nodes.values()), key=lambda x: x['start'])
    style = dict()
    style["stroke"] = "Gray"
    style["stroke-width"] = 0.2
    style["stroke-opacity"] = 1

    # each n connects to n+1 (named n1)
    # we have 4 combinations
    # n is + and n1 is + (connect end of n to start of n1)
    # n is + and n1 is - (connect end of n to end of n1)
    # n is - and n1 is + (connect start of n to start of n1)
    # n is - and n1 is - (connect start of n to end of n1)

    for idx, n in enumerate(ordered_nodes):
        if idx + 1 < len(ordered_nodes):
            n1 = ordered_nodes[idx + 1]
            n_len = n['x2'] - n['x1']
            n1_len = n1['x2'] - n1['x1']

            if n['orientation'] == "+" and n1['orientation'] == "+":
                output_svg.add_curve(n['x2'], n['y2'], n1['x1'], n1['y1'], n_len, n1_len, 1, 0, style)
            elif n['orientation'] == "+" and n1['orientation'] == "-":
                output_svg.add_curve(n['x2'], n['y2'], n1['x2'], n1['y2'], n_len, n1_len, 1, 1, style)
            elif n['orientation'] == "-" and n1['orientation'] == "+":
                output_svg.add_curve(n['x1'], n['y1'], n1['x1'], n1['y1'], n_len, n1_len, 0, 0, style)
            else:
                output_svg.add_curve(n['x1'], n['y1'], n1['x2'], n1['y2'], n_len, n1_len, 0, 1, style)


def draw_chained_contigs(in_paf, ctg_name, out_svg, max_dist, max_overlap, svg_width, svg_height, svg_legend):
    """
    produces the final SVG
	"""
    paf = get_alignment_obj(in_paf, ctg_name)
    if len(paf.alignments) == 0:
        logging.info(f"{ctg_name} was not found in input file {in_paf}")
        sys.exit()

    # split the alignments on chromosomes, if contig aligned to different chromosomes
    chromosomes = filter_chrom_paf(paf)

    # drawing alignments against the reference
    for chrom_name, chrom_paf in chromosomes.items():
        paf_chain = chrom_paf.chain(max_dist, max_overlap)
        nodes = paf_chain_to_nodes(paf_chain.alignments, contig_ref=False)
        reference = get_ref_info(nodes, chrom_name, svg_legend)
        output_svg = SVG(height=svg_height, width=svg_width, legend=svg_legend)
        reference['end_px'] = output_svg.width - svg_legend  # offset on the end by 10 pixels
        output_svg.add_ref(reference)
        ordered_nodes = order_nodes(nodes)

        calculate_nodes_coord(nodes, ordered_nodes, reference, output_svg)
        paf_chain_edges(nodes, output_svg)
        output_svg.out_svg(out_svg + "_" + chrom_name + ".svg")

    # drawing alignments against the contig itself as the reference
    paf_chain = paf.chain(max_dist, max_overlap)
    nodes = paf_chain_to_nodes(paf_chain.alignments, contig_ref=True)
    reference = get_ref_info(nodes, ctg_name, svg_legend)

    output_svg = SVG(height=svg_height, width=svg_width, legend=svg_legend)
    reference['end_px'] = output_svg.width - svg_legend
    output_svg.add_ref(reference)
    ordered_nodes = order_nodes(nodes)

    calculate_nodes_coord(nodes, ordered_nodes, reference, output_svg)
    covered = paf_chain.perc_query_covered()
    output_svg.add_txt(svg_legend, svg_height - svg_legend, f"{str(round(covered, 2))}% covered", color="Black",
                       size=14, rotation=0, anchor="start")
    paf_chain_edges(nodes, output_svg)
    output_svg.out_svg(out_svg + "_" + ctg_name + ".svg")


def draw_unchained_contigs(in_bam, ctg_name, out_svg, svg_width, svg_height, svg_legend):
    alignments = get_alignment_obj(in_bam, ctg_name)
    chromosomes = filter_chrom_bam(alignments)
    for chrom_name, alignments in chromosomes.items():
        alignments = sorted(alignments, key=lambda x: x.tBeg)

        nodes = paf_chain_to_nodes(alignments)
        reference = get_ref_info(nodes, chrom_name, svg_legend)

        output_svg = SVG(height=svg_height, width=svg_width, legend=svg_legend)
        reference['end_px'] = output_svg.width - svg_legend  # offset on the end by 10 pixels
        output_svg.add_ref(reference)
        ordered_nodes = order_nodes(nodes)

        calculate_nodes_coord(nodes, ordered_nodes, reference, output_svg)

        output_svg.out_svg(out_svg + "_" + chrom_name + ".svg")
