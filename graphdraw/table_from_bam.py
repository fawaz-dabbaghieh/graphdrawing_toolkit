import os
import sys
import pysam
import pickle
from graphdraw.utils import read_fasta_gen


def get_bam_info(in_bam):
    if not os.path.exists(in_bam):
        logging.error(f"The given input file {in_bam} does not exist!")
        sys.exit(1)

    if in_bam.endswith(".bam"):
        alignment_file = pysam.AlignmentFile(in_bam, "rb")
    else:
        alignment_file = pysam.AlignmentFile(in_bam, "r")
    # Data here is a dictionary with contig id as key and value is a dictionary of chromosomes
    # where that contig aligned. If a contig aligns to only one chromosome, then that dict is of length 1
    contig_info = dict()
    for a in alignment_file:
        if a.qname not in contig_info:
            contig_info[a.qname] = dict()
            contig_info[a.qname][a.reference_name] = [(a.reference_start, a.reference_end, a.mapq)]
        else:
            if a.reference_name in contig_info[a.qname]:
                contig_info[a.qname][a.reference_name].append((a.reference_start, a.reference_end, a.mapq))
            else:
                contig_info[a.qname][a.reference_name] = [(a.reference_start, a.reference_end, a.mapq)]
    return contig_info


def get_contig_lengths(in_fasta):
    """
    Get the contig lengths from a fasta file
    """
    return {key: len(value) for (key, value) in read_fasta_gen(in_fasta)}


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} bam_file contigs_fasta output_table")
        sys.exit(1)

    data = get_bam_info(sys.argv[1])
    contig_length = get_contig_lengths(sys.argv[2])
    with open(sys.argv[3], "w") as outfile:
        outfile.write("\t".join(["contig_name", "contig_length", "chromosome", "start", "end", "map_quality"]) + "\n")
        for ctg in data.keys():
            for chrom in data[ctg].keys():
                for alignment in data[ctg][chrom]:
                    outfile.write("\t".join([ctg, str(contig_length[ctg]), chrom,
                                             str(alignment[0]), str(alignment[1]), str(alignment[2])]) + "\n")
