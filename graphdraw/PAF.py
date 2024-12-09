"""
Module 'formats' - Contains classes for dealing with file formats such as fasta, bed, vcf, etc...
Original from Bernardo Rodriguez-martin, and modified by Fawaz Dabbaghie
"""

from graphdraw.gRanges import *
import pysam


class PAF:
    """
    Class for dealing with files in PAF format.
    """

    def __init__(self):
        """
        Initialize empty class instance
        """
        self.alignments = []

    def load_alignments(self, list_of_a, type):
        """
        Loads alignments from SAM, BAM, or PAF
        """
        # For line in the file
        for a in list_of_a:
            if type == "paf":
                self.alignments.append(PafAlignment(a))
            else:
                self.alignments.append(BamAlignment(a))
        # when calling infer_query_length with pysam, I don't always get the complete length (probably from
        # hard clipping), so I'll just take the maximum reported number and use that for all alignments
        max_q_length = max([a.qLen for a in self.alignments])
        for a in self.alignments:
            a.qLen = max_q_length

    def __iter__(self):
        return iter(self.alignments)

    def read(self, file_path):
        """
        PAF file reader. Read and store data line objects into a list:
        """
        paf_file = open(file_path)

        # For line in the file
        for line in paf_file:

            # Skip comments and blank lines
            if line.startswith('#') or not line:
                continue

            fields = line.split()
            line = PafAlignment(fields)
            self.alignments.append(line)

    def sort_by_len(self):
        """
        Sort alignments by query alignment length in descending order
        """
        sorted_alignments = sorted(self.alignments, key=lambda x: x.alignment_len(), reverse=True)
        return sorted_alignments

    def chain(self, max_dist, max_perc_overlap):
        """
        Chain PAF alignments based on alignment complementarity

        Input:
            1. maxDist: maximum distance between both ranges
            2. maxPercOverlap: maximum percentage of overlap between ranges

        Output:
            1. chain: PAF_chain object instance
        """
        # 1. Sort alignments by decreasing query alignment length
        sorted_alignments = self.sort_by_len()

        # 2. Pick the longest alignment and initiate chain
        longest = sorted_alignments[0]
        chain = PafChain([longest])

        # remove alignment
        del sorted_alignments[0]

        round_counter = 1

        # 3. Attemp to extend the chain with complementary alignments
        while True:

            # START ALIGNMENT CHAIN EXTENSION ROUND
            # Initialize boolean as not complementary alignment found
            compl_bool = False

            # Go through all the available alignments
            for index, alignment in enumerate(sorted_alignments):

                # Assess if alignment complementary to the chain
                chain_beg, chain_end = chain.interval()
                compl_bool, orientation = complementary(chain_beg, chain_end, alignment.qBeg, alignment.qEnd, max_dist,
                                                        max_perc_overlap)

                # Complementary alignment found
                if compl_bool:

                    # Add alignment to the chain
                    # a) Add to the chain begin
                    if orientation == "LEFT":
                        chain.alignments.insert(0, alignment)
                    # b) Add to the chain end
                    else:
                        chain.alignments.append(alignment)
                    # Remove from list
                    del sorted_alignments[index]

                    # Stop once complementary found
                    break

            round_counter += 1

            # STOP CHAIN EXTENSION IF:
            # a) No complementary alignment found in the last round OR
            # b) Mo alignments left
            if compl_bool is False or not sorted_alignments:
                break

        return chain

    def hits2dict(self):
        """
        Reorganize hits into a dictionary
        """

        hits_dict = {}
        # For each hit
        for hit in self.alignments:

            # Initialize list
            if hit.qName not in hits_dict:
                hits_dict[hit.qName] = []

            # Add hit to list
            hits_dict[hit.qName].append(hit)

        return hits_dict


class BamAlignment:
    """
    BAM entry class that PAF class can use as well
    """

    number = 0

    def __init__(self, pysam_alignment=None):
        BamAlignment.number += 1
        self.id = "BAM_alignment_" + str(BamAlignment.number)
        if not pysam_alignment:
            pass
        else:
            self.qName = pysam_alignment.query_name
            self.qLen = pysam_alignment.infer_query_length()
            self.qBeg = pysam_alignment.query_alignment_start
            self.qEnd = pysam_alignment.query_alignment_end
            if pysam_alignment.is_reverse:
                self.strand = "-"
            else:
                self.strand = "+"
            self.tName = pysam_alignment.reference_name
            self.tLen = pysam_alignment.reference_length
            self.tBeg = pysam_alignment.reference_start
            self.tEnd = pysam_alignment.reference_end
            self.MAPQ = pysam_alignment.mapq

    def alignment_len(self):
        """
        Compute the query alignment length
        """
        return self.qEnd - self.qBeg

    def alignment_perc(self):
        """
        Compute the query alignment length percentage
        """
        perc_len = float(self.alignment_len()) / self.qLen * 100
        return perc_len


class PafAlignment:
    """
    PAF entry class
    """
    number = 0  # Number of instances

    def __init__(self, fields=None):
        """
        Initialize paf line
        """
        PafAlignment.number += 1  # Update instances counter
        self.id = 'PAF_alignment_' + str(PafAlignment.number)
        if not fields:  # now can initiate an empty alignment
            pass
        else:
            self.qName = str(fields[0])
            self.qLen = int(fields[1])
            self.qBeg = int(fields[2])
            self.qEnd = int(fields[3])
            self.strand = str(fields[4])
            self.tName = str(fields[5])
            self.tLen = int(fields[6])
            self.tBeg = int(fields[7])
            self.tEnd = int(fields[8])
            self.MAPQ = int(fields[11])

    def alignment_len(self):
        """
        Compute the query alignment length
        """
        return self.qEnd - self.qBeg

    def alignment_perc(self):
        """
        Compute the query alignment length percentage
        """
        perc_len = float(self.alignment_len()) / self.qLen * 100
        return perc_len


class PafChain:
    """
    Chain of complementary PAF alignments
    """

    def __init__(self, alignments):
        """
        Initialize chain instance.

        Input:
            1. alignments. List of PAF_alignment instances
        """
        self.alignments = alignments

    def interval(self):
        """
        Return query interval covered by the chain
        """
        first_alignment = self.alignments[0]
        last_alignment = self.alignments[-1]

        return first_alignment.qBeg, last_alignment.qEnd

    def interval_template(self):
        """
        Return query interval covered by the chain
        """
        first_alignment = self.alignments[0]
        last_alignment = self.alignments[-1]

        return first_alignment.tBeg, last_alignment.tEnd

    def perc_query_covered(self):
        """
        Compute the percentage of the query sequence covered by the chain of alignments
        """
        # a) No alignments available
        if len(self.alignments) == 0:
            perc_covered = 0

        # b) Alignments available
        else:
            # Compute the number of bases covered
            beg, end = self.interval()
            alignment_len = end - beg

            # Compute the percentage of bases covered
            perc_covered = float(alignment_len) / self.alignments[0].qLen * 100

        return perc_covered
