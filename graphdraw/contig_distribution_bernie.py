##DEPENDENCIES ##
# External
import os
import sys
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statistics
import math

# Internal
import formats
import alignment 
import log
import unix 
import structures
import gRanges


###############
## Functions ##
###############
def read_chr_sizes(chrSizes):
    '''
    '''
    sizes_dict = {}

    ##Open file
    with open(chrSizes) as file:

        ##Read lines and format 
        lines = file.readlines()
        lines = [line.rstrip('\n').split('\t') for line in lines]

        ##Store into dictionary
        for line in lines:
            sizes_dict[line[0]] = int(line[1])
    
    return sizes_dict


def organize_hits_by_chrom(PAF_path):
    '''
    '''
    ## 1. Read PAF file
    PAF = formats.PAF()
    PAF.read(PAF_path)

    ## 2. Organize hits by query name into the dictionary
    hits_dict = {}

    for alignment in PAF.alignments:

        if alignment.tName not in hits_dict:
            hits_dict[alignment.tName] = []

        hits_dict[alignment.tName].append(alignment)

    return hits_dict
        

def nbHits_per_chr(hits_dict, outDir):
    '''
    '''
    chroms = []
    counts = []

    for chrId, hits in hits_dict.items():
        nbHits = len(hits) 

        chroms.append(chrId)
        counts.append(nbHits)

    ##Do bar plot
    fig = plt.figure(figsize=(15,5))
    plt.ylabel('# Hits', fontsize=12)    
    g = sns.barplot(x=chroms, y=counts, color=(0.2, 0.4, 0.6, 0.6))
    
    ##Save figure
    filePath = outDir + '/hits_per_chrom.pdf'
    plt.savefig(filePath)


def nbContigs_per_chr(hits_dict, outDir):
    '''
    '''
    chroms = []
    counts = []

    for chrId, hits in hits_dict.items():
        nbContigs = len(set([hit.qName for hit in hits]))

        chroms.append(chrId)
        counts.append(nbContigs)

    ##Do bar plot
    fig = plt.figure(figsize=(15,5))
    plt.ylabel('# Contigs', fontsize=12)    
    g = sns.barplot(x=chroms, y=counts, color=(0.2, 0.4, 0.6, 0.6))
    
    ##Save figure
    filePath = outDir + '/contigs_per_chrom.pdf'
    plt.savefig(filePath)


def hits2binDb(sizes_dict, hits_dict):
    '''
    '''

    ##Assign tBeg and tEnd as beg and end, create nested dictionary
    nested_dict = {}

    for chrId, hits in hits_dict.items():

        # Initialize 
        if chrId not in nested_dict:
            nested_dict[chrId] = {}
            nested_dict[chrId]['HIT'] = []

        # Set beg and end, add to the list
        for hit in hits:
            hit.beg = int(hit.tBeg)
            hit.end = int(hit.tEnd)

            nested_dict[chrId]['HIT'].append(hit)

    ##Create bin database
    hitsBinDb = structures.create_bin_database(sizes_dict, nested_dict)

    return hitsBinDb


def hits_per_bin(sizes_dict, hitsBinDb):
    '''
    '''
    binSize = 1000000 #1Mb
    counts = [['chr', 'bin', 'counts']]

    for chrom in sizes_dict.keys():

        ## Binning the reference
        targetRefs = [chrom]
        bins = gRanges.makeGenomicBins(sizes_dict, binSize, targetRefs)

        for index, interval in enumerate(bins):
            chrom, beg, end = interval  
        
            if chrom in hitsBinDb:
                hits = hitsBinDb[chrom].collect_interval(beg, end, ['HIT'])
                nbHits = len(hits)
            else:
                nbHits = 0
                
            counts.append([chrom, index, nbHits])

    ##Convert nested list into dataframe
    countsDf = pd.DataFrame(counts[1:],columns=counts[0])

    return countsDf


def contigs_per_bin(sizes_dict, hitsBinDb):
    '''
    '''
    binSize = 1000000 #1Mb
    counts = [['chr', 'bin', 'counts', 'contigIds']]

    for chrom in sizes_dict.keys():

        ## Binning the reference
        targetRefs = [chrom]
        bins = gRanges.makeGenomicBins(sizes_dict, binSize, targetRefs)

        for index, interval in enumerate(bins):
            chrom, beg, end = interval  
        
            if chrom in hitsBinDb:
                hits = hitsBinDb[chrom].collect_interval(beg, end, ['HIT'])
                nbContigs = len(set([hit[0].qName for hit in hits]))
                contigIds = ','.join(set([hit[0].qName for hit in hits]))

                if nbContigs == 0:
                    contigIds = 'NA'

            counts.append([chrom, index, nbContigs, contigIds])

    ##Convert nested list into dataframe
    countsDf = pd.DataFrame(counts[1:],columns=counts[0])
    return countsDf


def plot_counts_per_bin(counts, name, sizes_dict, outDir):
    '''
    '''
    chr_list = set(counts['chr'].tolist())

    for chr in chr_list:
        counts_chr = counts[counts['chr'] == chr]

        fig = plt.figure(figsize=(10,2))
        plt.ylabel('# ' + name, fontsize=12)
        ax = fig.add_subplot(111)
        sns.lineplot(data=counts_chr, x='bin', y='counts')
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.xlim(0, sizes_dict[chr]/1000000)
        
        ##Save figure
        filePath = outDir + '/' + name + '_counts_' + chr + '.pdf'
        plt.savefig(filePath)

def plot_contig_alignment_intervals(hits_dict, outDir):
    '''
    '''
    for chrId, hits in hits_dict.items():

        ##Resort hits by start genomic position
        hits.sort(key=lambda hit: hit.tBeg, reverse=False)

        ##Collect all the intervals:
        intervals = []

        for hit in hits:    
            intervals.append(hit.qName + '_' + str(hit.tBeg) + '_' + str(hit.tEnd))

        ##Remove redundancies
        intervals = list(set(intervals))

        ##Convert into list again
        intervals = [i.split('_') for i in intervals]
        intervals = [[interval[0], int(interval[1]), int(interval[2])]for interval in intervals]

        ##Resort
        intervals.sort(key=lambda x: x[1], reverse=False)  

        ##Plotting
        fig = plt.figure(figsize=(10,2))
        ax = plt.subplot(1, 1, 1)

        counter = 1
        for interval in intervals:
            contigId, tBeg, tEnd = interval
            coords = [tBeg, tEnd]

            sns.lineplot(x=coords, y=[counter, counter], markers=False, ci=None, dashes=False)
            counter += 1

        plt.xlim(0, hits[0].tLen)
        plt.ylim(0, counter)

        plt.xlabel('Coord(Mb)', fontsize=12)
        plt.ylabel('Contigs', fontsize=12)
        
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ## Save figure
        outFile = outDir + '/contig_alignments_' + chrId + '.pdf'
        fig.savefig(outFile)

def contigs_enrichment_per_bin(contig_counts):
    '''
    '''
    counts = contig_counts['counts'].tolist()
    mean_counts = statistics.mean(counts)
    enrichments = [i/mean_counts for i in counts]
 
    pseudocount = 0.0000001
    log_enrichments = [math.log(i+pseudocount, 2) for i in enrichments]
 
    ## Add to dataframe
    contig_counts['enrichments'] = enrichments
    contig_counts['log2_enrichments'] = log_enrichments

    ## Plot enrichments
    fig = plt.figure(figsize=(5,5))
    sns.histplot(data=contig_counts, x="enrichments")
    outFile = outDir + '/plots/contig_counts_enrichments.pdf'
    fig.savefig(outFile)

    ## Plot log enrichments
    fig = plt.figure(figsize=(5,5))
    sns.histplot(data=contig_counts, x="log2_enrichments")
    outFile = outDir + '/plots/contig_counts_log2_enrichments.pdf'
    fig.savefig(outFile)

    return contig_counts


def select_merge_enriched_bins(contig_counts, outDir):
    '''
    '''
    ##Select bins with an enrichment >= 2
    contig_counts_filt = contig_counts[contig_counts['log2_enrichments']>2]

    ##Merge bins into clusters
    clusters = []

    for index, row in contig_counts_filt.iterrows():

        ## Initialize list of clusters
        if not clusters:
            cluster = [row['chr'], row['bin'], (row['bin'] + 1), row['counts'], row['enrichments'], row['log2_enrichments'], [str(row['bin'])], row['contigIds']]
            clusters.append(cluster)
            continue

        ##Previous cluster
        prev_cluster = clusters[-1]
        prev_chr = prev_cluster[0]
        prev_bin = prev_cluster[2]

        ##Add bin to cluster
        if (row['chr'] == prev_chr) and (row['bin'] == prev_bin):
            clusters[-1]  = [prev_chr, prev_cluster[1], (row['bin'] + 1), (prev_cluster[3] + row['counts']), ((prev_cluster[4] + row['enrichments'])/2), ((prev_cluster[5] + row['log2_enrichments'])/2), (prev_cluster[6] + [str(row['bin'])]), (prev_cluster[7] + row['contigIds'])]

        ## Create new cluster containing the bin
        else:
            cluster = [row['chr'], row['bin'], (row['bin'] + 1), row['counts'], row['enrichments'], row['log2_enrichments'], [str(row['bin'])], row['contigIds']]
            clusters.append(cluster)

    ##Create dataframe with the clusters info
    binSize = 1000000 #1Mb

    for cluster in clusters:
        cluster[1] = cluster[1] * binSize
        cluster[2] = cluster[2] * binSize
        cluster[6] = ','.join(cluster[6])

    clustersDf = pd.DataFrame(clusters, columns = ['chr', 'start', 'end', 'counts', 'enrichments', 'log2_enrichments', 'binIds', 'contigIds'])
    
    outFile = outDir + '/contig_clusters.tsv'
    clustersDf.to_csv(outFile, sep="\t", index=False)

    
######################
## Get user's input ##
######################

## 1. Define parser ##
parser = argparse.ArgumentParser(description='')
parser.add_argument('PAF', help='Path to paf file containing contig alignments on the reference')
parser.add_argument('chrSizes', help='Path to tsv file containing two columns: chromId chromLen')
parser.add_argument('outDir', help='Output directory')

## 2. Parse user input ##
args = parser.parse_args()
PAF = args.PAF
chrSizes = args.chrSizes
outDir = args.outDir

## 3. Display configuration to standard output ##
scriptName = os.path.basename(sys.argv[0])
scriptName = os.path.splitext(scriptName)[0]

print()
print('***** ', scriptName, 'configuration *****')
print('PAF: ', PAF)
print('chrSizes: ', chrSizes)
print('outDir: ', outDir, "\n")


##########
## MAIN ##
##########

## 0. Create directory for plots
#################################
plotDir = outDir + '/plots/'
unix.mkdir(plotDir)


## 1. Create dictionary with chromosome lengths
####################################################
msg = '1. Create dictionary with chromosome lengths'
log.header(msg)
sizes_dict = read_chr_sizes(chrSizes)


## 2. Read paf file and organize hits by chromosome
####################################################
msg = '2. Read paf file and organize hits by chromosome'
log.header(msg)

##Organize hits into a dictionary by chromosome
hits_dict = organize_hits_by_chrom(PAF)

##filter dict
hits_dict = { chrId: hits_dict[chrId] for chrId in sizes_dict.keys()}


## 3. Count number of contigs/hits per chromosome
###################################################
msg = '3. Count number of contigs and hits per chromosome'
log.header(msg)

## Hits
nbHits_per_chr(hits_dict, plotDir)

##Contigs
nbContigs_per_chr(hits_dict, plotDir)


## 4. Plot number of contigs/hits per genomic bin along each chromosome
########################################################################
msg = '4. Plot number of contigs/hits per genomic bin along each chromosome'
log.header(msg)

##Organize hits in a bin database
hitsBinDb = hits2binDb(sizes_dict, hits_dict)

'''
##Count number of hits per 1Mb genomic bins along each chromosome
hit_counts = hits_per_bin(sizes_dict, hitsBinDb)

##Plot hits per bin for each chr
plot_counts_per_bin(hit_counts, 'hits', sizes_dict, plotDir)
'''

##Count number of contigs per 1Mb genomic bins along each chromosome
contig_counts = contigs_per_bin(sizes_dict, hitsBinDb)

'''
##Plot contigs per bin for each chr
plot_counts_per_bin(contig_counts, 'contigs', sizes_dict, plotDir)
'''

## 5. For each genomic bin assess enrichment on the number of contigs
#######################################################################
msg = '5. For each genomic bin assess enrichment on the number of contigs'
log.header(msg)

contig_counts = contigs_enrichment_per_bin(contig_counts)

###Select enriched genomic bins and merge
select_merge_enriched_bins(contig_counts, outDir)


'''
## 6. Plot contig alignments along each chromosome
###################################################
msg = '6. Plot contig alignments along each chromosome'
log.header(msg)

plot_contig_alignment_intervals(hits_dict, plotDir)
'''

