# cancergraph
This is a toolkit that can work on genomic graphs in GAF format, mostly designed to investigate and manipulate graphs for further analysis.

## Installation

After cloning or copying the code, you can simply run:
```
pip install -e .
```

the `-e` makes the code editable without having to recompile every time you change something in the code.


## Subcommands
So far, `cancergraph` consists of 6 subcommands:
* `color_graph` which produces a csv file that colors nodes based on chromosomes
* `plot_contigs` generates a PDF file with contig information plots for each chromosome
* `two_colors` if there were two sets of reads aligned back to the graph, it will color nodes depending on which of the two samples was mostly represented in that node, i.e. used to color the two different subclones.
* `subgraph` extracts subgraph using the chromosome coordinates, the graph nodes to be aligned to a reference first
* `draw_contig` draws the alignments of a contig against a reference and against the contig itself
* `draw_subgraph` draws a subgraph against the reference


## Example graph

The graph above represents a bluntified (overlaps removed) of a de Bruijn graph with *k*-mer size of 9.
These 4 sequences are traced with the different dotted lines in the graph, and they construct a bubble chains of 3 bubbles.
2 simple bubbles and 1 superbubble with a nested simple bubble inside.

The following sections will demonstrate some examples of using this tool.

### color_graph
Inputs:
- a SAM or BAM file of the nodes alignments against a reference

Output:
- TSV or CSV with basic information for each node, each node ends up with one alignment

Usage:
```
cancergraph color_graph -b in_alignment.bam --out_csv outfile.csv
```

### plot_contigs
Inputs:
- a SAM or BAM file of the contig alignments

Output:
- a PDF file with all plots

Usage:
```
cancergraph plot_contigs -b in_bam.bam --out_prefix output_plots
```

### two_colors
Inputs:
- GAF alignments of the first sample
- GAF alignments of the second sample
- What color to use for the first sample
- What color to use for the second sample

Output:
- Output prefix, the output is a CSV compatible with Bandage

Usage:
```
cancergraph two_colors --gaf1 in_sample1.gaf --gaf2 in_sample2.gaf --color1 red  --color2 blue --out_prefix output
```

### subgraph
Inputs:
- In GFA graph
- In BAM or SAM of the node alignments of the input graph
- coordinates similar to the coordinate scheme used in samtools
- neighborhood size around the area of extraction, an integer, when 0, then only the nodes in the coordinates given are returned and the graph is not explored further

Output:
- output GFA file of the subgraph, also a CSV with some information is outputted

Usage:
```
cancergraph subgraph -g in_graph.gfa -b in_graph_alignments.bam --coord chr1:1000000-2000000 -o subgraph.gfa
```

### draw_contig
Inputs:
- PAF alignments of contigs
- Name of the contig interested in plotting
- maximum distance between alignments
- maximum overlap allowed between alignments (0-100)

Optional inputs:
- SVG width
- SVG height
- SVG legend
- no chaining (this for now only works when giving BAM or SAM as input alignment and not PAF), it will draw all alignments

Output:
- output prefix

Usage:
```
cancergraph draw_contig -a contig_alignments.paf -c h1tg000003l_cluster20 --out_prefix h1tg000003l_cluster20 --max_dist 3000 --max_overlap 30
```

This will draw all alignments against the reference without chaining or adding edges
```
cancergraph draw_contig -a contig_alignments.paf -c h1tg000003l_cluster20 --out_prefix h1tg000003l_cluster20 --max_dist 3000 --max_overlap 30 --no_chaining
```

### draw_subgraph
Inputs:
- Alignment file in BAM or SAM format of the graph nodes
- The graph file in GFA format
- coordinates similar to the coordinate scheme used in samtools

Optional inputs:
- SVG width
- SVG height
- SVG legend

Output:
- Output SVG file

Usage:
```
draw_subgraph -a graph_node_alignments.bam -g graph.gfa --coord chr1:1000000-2000000 --out_svg testing_subgraph_drawing
```