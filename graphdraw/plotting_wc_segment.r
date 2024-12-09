library(ggplot2)
library(data.table)
library(scales)
library(assertthat)
library(glue)


in_range <- function(x, start, end) return(x >= start & x <= end)
format_Mb = function(x) {paste(comma(round(x/1e6,1)), "Mb")}


print_usage <- function() {
  message("                                                                                ")
  message("Usage: Rscript plot_segments.R count-file chrom start end out.pdf            ")
  message("                                                                                ")
  message("       Plot different segmentations for a single chromosome and coordinates.                    ")
  message("                                                                                ")
}


###########
# Arguments
args = commandArgs(trailingOnly = T)
if (length(args) != 5) {
  print_usage()
  stop("Stopped execution on purpose.")
}

# getting inputs
f_counts = args[1]
CHROM = args[2]
query_start = as.integer(args[3])
query_end = as.integer(args[4])
f_out = args[5]
zcat_command = "zcat"
invisible(assert_that(grepl('\\.pdf$', f_out)))


# reading and checking counts table
message(" * Reading count data ", f_counts, "...")
if (grepl('\\.gz$',f_counts)) {
  counts = fread(cmd=paste(zcat_command, f_counts))
} else {
  counts = fread(f_counts)
}
invisible(assert_that(all(c("chrom","start","end","class","sample",
                            "cell","w","c") %in% colnames(counts))))
setkey(counts, chrom)

# filter for the given chromosome
invisible(assert_that(CHROM %in% unique(counts$chrom)))
counts <- counts[chrom == CHROM]

# make sure the segment coordinates fit the chromosome
chrom_start = min(counts[[2]])
chrom_end = max(counts[[3]])
if (query_start < chrom_start){
  query_start <- chrom-start
}
if (query_end > chrom_end){
  query_end <- chrom_end
}


# filter counts based on interval given
filtered_counts <- counts[start >= query_start & end <= query_end]
invisible(assert_that(length(filtered_counts[[1]]) > 0))


# getting the cell ids
cells <- unique(filtered_counts[[5]])

cairo_pdf(f_out, width=14, height=10, onefile = T)

for (current_cell in cells){

  current_sample <- filtered_counts[cell == current_cell]
  plt <- ggplot(current_sample) + ggtitle(glue("Chrom {CHROM} coordinates {query_start} - {query_end}, cell {current_cell}")) +
    xlab("Chromosome Location") + ylab("SSeq counts")
  print(glue("Plotting cell {current_cell} with length {length(current_sample[[1]])}"))
  
  i = 0
  while (i <= length(current_sample[[1]])){
    start <- current_sample[i][[2]]
    end <- current_sample[i][[3]]
    
    w <- current_sample[i][[7]]
    c <- current_sample[i][[8]]
    
    plt <- plt + geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
      geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = c), fill='paleturquoise4')
    i <- i + 1
  }
  print(plt)
  
}
dev.off()


if (FALSE){
  
  plt <- ggplot(current_sample) + ggtitle(glue("Chrom {CHROM} coordinates {query_start} - {query_end}")) +
    xlab("Chromosome Location") + ylab("SSeq counts")
  i = 0
  while (i <= length(filtered_counts[[1]])){
    
    start <- filtered_counts[i][[2]]
    end <- filtered_counts[i][[3]]
    
    w <- filtered_counts[i][[7]]
    c <- filtered_counts[i][[8]]
    
    plt <- plt + geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = -w), fill='sandybrown') +
      geom_rect(aes(xmin = start, xmax=end, ymin=0, ymax = c), fill='paleturquoise4')
    i <- i + 1
  }
}