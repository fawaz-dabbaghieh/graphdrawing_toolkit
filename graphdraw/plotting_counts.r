library(ggplot2)
library(plotly)

setwd("/home/fawaz/projects/cancer_subclones/2d_histogram/new_counts/")
########### node counts
node_counts <- read.csv("filtered_99id_99len_counts_node_counts.csv", sep=",")
node_counts$normalized_cancer_counts = node_counts$cancer_counts/node_counts$node_length
node_counts$normalized_healthy_counts = node_counts$healthy_counts/node_counts$node_length


filtered_node_counts <- filter(node_counts, cancer_counts < 300 & healthy_counts < 300)
filtered_node_counts <- filter(filtered_node_counts, cancer_counts > 25 & healthy_counts > 25)
filtered_node_counts$normalized_cancer_counts = filtered_node_counts$cancer_counts/filtered_node_counts$node_length
filtered_node_counts$normalized_healthy_counts = filtered_node_counts$healthy_counts/filtered_node_counts$node_leng


p<- ggplot(node_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts aligned to MBG mixed graph") +
  theme_bw()
p

p<- ggplot(node_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts aligned to MBG mixed graph") +
  theme_bw()
p

p<- ggplot(filtered_node_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts aligned to MBG mixed graph (filtered)") +
  theme_bw()
p
ggplotly(p)

###### normalized counts
p<- ggplot(node_counts, aes(x=normalized_cancer_counts, y=normalized_healthy_counts) ) +
  geom_bin2d(bins = 300) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node normalized counts aligned to MBG mixed graph") +
  theme_bw()
p
ggplotly(p)

p<- ggplot(filtered_node_counts, aes(x=normalized_cancer_counts, y=normalized_healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node normalized and filtered counts aligned to MBG mixed graph") +
  theme_bw()
p
ggplotly(p)

########### edge counts
edge_counts <- read.csv("filtered_99id_99len_counts_edge_counts.csv", sep=",")
filtered_edge_counts <- filter(edge_counts, cancer_counts < 300 & healthy_counts < 300)


p<- ggplot(edge_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi edge counts aligned to MBG mixed graph") +
  theme_bw()
p

ggplot(filtered_edge_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi edge counts aligned to MBG mixed graph (filtered)") +
  theme_bw()



#########################################################


setwd("/home/fawaz/projects/cancer_subclones/2d_histogram/mbg_read_paths/")
########### node counts
node_counts <- read.csv("read_paths_node_counts.csv", sep=",")
node_counts$normalized_cancer_counts = node_counts$cancer_counts/node_counts$node_length
node_counts$normalized_healthy_counts = node_counts$healthy_counts/node_counts$node_length

filtered_node_counts <- filter(node_counts, cancer_counts < 300 & healthy_counts < 300)
filtered_node_counts$normalized_cancer_counts = filtered_node_counts$cancer_counts/filtered_node_counts$node_length
filtered_node_counts$normalized_healthy_counts = filtered_node_counts$healthy_counts/filtered_node_counts$node_length


p<- ggplot(node_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts MBG paths") +
  theme_bw()
p

p<- ggplot(filtered_node_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts MBG paths (filtered)") +
  theme_bw()
p


p<- ggplot(node_counts, aes(x=normalized_cancer_counts, y=normalized_healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts MBG paths") +
  theme_bw()
p

p<- ggplot(filtered_node_counts, aes(x=normalized_cancer_counts, y=normalized_healthy_counts) ) +
  geom_bin2d(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi node counts MBG paths (filtered)") +
  theme_bw()
p

########### edge counts
edge_counts <- read.csv("read_paths_edge_counts.csv", sep=",")
filtered_edge_counts <- filter(edge_counts, cancer_counts < 100 & healthy_counts < 100)

p<- ggplot(edge_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi edge counts MBG paths") +
  theme_bw()
p

ggplot(filtered_edge_counts, aes(x=cancer_counts, y=healthy_counts) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  ggtitle("Cancer/healthy Hifi edge counts MBG paths (filtered)") +
  theme_bw()
