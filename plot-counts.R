counts_filename <- '/Users/davids/OneDrive - The University of Melbourne/WEHI Rory CoVid19 RNA/DB CloudStor 20200812/Standard_R1.fastq.gz_UMIHashTrie.pkl_UMI-counts.csv'
counts <- read.csv(file = counts_filename, header = FALSE)
colnames(counts) <- c('UMI', 'count')

library(ggplot2)
rows.to.use <- nrow(counts)
# rows.to.use <- 20000
ggplot(counts[1:rows.to.use,], aes(x = 1:rows.to.use, y = count)) +
  geom_line() +
  scale_x_log10() + 
  scale_y_log10()