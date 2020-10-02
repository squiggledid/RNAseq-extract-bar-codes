library(RcppCNPy)

# data location
project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
hamming.npy.filename = file.path(project_dir, "Standard_R1.fastq.gz_UMIHashTrie.pkl_hamming.npy")
levenshtein.npy.filename = file.path(project_dir, "Standard_R1.fastq.gz_UMIHashTrie.pkl_levenshtein.npy")

# formula to recover original dm dimensions
num.dm.rows <- function(length.dm.raw) {
  return((1 + sqrt(1 + 8*length.dm.raw))/2)
}

# # load the raw distance data saved from Python
# hamming.dm.raw <- npyLoad(hamming.npy.filename)
# levenshtein.dm.raw <- npyLoad(levenshtein.npy.filename)
# 
# # make square distance matrices
# # TODO go directly to the tidy format from this
# n.dm.rows <- num.dm.rows(length(hamming.dm.raw))
# hamming.dm <- matrix(nrow = n.dm.rows, ncol = n.dm.rows, data = 0)
# for (i in 1:(n.dm.rows - 1)) {
#   for (j in (i + 1):n.dm.rows) {
#     hamming.dm[i, j] <- hamming.dm.raw[(i - 1)*n.dm.rows + j - 1]
#     hamming.dm[j, i] <- hamming.dm[i, j]
#   }
# }

hamming.dm <- npyLoad(hamming.npy.filename)
levenshtein.dm <- npyLoad(levenshtein.npy.filename)

### create the heapmap plots
# appearance parameters
reverse.heat = FALSE
direction <- ifelse(reverse.heat, 1, -1) # We actually want viridis reverse as the default

plot.heatmap <- function(melted.distance.matrix, pdf.file.name) {
  # plot the heatmap
  require(ggplot2)
  require(viridis)
  heatmap <- ggplot(data = melted.distance.matrix, aes(x=Var1, y=Var2)) +
    geom_raster(aes(fill=value)) +
    theme_bw() + # get rid of the grey background
    coord_equal() + # make it square
    scale_fill_viridis(option = "magma", direction = direction)
  pdf(pdf.file.name)
  cat("Saving heatmap to", pdf.file.name, "\n")
  print(heatmap)
  dev.off()
}

require(reshape2) # now deprecated. Apparently replaced by a generic in data.table?
nrows.to.use <- 1000

melted.distance.matrix <- melt(hamming.dm[1:nrows.to.use,1:nrows.to.use]) # transform into "tidy" format for ggplot
pdf.file.name <- paste0(hamming.npy.filename, ".pdf")
plot.heatmap(melted.distance.matrix, pdf.file.name)

melted.distance.matrix <- melt(levenshtein.dm[1:nrows.to.use,1:nrows.to.use]) # transform into "tidy" format for ggplot
pdf.file.name <- paste0(levenshtein.npy.filename, ".pdf")
plot.heatmap(melted.distance.matrix, pdf.file.name)

