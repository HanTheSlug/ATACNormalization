library(ggplot2)
library(PCAForQTL)
library(data.table)
library(optparse)

# Optparse
option_list = list(
  make_option(c("-i", "--infile"), type = "character", default = NULL, 
              help = "Path to CSV file containing TMM normalized TN5 counts", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = NULL, 
              help = "Path to output directory that contains a CSV of Principle Components, a CSV of Principle Component Loadings, and Elbow Plot", metavar = "character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$infile)) {
    stop("Missing TN5 Count file", call. = FALSE)
}
if(is.null(opt$out)) {
    stop("Missing output directory path", call. = FALSE)
}

# BE algorithm
BEAlgorithm <- function(data) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  resultRunBE <- PCAForQTL::runBE(data, B=20, alpha=0.05)
  return(resultRunBE$numOfPCsChosen)
}

# Create elbow plot of PCs
createElbowPlot <- function(pca_result, filepath) {
  pc_stddev <- pca_result$sdev
  elbow_data <- data.frame(PC = 1:length(pc_stddev), StdDev = pc_stddev)
  elbow_plot <- ggplot(elbow_data, aes(x = PC, y = StdDev)) +
    geom_line() +
    geom_point() +
    labs(title = "Elbow Plot for Principal Components",
         x = "Principal Component",
         y = "Standard Deviation") +
    theme_minimal()
  
  ggsave(filepath, elbow_plot, width = 6, height = 4, dpi = 300)
}


tmm_tn5_counts_df <- data.frame(fread(opt$infile), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
tmm_tn5_counts_df <- tmm_tn5_counts_df[,0:(ncol(tmm_tn5_counts_df)-1)]
#last_column_name <- colnames(tmm_tn5_counts_df)[ncol(tmm_tn5_counts_df)]
tmm_tn5_counts_df <- apply(tmm_tn5_counts_df, 2, function(x) as.numeric(as.character(x)))
peak_names <- rownames(tmm_tn5_counts_df)
tmm_tn5_counts_df <- t(tmm_tn5_counts_df)
tmm_tn5_counts_df <- tmm_tn5_counts_df[ , which(apply(tmm_tn5_counts_df, 2, var) != 0)]
row_names <- rownames(tmm_tn5_counts_df)

# Initialize output files
pc_file_path <- paste0(opt$out,"/pc_table.csv")
pc_loadings_file_path <- paste0(opt$out,"/pc_loadings.csv")
elbow_plot_path <- paste0(opt$out,"/elbow_plot.png")

# PCA
pca_result <- prcomp(tmm_tn5_counts_df[, -1], scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)
pc_stddev <- pca_result$sdev
pc_loadings <- as.data.frame(pca_result$rotation)
rownames(pca_data) <- row_names
rownames(pc_loadings) <- peak_names

write.csv(pca_data, file = pc_file_path)
write.csv(pc_loadings, file = pc_loadings_file_path)

# Run BE algorithm
opt_num_pcs <- BEAlgorithm(tmm_tn5_counts_df[, -1])
print(paste("The Optimal Number of PCs determined by the BE Algorithm is:",opt_num_pcs))

createElbowPlot(pca_result, elbow_plot_path)