library(edgeR)
library(data.table)
library(optparse)
library(rhdf5)

# Assuming the input dataset file is a h5 dataset that follows the following format for train set:
# train_x : (index_in_train, haplotype_index, window_size, OHE_nucleotide_for_each_base_in_window) example dimensions: (123456,2,600,4)
# train_y : (index_in_train, class_index (1=tn5, 2=asoc), value) example dimensions: (123456,1,1)
# train_z : (index_in_train,3) (peaks, (sample_index,chromosome_index,peak_index)) example dimensions: (123456,3)

# Optparse
option_list = list(
  make_option(c("-s", "--samplefile"), type = "character", default = NULL, 
              help = "Path to TSV file containing sample names", metavar = "character"),
  make_option(c("-i", "--infile"), type = "character", default = NULL, 
              help = "Path to input h5 dataset file", metavar = "character"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL, 
              help = "Output CSV file containing TMM normalized TN5 counts", metavar = "character")
)

parser = OptionParser(option_list = option_list)
args = parse_args(parser)

# Check if input and output files exist
if (is.null(args$samplefile)) {
    stop("Missing sample file file", call. = FALSE)
}
if (is.null(args$infile)) {
    stop("Missing h5 dataset", call. = FALSE)
}
if(is.null(args$outfile)) {
    stop("Missing output file path", call. = FALSE)
}

# Read the sample file as a data table
sample_table <- fread(args$samplefile, header = FALSE)
setnames(sample_table, c("SRX", "SRR"))  # Assuming the first column is the sample identifier

train_z <- h5read(args$infile, 'train_z')
train_y <- h5read(args$infile, 'train_y')

# Convert to data tables for efficient merging and manipulation
dt_z <- as.data.table(train_z)
dt_y <- as.data.table(train_y)

# Transpose dt_z and rename columns
dt_z_t <- transpose(dt_z)
setnames(dt_z_t, c("sample", "chromosome", "peak"))
# Make sure sample names are characters
dt_z_t[, sample := as.character(sample)]

dt_merged <- cbind(dt_z_t, value = dt_y[["value"]])
# Rename the column with the peak positions in bed file.
dt_merged[, peak_column_name := paste(chromosome, peak, sep = "_")]
# Create a wide data table where columns are peaks and rows are samples. Each entry is the tn5 value
# If a peak does not exist in one sample but exists in another, create a peak and set it to 0. 
#final_df <- dcast(dt_merged, sample ~ peak_column_name, value.var = "value", fill = 0)
count_df <- dcast(dt_merged, sample ~ peak_column_name, value.var = "value", fun.aggregate = sum, fill = 0)
# rename rownames to actual sample names
indices <- as.numeric(count_df$sample)
srx_values <- sample_table$SRX[indices]
rownames(count_df) <- srx_values
count_df$sample <- NULL
count_df <- t(count_df)
average_per_sample <- colMeans(count_df)
print(average_per_row)
# Create DGEList object. Object will be used for TMM normalization
dge <- DGEList(counts = count_df)

# Perform TMM normalization. Generates TMM normalizing factors 
dge <- calcNormFactors(dge, method = "TMM")

# Calculate the normalized counts
normalized_counts <- dge$counts * dge$samples$norm.factors
print(length(dge$samples$norm.factors))
# #normalized_counts <- t(normalized_counts)
# # Save the TMM normalized TN5 counts to output file.
# write.csv(normalized_counts, file = args$outfile)