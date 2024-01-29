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
setnames(sample_table, c("SRX", "SRR"))

# Read in train,test,valid h5 files
train_z <- h5read(args$infile, 'train_z')
train_y <- h5read(args$infile, 'train_y')
test_z <- h5read(args$infile, 'test_z')
test_y <- h5read(args$infile, 'test_y')
valid_z <- h5read(args$infile, 'valid_z')
valid_y <- h5read(args$infile, 'valid_y')

# Convert to data tables for efficient merging and manipulation
dt_z_train <- as.data.table(train_z)
dt_y_train <- as.data.table(train_y)
dt_z_test <- as.data.table(test_z)
dt_y_test <- as.data.table(test_y)
dt_z_valid <- as.data.table(valid_z)
dt_y_valid <- as.data.table(valid_y)

# Function that manipulates h5 data into a data frame with rows: samples columns: peaks
# and each entry is the tn5 count. 
# Also computes the average tn5 count value for each sample
makeCountDF <- function(dt_z, dt_y,sample_table) {
  
  # Transpose dt_z and rename columns
  dt_z_t <- transpose(dt_z)
  setnames(dt_z_t, c("sample", "chromosome", "peak"))
  
  # Add one to indicies in sample since R indexes from 1
  dt_z_t$sample = dt_z_t$sample + 1
  # Make sure sample names are characters
  dt_z_t[, sample := as.character(sample)]
  
  # Merge dt_z_t and dt_y by adding a 'value' column
  dt_merged <- cbind(dt_z_t, value = dt_y[["value"]])
  
  # Rename the column with the peak positions in bed file.
  dt_merged[, peak_column_name := paste(chromosome, peak, sep = "_")]
  # Create a wide data table where columns are peaks and rows are samples.
  # Each entry is the tn5 value. If a peak does not exist in one sample but exists in another, create a peak and set it to 0.
  count_df <- dcast(dt_merged, sample ~ peak_column_name, value.var = "value", fun.aggregate = sum, fill = 0)
  # Re order by how samples appear in input samples.txt file
  count_df[, sample := as.numeric(as.character(sample))]
  count_df <- count_df[order(sample)]
  # Rename rownames to actual sample names
  indices <- as.numeric(count_df$sample)
  srx_values <- sample_table$SRX[indices]
  rownames(count_df) <- srx_values
  count_df$sample <- NULL
  count_df <- t(count_df)
  
  # Calculate the average per sample
  average_per_sample <- colMeans(count_df)
  return(list(count_df = count_df, average_per_sample = average_per_sample))
}

# Get TMM normalizing factors
getTMMNormalizationFactors <- function(count_df) {
    # Create a DGE object using training data counts
    dge <- DGEList(counts = count_df)
    # Perform TMM normalization. Generates TMM normalizing factors 
    dge <- calcNormFactors(dge, method = "TMM")

    # Get the scaling factors
    return(dge$samples$norm.factors)
}

# For train dataset
count_df_and_average_train <- makeCountDF(dt_z_train,dt_y_train,sample_table)
average_per_sample_train <- count_df_and_average_train$average_per_sample
count_df_train <- count_df_and_average_train$count_df

# Normalizing factor
scaling_factors <- getTMMNormalizationFactors(count_df_train)

# For test dataset
count_df_and_average_test <- makeCountDF(dt_z_test,dt_y_test,sample_table)
average_per_sample_test <- count_df_and_average_test$average_per_sample
count_df_test <- count_df_and_average_test$count_df

# For valid dataset
count_df_and_average_valid <- makeCountDF(dt_z_valid,dt_y_valid,sample_table)
average_per_sample_valid <- count_df_and_average_valid$average_per_sample
count_df_valid <- count_df_and_average_valid$count_df

# Find indices in average_per_sample_train that correspond to the closest values in average_per_sample_test
# These indices are used to determine the closest scaling factor from the train set to apply to the test set to minimize loss

# For test
closest_scaling_factor_test <- sapply(average_per_sample_test, function(x) {
    which.min(abs(average_per_sample_train - x))
})

# For valid
closest_scaling_factor_valid <- sapply(average_per_sample_valid, function(x) {
    which.min(abs(average_per_sample_train - x))
})

# Normalizing Factors
print(scaling_factors)
print(scaling_factors[closest_scaling_factor_test])
print(scaling_factors[closest_scaling_factor_valid])
