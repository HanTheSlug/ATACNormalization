library(edgeR)
library(data.table)
library(optparse)

# Optparse
option_list = list(
  make_option(c("-i", "--infile"), type = "character", default = NULL, 
              help = "Path to CSV file containing raw TN5 counts", metavar = "character"),
  make_option(c("-o", "--outfile"), type = "character", default = NULL, 
              help = "Output CSV file containing TMM normalized TN5 counts", metavar = "character")
)

parser = OptionParser(option_list = option_list)
args = parse_args(parser)

# Check if input and output files exist
if (is.null(args$infile)) {
    stop("Missing raw TN5 Count file", call. = FALSE)
}
if(is.null(args$outfile)) {
    stop("Missing output file path", call. = FALSE)
}

# Read the input file
count_data <- data.frame(fread(args$infile), row.names = 1)

# Create DGEList object. Object will be used for TMM normalization
dge <- DGEList(counts = count_data)

# Perform TMM normalization. Generates TMM normalizing factors 
dge <- calcNormFactors(dge, method = "TMM")

# Calculate the normalized counts
normalized_counts <- dge$counts * dge$samples$norm.factors

# Save the TMM normalized TN5 counts to output file.
write.csv(normalized_counts, file = args$outfile)