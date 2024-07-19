# 
# Differential expression analysis with the DESeq2 package.
# 
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
#

# Load the library.
suppressPackageStartupMessages(library(DESeq2))

# Accept positional arguments.
# First argument: path to input count file.
# Second argument: path to output heatmap.
args <- commandArgs(trailingOnly = TRUE)

# The name of the file that contains the counts.
counts_file <- args[1]

# The sample file is in CSV format and must have the headers "sample" and "condition".
design_file = "data/design.csv"

# The final result file.
output_file <- args[2]

# Read the sample file.
colData <- read.csv(design_file, stringsAsFactors=F)

# Turn conditions into factors.
colData$condition = factor(colData$condition)

# The first level should correspond to the first entry in the file!
# Required later when building a model.
colData$condition = relevel(colData$condition, toString(colData$condition[1]))

# Isolate the sample names.
sample_names <- colData$sample

# Read the data from the standard input.
df = read.csv(counts_file, header=TRUE, row.names=1 )

# Created rounded integers for the count data
countData = round(df[, sample_names])

# Other columns in the dataframe that are not sample information. 
otherCols = df[!(names(df) %in% sample_names)]

#
# Running DESeq2
#

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)

# Run deseq
dse = DESeq(dds)

# Format the results.
res = results(dse)

#
# The rest of the code is about formatting the output dataframe.
#

# Turn the DESeq2 results into a data frame.
data = cbind(otherCols, data.frame(res))

# Create the foldChange column.
data$foldChange = 2 ^ data$log2FoldChange

# Rename columns to better reflect reality.
names(data)[names(data)=="pvalue"] <-"PValue"
names(data)[names(data)=="padj"] <- "FDR"

# Create a real adjusted pvalue
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Sort the data by PValue to compute false discovery counts.
data = data[with(data, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table.
data$falsePos = 1:nrow(data) * data$FDR

# Create the additional columns that we wish to present.
data$baseMeanA = 1
data$baseMeanB = 1

# Get the normalized counts.
normed = counts(dse, normalized=TRUE)

# Round normalized counts to a single digit.
normed = round(normed, 1)

# Merge the two datasets by row names.
total <- merge(data, normed, by=0)

# Sort again for output.
total = total[with(total, order(PValue, -foldChange)), ]

# Sample names for condition A
col_names_A = data.frame(split(colData, colData$condition)[1])[,1]

# Sample names for condition B
col_names_B = data.frame(split(colData, colData$condition)[2])[,1]

# Create the individual baseMean columns.
total$baseMeanA = rowMeans(total[, col_names_A])
total$baseMeanB = rowMeans(total[, col_names_B])

# Bringing some sanity to numbers. Round columns to fewer digits.
total$foldChange = round(total$foldChange, 3)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean  = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB =  round(total$baseMeanB, 1)
total$lfcSE = round(total$lfcSE, 2)
total$stat = round(total$stat, 2)
total$FDR = round(total$FDR, 4)
total$falsePos = round(total$falsePos, 0)

# Reformat these columns as string.
total$PAdj = formatC(total$PAdj, format = "e", digits = 1)
total$PValue = formatC(total$PValue, format = "e", digits = 1)

# Rename the first column.
colnames(total)[1] <- "name"

# Reorganize columns names to make more sense.
new_cols = c("name", names(otherCols), "baseMean","baseMeanA","baseMeanB","foldChange",
             "log2FoldChange","lfcSE","stat","PValue","PAdj", "FDR","falsePos",col_names_A, col_names_B)

# Slice the dataframe with new columns.
total = total[, new_cols]

# Write the results to the standard output.
write.csv(total, file=output_file, row.names=FALSE, quote=FALSE)

# Inform the user.
print("# Tool: DESeq2")
print(paste("# Design: ", design_file))
print(paste("# Input: ", counts_file))
print(paste("# Output: ", output_file))


