# Initial stuff
library(Biostrings)
library(optparse)

# Optparse
option_list <- list(
  make_option("--file1", action = "store", dest = "file1",
              help = "Input stringset ['name.fa'] for the positive condition."),
  make_option("--file0", action = "store", dest = "file0",
              help = "Input stringset ['name.fa'] for the negative condition."),
  make_option(c("-d", "--dframe"), action = "store", dest = "dframe",
              help = "Give a name to the data.frame with sequences. [Example: DF]")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Function to make dataframe from FASTA
makeProtSeqDF <- function(stringset) {
  first_DF <- readAAStringSet(stringset)
  seq_name <- names(first_DF)
  sequence <- paste(first_DF)
  protSeqDF <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  return(protSeqDF)
}

#######################################################################################
# Part 1: Load raw-data
algDF = makeProtSeqDF(stringset = file.path(getwd(), opt$file1))
nalgDF = makeProtSeqDF(stringset = file.path(getwd(), opt$file0))

# Part 2: Obtain sequence data.frame
algDF_seq = data.frame("binCondition" = 1, "sequence" = algDF[, 2])
nalgDF_seq = data.frame("binCondition" = 0, "sequence" = nalgDF[, 2])

write.table(algDF_seq, paste("algDFseq_", opt$dframe, ".txt", sep = ""),
            row.names = FALSE, quote = FALSE, sep = "\t")
write.table(nalgDF_seq, paste("nlgDFseq_", opt$dframe, ".txt", sep = ""), 
            row.names = FALSE, quote = FALSE, sep = "\t")
