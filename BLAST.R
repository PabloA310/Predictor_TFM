# Initial stuff
library(Biostrings)
library(msaR)
library(caTools)
library(seqinr)
library(dplyr)

# Function to load Allergen data, split it and write fasta files
prepare_data_alg <- function(stringSet){
  data <- readAAStringSet(stringSet)
  split <- caTools::sample.split(data, SplitRatio = 0.90)
  train <- data[split]
  train_fasta <- msaR::as.fasta(train)
  seqinr::write.fasta(train_fasta, names(train_fasta), "train.fasta")
  test_alg <- data[!split]
  test_alg_fasta <- msaR::as.fasta(test_alg)
  seqinr::write.fasta(test_alg_fasta, names(test_alg_fasta), "test_alg.fasta")
}

# Function to load Non allergen data
prepare_data_nlg <- function(stringSet){
  data <- readAAStringSet(stringSet)
  split <- caTools::sample.split(data, SplitRatio = 0.90)
  test_nlg <- data[!split]
  test_nlg_fasta <- msaR::as.fasta(test_nlg)
  seqinr::write.fasta(test_nlg_fasta, names(test_nlg_fasta), "test_nlg.fasta")
}

# Function to select best e-values for each query
select_eValue <- function(alg, nlg){
  eValue_alg <- summarise(V2 = min(V2), group_by(alg, V1))
  eValue_alg <- data.frame("condition" = "Allergen", "condBin" = 1, eValue_alg)
  eValue_nlg <- summarise(V2 = min(V2), group_by(nlg, V1))
  eValue_nlg <- data.frame("condition" = "Non allergen", "condBin" = 0, eValue_nlg)
  eValues <- rbind(eValue_alg, eValue_nlg)
  colnames(eValues) <- c("condition", "condBin", "Name", "e-value")
  return(eValues)
}
############################################################################
# Part 1: Prepare data
setwd("~/TFM/R/Data")
prepare_data_alg("alg_ltp.fa")
prepare_data_nlg("nlg_ltp.fa")

############################################################################
# Part 2: Run BLAST
system("makeblastdb -in train.fasta -dbtype 'prot'")
system("blastp -db train.fasta -query test_alg.fasta -out blast_alg_ltp -outfmt '6 qseqid evalue'")
system("blastp -db train.fasta -query test_nlg.fasta -out blast_nlgC_ltp -outfmt '6 qseqid evalue'")

############################################################################
# Part 3: Select best e-values for each query
alg <- data.frame(read.table("blast_alg_ltp"))
nlg <- data.frame(read.table("blast_nlg_ltp"))

eValues <- select_eValue(alg, nlg)
write.table(eValues, "eValues_ltp.txt", row.names = F, quote = F, col.names = T, sep = '\t')
