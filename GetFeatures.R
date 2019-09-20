# Initial stuff
library(Biostrings)
library(Peptides)
library(ggplot2)
library(optparse)

# Optparse
option_list <- list(
  make_option("--file1", action = "store", dest = "file1",
              help = "Input stringset ['name.fa'] for the positive condition."),
  make_option("--file0", action = "store", dest = "file0",
              help = "Input stringset ['name.fa'] for the negative condition."),
  make_option(c("-d", "--dframe"), action = "store", dest = "dframe",
              help = "Give a name to the data.frame with features. [Example: DF]"),
  make_option(c("-s", "--stdf"), action = "store_true", dest = "stdf",
              default = FALSE, help = "Make a stadistic data.frame.
              [default: %default]"),
  make_option(c("-b", "--boxp"), action = "store_true", dest = "boxplot",
              default = FALSE, help = "Make boxplot to compare features.
              [default: %default]"),
  make_option(c("-p", "--dplot"), action = "store_true", dest = "denplot",
              default = FALSE, help = "Make density to plots to compare features.
              [default: %default")
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

# Helper function to be called from within getProtFeatures
getAAComposition <- function(protSeqs) {
    AADataFrame <- data.frame(Peptides::aaComp(protSeqs))
    AADataFrame <- AADataFrame[, c(FALSE, TRUE)]
    AADataFrame <- data.frame(t(AADataFrame))
    rownames(AADataFrame) <- NULL
    colnames(AADataFrame) <- c("aaTiny_pc", "aaSmall_pc", "aaAliphatic_pc", "aaAromatic_pc",
                               "aaNonPolar_pc", "aaPolar_pc", "aaCharged_pc", "aaBasic_pc",
                               "aaAcidic_pc")
    return(AADataFrame)
}

# Function to obtain annotation
getProtFeatures <- function(protSeqDF) {
    protSeqs <- protSeqDF[, "sequence"]
    molecular_weight <- Peptides::mw(protSeqs, monoisotopic = FALSE)
    pepLength <- Peptides::lengthpep(protSeqs)
    isoelectric_point <- Peptides::pI(protSeqs, pKscale = "EMBOSS")
    instability <- Peptides::instaIndex(protSeqs)
    aliphaticIndex <- Peptides::aIndex(protSeqs)
    bomanIndex <- Peptides::boman(protSeqs)
    AADataFrame <- getAAComposition(protSeqs)
    protNames <- protSeqDF[, 1]
    protFeatsDF <- data.frame(protNames, molecular_weight, pepLength, 
                              isoelectric_point, instability, aliphaticIndex, bomanIndex, 
                              AADataFrame)
    return(protFeatsDF)
}

finalDF_seqs <- function(seqs, varName, condBin) {
    condition_seqs <- rep(varName, length(seqs[, 1]))
    binary_condition <- rep(condBin, length(seqs[, 1]))
    data.frame("condition" = condition_seqs, "binCondition" = binary_condition, seqs)
}
# Function to make boxplots for each feature
makeBoxplot <- function(DF, feature, yName){
  ggplot(DF, aes_string(x = DF$condition, y = feature)) + geom_boxplot() + ylab(yName)
}
getAllBoxplots <- function(DF){
  col_names_to_plot <- names(DF)[-c(1:3)]
  boxplot_list <- list()
  for (col_name in col_names_to_plot) {
    boxplot_list[[col_name]] <- makeBoxplot(DF, feature = DF[[col_name]], 
                                            yName = col_name)
  }
    pdf(paste(opt$dframe, "_boxplots.pdf", sep = ""))
    for (col_name in col_names_to_plot) {
      print(boxplot_list[[col_name]])
    }
}
# Function to make T-test and Mann Whitney test for each feature
makeAllTest <- function(DF) {
  col_names_to_test <- names(DF)[-c(1:3)]
    
  tTest_list <- list()
  for (col_name in col_names_to_test) {
    tTest_list[[col_name]] <- t.test(DF[[col_name]] ~ condition, data = DF)
  }
  mwTest_list <- list()
  for (col_name in col_names_to_test) {
    mwTest_list[[col_name]] <- wilcox.test(DF[[col_name]] ~ condition, data = DF)
  }
  test_list <- list("tTest_list" = tTest_list, "mwTest_list" = mwTest_list)
  return(test_list)
}
# Function to make a statistic data.frame
presentStatistic <- function(DF, statList) {
    # Define feature matrix - for each feature (row) we will put the results of
    # the statistical tests
    featureM <- matrix(nrow=dim(DF)[2]-3, ncol = 8)
    colnames(featureM) <- c("Allergen mean", "Non allergen mean", 
                            "Log2 means", "P value mean", "Allergen median",
                            "Non allergen median", "Log2 medians",
                            "P value median")
    i <- 1
    for (protType in unique(DF$condition)) {
      featuresDF <- DF[DF$condition == protType, -c(1:3)]
        meanV <- colMeans(featuresDF)
        featureM[, i] <- meanV
        i <- i + 1
        j <- 1
        medianV <- vector()
        for (feat in names(featuresDF)) {
            medianV[j] <- median(featuresDF[, feat])
            j <- j + 1
          }
        featureM[, i] <- medianV
        i <- i + 3
    }
    pValueMean <- vector(length=dim(featuresDF)[2])
    names(pValueMean) <- names(featuresDF)
    for (feat in names(featuresDF)) {
      pValueMean[feat] <- statList[["tTest_list"]][[feat]][[3]]
    }
    featureM[, "P value mean"] <- pValueMean
    pValueMedian <- vector(length=dim(featuresDF)[2])
    names(pValueMedian) <- names(featuresDF)
    for (feat in names(featuresDF)) {
      pValueMedian[[feat]] <- c(statList[["mwTest_list"]][[feat]][[3]])
    }
    featureM[, "P value median"] <- pValueMedian
    log2Means <- log2(featureM[, "Allergen mean"]/featureM[, "Non allergen mean"])
    featureM[, "Log2 means"] <- log2Means
    log2Medians <- log2(featureM[, "Allergen median"]/featureM[, "Non allergen median"])
    featureM[, "Log2 medians"] <- log2Medians
    row.names(featureM) <- names(featuresDF)
    as.data.frame(featureM)
}
# Function to make density curves for each feature
makeDensityPlot <- function(DF, feature, xName) {
  ggplot(DF, aes_string(x = feature,
                        color = DF$condition)) + geom_density()+ xlab(xName)
}
getAllDensityPlots <- function(DF) {
  col_names_to_plot <- names(DF[-c(1:3)])
  density_plot_list <- list()
  for (col_name in col_names_to_plot) {
    if (col_name == "molecular_weight" | col_name == "pepLength") {
      density_plot_list[[col_name]] <- makeDensityPlot(DF,
                                                       feature = log2(DF[[col_name]]), 
                                                       xName = col_name)
      
    } else { density_plot_list[[col_name]] <- makeDensityPlot(DF, 
                                                              feature = DF[[col_name]],
                                                              xName = col_name)
    } 
  }
  pdf(paste(opt$dframe, '_density_plots.pdf', sep = ""))
  for (col_name in col_names_to_plot) {
    print(density_plot_list[[col_name]])
  }
}

################################################################################
# Part 1: Load raw-data
algDF <- makeProtSeqDF(stringset = file.path(getwd(), opt$file1))
nalgDF <- makeProtSeqDF(stringset = file.path(getwd(), opt$file0))

################################################################################
# Part 2: Obtain feature annotation for each sequence using Peptides functions
algDF_features <- getProtFeatures(protSeqDF = algDF)
nalgDF_features <- getProtFeatures(protSeqDF = nalgDF)

################################################################################
# Part 3: Produce final data.frame
algsDF <- finalDF_seqs(seqs = algDF_features, varName = "Allergen",
                        condBin = 1)
nalgsDF <- finalDF_seqs(seqs = nalgDF_features, varName = "Non allergen",
                         condBin = 0)
DF <- rbind(algsDF, nalgsDF)

write.table(DF, paste(opt$dframe, ".txt", sep = ""), row.names = F,
            quote = F, sep = '\t')

#################################################################################
# Part 5: Produce boxplots and density plots for each feature and statistical tests
# Make Boxplot for features
if (opt$boxplot == TRUE) {
  getAllBoxplots(DF = DF)
}
# Represent densitiy for each feature and condition
if (opt$denplot == TRUE) {
  densityPlot_list <- getAllDensityPlots(DF = DF)
}
# Make T-test and Mann Whitney test for features
stat <- makeAllTest(DF = DF)
##############################################################################
# Part 6: Make a statistic dataframe
if (opt$stdf == TRUE) { 
  statDF <- presentStatistic(DF = DF, statList = stat)
  write.table(statDF, paste("stat", opt$dframe, ".tsv", sep = ""), row.names = T,
              quote = F,col.names = NA, sep = '\t')
}