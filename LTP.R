# Initial stuff
library(rvest)
library(magrittr)
library(Biostrings)
library(seqRFLP)
library(xml2)

# Function to make data.frame from FASTA
makeProtSeqDF <- function(stringset) {
  first_DF <- readAAStringSet(stringset)
  seq_name <- names(first_DF)
  sequence <- paste(first_DF)
  protSeqDF <- data.frame(seq_name, sequence, stringsAsFactors = FALSE)
  return(protSeqDF)
}

######################################################################################
# Part 1: Get names of allergens
html_lipid <- read_html("http://www.allergen.org/search.php?allergenname=&allergensource=&TaxSource=&TaxOrder=&foodallerg=all&bioname=lipid")
html_ltp <- read_html("http://www.allergen.org/search.php?allergenname=&allergensource=&TaxSource=&TaxOrder=&foodallerg=all&bioname=ltp")

# Get the names for the different allergens from the table
population_lipid <- html_lipid %>% html_nodes(xpath='//*[@id="example"]') %>%
  html_table %>% extract2(1)
population_ltp <- html_ltp %>% html_nodes(xpath='//*[@id="example"]') %>% 
  html_table %>% extract2(1)
alrgNames <- unique(c(population_lipid[population_lipid[,1] == "", 2], 
                      population_ltp[population_ltp[,1] == "", 2]))

# Get the links from the allergens from the table and them full links as input 
# for the next part
links <- union(
  paste0("http://www.allergen.org/", html_lipid %>% html_nodes(xpath = "//td/a") %>%
           html_attr("href")),
  paste0("http://www.allergen.org/", html_ltp %>% html_nodes(xpath = "//td/a") %>%
           html_attr("href"))
)
# Follow these links to get the Uniprot IDs
alr_up_list <- sapply(links,
                      function(u) u %>% read_html %>%
                        html_nodes(xpath = '//*[@id="isotable"]') %>% html_table() %>%
                        extract2(1) %>% extract("UniProt")
)
# Convert to vector and rename
alr_up <- setNames(unlist(alr_up_list, use.names=F),rep(alrgNames, lengths(alr_up_list)))
# Clean up and rename funny Ids
alr_up <- alr_up[! is.na(alr_up)]
alr_up[! alr_up == ""]
alr_up[alr_up == "P19656-1"] <- "P19656"
alr_up[alr_up == "P19656-2"] <- "P19656"
alr_up[alr_up == "C0HJR7"] <- "A0A161AT60"
alr_up[alr_up == "Q7X9Q5"] <- "P82007"
save(alr_up, file="alrgn_uniprot.RData")

######################################################################################
# Part 2: Load raw-data and subset allergens
ltp1 <- makeProtSeqDF("http://www.ebi.ac.uk/interpro/entry/IPR000528/proteins-matched?export=fasta")
ltp2 <- makeProtSeqDF("http://www.ebi.ac.uk/interpro/entry/IPR033872/proteins-matched?export=fasta")
ltp <- rbind(ltp1, ltp2)

alg_ltp <- subset(ltp, ltp[,1] %in% alr_up)
nlg_ltp <- ltp[ !(ltp$seq_name %in% alg_ltp$seq_name), ]
nlg_ltp <- nlg_ltp[1:540,]

alg_fa <- dataframe2fas(alg_ltp)
write.fasta(alg_fa, "alg_ltp.fa")
nlg_fa <- dataframe2fas(nlg_ltp)
write.fasta(nlg_fa, "nlg_ltp.fa")

