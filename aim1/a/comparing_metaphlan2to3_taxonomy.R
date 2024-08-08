# loading and installing required packages
list.of.packages <- c("tidyverse", "ggrepel", "plotly", "ggpubr", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# Reading the data
BASE_FOLDER <- '/pita/users/hila/biobakery_workflows/output_data'
path.metaphlan2 <- '/22_Aug_2021_HMP2_metagenome_samples_biobakery2/metaphlan/main/'
path.metaphlan3 <- '/16_Aug_2021_HMP2_metagenome_samples/metaphlan/main/'
#SAMPLE_ID <- 'CSM67UB1'
#SAMPLE_ID <- 'CSM67UF1_P'
#SAMPLE_ID <- 'CSM67UF1'
#SAMPLE_ID <- 'CSM7KOK1'
#SAMPLE_ID <- 'CSM7KOQ1'
SAMPLE_ID <- 'CSM9X211'
#SAMPLE_ID <- 'CSM5MCV1_P'
FILE_NAME <- paste(SAMPLE_ID, 'taxonomic_profile.tsv', sep = '_')
file.metaphlan2 <- paste(BASE_FOLDER, path.metaphlan2, FILE_NAME, sep="")
file.metaphlan3 <- paste(BASE_FOLDER, path.metaphlan3, FILE_NAME, sep="")

metaphlan2.data <- read.table(file.metaphlan2, sep = "\t", header = F, check.names = F)
colnames(metaphlan2.data) <- c("taxonomy", "abundance")

metaphlan3.data <- read.table(file.metaphlan3, sep = "\t", check.names = F)
colnames(metaphlan3.data) <- c("taxonomy", "ncbi", "abundance", "other.names")

taxonomy.data <- merge(x = metaphlan2.data, y = metaphlan3.data, by = "taxonomy", all = T)
names(taxonomy.data)[names(taxonomy.data) == "abundance.x"] <- 'abundance.mpln2'
names(taxonomy.data)[names(taxonomy.data) == "abundance.y"] <- 'abundance.mpln3'

# Omitting the kingdom row
taxonomy.data = taxonomy.data[-1,]

# Getting the phylum from the taxonomy column
taxonomy.data$phyla <- str_extract(taxonomy.data$taxonomy, 'p__\\w+[^|]')

# Getting the last taxonomic rank
taxonomy.data$lowest.rank.taxon <- str_extract(taxonomy.data$taxonomy, '\\w+$')

# Omit strains
#taxonomy.data <- taxonomy.data[!str_detect(taxonomy.data$lowest.rank.taxon, 't__'),]

# Set any NA abundance as 0
taxonomy.data[is.na(taxonomy.data$abundance.mpln2), "abundance.mpln2"] <- 0
taxonomy.data[is.na(taxonomy.data$abundance.mpln3), "abundance.mpln3"] <- 0


####################################################################################################
# Filtering only the species and removing the tag: s__
species.data <- taxonomy.data %>% filter(grepl('s__', lowest.rank.taxon)) 
species.data$lowest.rank.taxon <- gsub('s__', '', species.data$lowest.rank.taxon)
species.data$phyla <- str_remove(species.data$phyla, "p__")

# Find duplicated species' indices (provides only second duplicated entry)
dup.idx = duplicated(species.data$lowest.rank.taxon)

species.data.unique = species.data %>% select(abundance.mpln2, abundance.mpln3, phyla, lowest.rank.taxon)

# # Merging species in sample 'CSM5MCV1_P' 
# # Clostridium_bolteae, Escherichia_coli, Klebsiella_pneumoniae
# if(SAMPLE_ID == 'CSM5MCV1_P') {
#   mpln2.idx = c(7, 20, 21)
#   mpln3.idx = c(22, 31, 32)
#   species.data.unique$abundance.mpln2[mpln3.idx] <- species.data.unique$abundance.mpln2[mpln2.idx]
#   species.data.unique <- species.data.unique[-mpln2.idx,]
# }

  

# Removing the duplicated rows (by merging into one row)
# Go over all duplicated species
for(species in species.data$lowest.rank.taxon[dup.idx]) {
  # Indices (length(idx) == 2) of the rows which have said species.
  idx = which(species.data$lowest.rank.taxon == species)
  # If metaphlan2 abundace is 0 in row idx[1] then assign the abundance from the other row in idx[2]
  if(species.data$abundance.mpln2[idx[1]] == 0.0 & species.data$abundance.mpln3[idx[1]] != 0.0) {
    species.data.unique$abundance.mpln2[idx[1]] = species.data$abundance.mpln2[idx[2]]
  # If metaphlan3 abundace is 0 in row idx[1] then assign the abundance from the other row in idx[2]
  } else if(species.data$abundance.mpln2[idx[1]] != 0.0 & species.data$abundance.mpln3[idx[1]] == 0.0) {
    species.data.unique$abundance.mpln3[idx[1]] = species.data$abundance.mpln3[idx[2]]
  } else if(species.data$abundance.mpln2[idx[1]] == 0.0 & species.data$abundance.mpln3[idx[1]] == 0.0) {
    sprintf("Both abundances of %s are zero", species)
  } else {
    print("Both abundances of %s are non zero and that's a mistake", species)
  }
}
# remove row in idx[2]
species.data.unique = species.data.unique[!dup.idx,]
  

# Calculation of correlation coefficient
correlation_test <- function(df) {
  chosen_method = ""
  # Is the data normally distributed?
  
  # Sample size in Shapiro test must be between 3 and 5000
  if(dim(df)[1] <= 5000) {
    # I am including the zero abundance values 
    meta2.abund <- df$abundance.mpln2 # %>% replace(., . == 0, NA)
    meta3.abund <- df$abundance.mpln3 # %>% replace(., . == 0, NA)
    p.value2 <- shapiro.test(meta2.abund)$p.value
    p.value3 <- shapiro.test(meta3.abund)$p.value
    if(p.value2 > 0.05 & p.value3 > 0.05) {
      chosen_method = "pearson"
    } else {
      chosen_method = "spearman"
    }
  } else {
    chosen_method = "spearman"
  }
  cor.result = cor.test(df$abundance.mpln2, df$abundance.mpln3, 
                        method = chosen_method,
                        na.action=na.omit,
                        exact = FALSE)
  # if the p-value is smaller than 2e-16 then it results in 0, but that's the wrong output.
  if (cor.result$p.value != 0) {
    str <- sprintf("r = %.3f\np = %g\n%s", cor.result$estimate, 
                   cor.result$p.value, chosen_method)
  } else {
    str <- sprintf("r = %.3f\np < 2e-16\n%s", cor.result$estimate, 
                   chosen_method)
  }
  print(str)
  return(str)
}

str <- correlation_test(species.data.unique)


p <- species.data.unique %>% ggplot(aes(x = abundance.mpln2, y = abundance.mpln3)) +
  geom_smooth(aes(x = abundance.mpln2, y = abundance.mpln3), color = "steelblue3", method = lm, formula = y ~ x) +
  geom_abline(aes(intercept = 0, slope = 1), color = 'black', alpha = 0.5) + 
  geom_point(aes(color = phyla), size = 3, alpha = 0.7) +
  labs(x = "MetaPhlAn2", y = "MetaPhlAn3", 
       title = "Species Abundance", 
       subtitle = paste("Sample ID:", SAMPLE_ID)) +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5)) +
  annotate("text", x = 7, y = 40, label = str)

plot(p)

ggplotly(p)

###########################################################################################################################

# Some species are on different branches in metaphlan2 and in metaphlan 3. 
# Sometimes the branches are only different by a different spelling. 
# Here I am trying to manually find which species to merge. 

# Zero abundance in either metaphlan 2 or 3 means they are actually non existent in either one of them. 
# The zero abundance was added in order to deal with NA in the scatter plot.
zero.abund.df <- species.data %>% filter(abundance.mpln2 == 0 | abundance.mpln3 == 0)
duplicated.species.df <- species.data[dup.idx,]
duplicated.species <- duplicated.species.df$lowest.rank.taxon # names of duplicated species
only.dups.df <- species.data[species.data$lowest.rank.taxon %in% duplicated.species,] # df with only the duplicated species. 
# other species with zero abundance but not duplicated
zero.non.dups.df <- zero.abund.df[!(zero.abund.df$lowest.rank.taxon %in% duplicated.species),] 

# Species shared in both versions
joint.species <- species.data.unique %>% filter(abundance.mpln2 != 0 & abundance.mpln3 != 0)
colSums(joint.species[, c(1,2)])
