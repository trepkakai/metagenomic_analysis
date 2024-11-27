
#This is a document to run PRIOR to any data analysis of the MGS data. Specifically, this code will: 
# 1. Summarise all taxa abundances into tables
# 2. Summarise all normalizations into a single table
# 3. Convert UniRef90 gene family abundances into units of RPKG and output tables

# Libraries and directories
library(dplyr)
library(data.table)
pattern = "NE_" # grepl pattern that is present in all samples of interest
seqdir <- "/mnt/tank/labmainshare/qb3share/ktrepka/Sequencing_NE_Winter24/" # directory with outputs from MGS pipeline
tabledir <- "/mnt/tank/labmainshare/qb3share/ktrepka/NetherlandsStudy/Tables/" # directory to save outputs from this script

# Summarise taxa abundances
taxtabledir <- paste0(tabledir, "TaxaTables/")
suppressWarnings(dir.create(taxtabledir))

files <- list.files(paste0(seqdir, "metaphlan/"))
files <- files[grepl(pattern, files)]
Sample <- gsub("_S.*", "", files)
files_full <- paste0(seqdir, "metaphlan/", files)

taxa_df <- NULL
for (k in 1:length(files_full)){
  file <- files_full[k]
  samp <- fread(file, sep = "\t") %>% as.data.frame()
  samp <- samp %>% select(`#clade_name`, relative_abundance)
  colnames(samp) <- c("Taxa", Sample[k])
  if (length(taxa_df) == 0){
    taxa_df <- samp
  } else {
    taxa_df <- full_join(taxa_df, samp, by = "Taxa")
  }
}
taxa_df[is.na(taxa_df)] <- 0 # replace NA with 0

taxa_df$Level <- case_when(grepl("t__", taxa_df$Taxa) ~ "Strain",
                           grepl("s__", taxa_df$Taxa) ~ "Species",
                           grepl("g__", taxa_df$Taxa) ~ "Genus",
                           grepl("f__", taxa_df$Taxa) ~ "Family",
                           grepl("o__", taxa_df$Taxa) ~ "Order",
                           grepl("c__", taxa_df$Taxa) ~ "Class",
                           grepl("p__", taxa_df$Taxa) ~ "Phylum",
                           grepl("k__", taxa_df$Taxa) ~ "Kingdom",
                           TRUE ~ "Unclassified")
for (level in unique(taxa_df$Level)){ # save each taxa level
  taxa_export <- taxa_df %>% filter(Level == level)
  taxa_export <- taxa_export %>% select(-Level)
  fn <- paste0(taxtabledir, "Taxa", level, ".csv")
  write.csv(file = fn, taxa_export, row.names = FALSE)
}

# Summarise sequencing statistics and normalizations
normdir <- paste0(tabledir, "Normalizations/")
suppressWarnings(dir.create(normdir))

files <- list.files(paste0(seqdir, "normalizations/"))
files <- files[grepl(pattern, files)]
files_full <- paste0(seqdir, "normalizations/", files)

norm_df = NULL
for (k in 1:length(files_full)){
  file <- files_full[k]
  norm <- fread(file, sep = "\t") %>% as.data.frame()
  if (length(norm_df) == 0){
    norm_df <- norm
  } else {
    norm_df <- rbind(norm_df, norm)
  }
}
norm_df$Sample <- gsub("_S.*", "", norm_df$SampleID)
fn <- paste0(normdir, "Normalizations.csv")
write.csv(file = fn, norm_df, row.names = FALSE)

# Summarise gene family abundance in units of RPKG (i.e. divide counts by the genome equivalent)
genedir <- paste0(tabledir, "GeneTables/")
suppressWarnings(dir.create(genedir))

files <- list.files(paste0(seqdir, "humann3/"))
files <- files[grepl(pattern, files)]
files <- files[grepl("_genefamilies_unstratified", files)]
Sample <- gsub("_S.*", "", files)
files_full <- paste0(seqdir, "humann3/", files)

gene_df <- NULL
for (k in 1:length(files_full)){
  file <- files_full[k]
  samp <- fread(file, sep = "\t") %>% as.data.frame()
  colnames(samp) <- c("GeneFamily", Sample[k])
  geq <- norm_df %>% filter(Sample == Sample[k]) %>% select(GenomeEquivalents) %>% pull()
  samp[,Sample[k]] <- samp[,Sample[k]]/geq
  if (length(gene_df) == 0){
    gene_df <- samp
  } else {
    gene_df <- full_join(gene_df, samp, by = "GeneFamily")
  }  
}
gene_df[is.na(gene_df)] <- 0 # replace NA with 0

fn <- paste0(genedir, "GeneFamiliesUnstratified.csv")
write.csv(file = fn, gene_df, row.names = FALSE, quote = FALSE)

# Convert to Kegg Orthologs (KO)
fn <- "/mnt/tank/labmainshare/qb3share/ktrepka/GOStudy/SequencingAnalysis/Tables/humann/Resources/map_ko_uniref90_keyvalue.csv"
dic <- read.csv(fn)
colnames(dic) <- c("KO", "GeneFamily")

merged <- inner_join(gene_df, dic, by = "GeneFamily") 
merged <- merged %>% select(-"GeneFamily")
y_columns <- colnames(merged %>% select(-KO))
summed <- merged %>%
  group_by(KO) %>%
  summarize(across(all_of(y_columns), \(x) sum(x, na.rm = TRUE)))

fn <- paste0(genedir, "KO.csv")
write.csv(summed, file = fn, quote = FALSE)

# Summarise pathway abundance in units of RPKG (i.e. divide counts by the genome equivalent)
files <- list.files(paste0(seqdir, "humann3/"))
files <- files[grepl(pattern, files)]
files <- files[grepl("_pathabundance_unstratified", files)]
Sample <- gsub("_S.*", "", files)
files_full <- paste0(seqdir, "humann3/", files)

path_df <- NULL
for (k in 1:length(files_full)){
  file <- files_full[k]
  samp <- fread(file, sep = "\t") %>% as.data.frame()
  colnames(samp) <- c("Pathway", Sample[k])
  samp$Pathway <- gsub(",", "", samp$Pathway)
  geq <- norm_df %>% filter(Sample == Sample[k]) %>% select(GenomeEquivalents) %>% pull()
  samp[,Sample[k]] <- samp[,Sample[k]]/geq
  if (length(path_df) == 0){
    path_df <- samp
  } else {
    path_df <- full_join(path_df, samp, by = "Pathway")
  }  
}
path_df[is.na(path_df)] <- 0 # replace NA with 0

fn <- paste0(genedir, "PathAbundanceUnstratified.csv")
write.csv(file = fn, path_df, row.names = FALSE, quote = FALSE)

# Summarise pathway abundance in units of RPKG (i.e. divide counts by the genome equivalent)
files <- list.files(paste0(seqdir, "humann3/"))
files <- files[grepl(pattern, files)]
files <- files[grepl("_pathabundance_stratified", files)]
Sample <- gsub("_S.*", "", files)
files_full <- paste0(seqdir, "humann3/", files)

path_df <- NULL
for (k in 1:length(files_full)){
  file <- files_full[k]
  samp <- fread(file, sep = "\t") %>% as.data.frame()
  colnames(samp) <- c("Pathway", Sample[k])
  samp$Pathway <- gsub(",", "", samp$Pathway)
  geq <- norm_df %>% filter(Sample == Sample[k]) %>% select(GenomeEquivalents) %>% pull()
  samp[,Sample[k]] <- samp[,Sample[k]]/geq
  if (length(path_df) == 0){
    path_df <- samp
  } else {
    path_df <- full_join(path_df, samp, by = "Pathway")
  }  
}
path_df[is.na(path_df)] <- 0 # replace NA with 0

fn <- paste0(genedir, "PathAbundanceStratified.csv")
write.csv(file = fn, path_df, row.names = FALSE, quote = FALSE)
