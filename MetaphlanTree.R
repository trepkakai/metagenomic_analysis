# This script is for plotting trees of select MetaPhlan3 organisms.

# Packages
library(ape)

# Location to nwk tree
fn_metaphlan_tree <- "/mnt/tank/labmainshare/qb3share/shared_resources/databases/Metaphlan3_Tree/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk"

# Read in file and get tips
metaphlan_tree <- ape::read.tree(fn_metaphlan_tree)
meta_tips <- metaphlan_tree$tip.label
df_tips <- data.frame(TipLabel = meta_tips)
df_tips$Species <- gsub(".*s__", "", df_tips$TipLabel)

# Select species of interest, meta_soi. Some examples are included below
meta_soi <- c("Eggerthella_lenta", "Adlercreutzia_equolifaciens", "Asaccharobacter_celatus", "Gordonibacter_pamelaeae", "Slackia_isoflavoniconvertens", "Collinsella_aerofaciens", "Enterorhabdus_caecimuris", "Roseburia_inulinivorans", "Ruminococcus_gnavus", "Roseburia_intestinalis")
meta_soi <- stringr::str_extract(meta_soi, "[^_]*_[^_]*")
df_tips_soi <- df_tips %>% filter(Species %in% meta_soi)

# Make sure all species present
length(df_tips_soi$Species) == length(meta_soi)

# Trim metaphlan tree and reformat labels
meta_tree <- keep.tip(metaphlan_tree, df_tips_soi$TipLabel)
meta_tree$tip.label <- gsub(".*s__", "", meta_tree$tip.label)
meta_tree$tip.label <- gsub(" ", "_", meta_tree$tip.label)

# Save plot
fn <- "savename.pdf"
pdf(file = fn, height = 4, width = 4, pointsize = 1/300)
plot(meta_tree, type = 'tidy', use.edge.length = TRUE)
dev.off()
