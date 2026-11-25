library(clusterProfiler)
library(dplyr)

# 1. Your regulated UniProt accessions (from column 1)
regulated_proteins <- c(
  "Q2G2H5", "Q2G2H4", "Q2G276", "Q2G275", "Q2G274", "Q2G2Q0",
  "Q2G2P8", "Q2G2P7", "Q2FXZ4", "P95689", "Q2G2P6", "Q2G2T5",
  "Q2G2T4", "Q2G2T2", "Q2G2T6", "Q2G2T3", "Q2G1S4", "Q2G1S3",
  "Q2G2U6", "Q2G2U4", "Q2G2U3", "Q2G2U2", "Q2G2U5", "Q2G253", "Q2G252"
)

# 2. Run KEGG pathway enrichment (organism = "sao" for NCTC 8325)
# Tip: pass your full quantified proteome as `universe` for accurate proteomics background!
kk <- enrichKEGG(
  gene          = regulated_proteins,
  organism      = "sao",
  keyType       = "uniprot",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# View results
res_df <- as.data.frame(kk)
head(res_df[, c("ID", "Description", "GeneRatio", "p.adjust", "geneID")])

# Quick visualization
dotplot(kk, showCategory = 15)





library(tidyverse)

# 1. Your regulated UniProt accessions (from column 1)
regulated_proteins <- c(
  "Q2G2H5", "Q2G2H4", "Q2G276", "Q2G275", "Q2G274", "Q2G2Q0",
  "Q2G2P8", "Q2G2P7", "Q2FXZ4", "P95689", "Q2G2P6", "Q2G2T5",
  "Q2G2T4", "Q2G2T2", "Q2G2T6", "Q2G2T3", "Q2G1S4", "Q2G1S3",
  "Q2G2U6", "Q2G2U4", "Q2G2U3", "Q2G2U2", "Q2G2U5", "Q2G253", "Q2G252"
)


# 1. Fetch official UniProt <-> KEGG (sao) mapping and Pathway definitions
conv_df <- read_tsv("https://rest.kegg.jp/conv/sao/uniprot", col_names = c("uniprot", "kegg_gene")) %>%
  mutate(uniprot = str_remove(uniprot, "^up:"))

link_df <- read_tsv("https://rest.kegg.jp/link/pathway/sao", col_names = c("kegg_gene", "pathway_id"))

path_names <- read_tsv("https://rest.kegg.jp/list/pathway/sao", col_names = c("pathway_id", "pathway_name")) %>%
  mutate(pathway_name = str_remove(pathway_name, " - Staphylococcus aureus.*"))

# 2. Join with your regulated list
my_proteins <- tibble(uniprot = regulated_proteins)

classified_proteins <- my_proteins %>%
  left_join(conv_df, by = "uniprot") %>%
  left_join(link_df, by = "kegg_gene") %>%
  left_join(path_names, by = "pathway_id")

# View proteins grouped by pathway
pathway_summary <- classified_proteins %>%
  filter(!is.na(pathway_id)) %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(
    count = n(),
    proteins = paste(uniprot, collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

