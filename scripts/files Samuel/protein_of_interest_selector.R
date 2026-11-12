##make a list of proteins of interest
library(tidyverse)
library(readr)
library(Hmisc)


#load data
proteome <- readxl::read_excel("data/raw/total/DE_total.xlsx", sheet = "diff_exp_analysis")


#add logs of data
proteome <- proteome |> mutate(neg_log10_FDR = -log10(FDR))

#add identifier of up/down
proteome <- proteome |> mutate(regulation = case_when(
  diff <= log2(1/2) & neg_log10_FDR >= 1 ~ "down",
  diff >= log2(2) & neg_log10_FDR >= 1~ "up",
  .default = "n.s."
))

#reduce to TSB vs SN
proteome_TSB_SN <- proteome |> filter(contrast == "TSB_vs_SN") |> filter(CON != TRUE) |> filter(regulation != "n.s.")


#threshold regulation
 
th_up_table <- proteome_TSB_SN |> filter(regulation == "up") |> arrange(desc(diff))
th_up_ids <- data.frame(th_up_table$protein_Id[1:30])
colnames(th_up_ids) <- "ID"


th_down_table <- proteome_TSB_SN |> filter(regulation == "down") |> arrange(desc(diff))
th_down_ids <- data.frame(th_down_table$protein_Id[1:30])
colnames(th_down_ids) <- "ID"

protein_of_interest <- data.frame(rbind(th_down_ids, th_up_ids))

#get uniprot labels of proteins
uniprot_conversion <- read_tsv(file = "uniprot reference/JE2_protein_metadata.tsv")

proteins_of_interest_uniprot <- uniprot_conversion |> filter(fasta_id %in% protein_of_interest$ID) |> select(uniprot_entry_name)


write_csv(proteins_of_interest_uniprot, "proteins_of_interest.csv")







#isolate for down / up
proteins_of_interest_uniprot_up <- uniprot_conversion |> filter(fasta_id %in% th_up_ids$ID) |> select(uniprot_entry_name)

write_csv(proteins_of_interest_uniprot_up, "proteins_of_interest_up.csv")



proteins_of_interest_uniprot_down <- uniprot_conversion |> filter(fasta_id %in% th_down_ids$ID) |> select(uniprot_entry_name)

write_csv(proteins_of_interest_uniprot_down, "proteins_of_interest_down.csv")



#take all up/downregulated proteins 

proteome_TSB_SN_uniprot <- uniprot_conversion |> filter(fasta_id %in% proteome_TSB_SN$protein_Id) |> select(uniprot_entry_name)
write_csv(proteome_TSB_SN_uniprot, "proteome_TSB_SN_uniprot.csv")



