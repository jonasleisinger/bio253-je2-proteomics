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
  diff <= log2(1/2) & neg_log10_FDR >= 1 ~ "up",
  diff >= log2(2) & neg_log10_FDR >= 1~ "down",
  .default = "n.s."
))

#reduce to TSB vs SN
proteome_TSB_SN <- proteome |> filter(contrast == "TSB_vs_SN") |> filter(CON != TRUE) |> filter(regulation != "n.s.")


#threshold regulation

thresholds <- cutGn(proteome_TSB_SN$diff, m = length(proteome_TSB_SN$diff) / 5, what = "summary")
th_down_ids <- data.frame(proteome_TSB_SN$protein_Id[proteome_TSB_SN$diff > thresholds[5,1]])
th_up_ids <- data.frame(proteome_TSB_SN$protein_Id[proteome_TSB_SN$diff < thresholds[1,2]])
colnames(th_up_ids) <- "ID"
colnames(th_down_ids) <- "ID"


protein_of_interest <- data.frame(rbind(th_down_ids, th_up_ids))

#get uniprot labels of proteins
uniprot_conversion <- read_tsv(file = "data/reference/uniprot reference/JE2_protein_metadata.tsv")

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



