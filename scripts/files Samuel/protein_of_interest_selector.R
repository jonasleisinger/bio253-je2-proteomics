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
th_down_ids <- data.frame(proteome_TSB_SN$protein_Id[proteome_TSB_SN$diff > 0])
th_up_ids <- data.frame(proteome_TSB_SN$protein_Id[proteome_TSB_SN$diff < 0 ])
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

a <- merge(proteome_TSB_SN,
           uniprot_conversion,
           by.x= "protein_Id" ,
           by.y = "fasta_id",
           all=F)

proteome_TSB_SN_uniprot_diff <- a |> select(uniprot_entry_name, diff)


write_csv(proteome_TSB_SN_uniprot_diff, "proteome_TSB_SN_uniprot_diff.csv")

























##for schwyter et al

#load data 
data_schwyter <- readxl::read_excel("data/raw/schwyter/Original_Data_Schwyter_etal_2025_mSystems.xlsx", sheet = "Figure3(A)_S3")

#select relevant variables
data_schwyter_red <- data_schwyter |> select(UniProt,contrast, diff, FDR)
data_schwyter_red <- data_schwyter_red |> mutate(neg_log10_FDR = -log10(FDR))

data_schwyter_sig <- data_schwyter_red |> mutate(regulation = case_when(
  diff <= log2(1/2) & neg_log10_FDR >= 1 ~ "up",
  diff >= log2(2) & neg_log10_FDR >= 1~ "down",
  .default = "n.s."
))


#take out not significants
data_schwyter_sig <- data_schwyter_sig |> filter(regulation != "n.s.")



#select proteins of interst:
th_down_ids_schwyter <- data.frame(data_schwyter_sig$UniProt[data_schwyter_sig$diff > 0])
th_up_ids_schwyter <- data.frame(data_schwyter_sig$UniProt[data_schwyter_sig$diff < 0 ])
colnames(th_up_ids_schwyter) <- "ID"
colnames(th_down_ids_schwyter) <- "ID"


protein_of_interest_schwyter <- data.frame(rbind(th_down_ids_schwyter, th_up_ids_schwyter))


protein_up_scwyter <- data_schwyter_sig |> filter(UniProt %in% th_up_ids_schwyter$ID) |> filter(regulation == "up")
STRG_up_schwyter <- protein_up_scwyter |> select(UniProt,diff)
protein_down_schwyter <- data_schwyter_sig |> filter(UniProt %in% th_down_ids_schwyter$ID) |> filter(regulation == "down")
STRG_dwon_schwyter <- protein_down_schwyter |> select(UniProt, diff)

write_csv(STRG_up_schwyter, "STRG_up_schwyter.csv")
write_csv(STRG_dwon_schwyter, "STRG_down_schwyter.csv")




sink(file = "protein_new_down.txt")

df <- proteins_of_interest_uniprot_down
print(df, n = 154)

sink()


sink(file = "protein_new_up.txt")

df1 <- proteins_of_interest_uniprot_up

print(df1, n = 287)
sink()


sink(file = "protein_sch_down.txt")

df2 <- STRG_dwon_schwyter |> select(UniProt)
print(df2, n = 219)
sink()

sink(file = "protein_sch_up.txt")
df3 <- STRG_up_schwyter |> select(UniProt)
print(df3, n = 155)
sink()
