library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

#data import
data_raw <- read_excel("data/raw/total/DE_total.xlsx", sheet="diff_exp_analysis")
mapping <- read.delim("data/raw/resources/mapping_uniprot/JE2_protein_metadata.tsv")

data <- data_raw %>% filter(contrast=="TSB_vs_SN_at_anc", CON=="FALSE") %>%
  mutate(diff=-diff)

data <- data %>%
  left_join(select(mapping, protein_key, uniprot_acc),
            by = c("proteinname" = "protein_key"))

#export proteinnames for string
tot <- data %>%
  select(uniprot_acc)

newup <- data %>%
  filter(diff>1.5) %>%
  select(uniprot_acc)

newdwn <- data %>%
  filter(diff< -1.5) %>%
  select(uniprot_acc)

write_csv(tot, "data/outputs/only_new/newbackground.csv")
write_csv(newup, "data/outputs/only_new/newup.csv")
write_csv(newdwn, "data/outputs/only_new/newdwn.csv")
