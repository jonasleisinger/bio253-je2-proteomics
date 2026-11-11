library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

#data import
data_schwyter_raw <- read_excel("data/schwyter_etal/Original_Data_Schwyter_etal_2025_mSystems.xlsx", sheet="Figure3(A)_S3")
data_new_raw <- read_excel("data/raw/total/DE_total.xlsx", sheet="diff_exp_analysis")
mapping <- read.delim("~/bio253-je2-proteomics/data/raw/resources/mapping_uniprot/JE2_protein_metadata.tsv")

#data tidying
data_schwyter <- data_schwyter_raw %>% rename(proteinname = Protein_ID)
data_new <- data_new_raw %>% filter(contrast=="TSB_vs_SN_at_anc")

data_schwyter <- data_schwyter %>% mutate(in_schwyter = TRUE)
data_new <- data_new %>% mutate(in_new = TRUE)

merged <- full_join(data_schwyter, data_new, by = "proteinname", suffix = c("_schwyter", "_new")) %>%
  mutate(
    in_schwyter = coalesce(in_schwyter, FALSE),
    in_new = coalesce(in_new, FALSE),
    origin = case_when(
      in_schwyter & in_new ~ "both",
      in_schwyter          ~ "only schwyter",
      in_new          ~ "only new"
    )
  ) %>%
  select(-in_schwyter, -in_new)

merged <- merged %>%
  left_join(select(mapping, protein_key, uniprot_acc),
            by = c("proteinname" = "protein_key"))

#export datasets for string analysis
str_new <- merged %>%
  filter(origin=="only new") %>%
  select(uniprot_acc)

str_all <- merged %>%
  select(uniprot_acc)

write_csv(str_new, "data/outputs/newprot.csv")
write_csv(str_all, "data/outputs/allprot.csv")

#plot differences between datasets

# p1<- ggplot(merged, aes(x = origin, y = avgAbd_new)) +
#   geom_boxplot() +
#   theme_bw()
# 
# p2<- ggplot(merged, aes(x = origin, y = avgAbd_schwyter)) +
#   geom_boxplot() +
#   theme_bw()
# 
# p1+p2


p3 <- ggplot(merged, aes(x = avgAbd_new, color=origin)) +
  geom_density(linewidth = 1) +
  labs(title = "Histogram of variable", x = "variable", y = "Count") +
  theme_minimal()

p3
