library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

data_schwyter_raw <- read_excel("data/schwyter_etal/Original_Data_Schwyter_etal_2025_mSystems.xlsx", sheet="Figure3(A)_S3")
data_new_raw <- read_excel("data/raw/total/DE_total.xlsx", sheet="diff_exp_analysis")

data_schwyter <- data_schwyter_raw %>% rename(proteinname = Protein_ID)
data_new <- data_new_raw %>% filter(contrast=="TSB_vs_SN_at_anc")


data_schwyter <- data_schwyter %>% mutate(in_schwyter = TRUE)
data_new <- data_new %>% mutate(in_new = TRUE)

merged <- full_join(data_schwyter, data_new, by = "proteinname", suffix = c("_schwyter", "_new")) %>%
  mutate(
    in_schwyter = coalesce(in_schwyter, FALSE),
    in_new = coalesce(in_new, FALSE),
    herkunft = case_when(
      in_schwyter & in_new ~ "beide",
      in_schwyter          ~ "nur schwyter",
      in_new          ~ "nur new"
    )
  ) %>%
  select(-in_schwyter, -in_new)



p1<- ggplot(merged, aes(x = herkunft, y = avgAbd_new)) +
  geom_boxplot() +
  theme_bw()

p2<- ggplot(merged, aes(x = herkunft, y = avgAbd_schwyter)) +
  geom_boxplot() +
  theme_bw()

p1+p2
