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
data_new <- data_new_raw %>% filter(contrast=="TSB_vs_SN_at_anc", CON=="FALSE")

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

#*IMPROVEMENT OF MEASUREMENT: NEWLY DISCOVERED PROTEINS IN OUR DATASET*

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
  labs(x = "Average abundance per protein", y = "Count") +
  theme_minimal()

p3


#*DIFF CORRELATION BETWEEN DATASETS*
merged_corr <- merged %>%
  filter(origin=="both")%>%
  mutate(sig = if_else(FDR_schwyter < 0.01 & FDR_new < 0.01,
                       "FDR < 0.01 in both", "otherwise"))

# correlation value
r  <- cor(merged_corr$diff_schwyter, merged_corr$diff_new, use = "complete.obs", method = "pearson")
rho <- cor(merged_corr$diff_schwyter, merged_corr$diff_new, use = "complete.obs", method = "spearman")

p4 <- ggplot(merged_corr, aes(x = diff_schwyter, y = diff_new)) +
  geom_point(alpha = 0.6,  colour="darkblue") +
  geom_smooth(method = "lm", se = TRUE,  colour="darkblue") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d",
                           r, rho, nrow(merged))) +
  labs(x = "diff schwyter", y = "diff new") +
  theme_bw()

p5 <- ggplot(merged_corr, aes(x = diff_schwyter, y = diff_new, colour=sig)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("FDR < 0.01 in both" = "magenta", "otherwise" = "darkblue"),
                     name = NULL) +
  # annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
  #          label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d",
  #                          r, rho, nrow(merged))) +
  labs(x = "diff schwyter", y = "diff new") +
  theme_bw()

p4

#correlated and anticorrelated points
anticorr_cut=-0.5
corr_cut=1

merged_corr <- merged_corr %>%
  mutate(contrib = (diff_schwyter - mean(diff_schwyter, na.rm = TRUE)) *
           (diff_new - mean(diff_new, na.rm = TRUE)))

anticorr   <- merged_corr %>% filter(contrib < anticorr_cut)
corr   <- merged_corr %>% filter(contrib > corr_cut)


crup <- corr %>% 
  filter(diff_new>0) %>% 
  select(uniprot_acc)
atcrup <- anticorr %>% 
  filter(diff_new>0) %>% 
  select(uniprot_acc)
anticorr_str <- anticorr %>%
  select(uniprot_acc)
corr_str <- corr %>%
  select(uniprot_acc)


write_csv(anticorr_str, "data/outputs/anticorr.csv")
write_csv(corr_str, "data/outputs/corr.csv")

write_csv(crup, "data/outputs/corrup.csv")
write_csv(atcrup, "data/outputs/anticorrup.csv")


p6 <- ggplot(merged_corr, aes(x = diff_schwyter, y = diff_new, color = contrib < anticorr_cut)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "magenta"),
                     labels = c("FALSE" = "Other", "TRUE" = "anticorrelated"),
                     name = NULL) +
  labs(x = "diff schwyter", y = "diff new") +
  theme_bw()

p7 <- ggplot(merged_corr, aes(x = diff_schwyter, y = diff_new, color = contrib > corr_cut)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "magenta"),
                     labels = c("FALSE" = "Other", "TRUE" = "correlated"),
                     name = NULL) +
  labs(x = "diff schwyter", y = "diff new") +
  theme_bw()


p7

#do only new proteins also show a big/small diff

p8<-ggplot(merged, aes(x = diff_new, y = -log10(FDR_new), color = origin)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "Difference (log2 fold change)", y = "-log10(FDR)") +
  theme_minimal()

p8
