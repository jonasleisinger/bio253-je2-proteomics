library(readxl)
library(dplyr)
library(ggplot2)

data_schwyter_raw <- read_excel("data/schwyter_etal/Original_Data_Schwyter_etal_2025_mSystems.xlsx", sheet="Figure3(A)_S3")
data_new_raw <- read_excel("data/raw/total/DE_total.xlsx", sheet="diff_exp_analysis")

data_schwyter <- data_schwyter_raw %>% rename(proteinname = Protein_ID)
data_new <- data_new_raw %>% filter(contrast=="TSB_vs_SN_at_anc", CON=="FALSE")


# df1_mean and df2_mean = your averaged data frames
merged <- inner_join(data_schwyter, data_new,
                     by = c("proteinname"),
                     suffix = c("_schwyter", "_new"))

merged <- merged %>%
  mutate(sig = if_else(FDR_schwyter < 0.01 & FDR_new < 0.01,
                       "FDR < 0.01 in both", "otherwise"))

# correlation value
r  <- cor(merged$diff_schwyter, merged$diff_new, use = "complete.obs", method = "pearson")
rho <- cor(merged$diff_schwyter, merged$diff_new, use = "complete.obs", method = "spearman")

ggplot(merged, aes(x = diff_schwyter, y = diff_new, colour=sig)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  # annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
  #          label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d",
  #                          r, rho, nrow(merged))) +
  labs(x = "diff schwyter", y = "diff new") +
  theme_bw()
