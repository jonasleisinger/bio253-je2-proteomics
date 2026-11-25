library(readxl)
library(dplyr)
library(ggplot2)

data_phospho <- read_excel("data/raw/phospho_STY/DE_phospho_STY.xlsx", sheet="diff_exp_analysis")
data_total <- read_excel("data/raw/total/DE_total.xlsx", sheet="diff_exp_analysis")
mapping <- read.delim("data/raw/resources/mapping_uniprot/JE2_protein_metadata.tsv")

phospho <- data_phospho %>%
  left_join(select(mapping, protein_key, uniprot_acc),
            by = c("proteinname" = "protein_key"))
#phospho <-data_phospho %>% select(proteinname, contrast, diff)
phospho_mean <- phospho %>%
  group_by(proteinname, contrast) %>%
  summarise(diff = mean(diff, na.rm = TRUE), .groups = "drop")

# phospho_str <- phospho_mean %>%
#   filter(contrast=="TSB_vs_SN_at_anc", diff>1) %>%
#   select(uniprot_acc)
# 
# write_csv(phospho_str, "data/outputs/phospho/phoshpoup.csv")


# df1_mean and df2_mean = your averaged data frames
merged <- inner_join(phospho_mean, data_total,
                     by = c("proteinname", "contrast"),
                     suffix = c("_phospho", "_total"))

# correlation value
r  <- cor(merged$diff_phospho, merged$diff_total, use = "complete.obs", method = "pearson")
rho <- cor(merged$diff_phospho, merged$diff_total, use = "complete.obs", method = "spearman")

ggplot(merged, aes(x = diff_phospho, y = diff_total)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "red") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Pearson r = %.2f\nSpearman rho = %.2f\nn = %d",
                           r, rho, nrow(merged))) +
  labs(x = "diff (phospho)", y = "diff (total)") +
  theme_bw()
