library(readxl)
library(dplyr)
library(ggplot2)

data_phospho <- read_excel("data/raw/phospho_STY/DE_phospho_STY.xlsx", sheet="diff_exp_analysis")
data_total <- read_excel("data/raw/total/DE_total.xlsx", sheet="diff_exp_analysis")

df<-data_phospho %>% select(proteinname, contrast, diff)
phospho_mean <- df %>%
  group_by(proteinname, contrast) %>%
  summarise(diff = mean(diff, na.rm = TRUE), .groups = "drop")


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
