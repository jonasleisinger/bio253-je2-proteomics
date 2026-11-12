#import needed libraries
library(tidyverse)


#make volcano plot for our data



#get our data in

protein_abundances <- readxl::read_excel("data/raw/total/DE_total.xlsx", sheet = "diff_exp_analysis")

#get relevant data seperated

protein_abundances_clean <- protein_abundances |> select(protein_Id
                                                        ,contrast,
                                                        p.value,
                                                        diff,
                                                        FDR, proteinname,
                                                        avgAbd)


#add axis of vulcan plot to dataset (-log10_FDR, log2_FC)
protein_abundances <- protein_abundances |> mutate(neg__log10_FDR = -log10(FDR),
                                                               log2_FC = log2(diff))




#take only data of contrast TSB vs SN:
protein_abundances_TSB_SN <- protein_abundances |> filter(contrast == "TSB_vs_SN")

#add classification of cutoff 
protein_abundances_TSB_SN <- protein_abundances_TSB_SN |> mutate(regulation = case_when(
  diff <= log2(1/2) & neg__log10_FDR >= 1 ~ "up",
  diff >= log2(2) & neg__log10_FDR >= 1~ "down",
  .default = "n.s."
))


#make color labels
manual_colors = c(
  "down" = "tan",
  "up" = "turquoise",
  "n.s" = "lightgrey"
)
#make vulcan plot for this case alone:
a <- ggplot(data = protein_abundances_TSB_SN, aes(x = diff, y = neg__log10_FDR, col = regulation)) +
  geom_point() +
  geom_hline(yintercept = 1, col = "gray", linetype = "dashed") +
  geom_vline(xintercept = log2(1/2), col = "gray", linetype = "dashed") +
  geom_vline(xintercept = log2(2), col = "gray", linetype = "dashed") +
  ggtitle("Change in protein regulation TSB vs. SN") +
  xlab("log2_FC (TSB / PASN)") +
  ylab("-log10_FDR") +
  scale_color_manual(values = manual_colors)


plotly::ggplotly(a)




protein_abundances <- protein_abundances |> mutate(regulation = case_when(
  diff <= log2(1/2) & neg__log10_FDR >= 1 ~ "down",
  diff >= log2(2) & neg__log10_FDR >= 1~ "up",
  .default = "n.s."
))


#make vulcan plot for each set, with contrast as 
ggplot(data = protein_abundances, aes(x = diff, y = neg__log10_FDR, col = regulation)) +
  geom_point() +
  facet_wrap(vars(contrast)) +
  geom_hline(yintercept = 1, col = "gray", linetype = "dashed") +
  geom_vline(xintercept = log2(1/2), col = "gray", linetype = "dashed") +
  geom_vline(xintercept = log2(2), col = "gray", linetype = "dashed") +
  scale_color_manual(values = manual_colors)
