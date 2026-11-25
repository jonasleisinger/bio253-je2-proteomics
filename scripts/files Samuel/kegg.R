#import packages
library(ggplot2)
library(readr)
library(dplyr)


#import KEGG file
KEGG_new <- read_tsv("data/reference/KEGG analysis/enrichment.KEGG.new.tsv")
KEGG_new_down <- read_tsv("data/reference/KEGG analysis/enrichment.KEGG.down.tsv")
KEGG_new_up <- read_tsv("data/reference/KEGG analysis/enrichment.KEGG.new.tsv")


KEGG_sch_down <- read_tsv("data/reference/KEGG analysis/enrichment.KEGG.schwyter.down.tsv")
KEGG_sch_up <- read_tsv("data/reference/KEGG analysis/enrichment.KEGG.schwyter.up.tsv")


KEGG_pathways <- unique(c(KEGG_new$`term description`, KEGG_new_down$`term description`, KEGG_new_up$`term description`))
limits <- c(0:(length(KEGG_pathways)-1))
df_kegg <- data.frame(KEGG_pathways, limits)
#mark the up and down in the total KEGG data

KEGG_new_down <- KEGG_new_down |> mutate(reg = "down")
KEGG_new_up <- KEGG_new_up |> mutate(reg = "up")
KEGG_total <- merge(KEGG_new_down, KEGG_new_up, by = "term description", all = TRUE)


KEGG_sch_down <- KEGG_sch_down |> mutate(reg = "down")
KEGG_sch_up <- KEGG_sch_up |> mutate(reg = "up")
KEGG_total_sch <- merge(KEGG_sch_down, KEGG_sch_up, by = "term description", all = TRUE)

#map term describtion to gene count
ggplot(KEGG_new_down, aes(x = `observed gene count`, y = `term description`)) +
  geom_col()

ggplot(KEGG_total) +
  geom_col(
    aes(x = `observed gene count.x`, y = `term description`, fill = reg.x),
    width = 0.3,
    position = position_nudge(y = -0.15)
  ) +
  geom_col(
    aes(x = `observed gene count.y`, y = `term description`, fill = reg.y),
    width = 0.3,
    position = position_nudge(y = 0.15)
  )



ggplot(KEGG_total_sch) +
  geom_col(
    aes(x = `observed gene count.x`, y = `term description`, fill = reg.x),
    width = 0.3,
    position = position_nudge(y = -0.15)
  ) +
  geom_col(
    aes(x = `observed gene count.y`, y = `term description`, fill = reg.y),
    width = 0.3,
    position = position_nudge(y = 0.15))

