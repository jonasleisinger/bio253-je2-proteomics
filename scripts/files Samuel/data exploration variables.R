#import used packages
library(tidyverse)
library(GGally)
library(stringr)



#load data 
data_proteome <- readxl::read_excel("data/raw/total/DE_total.xlsx", sheet = "diff_exp_analysis")


#make reduced data set for pairs plot
data_proteome_red_pairs <- data_proteome |> select(nrPeptides, protein_length, nr_tryptic_peptides,CON,modelName,estimate_type,contrast,diff,avgAbd,std.error.unmoderated,df.unmoderated,std.error,statistic,df,p.value,sigma,FDR)

#make all non numerics into factors
for(i in 1:ncol(data_proteome_red_pairs)){
  if(!is.numeric(data_proteome_red_pairs[[i]])){
    data_proteome_red_pairs[[i]] <- as.factor(data_proteome_red_pairs[[i]])
  }
}

#make pairs plot
pairs(head(data_proteome_red_pairs, 800), lower.panel = NULL)



#check what values number of peptide can take
data_proteome_red_pairs |> group_by(nrPeptides) |> summarise(n = n())
#check distribution of nrPeptides
ggplot(data_proteome_red_pairs, aes(x = nrPeptides)) +
  geom_histogram(bins = 200)
#nrPeptides is skewed towards left side, why?

#add marker for int.contrast
data_proteome <- data_proteome |> mutate(interaction = case_when(str_detect(contrast, "interaction")== TRUE ~ 1 ,
                                                                 .default = 0))

#check what values protein_length can take
data_proteome_red_pairs |> group_by(protein_length) |> summarise(n = n())
#check distribution
ggplot(data_proteome_red_pairs, aes(x = protein_length)) +
  geom_histogram(bins = 300)

#check correlation protein_length, nrPeptides
ggplot(data_proteome_red_pairs, aes(x = protein_length, y = nrPeptides)) +
  geom_point()



#contrast vs.number peptides
ggplot(data_proteome_red_pairs, aes(x = contrast, y = nrPeptides)) +
  geom_point()
#diff vs.number of peptides
ggplot(data_proteome_red_pairs, aes(x = diff, y = nrPeptides)) +
    geom_s
ggplot(data_proteome_red_pairs, aes(x = nrPeptides, y = diff)) +
  geom_point()



with(data_proteome_red_pairs, scatter.smooth(x = diff, y = nrPeptides))

with(data_proteome_red_pairs, smoothScatter(x = diff, y = nrPeptides))


#diff vs. contrast
ggplot(data_proteome_red_pairs, aes(x= diff, y = contrast)) +
  geom_point()

#avgabund vs. nrPept  
ggplot(data_proteome, aes(x = avgAbd, y = nrPeptides, col =interaction)) +
  geom_point()


#std.error vs.nrPept
ggplot(data_proteome, aes(x = nrPeptides, y = std.error.unmoderated)) +
  geom_point()
ggplot(data_proteome, aes(x = nrPeptides, y = std.error)) +
  geom_point()



#diff vs. std.error
ggplot(data_proteome, aes(x = diff, y = std.error)) +
  geom_point()
ggplot(data_proteome, aes(x = avgAbd, y = std.error, col =interaction)) +
  geom_point()
ggplot(data_proteome, aes(x = df.unmoderated, y = std.error, col =interaction)) +
  geom_point()

#statistic vs. nrPeptides
ggplot(data_proteome, aes(x = nrPeptides, y = statistic, col =interaction)) +
  geom_point()


ggplot(data_proteome, aes(x = diff, y = statistic, col =interaction)) +
  geom_point()



#diff vs FDR
ggplot(data_proteome, aes(x = diff, y = FDR)) +
  geom_point()

