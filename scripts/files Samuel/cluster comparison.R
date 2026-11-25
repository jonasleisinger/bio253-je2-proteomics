##comparing leiden and mcl algorithm in string db

#load packages
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)



#load data 
clusters_mcl <- read_tsv("data/string clusters/string_MCL_cluster_descriptions.tsv")
clusters_leiden <- read_tsv("data/string clusters/string_leiden_cluster_descriptions.tsv")

#make function to see, what each protein is assigned to

#define function
protein_name_split <- function(cluster_db){
  #make data frame
  describtion_by_protein <- data.frame(primary = 0, secondary = 0, tertiary = 0, protein = 0, cluster_Id = 0)
  
  #split the protein names
  a <- cluster_db["protein names"]
  a <- strsplit(a[[1]], ",")
  
  #go through all of the clusterdb
  for(i in 1:dim(cluster_db)[1]){
    primary <- cluster_db["primary description"][[1]][i]
    secondary <- cluster_db["secondary description"][[1]][i]
    tertiary <- cluster_db["tertiary description"][[1]][i]
    cluster_Id <- cluster_db["cluster ID"][[1]][i]
    for(j in 1:length(a[[i]])){
      df <- data.frame(primary = primary, secondary = secondary, tertiary = tertiary, protein = a[[i]][j], cluster_Id = cluster_Id)
      describtion_by_protein <- rbind(describtion_by_protein, df)
      
    }
    
    
  }
  describtion_by_protein <- describtion_by_protein |> filter(protein != 0)
  return(describtion_by_protein)
}


mcl_protein_names <- protein_name_split(clusters_mcl)
leiden_protein_names <- protein_name_split(clusters_leiden)

total_protein_cl <- merge(mcl_protein_names, leiden_protein_names, by = "protein", suffixes = c("_mcl", "_leiden") ) 

#make a list with all clusters in mcl
cl_mcl <- unique(total_protein_cl$cluster_Id_mcl)
cl_leiden <- unique(total_protein_cl$cluster_Id_leiden)


ploting <- list()


for(i in 1:length(cl_mcl)){
  
  a <- data.frame(table(total_protein_cl[total_protein_cl$cluster_Id_leiden == cl_leiden[i], 5]))
  colnames(a) <- c("protein", "Freq")
  a$Freq <- a$Freq / sum(a$Freq)
  p <- ggplot(a, aes(x = "", y = Freq, fill = protein)) + 
    geom_col() +
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = paste("Cluster", cl_leiden[i])) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      
      # Legend formatting
      legend.position = "bottom",
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.2, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.spacing.x = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    ) +
    labs(fill = NULL) +
    theme(
      plot.margin = margin(
        t = 5,
        r = 20,
        b = 5,
        l = 20,
        unit = "pt"
      )
    )
  
    
  ploting[[i]] <- p
}
patchwork::wrap_plots(ploting, ncols = 7)
gridExtra::grid.arrange(grobs = ploting)


