library(bayesplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggExtra)
library(factoextra)
library(patchwork)
library(posterior)
library(pracma)
library(ggpmthemes)
library(ggdist)
library(mvtnorm)
library(bipartite)
library(plyr)


df <- list.files(path = "Data/Bernoulli_fruitset", pattern = "nseeds_fruitset_", full.names = T, recursive = TRUE) 
dd.list <- lapply(df, read.csv, sep=' ') 
names(dd.list) <- sapply(
	df,
	function(x){
		nn <- strsplit(x, '_')[[1]]
		nn[4]
	}
)

dd.fl <- read.csv("Data/flowering_weeks_plants.csv", sep=';') %>% select(Plant_id, Flowering_min)

for (i in seq_along(dd.list)){
  temp <- dd.list[[i]] 
  temp <- temp %>% mutate(Plant_id=as.character(Plant_id)) %>% left_join(dd.fl, by="Plant_id")
  
  cn <- colnames(temp)
  pollinators <- cn[which(!cn %in% c('Plant_id','Plant_sp','Response',
                                   'Plot','N_flowers','Type', 'Flowering_min') )]

  temp <- temp %>% mutate(degree= rowSums(select(.,all_of(pollinators))!=0))
  temp <- temp %>% filter(degree !=0)
  dd.list[[i]] <- temp
}

# perform the model fitting process for each focal species
dd.list <- dd.list[c("CLIB", "HCOM", "HHAL")]
  
for(i in names(dd.list)){
	focal <- i

	# prep the data for model fitting
	source('Code/setup_data.R')

	# fit things with cmdstanr and generate some pretty figures
	source('Code/fit_joint_model_stan.R')
}

