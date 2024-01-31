## Basics and Wrangle
require(magrittr)
library(plotly)
require(tidyverse)
require(data.table)
require(dplyr)
require(tidyr)

## NetSci
require(igraph)
require(NetSci)

`%ni%` <- Negate(`%in%`)

proximity_sign = function(N = N_reps, 
                          disease_genes = disease_genes_RA, 
                          Drugs = Targets,
                          distance = distance_PPI, 
                          disease_name = "arthritis rheumatoid", 
                          possible_targets = targets_DB){
  
  index = data.frame(gene = colnames(distance))
  index$index = 1:nrow(index)
  distance = distance %>% as.data.frame()
  
  id = Drugs$ID %>% unique()
  o = list()
  
  dis = index %>%
    filter(gene %in% disease_genes) %>%
    pull(index)
  
  for(i in 1:length(id)){
    # message(i)
    
    t_aux = Drugs %>%
      filter(ID == id[i]) %>%
      pull(Target)
    
    t = index %>%
      filter(gene %in% t_aux) %>%
      pull(index)
    
    
    d = distance[t,dis]
    
    if(length(t) >=1){
      if(length(t) == 1){
        out = min(d)
      } else if (length(t) > 1) {
        out = apply(d, 2, min) %>% 
          mean()
      }
      o[[i]] = data.frame(Proximity = out, 
                          Targets = length(t),
                          Drug = id[i], 
                          Disease = disease_name)
    }
  }
  o %<>% bind_rows()
  
  #### 
  #### Resample
  #### 
  
  target_amount = o$Targets %>% unique() %>% sort
  
  t = index %>%
    filter(gene %in% possible_targets) %>%
    pull(index)
  
  D_cl = distance[t, dis]
  
  distributions = list()
  
  for(i in target_amount){
    distributions[[i]] = ecdf(replicate(N, prox_helper(D_cl, i)))
  }
  
  for(i in 1:nrow(o)){
    o$p[i] = distributions[[o$Targets[i]]](o$Proximity[i])
  }
  
  return(o)
}

prox_helper = function(D_cl, no_targets){
  if(no_targets > 1){
    p_star = D_cl[sample(1:nrow(D_cl), size = no_targets),] %>% 
      apply(., 2, min) %>% 
      mean()
  } else if(no_targets == 1){
    p_star = D_cl[sample(1:nrow(D_cl), size = no_targets),] %>% 
      min()
  } 
  return(p_star)
}
