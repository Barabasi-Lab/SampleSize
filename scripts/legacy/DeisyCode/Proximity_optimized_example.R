setwd("~/Desktop/PostDoc/00_Projects/GDA_NC/01_code")

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

## ProgressBar
require(progress)

## Graphics
require(Cairo)
require(ggplot2)
require(ggrepel)
require(patchwork)
require(scales)
library(UpSetR)

`%ni%` <- Negate(`%in%`)

proximity_sign = function(N = N_reps, 
                          disease_genes = disease_genes_RA, 
                          Drugs = Targets,
                          distance = distance_PPI, 
                          disease_name = "arthritis rheumatoid", 
                          possible_targets = targets_DB){
  
  index = data.frame(gene = colnames(distance))
  index$index = 1:nrow(index)
  
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

load("../00_data/data/graphs_new.Rdata")

DT = fread("~/Dropbox (CCNR)/Biology/99_Toolbox/data/DrugBank/data/out/DB_Drug_Targets_v5.1.9.csv")
head(DT)

DT %<>% 
  select(DB_ID, Name, Gene_Target, Indication) %>%
  unique()

Targets = DT %>%
  select(ID = DB_ID, Target = Gene_Target) %>%
  filter(Target %in% V(gPPInc)$name)

Drugs = DT$DB_ID %>% unique()


diseases = c("arthritis rheumatoid", 
             "pre eclampsia",
             "inflammatory bowel diseases",
             "inflammation",
             "crohn disease")

targets_DB = Targets$Target %>% unique()
distance_PPI = distances(gPPI)  
distance_PPINC = distances(gPPInc)  

prox_dis_all2 = prox_dis_all = prox_dis_lcc = prox_dis_lcc2 = list()
# Targets_fake = Targets %>%
#   filter(ID %in% c("DB00001", "DB00002", "DB12010" ))

pb <- progress_bar$new(
  format = "  :spin Running [:bar] :percent eta: :eta",
  total = 100, 
  clear = FALSE, 
  width= 60)

N_ = 5000
for(j in 1:length(diseases)){
  pb$tick()
  
  disease_genes_x = GDA %>%
    filter(hgnc_symbol %in% V(gPPI)$name) %>% 
    filter(DiseaseName == diseases[j]) %>%
    pull(hgnc_symbol)
  
  prox_dis_all[[j]] = proximity_sign(N = N_, 
                                     disease_genes = disease_genes_x, 
                                     Drugs = Targets, 
                                     distance = distance_PPI, 
                                     disease_name = diseases[j], 
                                     possible_targets = targets_DB) %>%
    mutate(graph = "PPI") %>%
    mutate(type = "Complete")
  
  disease_genes_lcc = gPPI %>% 
    induced_subgraph(., disease_genes_x) %>% 
    extract_LCC() %>% 
    V() %>%
    unclass() %>%
    names()
  pb$tick()
  prox_dis_lcc[[j]] = proximity_sign(N = N_, 
                                     disease_genes = disease_genes_x, 
                                     Drugs = Targets, 
                                     distance = distance_PPI, 
                                     disease_name = diseases[j], 
                                     possible_targets = targets_DB) %>%
    mutate(graph = "PPI") %>%
    mutate(type = "LCC")
  
  disease_genes_x = GDA %>%
    filter(hgnc_symbol %in% V(gPPInc)$name) %>% 
    filter(DiseaseName == diseases[j]) %>%
    pull(hgnc_symbol)
  pb$tick()
  prox_dis_all2[[j]] = proximity_sign(N = N_, 
                                      disease_genes = disease_genes_x, 
                                      Drugs = Targets, 
                                      distance = distance_PPINC, 
                                      disease_name = diseases[j], 
                                      possible_targets = targets_DB)%>%
    mutate(graph = "PPI & NCI") %>%
    mutate(type = "Complete")
  
  disease_genes_lcc = gPPInc %>% 
    induced_subgraph(., disease_genes_x) %>% 
    extract_LCC() %>% 
    V() %>%
    unclass() %>%
    names()
  pb$tick()
  prox_dis_lcc2[[j]] = proximity_sign(N = N_, 
                                      disease_genes = disease_genes_lcc, 
                                      Drugs = Targets, 
                                      distance = distance_PPINC, 
                                      disease_name = diseases[j], 
                                      possible_targets = targets_DB) %>%
    mutate(graph = "PPI & NCI") %>%
    mutate(type = "LCC")
}



prox_dis = bind_rows(prox_dis_all, prox_dis_all2, prox_dis_lcc, prox_dis_lcc2) 

prox_dis %>%
  filter(Disease %in% "crohn disease") %>%
  ggplot() +
  aes(x = Proximity) +
  geom_histogram(bins = 100L, fill = "#112446") +
  theme_minimal() +
  facet_grid(vars(graph), vars(type))

prox_dis %>%
  filter(Disease %in% "crohn disease") %>%
  ggplot() +
  aes(x = Proximity, y = p) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal() +
  facet_grid(vars(graph), vars(type))

prox_dis %<>%
  ungroup() %>% 
  group_by(Disease, graph, type, Targets) %>%
  mutate(padj = p.adjust(p, method = "fdr"))

fwrite(prox_dis, "../03_out/proximity.csv")



RA_drugs = DT$DB_ID[grep(pattern = "rheumatoid arthritis", ignore.case = T, x = DT$Indication)] %>%
  unique()

prox_RA = prox_dis  %>% 
  filter(Disease %in% "arthritis rheumatoid") %>%
  filter(Drug %in% RA_drugs)

PPI_colors = c(#NCI = "#86C188",
  PPI = "#8D3B72",
  `PPI & NCI` = "#72E1D1")



prox_RA %>%
  mutate(sign = ifelse(p < 0.05, "p < 0.05", "p > 0.05")) %>% 
  ggplot() +
  aes(y = Proximity, 
      x = graph, 
      colour = graph, 
      group = Drug, 
      shape = sign) +
  geom_line(size = 0.5, color = "gray75") +
  geom_point(size = 3, color = "white", pch = 16) +
  geom_point(size = 3, alpha = 0.2) +
  scale_color_manual(values = PPI_colors ) +
  scale_shape_manual(values = c("p < 0.05" = 16, "p > 0.05" = 1))+
  theme_minimal() +
  facet_wrap(vars(type))+
  theme(legend.position = "bottom")


ggplot(prox_dis) +
  aes(x = p, fill = graph, colour = graph) +
  geom_density(adjust = 3, alpha = 0.2) +
  scale_color_manual(values = PPI_colors ) +
  scale_fill_manual(values = PPI_colors ) +
  theme_minimal()+
  facet_wrap(vars(type))

GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, 
                      xminv = x - violinwidth * (x - xmin), 
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}

prox_dis %>%
  filter(type == "LCC") %>% 
  ggplot() +
  aes(y = Proximity,  fill = graph, colour = graph) +
  geom_half_violin(alpha = 0.2, side = "l", scale = "count") +
  scale_color_manual(values = PPI_colors ) +
  scale_fill_manual(values = PPI_colors ) +
  theme_minimal() +
  facet_grid( vars(Disease))


prox_dis %>%
  # filter(type == "LCC") %>% 
  ggplot() +
  aes(x = type, y = Proximity,  fill = graph, colour = graph) +
  geom_split_violin(trim = FALSE, alpha = .4, scale = "count")+
  scale_color_manual(values = PPI_colors ) +
  scale_fill_manual(values = PPI_colors ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_grid( vars(Disease))
