library(plotly)
library(UpSetR)
require(tidyverse)
require(data.table)
require(scales)
require(dplyr)
require(ggplot2)
require(ggrepel)
require(magrittr)
require(tidyr)
require(igraph)
require(patchwork)
library(showtext)
require(Cairo)
require(progress)
require(ROCR)
`%ni%` <- Negate(`%in%`)

load("~/Desktop/PostDoc/00_Projects/GDA_NC/00_data/data/graphs_new.Rdata")
load("~/Desktop/PostDoc/00_Projects/GDA_NC/00_data/Correlation.RData")

# head(MM2)
# 
MM2$bip = ifelse(is.na(MM2$bip), 0, MM2$bip)
genes_cor = unique(unique(MM2$from), unique(MM2$to))

protein_coding = n %>%
  filter(Category_cl %in% c("Protein Coding", "TF"))


PC_interactions = MM2 %>%
  filter( from %in% protein_coding$Symbol & 
            to %in% protein_coding$Symbol)

rm(MM2)
PC_interactions$bip = ifelse(PC_interactions$bip > 0, 1, 0)
PC_interactions$bipppi = ifelse(PC_interactions$bip > 0 | PC_interactions$PCPC > 0, 1, 0)


# pc_l = PC_interactions %>% 
#   select(- PPINC) %>% 
#   pivot_longer(cols = c("PCPC", "bip", "bipppi"))
# 
# teste = pc_l[sample(1:nrow(pc_l), size = 200000), ] %>%
#   rename(bind = value) %>% 
#   rename(base = name) %>% 
#   mutate(bind = factor(bind, levels = c(0, 1), labels = c("no bind", "bind"))) %>% 
#   pivot_longer(cols = c(wto, cor))
# 
# pc_l %>% 
#   rename(bind = value) %>% 
#   rename(base = name) %>% 
#   mutate(bind = factor(bind, levels = c(0, 1), labels = c("no bind", "bind"))) %>% 
#   pivot_longer(cols = c(wto, cor)) %>% 
#   ggplot() +
#   aes(x = base, 
#       y = abs(value), 
#       fill = bind, 
#       colour = bind) +
#   geom_boxplot(outlier.size = 0, 
#                outlier.alpha = 0, 
#                alpha = 0.5) +
#   # scale_fill_manual(direction = 1) +
#   # scale_color_manual(direction = 1) +
#   theme_minimal() +
#   facet_wrap(vars(name))+
#   theme(legend.position = "bottom")
# 
# 
# summary_data = pc_l  %>%
#   group_by(name, value) %>%
#   summarise(mean_cor_abs = mean(abs(cor)), 
#             mean_cor = mean(cor),
#             
#             median_cor_abs = median(abs(cor)), 
#             median_cor = median(cor),
#             
#             Q1_cor_abs = quantile(abs(cor), 0.25), 
#             Q1_cor = quantile(cor, 0.25),
#             
#             Q3_cor_abs = quantile(abs(cor), 0.75), 
#             Q3_cor = quantile(cor, 0.75),
#             
#             wto_abs = mean(abs(wto)), 
#             wto = mean(wto),
#             
#             median_wto_abs = median(abs(wto)), 
#             median_wto = median(wto),
#             
#             Q1_wto_abs = quantile(abs(wto), 0.25), 
#             Q1_wto = quantile(wto, 0.25),
#             
#             Q3_wto_abs = quantile(abs(wto), 0.75), 
#             Q3_wto = quantile(wto, 0.75))
# 
# summary_data

#### Opt2 
binds = list()

binds[[1]] = PC_interactions %>%
  filter(PCPC > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Direct") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds[[2]] = PC_interactions %>%
  filter(bip > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Indirect") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds[[3]] = PC_interactions %>%
  filter(bip > 0 | PCPC > 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "Indirect + Direct") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds[[4]] = PC_interactions %>%
  filter(bip == 0 | PCPC == 0) %>%
  dplyr::select(from,to,wto,cor) %>%
  mutate(type = "None") %>%
  pivot_longer(cols = c( 'wto','cor' ))

binds %<>% bind_rows()
binds$value = round(binds$value, 3)
# rm(PC_interactions)
# binds2 = binds[sample(1:nrow(binds), size = 10000),]

KW_wto = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "wto") %>% 
  kruskal.test(abs(value) ~ type, .)

Dunn_wto = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "wto") %>% 
  FSA::dunnTest(abs(value) ~ type, .)

Dunn_wto$res %>%
  as.data.frame() %>% 
  # mutate(padj = p.adjust(None, method = "fdr")) %>% 
  mutate(KW = KW_wto$p.value) %>% 
  fwrite("~/Desktop/PostDoc/00_Projects/GDA_NC/00_data/Dunn_test_wto.csv")

rm(KW_wto)
rm(Dunn_wto)

KW_cor = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "cor") %>% 
  kruskal.test(abs(value) ~ type, .)

Dunn_cor = binds %>% 
  mutate(type = factor(type, levels = c("None","Direct", "Indirect","Indirect + Direct")))%>% 
  filter(name == "cor") %>% 
  FSA::dunnTest(abs(value) ~ type, .)

Dunn_cor$res %>%
  as.data.frame() %>% 
  # mutate(padj = p.adjust(None, method = "fdr")) %>% 
  mutate(KW = KW_cor$p.value) %>% 
  fwrite("~/Desktop/PostDoc/00_Projects/GDA_NC/00_data/Dunn_test_cor.csv")

rm(KW_cor)
rm(Dunn_cor)


colors_defined = c("Direct" = "#DBA9CA",
                   
                   "Indirect" = "#BCF1E9",
                   
                   "Indirect + Direct" = "#FAF69E",
                   
                   "None" = "#B8BECC")

Cairo::CairoPDF("~/Desktop/PostDoc/00_Projects/GDA_NC/02_figs/V5_newdata/Correlation_Differences.pdf",
                width = 10, height = 8)
p0 = binds %>%
  ggplot() +
  aes(x = abs(value),
      y = type,
      fill = type,
      colour = type) +
  geom_boxplot(alpha = 0.8,
               outlier.size = 0,
               outlier.alpha = 0) +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "bottom",
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5),
        axis.title.x.top = element_text(face = "bold",
                                        hjust = 0.5),
        strip.text = element_text(face = "bold",
                                  hjust = 0.5)) +
  facet_wrap(vars(name)) +
  labs (x = "Absolute co-expression weight",
        y = "Interaction type",
        color = "",
        fill = "")
p0
dev.off()




####
pred_PPI_cor = ROCR::prediction(predictions = abs(PC_interactions$cor),
                                labels = PC_interactions$PCPC)

pred_PPI_wto = ROCR::prediction(predictions = abs(PC_interactions$wto),
                                labels = PC_interactions$PCPC)

pred_BIP_cor = ROCR::prediction(predictions = abs(PC_interactions$cor),
                                labels = PC_interactions$bip)

pred_BIP_wto = ROCR::prediction(predictions = abs(PC_interactions$wto),
                                labels = PC_interactions$bip)

pred_BIPPPI_cor = ROCR::prediction(predictions = abs(PC_interactions$cor),
                                   labels = PC_interactions$bipppi)

pred_BIPPPI_wto = ROCR::prediction(predictions = abs(PC_interactions$wto),
                                   labels = PC_interactions$bipppi)

require(pROC)
PC_interactions$abscor = abs(PC_interactions$cor)
PC_interactions$abswto = abs(PC_interactions$wto)

roc_ppi_wto = roc(PC_interactions, PCPC, abswto, ci = T); roc_ppi_wto
roc_ppi_cor = roc(PC_interactions, PCPC, abscor, ci = T); roc_ppi_cor

roc_bip_wto = roc(PC_interactions, bip, abswto, ci = T); roc_bip_wto
roc_bip_cor = roc(PC_interactions, bip, abscor, ci = T); roc_bip_cor

roc_bipppi_wto = roc(PC_interactions, bipppi, abswto, ci = T); roc_bipppi_wto
roc_bipppi_cor = roc(PC_interactions, bipppi, abscor, ci = T); roc_bipppi_cor

get_auc = function(rocs, base = "xxx", measure = 'wto'){
  o = rocs$ci %>% 
    as.vector() %>% 
    t %>%
    as.data.frame()  %>%
    mutate(base = base, 
           measure = measure) 
  names(o)[1:3] = c("CI_l", "AUC", "CI_h")
  
  return(o)
}

perf_sum_rocs = list()
perf_sum_rocs[[1]] = roc_ppi_wto %>%
  get_auc(., base = "Direct", 
          measure = "wTO")

perf_sum_rocs[[2]] = roc_ppi_cor %>%
  get_auc(., base = "Direct", 
          measure = "cor")

perf_sum_rocs[[3]] = roc_bip_wto %>%
  get_auc(., base = "Indirect", 
          measure = "wTO")

perf_sum_rocs[[4]] = roc_bip_cor %>%
  get_auc(., base = "Indirect", 
          measure = "wTO")

perf_sum_rocs[[5]] = roc_bipppi_wto %>%
  get_auc(., base = "Direct + Indirect", 
          measure = "wTO")

perf_sum_rocs[[6]] = roc_bipppi_cor %>%
  get_auc(., base = "Direct + Indirect", 
          measure = "cor")

perf_sum_rocs %<>% bind_rows()



performances = function(data, label){
  roc <- performance(data,"tpr","fpr")
  roc = data.frame(x = as.numeric(roc@x.values[[1]]),
                   y = as.numeric(roc@y.values[[1]]),
                   type = "ROC")

  auc <- performance(data,"auc")
  auc = auc@y.values

  prec_rec <- performance(data, "prec", "rec")
  prec_rec = data.frame(x = as.numeric(prec_rec@x.values[[1]]),
                        y = as.numeric(prec_rec@y.values[[1]]),
                        type = "ROCPR")


  aucpr <- performance(data, "aucpr")
  aucpr = aucpr@y.values

  out =  bind_rows(roc, prec_rec) %>%
    mutate(base = label) %>%
    mutate(AUC = auc) %>%
    mutate(AUCPR = aucpr)

  return(out)
}
perfs = list()

perfs[[1]] = pred_PPI_cor  %>%
  performances(., label = "PPI cor")

perfs[[2]] = pred_PPI_wto  %>%
  performances(., label = "PPI wto")

perfs[[3]] = pred_BIP_cor  %>%
  performances(., label = "BIP cor")

perfs[[4]] = pred_BIP_wto  %>%
  performances(., label = "BIP wto")

perfs[[5]] = pred_BIPPPI_cor  %>%
  performances(., label = "BIPPPI cor")

perfs[[6]] = pred_BIPPPI_wto  %>%
  performances(., label = "BIPPPI wto")

rm(pred_BIP_cor)
rm(pred_BIP_wto)

rm(pred_BIPPPI_cor)
rm(pred_BIPPPI_wto)

rm(pred_PPI_cor)
rm(pred_PPI_wto)

perfs %<>% bind_rows()

fwrite(perfs, "~/Desktop/PostDoc/00_Projects/GDA_NC/00_data/Performance.csv")

perfs = fread("~/Desktop/PostDoc/00_Projects/GDA_NC/00_data/Performance.csv")

perf_sum = perfs %>%
  select(base, AUC,AUCPR) %>%
  unique() %>% 
  pivot_longer(cols = c("AUC", "AUCPR")) 

tmp = perf_sum$base %>% 
  stringr::str_split(., " ", simplify = T) %>%
  as.data.frame()
perf_sum$base = tmp$V1
perf_sum$measure = tmp$V2

p1 = perf_sum %>% 
  mutate(base = factor(base, 
                       levels = c("PPI", 
                                  "BIP", 
                                  "BIPPPI"), 
                       labels = c("Direct",
                                  "Indirect",
                                  "Indirect + Direct"))) %>% 
  ggplot() +
  aes(x = base, 
      fill = base, 
      colour = base, 
      weight = value) +
  geom_bar() +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  theme_minimal() +
  facet_grid(vars(name), vars(measure))+ 
  theme(legend.position = "bottom", 
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5), 
        axis.title.x.top = element_text(face = "bold", 
                                        hjust = 0.5), 
        strip.text = element_text(face = "bold", 
                                  hjust = 0.5)) +
  labs(y = NULL, #"AUC / AUCPR", 
       x = "Base", 
       fill = "Base", 
       color = "Base")

p1_1 = perf_sum %>% 
  filter(name == "AUC") %>% 
  mutate(base = factor(base, 
                       levels = c("PPI", 
                                  "BIP", 
                                  "BIPPPI"), 
                       labels = c("Direct",
                                  "Indirect",
                                  "Indirect + Direct"))) %>% 
  ggplot() +
  aes(x = base, 
      fill = base, 
      colour = base, 
      weight = value) +
  geom_bar() +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  theme_minimal() +
  # ylim(c(0.5,1)) +
  facet_grid(vars(name), vars(measure))+ 
  theme(legend.position = "bottom", 
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5), 
        axis.title.x.top = element_text(face = "bold", 
                                        hjust = 0.5), 
        strip.text = element_text(face = "bold", 
                                  hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, 0.7))+
  labs(y = NULL, #"AUC / AUCPR", 
       x = "Base", 
       fill = "Base", 
       color = "Base")

p1_2 = perf_sum %>% 
  filter(name == "AUCPR") %>% 
  mutate(base = factor(base, 
                       levels = c("PPI", 
                                  "BIP", 
                                  "BIPPPI"), 
                       labels = c("Direct",
                                  "Indirect",
                                  "Indirect + Direct"))) %>% 
  ggplot() +
  aes(x = base, 
      fill = base, 
      colour = base, 
      weight = value) +
  geom_bar() +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  theme_minimal() +
  # ylim(c(0.5,1)) +
  facet_grid(vars(name), vars(measure))+ 
  theme(legend.position = "bottom", 
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5), 
        axis.title.x.top = element_text(face = "bold", 
                                        hjust = 0.5), 
        strip.text = element_text(face = "bold", 
                                  hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1))+
  labs(y = NULL, #"AUC / AUCPR", 
       x = "Base", 
       fill = "Base", 
       color = "Base")

# 
# colors_defined = c("PPI cor" = "#DBA9CA",
#                    "PPI wto" = "#652A51",
#                    
#                    "BIP cor" = "#BCF1E9",
#                    "BIP wto" = "#197668",
#                    
#                    "BIPPPI cor" = "#FAF69E",
#                    "BIPPPI wto" = "#756F06")



p1 = perf_sum %>% 
  mutate(base = factor(base, 
                       levels = c("PPI", 
                                  "BIP", 
                                  "BIPPPI"), 
                       labels = c("Direct",
                                  "Indirect",
                                  "Indirect + Direct"))) %>% 
  ggplot() +
  aes(x = base, 
      fill = base, 
      colour = base, 
      weight = value) +
  geom_bar() +
  scale_fill_manual(values = colors_defined) +
  scale_color_manual(values = colors_defined) +
  theme_minimal() +
  facet_grid(vars(name), vars(measure))+ 
  theme(legend.position = "bottom", 
        text = element_text(size = 18),
        axis.title  = element_text(face = "bold",
                                   hjust = 0.5), 
        axis.title.x.top = element_text(face = "bold", 
                                        hjust = 0.5), 
        strip.text = element_text(face = "bold", 
                                  hjust = 0.5)) +
  labs(y = NULL, #"AUC / AUCPR", 
       x = "Base", 
       fill = "Base", 
       color = "Base")


Cairo::CairoPDF("~/Desktop/PostDoc/00_Projects/GDA_NC/02_figs/V5_newdata/AUC_AUPR_cor.pdf", 
                width = 10, height = 8)
p1
dev.off()

p0 = p0 + 
  labs(y = NULL) + 
  theme(legend.position = "none") + 
  xlim(c(0,1))

p1 = p1 + labs(x = "Interaction Type") + ylim(c(0,1))
require(patchwork)

Cairo::CairoPDF("~/Desktop/PostDoc/00_Projects/GDA_NC/02_figs/V5_newdata/AUC_AUPR_cors.pdf", 
                width = 18, height = 15)
(p0/p1) + plot_layout(heights = c(1.5, 3))
dev.off()

Cairo::CairoPDF("~/Desktop/PostDoc/00_Projects/GDA_NC/02_figs/V5_newdata/AUC_AUPR_cors2.pdf", 
                width = 18, height = 15)
(p0/(p1_1 + theme(legend.position = "none") + labs(x = "Network"))) + plot_layout(heights = c(1.5, 1.5))
dev.off()

##### Calculate Performance and 
##### cors for each disease
