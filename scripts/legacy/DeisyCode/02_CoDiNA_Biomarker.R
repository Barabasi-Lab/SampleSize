require(dplyr)
require(magrittr)
require(data.table)
require(ggplot2)
nR = fread(
  "~/Desktop/PostDoc/00_Projects/RA/00_Data/wTO/NonResponder_CN.csv")
R = fread(
  "~/Desktop/PostDoc/00_Projects/RA/00_Data/wTO/Responder_CN.csv")
Genes = c(nR$Node.1, nR$Node.2, R$Node.1, R$Node.2) %>% unique()
Genes = data.frame(Node.1 = Genes,
                   Node.2  = Genes,
                   CN = 0,
                   pval.fisher = 0)

nR = rbind(nR, Genes)
R = rbind(R, Genes)
PPI = fread("~/Desktop/PostDoc/02_Toolbox/03_Data/PPI_Symbol_Entrez.csv")
PPI = PPI[,-c(1,2)] %>% unique()
PPI %<>%
  filter(Symbol_A != "" & Symbol_B != "") %>%
  filter(Symbol_A != Symbol_B)

PPI$PPI = 1
PPI = OrderNames(PPI) %>% unique()

require(CoDiNA)
Diff = MakeDiffNet(Data = list(PPI, R, nR), Code = c('PPI', "R", "nR"), stretch = F)
Diff_cl = Diff %>% filter(Score_ratio > 1)
Diff_cl %>% group_by(Phi_tilde) %>% summarize(n = n())
Diff %>% group_by(Phi_tilde) %>% summarize(n = n())

Diff_cl_PPI = Diff_cl %>% filter(Phi_tilde %in% c("a", "g.PPI.nR", "g.PPI.R"))

fwrite(Diff_cl_PPI, '~/Desktop/PostDoc/00_Projects/RA/CN/out/CoDiNA_PPI_CL.csv')

Nodes = ClusterNodes(DiffNet = Diff_cl_PPI, cutoff.external = 0, cutoff.internal = 1)
table(Nodes$Phi_tilde)
Node_Cat = Nodes %>% mutate(p_adj = p.adjust(pval_Phi_Tilde, method = 'fdr')) %>%
  filter(p_adj < 0.01) %>%
  select(Node, Phi_tilde)


Nodes %>% 
  mutate(p_adj = p.adjust(pval_Phi_Tilde, method = 'fdr')) %>% 
  filter(p_adj < 0.01) %>% 
  group_by(Phi_tilde) %>%
  summarise(n = n())

Node_Cat = Nodes %>% 
  mutate(p_adj = p.adjust(pval_Phi_Tilde, method = 'fdr')) %>% 
  filter(p_adj < 0.01) %>%
  select(Node, Phi_tilde)

Node_Cat %>% 
  filter(Phi_tilde == "g.PPI.R") %>% 
  pull(Node) %>%
  write.table("~/Desktop/PostDoc/00_Projects/RA/CN/out/CN_R.txt", quote = F, 
              row.names = F, col.names = F)

Node_Cat %>% 
  filter(Phi_tilde == "g.PPI.nR") %>% 
  pull(Node) %>%
  write.table("~/Desktop/PostDoc/00_Projects/RA/CN/out/CN_nR.txt", quote = F, 
              row.names = F, col.names = F)

Node_Cat %>% 
  filter(Phi_tilde == "a") %>% 
  pull(Node) %>%
  write.table("~/Desktop/PostDoc/00_Projects/RA/CN/out/CN_a.txt", quote = F, 
              row.names = F, col.names = F)


Node_Cat %>% 
  filter(Phi_tilde != "U") %>% 
  # pull(Node) %>%
  write.table("~/Desktop/PostDoc/00_Projects/RA/CN/out/CN_All.txt", quote = F, 
              row.names = F, col.names = F)

source("~/Desktop/PostDoc/02_Toolbox/01_NetSci/LCC.R")
require(igraph)
`%ni%` <- Negate(`%in%`)

histLCC = function(LCC_L, Name){
  Fn = ecdf(LCC_L$LCCZ)
  p = 1 - Fn(LCC_L$LCC)
  
  lim = c(LCC_L$LCC, LCC_L$LCCZ)
  LCC_L[[1]]%>% hist(., las = 1, 
                     main = "",
                     # ylim = c(0,350),
                     xlim = c(min(lim -10), max(lim + 10)),
                     col = 'gray75', ylab = "")
  abline(v = LCC_L$LCC, col = "red")
  title(main = Name, sub = paste0("LCC: ",
                                  round(LCC_L$LCC,0),
                                  " (", 
                                  round(LCC_L$mean,2),
                                  " ± ", 
                                  round(LCC_L$sd,2),"; ",
                                  "p: " ,round(p,2),
                                  ")"))
  
  
}


require(igraph)
Biom_R =   Node_Cat %>% filter(Phi_tilde == "g.PPI.R") %>% 
  pull(Node) 
Biom_nR =   Node_Cat %>% filter(Phi_tilde == "g.PPI.nR") %>% 
  pull(Node) 
Biom_a =   Node_Cat %>% filter(Phi_tilde == "a") %>% 
  pull(Node) 

gPPI = graph_from_data_frame(PPI, directed = FALSE) %>% simplify()
Biom_R = subset(Biom_R, Biom_R %in% V(gPPI)$name)
Biom_nR = subset(Biom_nR, Biom_nR %in% V(gPPI)$name)
Biom_a = subset(Biom_a, Biom_a %in% V(gPPI)$name)
Codina_Biomarker_R = Biom_R
Codina_Biomarker_nR = Biom_nR
Codina_Biomarker_a = Biom_a

LCC_R = LCC_Significance(N = 1000, 
                         Targets = Biom_R, 
                         PPIg = gPPI, bins = 100)

LCC_nR = LCC_Significance(N = 1000, 
                          Targets = Biom_nR, 
                          PPIg = gPPI, bins = 100)

LCC_a= LCC_Significance(N = 1000, 
                        Targets = Biom_a, 
                        PPIg = gPPI, bins = 100)
par(mfrow=c(1,3))
histLCC(LCC_L = LCC_R, "Biomarker - Responder")
histLCC(LCC_L = LCC_nR, "Biomarker - non Responder")
histLCC(LCC_L = LCC_a, "Biomarker - both")

CODINA_OR_LCC = list(LCC_R = LCC_R, 
                     LCC_a = LCC_a, 
                     LCC_nR = LCC_nR)

################
################
################

GDA = fread("~/Desktop/PostDoc/00_Projects/GDDA/New/finaldata/GDA_18122020_ClassFromDisGeNet.csv")
GDA_S = GDA %>%
  mutate(SW = Strong + Weak) %>%
  filter(SW > 1) 


Dis = GDA_S %>% group_by(DiseaseName) %>%
  summarise(n = n()) %>% 
  filter(n > 5) %>% 
  pull(DiseaseName)

GDA_S %<>% filter(DiseaseName %in% Dis)
GDA_S_Dis_Gene = GDA_S %>% select(DiseaseName, hgnc_symbol) %>% unique()

Genes = Node_Cat
names(Genes ) = c('hgnc_symbol', 'DiseaseName')

data = rbind(GDA_S_Dis_Gene, Genes) %>% unique()




source("~/Desktop/PostDoc/02_Toolbox/01_NetSci/LCC.R")
source("~/Desktop/PostDoc/02_Toolbox/01_NetSci/S_ab.R")
source("~/Desktop/PostDoc/02_Toolbox/01_NetSci/01_Shapes.R")

SAB = new_sab(gPPI, data)
Sab = SAB$Sab %>% as.matrix %>% reshape2::melt() 

x = Sab %>% filter(Var2 %in% c("a","g.PPI.nR","g.PPI.R")) %>% na.exclude()

ggplot(x) +
  aes(x = Var2, y = Var1, fill = value) +
  geom_tile(size = 1L) +
  scale_fill_gradient() +
  theme_minimal()

require(superheat)
require(tidyr)
x %<>% pivot_wider(names_from = Var2, values_from = value) 
x[, -1] %>% superheat(., pretty.order.rows = TRUE, pretty.order.cols = T, row.dendrogram = T)

#################
#################
#################
RA = data %>% filter(DiseaseName %in% "rheumatoid arthritis") %>% pull(hgnc_symbol)
g = Diff_cl_PPI %>% 
  graph_from_data_frame(., directed = F) %>% 
  delete.vertices(degree(.) == 0)

colors = c("5dc5fe","92eefc","F7E6A1","fab2dc","f99fe6") %>% 
  paste0("#", .)
c1 = adjustcolor("gray80", alpha.f = 0.8)
c4 = adjustcolor(colors[1], alpha.f = 0.8)
c2 = adjustcolor(colors[3], alpha.f = 0.8)
c3 = adjustcolor(colors[5], alpha.f = 0.8)

cc1 = adjustcolor("gray80", alpha.f = 0.3)
cc4 = adjustcolor(colors[1], alpha.f = 0.3)
cc2 = adjustcolor(colors[3], alpha.f = 0.3)
cc3 = adjustcolor(colors[5], alpha.f = 0.3)


V(g)$color = ifelse(V(g)$name %in% Biom_a , c2, c1)
V(g)$color = ifelse(V(g)$name %in% Biom_nR, c3, V(g)$color)
V(g)$color = ifelse(V(g)$name %in% Biom_R , c4, V(g)$color)
V(g)$frame.color = ifelse(V(g)$name %in% RA, "red", V(g)$color)

E(g)$color = ifelse(E(g)$Phi_tilde %in% "a" , cc2, cc1)
E(g)$color = ifelse(E(g)$Phi_tilde %in% "g.PPI.nR" , cc3, E(g)$color)
E(g)$color = ifelse(E(g)$Phi_tilde %in% "g.PPI.R" , cc4, E(g)$color)


g %<>%   delete.edges(E(g)[E(g)$color == cc1])
g %<>%   delete.vertices(V(g)$frame.color == c1 )
E(g)$W = 1 - E(g)$Score_ratio %>% CoDiNA::normalize()
E(g)$width = E(g)$Score_ratio %>% CoDiNA::normalize() + 0.1
E(g)$weight = E(g)$Score_ratio # (E(g)$width + 0.1)^6


V(g)$size = strength(g) %>% CoDiNA::normalize() + 0.1
V(g)$size = V(g)$size*5
V(g)$size = ifelse(V(g)$name %in% RA, V(g)$size*3, V(g)$size)
V(g)$label = NA
par(mfrow = c(1,1))



# g %<>% delete.vertices(V(g)$frame.color == c2 )
# g %<>%   delete.vertices(strength(g) < 3)

V(g)$memb = components(g)$membership
g %<>%   delete.vertices(V(g)$memb != 1)

coord = layout_with_fr(g, weights = E(g)$weight ^ 4)
plot(g, layout = coord)
legend('bottom', legend = c("Common", "Non Responder", "Responder", "RA"), 
       fill = c(c2,c3,c4, "red"), bty ='n', ncol = 4)



A = g %>% get.adjacency(attr = 'W') %>% as.matrix()

B = c(Biom_a, Biom_nR, Biom_R)
A[B, B] %>% 
superheat::superheat(., pretty.order.rows = T, pretty.order.cols = T, scale = F)
d = A[B, B] %>% as.matrix() %>% dist(method = 'binary')%>% as.matrix() 
hc = d %>% hclust()
mem = hc %>% cutree(h = 2.5)
mem = mem %>% as.data.frame()
mem$gene = row.names(mem)
names(mem)[1]= "Member"
gps = mem %>% group_by(Member) %>% summarise(n = n()) %>% filter(n > 5) %>% pull(Member)

genes = mem %>% filter(Member %in% gps )
genes$gene %in% Biom_R


A[genes$gene, genes$gene] %>% 
  superheat::superheat(., pretty.order.rows = T, pretty.order.cols = T, scale = F)

g_biom = g %>% delete.vertices(., V(g)$name %ni% genes$gene)
plot(g_biom)
