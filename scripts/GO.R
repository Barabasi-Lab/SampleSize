####################################
### Author: Deisy Morselli Gysi 
### Date: 03/10/2022
####################################
### Keyworks: Tool; Gene Enrichment; GO
### Description: Calculates GO enrichment
### Input: 
### The GO function has three parameters:
### - `ID_type`: The gene identification system. One of the following:
###   - "symbol", "entrez", "ensembl"
### - `bg`: Vector containing the background set of genes (the ones expressed in the platform, the ones in the network)
### - `g`: vector containing the set of genes you want to define the enrichment. 
### - `ONTO`: The ontology you are using. One of the following:
###   - "BP" (biological process), "MF" (molecular function), "CC" (celular component)
### 
### Ouput: a named list
### - Res: All possible enrichment.
### - Sign: Significant enriched terms.
### - Gene2GO = Genes associated to each GO.
####################################

`%ni%` <- Negate(`%in%`)
GO = function(ID_type = "symbol", 
              bg, 
              g, 
              ONTO = 'BP'){
  require(topGO)
  require(magrittr)
  require(dplyr)
  if(ID_type %ni% c("symbol", "entrez", "ensembl")){
    stop('Please, use one of the following ID_type for your gene: "symbol", "entrez", "ensembl"')
  }
  if(ONTO %ni% c("BP", "MF", "CC")){
    stop('Please, use one of the following ontology: "BP", "MF", "CC"')
  }
  geneID2GO <- annFUN.org(ONTO, 
                          mapping = "org.Hs.eg.db",
                          ID = ID_type) %>% 
    inverseList()
  
  Background <- data.frame(V1 = bg)
  GenesOI <- data.frame(V1 = g)
  geneList <- factor(as.integer(Background$V1 %in% GenesOI$V1))
  names(geneList) <- Background$V1
  
  GOdata <- new("topGOdata",
                ontology = ONTO,
                allGenes = geneList,
                geneSel = GenesOI,
                annot = annFUN.gene2GO,  # the new annotation function
                gene2GO = geneID2GO)    ## the gene ID to GOs dataset
  test.statk <- new("classicScore", testStatistic = GOKSTest, name = "KS test")
  test.statf <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  test.weight<- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test.statf)
  resultKS <- getSigGroups(GOdata, test.statk)
  pvalFis <- score(resultFisher)
  pvalKS <- score(resultKS)
  resultWeight <- getSigGroups(GOdata, test.weight)
  pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
  cor(pvalFis, pvalWeight)
  geneData(resultWeight)
  allRes_BP <- GenTable(GOdata, 
                        classic = resultFisher, 
                        KS = resultKS, 
                        weight = resultWeight, 
                        orderBy = "weight", 
                        ranksOf = "classic", 
                        topNodes = 500, 
                        numChar = 400)
  
  
  goID <- allRes_BP[, "GO.ID"]
  ID2GO = list()
  for(i in 1:length(geneID2GO)){
    ID2GO[[i]] = data.frame(gene = names(geneID2GO)[i],
                            GOs = geneID2GO[[i]])
  }
  ID2GO %<>% bind_rows()
  
  ID2GO_sign = ID2GO %>% 
    filter(GOs %in% goID)
  
  allRes_BP_f = subset(allRes_BP, allRes_BP$weight < 0.01)
  return(list(Res = allRes_BP, 
              Sign = allRes_BP_f, 
              Gene2GO = ID2GO_sign))
}