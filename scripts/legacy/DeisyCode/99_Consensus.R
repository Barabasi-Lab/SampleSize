

wTO.Consensus = function(data){
  if (!is.list(data)){
    stop("data must be a list of data.frames.")
  }
  ### Weight
  weight = pval = nodes = list()
  
  for ( i in 1:length(data)){
    weight[[i]] = data[[i]][,1:3]
    pval[[i]] = data[[i]][,c(1:2,4)]
    ID = unique(c(as.character(data[[i]]$Node.1), as.character(data[[i]]$Node.2)))
    nodes[[i]] = data.frame(ID =  ID)
    names(weight[[i]])[3] = paste0(names(weight[[i]])[3], i)
    names(pval[[i]])[3] = paste0(names(pval[[i]])[3], i)
  }
  
  weight = plyr::join_all(weight, type = 'full')
  pval = plyr::join_all(pval, type = 'full')
  nodes = plyr::join_all(nodes, type = 'inner')
  
  message(paste('Total common nodes:', nrow(nodes)))
  weight = subset(weight, weight$Node.1 %in% nodes$ID & weight$Node.2 %in% nodes$ID)
  pval = subset(pval, pval$Node.1 %in% nodes$ID & pval$Node.2 %in% nodes$ID)
  pval[is.na(pval)] <- 1
  weight[is.na(weight)] <- 0.01
  
  wTOCN = CN_aux(weight[, -c(1:2)])
  pvalue_fisher = fishermethod(pval[, -c(1:2)])
  
  Out = data.frame(Node.1 = pval[,1], Node.2 = pval[,2],
                   CN = wTOCN, pval.fisher = pvalue_fisher)
  return(Out)
}


fishermethod = function(data_x){
  chi = rowSums(log(data_x))*-2
  pval = sapply(chi, function(x) stats::pchisq(x, 2*ncol(data_x), lower.tail = FALSE))
  return(pval)
}

CN_aux = function(data_x){
  abs_x = apply(data_x, 2, abs)
  sum_abs_x = apply(abs_x, 1, sum)
  div = (abs_x/sum_abs_x) * data_x
  wTO_cons = apply(div, 1, sum)
  return(wTO_cons)
}