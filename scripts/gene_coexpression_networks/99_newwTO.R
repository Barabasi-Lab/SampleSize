library(Rcpp)
library(Rfast)

initial.options <- commandArgs(trailingOnly = FALSE)
scripts.dir <- dirname(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
matrixmult.file <- paste(scripts.dir, "99_matrixmult.cpp", sep="/")
sourceCpp(matrixmult.file)

wTO2 = 
  function (A, sign = c("abs", "sign")) 
  {
    if (sign %in% c("abs", "absolute")) {
      A = abs(A)
    }
    A_TF = as.data.frame(subset(A, select = row.names(A)))
    C = eigenMapMatMult(A, Rfast::transpose(A))
    W = C + A_TF
    K = matrix(NA, nrow(A_TF), ncol(A_TF))
    KI = rowSums(abs(A), na.rm = T)
    for (ii in 1:nrow(A_TF)) {
      for (jj in 1:ncol(A_TF)) {
        K[ii, jj] = min(KI[ii], KI[jj])
      }
    }
    WTO = round(W/(K + 1 - abs(A_TF)), 3)
    return(WTO)
  }

wTO.faster = function(Data, Overlap = row.names(Data), 
                      method = 'p', sign = 'sign', 
                      delta = 0.2, n = 10,
                      method_resampling = 'Bootstrap', lag = NULL, ID = NULL){
  Overlap = unique(as.character(Overlap))
  `%ni%` <- Negate(`%in%`)
  ##### Messages
  
  if(is.numeric(n) == F){
    stop("n must be numeric.")
  }
  if(n <= 0){
    stop("n must be greater than 0.")
  }
  if(is.data.frame(Data) == F){
    stop("Data must be a data.frame.")
  }
  
  if(method %ni% c("s", "spearman", "p", "pearson")){
    stop('Method must be: "s", "spearman", "p" or "pearson".')
  }
  
  if(method_resampling %ni% c("Bootstrap", "BlockBootstrap")){
    stop('Method must be: "Bootstrap" or "BlockBootstrap".')
  }
  if(method_resampling %in% "BlockBootstrap"){
    if (is.null(lag)&is.null(ID)){
      stop('If you want to use the "BlockBootstrap" please give a lag or the indivuals ID.')
    }
    if(!is.null(lag)&!is.null(ID)){
      stop('If you want to use the "BlockBootstrap" please give a lag OR the indivuals ID.')
    }
  }
  
  DIM_Overlap = nrow(subset(Data, row.names(Data) %in% Overlap))
  if(DIM_Overlap == 0){
    stop('There is no overlapping nodes. Please check your input "Overlap"')
  }
  if(!is.null(DIM_Overlap)){
    message(paste('There are',DIM_Overlap, "overlapping nodes,",dim(Data)[1],
                  "total nodes and" , dim(Data)[2],"individuals." ))
  }
  
  message("This function might take a long time to run. Don't turn off the computer.")
  
  wtomelt0 =  wTO::CorrelationOverlap(Data = Data, Overlap = Overlap, method = method) %>% 
    wTO2(., sign)
  `%>%` <- magrittr::`%>%`
  . <- NULL
  for ( i in 1:n){
    message(' ',i,' ', appendLF = FALSE)
    
    if(method_resampling == 'BlockBootstrap'){
      if (!is.null(lag)){
        nsampl = ifelse (ncol(Data) %% lag == 0, ncol(Data) %/% lag, ncol(Data) %/% lag +1)
        Y = sample(1:nsampl, size = nsampl, replace =  T)
        Vect = Y*lag
        j = lag - 1
        while( j > 0){
          Vect = cbind(Vect , Y*lag - j)
          j = j - 1
        }
        
        SAMPLES = c(Vect)
        SAMPLES[SAMPLES > ncol(Data)] <- NA
        SAMPLE = stats::na.exclude(SAMPLES)
        Data_boot = Data[,SAMPLE]
      }
      
      if(!is.null(ID)){
        ID %<>% as.factor
        bootID = sample(levels(ID), replace = TRUE)
        
        
        Data_boot = subset(Data, select = ID %in% bootID[1])
        
        for (k in 2:length(bootID)){
          Data_boot = cbind(Data_boot,
                            subset(Data, select = ID %in% bootID[k]))
        }
      }
      
      res = wTO::CorrelationOverlap(Data = Data_boot, Overlap = Overlap, method = method) %>% 
        wTO2(., sign)
      
    }
    else if (method_resampling != 'BlockBootstrap'){
      res = wTO::CorrelationOverlap(Data = Data[,sample(1:ncol(Data), replace  = TRUE)], Overlap = Overlap, method = method) %>% 
        wTO2(., sign)
    }
    
    
    U  = (res < wtomelt0 - delta) + (res > wtomelt0 + delta)
    if ( i == 1){
      out = U}
    if (i != 1){
      out = out + U
    }
    rm(res)
    rm (U)
  }
  
  wtomelt0 = wTO.in.line(wtomelt0)
  cor      = wTO.in.line(out)
  adj.pval = p.adjust(cor$wTO/n, method = 'BH')
  
  pval = data.table::data.table(wtomelt0, pval = cor$wTO/n, pval.adj = adj.pval)
  
  message('Done!')
  return(pval)
}


require(dplyr)
require(magrittr)
require(data.table)
require(wTO)
require(parallel)
