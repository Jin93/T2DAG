SEM.fit = function(p.model, ntotal, prefix, Data, A){
  if (p.model == ncol(A)){
    Qest = matrix(0,p.model,p.model) # estimated coefficient matrix Q
    Epsi = matrix(0,ntotal,p.model) # estimated residuals
    Rdiag = numeric()
    for (i in 1:p.model){
      non0i = which(A[i,]!=0) # index for the xi's that are included in the i-th regression:
      if (length(non0i)>0){ # if there exist directed edges from xj's to this xi
        formulai=formula(paste0(paste(paste0(prefix,i), paste(paste0(prefix,non0i),collapse="+"), sep="~"), " -1"))
        fiti = lm(formulai,data = Data)
        Qest[i,non0i] = fiti$coefficients
        Epsi[,i] = Data[,i] - as.matrix(Data[,non0i])%*%t(t(fiti$coefficients))
      }
      # if there is no directed edge from any xj to this xi
      if (length(non0i)==0) Epsi[,i] = Data[,i]
      Rdiag[i] = var(Epsi[,i]) # estimated Ri's
    }
  }
  if (p.model < ncol(A)){
    # p: the # of regression! not the number of nodes in the gene expression data
    # A: p x p+n.confounder
    Qest = matrix(0,ncol(A),ncol(A)) # estimated coefficient matrix Q
    Epsi = matrix(0,ntotal,ncol(A)) # estimated residuals
    Rdiag = numeric()
    for (i in 1:p.model){
      non0i = which(A[i,]!=0) # indices of the xi's that are included in the i-th regression:
      if (length(non0i)>0){ # if there exist directed edges from xj's to this xi
        formulai=formula(paste0(paste(paste0(prefix,i), paste(paste0(prefix,non0i),collapse="+"), sep="~"), " -1"))
        fiti = lm(formulai,data = Data)
        Qest[i,non0i] = fiti$coefficients
        Epsi[,i] = Data[,i] - as.matrix(Data[,non0i])%*%t(t(fiti$coefficients))
      }
      # if there is no directed edge from any xj to this xi
      if (length(non0i)==0) Epsi[,i] = Data[,i]
      Rdiag[i] = var(Epsi[,i]) # estimated Ri's
    }
  }
  return(list(Qest=Qest, Rdiag=Rdiag))
}
