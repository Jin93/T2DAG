#' Performs the T2DAG test.
#' @param X The n1 x p sample matrix from group 1. Each row corresponds to a p-dimensional sample.
#' @param Y The n2 x p sample matrix from group 2. Each row corresponds to a p-dimensional sample.
#' @param A The binary adjacency matrix, where an entry aij=0 indicates no directed relationship
#' from variable i to variable j, and aij=1 indicates that there exists a directed relationship
#' from variable i to variable j.
#' @return A list that contains the following elements:
#' @param Qest The estimate of the coefficient matrix Q.
#' @param Rdiag The estimate of the covariance matrix of the residuals.
#' @param T2DAG.stat T2DAG test statistic.
#' @param T2DAG.pval P-value of the T2DAG test.
#' @param T2DAG.pval Conclusion of the T2DAG test. 1 if H0 rejected and 0 otherwise.
#' @export
T2DAG = function(X, Y, A){
  p = nrow(A) # number of genes in the pathway
  n = c(nrow(X),nrow(Y)) # sample sizes
  N = sum(n) # total sample size
  nmin = min(n[1],n[2]) # minimum of the two sample sizes

  meandiff = colMeans(X)-colMeans(Y)
  Z = rbind(scale(X, center = T, scale = F), scale(Y, center = T, scale = F))
  Z = as.data.frame(Z)
  colnames(Z) = c(paste0("Z", 1:p))

  SEM.results = SEM.fit(p, N, 'Z', Z, A[1:p,1:p])
  Qest = SEM.results$Qest
  Rdiag = SEM.results$Rdiag
  #### T= n1*n2/(n1+n2) (Xbar-Ybar)^T (I-Qhat) R^(_1)hat (I-Qhat) (Xbar-Ybar)
  Sigma.hat = solve(diag(1,p)-Qest)%*%diag(Rdiag)%*%solve(t(diag(1,p)-Qest))
  T.graph = (n[1]*n[2])/(N)*(t(colMeans(X)-colMeans(Y))%*% solve(Sigma.hat[1:p,1:p]) %*%t(t(colMeans(X)-colMeans(Y))))

  # T2DAG test statistic
  T2DAG.stat = (T.graph - p)/sqrt(2*p)
  T2DAG.pval = 2*pnorm(q=abs(T2DAG.stat),mean=0,sd=1,lower.tail=F)
  T2DAG.rej = ifelse(T2DAG.pval<alpha,1,0)
  return(list(Qest = Qest, Rdiag = Rdiag,
              T2DAG.stat = T2DAG.stat, T2DAG.pval = T2DAG.pval, T2DAG.rej = T2DAG.rej))
}


#' Fit the SEM models used for estimating the coefficient matrix Q.
#' @param p.model the # of regressions to run, which is equal to the number of variables (genes)
#' that are children nodes of at least one other variable in the directed acyclic graph
#' that describes the pathway.
#' @param ntotal The total sample size n1 + n2
#' @param prefix The prefix of the column names in the data matrix 'Data'
#' @param Data A ntotal x p data frame, the first n1 rows store the standardized data
#' of the n1 samples from group 1, the second n1 rows store the standardized data
#' of the n2 samples from group 2.
#' @param A The binary adjacency matrix, where an entry aij=0 indicates no directed relationship
#' from variable i to variable j, and aij=1 indicates that there exists a directed relationship
#' from variable i to variable j.
#' @return A list that contains the following elements:
#' @param Qest The estimate of the coefficient matrix Q.
#' @param Rdiag The estimate of the covariance matrix of the residuals.
#' @export
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
