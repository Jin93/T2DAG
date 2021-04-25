#' T2-DAG test
#'
#' @description DAG-informed two-sample test for mean difference in a vector of gene expression levels of a pathway.
#'
#' @param X The n1 x p sample matrix from group 1. Each row corresponds to a p-dimensional sample.
#' @param Y The n2 x p sample matrix from group 2. Each row corresponds to a p-dimensional sample.
#' @param A The binary adjacency matrix, where an entry aij=0 indicates no directed relationship
#' from variable i to variable j, and aij=1 indicates that there exists a directed relationship
#' from variable i to variable j.
#' @return A list that contains the following elements:
#' \item{Qest}{The estimate of the coefficient matrix Q.}
#' \item{Rdiag}{The estimate of the covariance matrix of the residuals.}
#' \item{T2DAG.stat}{T2DAG test statistic.}
#' \item{T2DAG.pval}{P-value of the T2DAG test.}
#' @references
#' @export
#' @examples
#' data('pathway.list',package='T2DAG')
#' idx = 8
#' pathwayID = pathway.list[[1]][idx]
#' kegg.names = pathway.list[[2]][idx]
#' getKGMLurl(pathwayID)
#' tmp <- paste0(pathway.dir,pathwayID,'.xml')
#' if (!file.exists(tmp)) retrieveKGML(pathwayid=substr(pathwayID,4,8), organism='hsa', destfile=tmp, method="curl")
#' pathway <- parseKGML(tmp)
#'
#' ### Gene names and indices in the DAG
#' nodes <- nodes(pathway)
#' types = sapply(1:length(nodes),function(x){getType(nodes[[x]])})
#' genes = sapply(1:length(nodes),function(x){substr(getName(nodes[[x]])[1],5,nchar(getName(nodes[[x]])[1]))})[types == 'gene'] # 279, only keep the first gene in each node
#' genes.index = names(nodes)[(types == 'gene')]
#' genes.index = genes.index[(genes %in% ge.genes)] # extract indices of the genes that have expression data available
#' genes = genes[genes %in% ge.genes]
#' names(genes) = genes.index
#'
#' edges <- edges(pathway)
#' edge.node = t(sapply(1:length(edges),function(x){getEntryID(edges[[x]])}))
#' keep.indx = which(((edge.node[,1] %in% genes.index)&(edge.node[,2] %in% genes.index)) == T)
#' edge.node = edge.node[keep.indx,] # remove nodes that are not genes
#'
#' ## get subtype name
#' # https://www.kegg.jp/kegg/xml/docs/
#' edge.subtype = lapply(1:length(edges),function(x){getSubtype(edges[[x]])})
#' edge.subtype.length = sapply(1:length(edge.subtype),function(x){length(edge.subtype[[x]])})
#' edge.types = vector('list',length(edge.subtype.length))
#' tem.indx = which(edge.subtype.length == 0)
#' if (length(tem.indx)>0){
#'   edge.types[-tem.indx] <- lapply(c(1:length(edges))[-tem.indx],function(x){getSubtype(edges[[x]])[[1]]})
#'   edge.type = rep(NA,length(edge.types))
#'   edge.type[-tem.indx] <- sapply(c(1:length(edges))[-tem.indx],function(x){getName(edge.types[[x]])})
#' }
#' if (length(tem.indx)==0){
#'    edge.types <- lapply(c(1:length(edges)),function(x){getSubtype(edges[[x]])[[1]]})
#'    edge.type = rep(NA,length(edge.types))
#'    edge.type <- sapply(c(1:length(edges)),function(x){getName(edge.types[[x]])})
#' }
#' edge.type = edge.type[keep.indx]
#'
#' keep.edge.types = c('activation','inhibition','expression','repression')
#' edge.info = cbind(matrix(edge.node,ncol=2),matrix(edge.type,ncol=1))
#' if (nrow(edge.info) > 0) edge.info = matrix(edge.info[complete.cases(edge.info),],ncol=3)
#' if (nrow(edge.info) > 0) edge.info = matrix(edge.info[edge.info[,3] %in% keep.edge.types,],ncol=3)
#' if (length(edge.info) == 3) edge.info = matrix(edge.info,nrow=1)
#'
#' #Remove duplicated genes
#' edge.info[,1] = genes[edge.info[,1]]
#' edge.info[,2] = genes[edge.info[,2]]
#' tem.idx = which(duplicated(edge.info))
#' if (length(tem.idx)>0) edge.info = edge.info[-tem.idx,]
#' genes = unique(genes)
#' n.activation = sum(edge.info[,3] == 'activation'); n.inhibition = sum(edge.info[,3] == 'inhibition')
#' n.expression = sum(edge.info[,3] == 'expression'); n.repression = sum(edge.info[,3] == 'repression')
#'
#' # ------- Construct A, the binary adjacency matrix used for the T2DAG test -------
#' A = matrix(0, length(genes), length(genes))
#' #### Z = QZ + E, Qij: j->i
#' for (i in 1:nrow(edge.info)) A[which(genes == edge.info[i,2]), which(genes == edge.info[i,1])] = 1
#' diag(A) = 0
#' #### remove circles
#' a=adj.remove.cycles(A,maxlength=length(genes))
#' A = a$adjmat.acyclic
#'
#' #### Load example gene expression data.
#' data('X',package='T2DAG')
#' data('Y',package='T2DAG')
#'
#' # T2DAG test
#' T2DAG.results = T2DAG(X, Y, A)
#'
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
  return(list(Qest = Qest, Rdiag = Rdiag,
              T2DAG.stat = T2DAG.stat, T2DAG.pval = T2DAG.pval))
}
