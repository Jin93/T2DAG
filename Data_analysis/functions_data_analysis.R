KEGG_testing = function(pathwayID){
  pathway.index = which(kegg.pathways == pathwayID)
  getKGMLurl(pathwayID)
  tmp <- paste0(pathway.dir,pathwayID,'.xml')
  if (!file.exists(tmp)) retrieveKGML(pathwayid=substr(pathwayID,4,8), organism='hsa', destfile=tmp, method="curl")
  pathway <- parseKGML(tmp)
  ### get gene names and indexes in the DAG
  nodes <- nodes(pathway)
  types = sapply(1:length(nodes),function(x){getType(nodes[[x]])})
  genes = sapply(1:length(nodes),function(x){substr(getName(nodes[[x]])[1],5,nchar(getName(nodes[[x]])[1]))})[types == 'gene'] # 279, only keep the first gene in each node
  genes.index = names(nodes)[(types == 'gene')]
  genes.index = genes.index[(genes %in% ge.genes)]
  genes = genes[genes %in% ge.genes]
  names(genes) = genes.index
  ###
  edges <- edges(pathway) # raw: 282 # if error, detach spatstat

  if (length(edges) > 0){
    edge.node = t(sapply(1:length(edges),function(x){getEntryID(edges[[x]])}))
    keep.indx = which(((edge.node[,1] %in% genes.index)&(edge.node[,2] %in% genes.index)) == T)
    edge.node = edge.node[keep.indx,] # filter out the nodes that are not genes

    ## get subtype name
    # https://www.kegg.jp/kegg/xml/docs/
    edge.subtype = lapply(1:length(edges),function(x){getSubtype(edges[[x]])})
    edge.subtype.length = sapply(1:length(edge.subtype),function(x){length(edge.subtype[[x]])})
    edge.types = vector('list',length(edge.subtype.length))
    tem.indx = which(edge.subtype.length == 0)
    if (length(tem.indx)>0){
      edge.types[-tem.indx] <- lapply(c(1:length(edges))[-tem.indx],function(x){getSubtype(edges[[x]])[[1]]})
      edge.type = rep(NA,length(edge.types))
      edge.type[-tem.indx] <- sapply(c(1:length(edges))[-tem.indx],function(x){getName(edge.types[[x]])})
    }
    if (length(tem.indx)==0){
      edge.types <- lapply(c(1:length(edges)),function(x){getSubtype(edges[[x]])[[1]]})
      edge.type = rep(NA,length(edge.types))
      edge.type <- sapply(c(1:length(edges)),function(x){getName(edge.types[[x]])})
    }
    edge.type = edge.type[keep.indx]

    #------------------------ AOAS
    keep.edge.types = c('activation','inhibition','expression','repression')
    edge.info = cbind(matrix(edge.node,ncol=2),matrix(edge.type,ncol=1))
    if (nrow(edge.info) > 0) edge.info = matrix(edge.info[complete.cases(edge.info),],ncol=3)
    if (nrow(edge.info) > 0) edge.info = matrix(edge.info[edge.info[,3] %in% keep.edge.types,],ncol=3)
    if (length(edge.info) == 3) edge.info = matrix(edge.info,nrow=1)
    if (nrow(edge.info) > 0){
      #### starting from here, we do not gene.index anymore
      #### combine duplicated genes:
      edge.info[,1] = genes[edge.info[,1]]
      edge.info[,2] = genes[edge.info[,2]]
      tem.idx = which(duplicated(edge.info))
      if (length(tem.idx)>0) edge.info = edge.info[-tem.idx,]
      genes = unique(genes)
      ## with self loops:
      A0 = matrix(0, length(genes), length(genes))
      ##### Z = QZ + E, Qij: j->i
      for (i in 1:nrow(edge.info)){
        A0[which(genes == edge.info[i,2]), which(genes == edge.info[i,1])] = ifelse(edge.info[i,3] %in% c('activation','expression'),1,-1)
      }

      #### remove loops:
      loops.index = which(edge.info[,1] == edge.info[,2])
      if (length(loops.index) > 1) edge.info = edge.info[-loops.index,]

      n.activation = sum(edge.info[,3] == 'activation')
      n.inhibition = sum(edge.info[,3] == 'inhibition')
      n.expression = sum(edge.info[,3] == 'expression')
      n.repression = sum(edge.info[,3] == 'repression')

      ##### Construct the adjacency matrix A for AOAS
      A = matrix(0, length(genes), length(genes))
      ##### Z = QZ + E, Qij: j->i
      for (i in 1:nrow(edge.info)){
        A[which(genes == edge.info[i,2]), which(genes == edge.info[i,1])] = 1
      }
      ### remove self loops
      diag(A) = 0
      abssum = function(x) sum(abs(x))
      nonzero.A = apply(A,1,abssum)
      d = max(nonzero.A)
      p = length(genes)
      p0 = sum(nonzero.A > 0)
      sparsity = sum(abs(A) > 0)/p^2

      ###
      a=adj.remove.cycles(A,maxlength=p)
      A = a$adjmat.acyclic
      n.circles = sum(a$adjmat.removed)


      ##### Gene expression data
      ge = GE[ge.genes %in% genes,]


      ##### Load patient info
      if (comparison == 'normal-tumor'){
        control.sample.ID = sampleinfo[sampleinfo[,'Sam_Tissue'] == 2, 'Sam_Name']
        ##### remove the patients that have both benigh tissue and tumor
        X = t(ge[,colnames(ge) %in% control.sample.ID]) # control
        Y = ge[,!(colnames(ge) %in% control.sample.ID)]
        Y = t(Y[,-1]) #remove the patient ID column
      }
      if (comparison == 'tumor-tumor'){
        #grp1
        #grp2
        ### only consider the 515 cases that are from separate patients
        case.patient.ID = sampleinfo[sampleinfo[,'Sam_Tissue'] == 1, 'Sam_Patient']
        case.sample.ID = sampleinfo[sampleinfo[,'Sam_Tissue'] == 1, 'Sam_Name']
        # remove the case of the dulicated patient IDs
        case.sample.ID = case.sample.ID[-which(duplicated(case.patient.ID))]
        # remove dulicated patient IDs
        case.patient.ID = case.patient.ID[-which(duplicated(case.patient.ID))]
        rownames(patientinfo) = patientinfo[,'Pat_ID']

        case.stage = patientinfo[case.patient.ID,'Pat_Stage']
        ###### no need to remove the patients that have missing tumor stage
        #missing.ind = which(is.na(case.stage))
        I.sample.ID = case.sample.ID[which(case.stage == grp1)]
        II.sample.ID = case.sample.ID[which(case.stage == grp2)]

        ##### remove the patients that have both benigh tissue and tumor
        X = t(ge[,colnames(ge) %in% I.sample.ID]) # control
        Y = t(ge[,(colnames(ge) %in% II.sample.ID)])
      }

      ##### Hypothesis testing
      n = list(X=nrow(X),Y=nrow(Y))
      nx = n[['X']]; ny = n[['Y']]
      N = n[['X']]+n[['Y']] ### in BS paper 'n' = N-2
      nmin = min(nx,ny) # minimum sample size
      n.confounder = 0
      graph.type = 'DAG'


      pval=rep(NA,length(methods)) # calculated p values
      rej=rep(NA,length(methods)) # 1: H0 rejected, 0: H0 not rejected
      names(pval) = names(rej) = methods

      Z = rbind(scale(X, center = T, scale = F),
                scale(Y, center = T, scale = F))
      Z0 = Z
      Z = as.data.frame(Z)
      colnames(Z) = c(paste0("Z", 1:(p+n.confounder)))#,paste0('M',1:n.confounder))
      meandiff = colMeans(X)-colMeans(Y)

      ####### 1. T2DAG
      SEM.results = SEM.fit(p+n.confounder,N,'Z', Z,A[1:(p+n.confounder),1:(p+n.confounder)])
      Qest = SEM.results$Qest
      Rdiag = SEM.results$Rdiag
      #### T= n1*n2/(n1+n2) (Xbar-Ybar)^T (I-Qhat) R^(_1)hat (I-Qhat) (Xbar-Ybar)
      Sigma.hat = solve(diag(1,p+n.confounder)-Qest)%*%diag(Rdiag)%*%solve(t(diag(1,p+n.confounder)-Qest))
      T.graph = (n[["X"]]*n[["Y"]])/(N)*(t(colMeans(X)-colMeans(Y))%*% solve(Sigma.hat[1:p,1:p]) %*%t(t(colMeans(X)-colMeans(Y))))

      ####### T2DAG test
      T.graph2 = (T.graph - p)/sqrt(2*p)
      methodname = 'T2DAG'
      pval[methodname] =  2*pnorm(q=abs(T.graph2),mean=0,sd=1,lower.tail=F)
      rej[methodname] = ifelse(pval[methodname]<alpha,1,0) # decision: 1 if reject, 0 o.w.

      ###################### Graph T2 test in Jacob et al. (2012): ######################
      Asym = A0[1:p,1:p]
      ncp <- 0.5
      sigma <- diag(p)/sqrt(p)
      ## Build graph, decompose laplacian
      lfA <- laplacianFromA(Asym)#,ltype="unnormalized"): note: this will create negative eigenvalues, but the results are the same
      k= 0.2 * p
      # Test statistic
      graph.T2.results <- graph.T2.test(X,Y,lfA=lfA,k=k,nmin=nmin)# T2.test(X,Y,k=k)
      pval["Graph.T2"] = graph.T2.results$p.value
      rej["Graph.T2"] = ifelse(graph.T2.results$p.value<alpha,1,0)

      ###################### Standard Hotelling's T2 test: ######################
      if (N > p){
        S.pool = ((n[["X"]]-1)*cov(X) + (n[["Y"]]-1)*cov(Y))/(N - 2)
        if (qr(S.pool)$rank >= p) #even if n1+n2>=p, sample covariance can still be singular
        {
          #### test statistic & p-value::
          T.h = (n[["X"]]*n[["Y"]])/(N)*(t(meandiff)%*%solve(S.pool, toler=1e-30)%*%t(t(meandiff)))
          pval["T2"] = pf(q=(N-p-1)/((N-2)*p)*T.h,df1=p,df2=(N-1-p),lower.tail=FALSE)
          rej["T2"] = ifelse(pval["T2"]<alpha,1,0)
        }
      }

      ###################### CHQ: ######################
      CH_Q = apval_Chen2010(sam1=X, sam2=Y, eq.cov = TRUE) # from aSPU package, faster than ChenQin.test(X,Y) from highD2pop
      pval["CH-Q"] = CH_Q$pval
      rej["CH-Q"] = ifelse(pval["CH-Q"]<alpha,1,0)

      ###################### SK: ######################
      SK = SK.test(X,Y)
      T.sk = SK$TSvalue
      pval["SK"] = SK$pvalue
      rej["SK"] = ifelse(SK$pvalue<alpha,1,0)

      ###################### CLX: ######################
      CLX = apval_Cai2014(sam1=X, sam2=Y) #from aSPU package, can also use CLX.test.equalcov(X,Y) from highD2pop
      pval["CLX"] = CLX$pval
      rej["CLX"] = ifelse(pval["CLX"]<alpha,1,0)

      ###################### GCT: ######################
      if (p >= 20){
        GCT = GCT.test(X,Y,r=ceiling(2/3*p^(1/2)), ntoorderminus = 0)#r=10) #in their paper r=L=10,15,20, smaller l leads to larger pvalue and smaller power
        # L = 2/3*p^(1/2), smaller L leads to more strict type I error control
        pval["GCT"] = GCT$pvalue
        rej["GCT"] = ifelse(GCT$pvalue<alpha,1,0)
      }
      if (p < 20){
        GCT = GCT.test(X,Y,r=ceiling(2/3*p^(1/2))+1, ntoorderminus = 0)#r=10) #in their paper r=L=10,15,20, smaller l leads to larger pvalue and smaller power
        # L = 2/3*p^(1/2), smaller L leads to more strict type I error control
        pval["GCT"] = GCT$pvalue
        rej["GCT"] = ifelse(GCT$pvalue<alpha,1,0)
      }

      ####################### aSPU: ######################
      if (p>=50){
        aspu = apval_aSPU(sam1=X,sam2=Y) #default: eq.cov = TRUE
        pval["aSPU"] = aspu$pval["aSPU"]
        rej["aSPU"] = ifelse(pval["aSPU"]<alpha,1,0)
        aspu.lambda = aspu$pow[which.min(aspu$pval)] # selected lambda
      }
      if (p<50){
        aspu = apval_aSPU(sam1=X,sam2=Y,bandwidth = 1) #default: eq.cov = TRUE
        pval["aSPU"] = aspu$pval["aSPU"]
        rej["aSPU"] = ifelse(pval["aSPU"]<alpha,1,0)
        aspu.lambda = aspu$pow[which.min(aspu$pval)] # selected lambda
      }
      basic.info = c(as.character(round(c(unlist(n),p,nrow(edge.info),length(loops.index),d,p0),0)),signif(sparsity,3),n.activation,n.inhibition,n.expression,n.repression,n.circles)
      names(basic.info) = c('n_x', 'n_y','p','n.edge','n.loop','d','p0','sparisty','n.activation','n.inhibition','n.expression','n.repression','n.circles')
      print(paste0('Pathway ', pathway.index, ': ', pathwayID,' Completed'))
      return(list(basic.info=basic.info,pval=pval))
    }
  }
}
