# T2DAG

R package for T2-DAG, a DAG-informed two-sample test for mean difference in the vector of gene expression levels of a pathway. In addition to gene expression data, the method efficiently leveraging axiliary pathway information on gene interactions through a linear structural equation model. 


## Installation


### Install and load the following R packages
```r
BiocManager::install("gage") # for KEGG pathway analysis
BiocManager::install("gageData") # for KEGG pathway analysis
BiocManager::install("RCy3")
BiocManager::install("KEGGlincs") 
BiocManager::install("KEGGprofile") 
BiocManager::install("KEGGgraph") 
BiocManager::install("DEGraph")
BiocManager::install("predictionet")
BiocManager::install("clipper")
# if installation of clipper is unsuccessful: install an old version of gRbase so as to install the package clipper
# install.packages('~/Downloads/gRbase_1.8-3.tar.gz', repos = NULL, type ='source')

library(gage)
library(gageData)
library(KEGGlincs)
library(KEGGprofile)
library(KEGGgraph)
library(DEGraph) # for graph.T2.test
library(predictionet)
library(clipper) # for removing an edge from graphNEL
```

### Install package "T2DAG"
``` r
devtools::install_github("Jin93/T2DAG")
library(T2DAG)
# download T2DAG github repository and set path to the folder "T2DAG-main"
setwd('~/T2DAG-main/') # set path to the Github directory
```

## Example data analysis:
Testing the mean difference in the expression levels of genes in a KEGG pathway between stage I and stage II lung cancer.


### Step 1: data preparation

#### Data sources (All GWAS samples are of European ancestry):
  1. Gene expression data collected from lung tissues of different lung cancer stages (normal, stage I, II, III, and IV) for lung cancer patients obtained from the Cancer Genome Atlas (TCGA) Program<sup>1,2</sup>.  
  2. KEGG pathways<sup>3,4,5</sup>.

Load additional packages needed for data analysis
```r
library(mvtnorm)
library(devtools)
library(MASS)
library(data.table)
library(readr)
library(stringr) # for str_detect
library(dplyr)
library(rrcov)
library(reshape2) # wide to long
library(corpcor) # campute partial correlation from correlation matrix
library(highD2pop)  # for GCT.test - R package from GCT paper (2015)
library(highmean)
library(expm) # for sqrtm
library(readxl)
```

Download the file "lung_cancer_gene_expression" from [this link](https://www.dropbox.com/scl/fi/bb8nco5y3dlked6dckjmo/lung_cancer_gene_expression.xlsx?dl=0&rlkey=kc9ec0bq60qjysxw6tfi40pwf). If the link does not work, please contact jjin31@jhu.edu to request for the file. 

Save the file to the directory "data/" and load lung-tissue gene expression data.

```r
ge.file = 'lung_cancer_gene_expression.xlsx'
ge.dir = 'data/'
output.dir = 'data/gene_expression/'
GE = read_xlsx(paste0(ge.dir,ge.file),sheet = 'Expression', col_names = T, col_types = c('text',rep('numeric',576)))
ge.genes = GE[,1][[1]]
GE = as.data.frame(GE)
colnames(GE)[2:ncol(GE)] = gsub('\\.','-',colnames(GE)[2:ncol(GE)])
# Patient info
patientinfo = read_xlsx(paste0(ge.dir,ge.file),sheet = 'Patients', col_names = T)
patientinfo = as.data.frame(patientinfo)
patientinfo = patientinfo[!(is.na(patientinfo$Pat_Stage)),]
rownames(patientinfo) = patientinfo[,'Pat_ID']

# Sample info
sampleinfo = read_xlsx(paste0(ge.dir,ge.file),sheet = 'Samples', col_names = T)
sampleinfo = as.data.frame(sampleinfo)
rownames(sampleinfo) = sampleinfo[,'Sam_Name']
sampleinfo = sampleinfo[sampleinfo[,'Sam_Patient'] %in% rownames(patientinfo),]
```

Load pathway information.
```r
data('pathway.list',package='T2DAG')
kegg.pathways = pathway.list[[1]]
kegg.names = pathway.list[[2]]
```


### Now we show an example of the analysis on the pathway "p53 signaling pathway"

Prepare the corresponding pathway information and gene expression data
```r
set.seed(1234)
pathway.index = 28
pathway.dir = 'data/'
pathwayID = kegg.pathways[pathway.index]
getKGMLurl(pathwayID)
tmp <- paste0(pathway.dir,pathwayID,'.xml')
if (!file.exists(tmp)) retrieveKGML(pathwayid=substr(pathwayID,4,8), organism='hsa', destfile=tmp, method="curl")
pathway <- parseKGML(tmp)
```


Summarize node and edge information of the pathway
```r
### Gene names and indices in the DAG
nodes <- nodes(pathway)
types = sapply(1:length(nodes),function(x){getType(nodes[[x]])})
genes = sapply(1:length(nodes),function(x){substr(getName(nodes[[x]])[1],5,nchar(getName(nodes[[x]])[1]))})[types == 'gene'] # 279, only keep the first gene in each node
genes.index = names(nodes)[(types == 'gene')]
genes.index = genes.index[(genes %in% ge.genes)] # extract indices of the genes that have expression data available
genes = genes[genes %in% ge.genes]
names(genes) = genes.index
###
edges <- edges(pathway)

edge.node = t(sapply(1:length(edges),function(x){getEntryID(edges[[x]])}))
keep.indx = which(((edge.node[,1] %in% genes.index)&(edge.node[,2] %in% genes.index)) == T)
edge.node = edge.node[keep.indx,] # remove nodes that are not genes

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

keep.edge.types = c('activation','inhibition','expression','repression')
edge.info = cbind(matrix(edge.node,ncol=2),matrix(edge.type,ncol=1))
if (nrow(edge.info) > 0) edge.info = matrix(edge.info[complete.cases(edge.info),],ncol=3)
if (nrow(edge.info) > 0) edge.info = matrix(edge.info[edge.info[,3] %in% keep.edge.types,],ncol=3)
if (length(edge.info) == 3) edge.info = matrix(edge.info,nrow=1)

#Remove duplicated genes
edge.info[,1] = genes[edge.info[,1]]
edge.info[,2] = genes[edge.info[,2]]
tem.idx = which(duplicated(edge.info))
if (length(tem.idx)>0){
  edge.info = edge.info[-tem.idx,]
}
genes = unique(genes)
n.activation = sum(edge.info[,3] == 'activation')
n.inhibition = sum(edge.info[,3] == 'inhibition')
n.expression = sum(edge.info[,3] == 'expression')
n.repression = sum(edge.info[,3] == 'repression')
```

Construct adjacency matrices that will be used in the Graph T2 test and T2DAG test.
```r
# ------- Construct A0, the adjacency matrix for Graph.T2 test -------
# which does not necessarily corresponds to a DAG 
# and can additionally incorporate the direction of gene interactions
## A0(i,j) = 1 indicates activation/expression of gene j by gene i; 
## A0(i,j) = -1 indicates inhibition/repression of gene j by gene i; 
## A0(i,j) = 0 indicates no effect of the expression of gene i on the expression of gene j.
A0 = matrix(0, length(genes), length(genes))
##### Z = QZ + E, Qij: j->i
for (i in 1:nrow(edge.info)){
  A0[which(genes == edge.info[i,2]), which(genes == edge.info[i,1])] = ifelse(edge.info[i,3] %in% c('activation','expression'),1,-1)
}

# ------- Construct A, the binary adjacency matrix used for the T2DAG test -------
## A(i,j) = 1 indicates the existence of effect of the expression of gene i on the expression of gene j; 
## A0(i,j) = 0 indicates no effect of the expression of gene i on the expression of gene j.
#### number of self loops:
loops.index = which(edge.info[,1] == edge.info[,2])
if (length(loops.index) > 1) edge.info = edge.info[-loops.index,]
n.loops = length(loops.index)

A = matrix(0, length(genes), length(genes))
#### Z = QZ + E, Qij: j->i
for (i in 1:nrow(edge.info)) A[which(genes == edge.info[i,2]), which(genes == edge.info[i,1])] = 1
diag(A) = 0

#### remove circles
a=adj.remove.cycles(A,maxlength=length(genes))
A = a$adjmat.acyclic

#### Finally, summaries of the DAG
n.circles = sum(a$adjmat.removed)
abssum = function(x) sum(abs(x))
nonzero.A = apply(A,1,abssum)
d = max(nonzero.A)
p = length(genes)
p0 = sum(nonzero.A > 0)
sparsity = sum(abs(A) > 0)/p^2
```

Extract gene expression data
```r
comparison = 'tumor-tumor'; 
grp1 = 'I'; grp2 = 'II' # stages of lung cancer to compare
# 'tumor-tumor': comparison between tumor tissues of two lung cancer stages
# 'normal-tumor': comparison between normal tissue and tumor tissue of one lung cancer stage

#### Gene expression data.
ge = GE[ge.genes %in% genes,]

##### Patient information
if (comparison == 'normal-tumor'){
  control.sample.ID = sampleinfo[sampleinfo[,'Sam_Tissue'] == 2, 'Sam_Name']
  ##### remove the patients that have both benigh tissue and tumor to avoid between-sample correlation
  X = t(ge[,colnames(ge) %in% control.sample.ID]) # control
  Y = ge[,!(colnames(ge) %in% control.sample.ID)]
  Y = t(Y[,-1]) #remove the patient ID column
}
if (comparison == 'tumor-tumor'){
  ### only consider the tumor tissues that are from separate patients
  case.patient.ID = sampleinfo[sampleinfo[,'Sam_Tissue'] == 1, 'Sam_Patient']
  case.sample.ID = sampleinfo[sampleinfo[,'Sam_Tissue'] == 1, 'Sam_Name']
  # remove the cases of the dulicated patient IDs
  case.sample.ID = case.sample.ID[-which(duplicated(case.patient.ID))]
  # remove dulicated patient IDs
  case.patient.ID = case.patient.ID[-which(duplicated(case.patient.ID))]
  rownames(patientinfo) = patientinfo[,'Pat_ID']
  
  case.stage = patientinfo[case.patient.ID,'Pat_Stage']
  I.sample.ID = case.sample.ID[which(case.stage == grp1)]
  II.sample.ID = case.sample.ID[which(case.stage == grp2)]
  
  ##### remove patients that have both benigh tissue and tumor
  X = t(ge[,colnames(ge) %in% I.sample.ID]) # control
  Y = t(ge[,(colnames(ge) %in% II.sample.ID)])
}
```

### Step 2. Conduct hypothesis tests

```r
alpha = 0.05 # significance level for the hypothesis test.
methods = c("T2DAG","Graph.T2","T2","CH-Q","SK","CLX","GCT","aSPU")
pval=rep(NA,length(methods)) # a vector that stores p values of the various tests
rej=rep(NA,length(methods)) # a vector that stores conclusion of the various tests
# 1: H0 rejected, 0: H0 not rejected
names(pval) = names(rej) = methods
```

#### 1. T2-DAG test<sup>6</sup>
```r
T2DAG.results = T2DAG(X, Y, A)
pval['T2DAG'] =  T2DAG.results$T2DAG.pval
rej['T2DAG'] = ifelse(T2DAG.results$T2DAG.pval<alpha,1,0)
```

#### 2. Graph T2 test<sup>7</sup>
```r
Asym = A0[1:p,1:p]
ncp <- 0.5
sigma <- diag(p)/sqrt(p)
## Build graph, decompose laplacian
lfA <- laplacianFromA(Asym)
k= 0.2 * p
# Test statistic
graph.T2.results <- graph.T2.test(X,Y,lfA=lfA,k=k,nmin=nmin)# T2.test(X,Y,k=k)
pval["Graph.T2"] = graph.T2.results$p.value
rej["Graph.T2"] = ifelse(graph.T2.results$p.value<alpha,1,0)
```

#### 3. Hotelling's T2 test<sup>8</sup>
```r
p = nrow(A) # number of genes in the pathway
n = c(nrow(X),nrow(Y)) # sample sizes
N = sum(n) # total sample size
nmin = min(n[1],n[2]) # minimum of the two sample sizes
meandiff = colMeans(X)-colMeans(Y)

if (N > p){
  S.pool = ((n[1]-1)*cov(X) + (n[2]-1)*cov(Y))/(N - 2)
  if (qr(S.pool)$rank >= p){ #even if n1+n2>=p, sample covariance can still be singular
    #### test statistic & p-value:
    T.h = (n[1]*n[2])/(N)*(t(meandiff)%*%solve(S.pool, toler=1e-30)%*%t(t(meandiff)))
    pval["T2"] = pf(q=(N-p-1)/((N-2)*p)*T.h,df1=p,df2=(N-1-p),lower.tail=FALSE)
    rej["T2"] = ifelse(pval["T2"]<alpha,1,0)
  }
}
```
#### 4. CH-Q test<sup>9</sup>
```r
CH_Q = apval_Chen2010(sam1=X, sam2=Y, eq.cov = TRUE) # This function is from the package "highmean". Can also use the function ChenQin.test() from the package "highD2pop".
pval["CH-Q"] = CH_Q$pval
rej["CH-Q"] = ifelse(pval["CH-Q"]<alpha,1,0)
```

#### 5. SK test<sup>10</sup>
```r
SK.results = SK.test(X,Y)
T.sk = SK.results$TSvalue
pval["SK"] = SK.results$pvalue
rej["SK"] = ifelse(SK.results$pvalue<alpha,1,0)
```

#### 6. CLX test<sup>11</sup>
```r
CLX.results = apval_Cai2014(sam1=X, sam2=Y) # This function is from the package "highmean". Can also use CLX.test.equalcov() from the package "highD2pop".
pval["CLX"] = CLX.results$pval
rej["CLX"] = ifelse(pval["CLX"]<alpha,1,0)
```

#### 7. GCT test<sup>12</sup>
```r
if (p >= 20){
  GCT.results = GCT.test(X,Y,r=ceiling(2/3*p^(1/2)), ntoorderminus = 0) # In their paper r=L=10,15,20, smaller l leads to larger pvalue and smaller power.
  pval["GCT"] = GCT.results$pvalue
  rej["GCT"] = ifelse(GCT.results$pvalue<alpha,1,0)
}
if (p < 20){
  GCT.results = GCT.test(X,Y,r=ceiling(2/3*p^(1/2))+1, ntoorderminus = 0)#r=10)
  pval["GCT"] = GCT.results$pvalue
  rej["GCT"] = ifelse(GCT.results$pvalue<alpha,1,0)
}
```

#### 8. aSPU test<sup>13</sup>
```r
if (p>=50){
  aspu = apval_aSPU(sam1=X,sam2=Y) # default: eq.cov = TRUE
  pval["aSPU"] = aspu$pval["aSPU"]
  rej["aSPU"] = ifelse(pval["aSPU"]<alpha,1,0)
  aspu.lambda = aspu$pow[which.min(aspu$pval)] # selected lambda
}
if (p<50){
  aspu = apval_aSPU(sam1=X,sam2=Y,bandwidth = 1)
  pval["aSPU"] = aspu$pval["aSPU"]
  rej["aSPU"] = ifelse(pval["aSPU"]<alpha,1,0)
  aspu.lambda = aspu$pow[which.min(aspu$pval)] # selected lambda
}
```



### Hypothesis testing results
```r
output = list(basic.info = data.frame(n1=n[1], n2=n[2], p=p, n.edges=nrow(edge.info),
                                      n.loops = n.loops, d=d, p0=p0, sparsity = signif(sparsity,2),
                                      n.activation = n.activation, n.inhibition = n.inhibition,
                                      n.expression = n.expression, n.repression = n.repression,
                                      n.circles = n.circles),
              test.results = signif(pval,2))
# output
$basic.info
   n1  n2  p n.edges n.loops d p0 sparsity n.activation n.inhibition
  278 124 60      56       0 6 44    0.015           12            7
  n.expression n.repression n.circles
            37            0         0

$test.results
   T2DAG Graph.T2       T2     CH-Q       SK      CLX      GCT     aSPU 
 4.2e-11  2.3e-03  6.7e-03  3.6e-04  3.4e-03  1.9e-02  8.9e-04  3.0e-03
```




##  References
  1. Cancer Genome Atlas Research Network, 2014. Comprehensive molecular profiling of lung adenocarcinoma. Nature, 511(7511), p.543.
  2. Cai, L., Lin, S., Girard, L., Zhou, Y., Yang, L., Ci, B., Zhou, Q., Luo, D., Yao, B., Tang, H. and Allen, J., 2019. LCE: an open web portal to explore gene expression and clinical associations in lung cancer. Oncogene, 38(14), pp.2551-2564.
  3. KEGG pathways. 
      [ftp://ftp.genome.ad.jp/pub/kegg/pathways](ftp://ftp.genome.ad.jp/pub/kegg/pathways)
  4. Kanehisa, M. and Goto, S., 2000. KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), pp.27-30.
  5. Kanehisa, M., Sato, Y., Furumichi, M., Morishima, K. and Tanabe, M., 2019. New approach for understanding genome variations in KEGG. Nucleic acids research, 47(D1), pp.D590-D595.
  6. A powerful test for differentially expressed gene pathways via graph-informed structural equation modeling.
  7. Jacob, L., Neuvial, P. and Dudoit, S., 2012. More power via graph-structured tests for differential expression of gene networks. The Annals of Applied Statistics, pp.561-600.
  8. Hotelling, H., 1992. The generalization of Student’s ratio. In Breakthroughs in statistics (pp. 54-65). Springer, New York, NY.
  9. Chen, S.X. and Qin, Y.L., 2010. A two-sample test for high-dimensional data with applications to gene-set testing. The Annals of Statistics, 38(2), pp.808-835.
  10. Srivastava, M.S. and Kubokawa, T., 2013. Tests for multivariate analysis of variance in high dimension under non-normality. Journal of Multivariate Analysis, 115, pp.204-216.
  11. Cai, T.T., Liu, W. and Xia, Y., 2014. Two-sample test of high dimensional means under dependence. Journal of the Royal Statistical Society: Series B: Statistical Methodology, pp.349-372.
  12. Gregory, K.B., Carroll, R.J., Baladandayuthapani, V. and Lahiri, S.N., 2015. A two-sample test for equality of means in high dimension. Journal of the American Statistical Association, 110(510), pp.837-849.
  13. Xu, G., Lin, L., Wei, P. and Pan, W., 2016. An adaptive two-sample test for high-dimensional means. Biometrika, 103(3), pp.609-624.


##  Author Information

* Jin Jin  Department of Biostatistics, Johns Hopkins Bloomberg School of Public Health jjin31@jhu.edu
* Yue Wang  School of Mathematical and Natural Sciences, Arisona State University Yue.Wang.Stat@asu.edu
