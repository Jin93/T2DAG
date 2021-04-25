setwd('~/T2DAG/')
# BiocManager::install("gage")
# BiocManager::install("gageData")
# BiocManager::install("KEGGprofile")
# BiocManager::install("DEGraph")
# BiocManager::install("predictionet")
# BiocManager::install("KEGGgraph")
library(T2DAG)
library(mvtnorm)
library(devtools)
library(MASS)
library(data.table)
library(readr)
library(stringr) # for str_detect
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(maptools)
library(ggExtra)
library(gtable)
library(gridExtra)
library(gbm)
library(gageData)
library(KEGGgraph)
library(readxl)
library(predictionet)
library(reshape2) # wide to long
library(corpcor) # campute partial correlation from correlation matrix
library(KEGGprofile)
library(gage)
library(KEGGlincs)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(highD2pop)  # for GCT.test - R package from GCT paper (2015)
library(highmean)
library(DEGraph) # for graph.T2.test
library(clipper) # for removing an edge from graphNEL
library(expm) # for sqrtm
library(gRbase)
library(rrcov)


graphT2.p = 0.8 # parameter used in graph T2 test
alpha = 0.05 # significance level for the hypothesis test.
methods = c("T2DAG","Graph.T2","T2","CH-Q","SK","CLX","GCT","aSPU")
#Load lung-tissue gene expression data.
ge.file = 'lung_cancer_gene_expression.xlsx'
ge.dir = 'data/gene_expression/'
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

# Load pathway information.
pathway.dir = 'data/pathways/'
data(kegg.sets.hs)
kegg.pathways1 = data.frame(id=substr(names(kegg.sets.hs),1,8), name=substr(names(kegg.sets.hs),10,1000)) # 229
#kegg.pathways2 = fread('data/pathways/kegg.hs.list.txt',header=F) # 337, 226 overlap with kegg.pathways1
data('kegg.hs.list',package='T2DAG')
kegg.pathways2 = kegg.hs.list # 337, 226 overlap with kegg.pathways1
kegg.pathways2 = as.data.frame(kegg.pathways2)
kegg.pathways2[,1] = substr(kegg.pathways2[,1], 6,13)
kegg.pathways2[,2] = sapply(1:nrow(kegg.pathways2),function(x){strsplit(kegg.pathways2[x,2], ' -')[[1]][1]})
colnames(kegg.pathways2) = c('id','name')
kegg.pathways.combine = kegg.pathways2
for (i in 1:nrow(kegg.pathways1)){
  if (!(kegg.pathways1[i,1] %in% kegg.pathways.combine[,1])){
    kegg.pathways.combine = rbind(kegg.pathways.combine,kegg.pathways1[i,1])
  }
}
kegg.pathways = kegg.pathways.combine[,1] # 340 pathways
kegg.names=kegg.pathways.combine[,2]
results = read.table(paste0(output.dir,'results_',dat,'_comparison_I_II.txt'), header = T)
kegg.names.subset = data.frame(pathways=rownames(results)) # 206
kegg.pathways.subset = kegg.pathways.combine[which(kegg.names %in% kegg.names.subset[,1]),]  # 210
kegg.pathways.subset = kegg.pathways.subset[-c(32,39,40,55),] # remove duplicated pathways
# remove duplicate
#fwrite(kegg.pathways.subset, 'data/pathway.list.txt', sep=' ', col.names=T)
save(kegg.pathways.subset, 'data/pathway.list.RData')



comparison = 'tumor-tumor'
grps = rbind(c('I','II'),c('I','III'),c('I','IV'),c('II','III'),c('II','IV'),c('III','IV'))
grps_ns = rbind('I','II','III','IV')
for (G in 1:nrow(grps)){
  if (comparison == 'tumor-tumor'){
    grp1 = grps[G,1]
    grp2 = grps[G,2]
  }
  if (comparison == 'normal-stage'){
    grp1 = 'normal'
    grp2 = grp_ns[G,1]
  }
  kegg.results = list()
  for (i in 1:length(kegg.pathways)) kegg.results[[i]] = KEGG_testing(kegg.pathways[i])

  ##### Extract results
  basic.info = pval = NULL
  l = 1
  for (i in 1:length(kegg.results)){
    if (length(kegg.results[[i]]) > 0){
      basic.info = rbind(basic.info,as.numeric(kegg.results[[i]]$basic.info))
      basic.info.colnames = names(kegg.results[[i]]$basic.info)
      pval = rbind(pval, kegg.results[[i]]$pval)
      rownames(basic.info)[l] = rownames(pval)[l] = kegg.names[i]
      l = l + 1
    }
  }
  colnames(basic.info) = basic.info.colnames
  pval = signif(pval,2)
  results = cbind(basic.info,pval)
  write.table(results,file=paste0(out.dir,'results_',dat,'_comparison_',grp1,'_',grp2,'.txt'),col.names = T,row.names = T)
}
