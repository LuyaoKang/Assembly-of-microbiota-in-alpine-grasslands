# using the normalized stochasticity ratio (NST) to quantified the relative importance of stochastic processes (Ning et al., 2019)
library(NST)
NST<-function(otu_table){
  group<-data.frame(c(rep('all',nrow(otu_table))))
  rownames(group)<-rownames(otu_table)
  tnst.dat=tNST(comm=otu_table, group=group, dist.method="jaccard",
                abundance.weighted=TRUE, rand=100,output.rand = T,
                nworker=4, null.model="PF", between.group=F,
                SES=T, RC=T)
  nst.bt=nst.boot(nst.result=tnst.dat, group=NULL, rand=99,
                  trace=TRUE, two.tail=FALSE, out.detail=T,
                  between.group=FALSE, nworker=1)
  NST<-nst.bt$detail$NST.boot$all
  return(NST)
}
