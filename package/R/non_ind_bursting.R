non_ind_bursting=function(alleleA, alleleB, gene.category, results.list){
  bursty.gene= which(gene.category=='Biallelic.bursty')
  pval.ind=non.ind.type=rep(NA,nrow(alleleA))
  prop=matrix(0,ncol=4,nrow=nrow(alleleA))
  colnames(prop)=c('A','B','AB','Off')
  for(i in 1:length(results.list)){
    results.temp=results.list[[i]]
    prop.temp=round(results.temp/sum(results.temp),4)
    prop[i,names(prop.temp)]=prop.temp
  }
  
  for(i in bursty.gene){
    n.cells=sum(results.list[[i]])
    obs=rep(0,4)
    names(obs)=c('Off','B','A','AB')
    obs[names(results.list[[i]])]=results.list[[i]]
    
    pA=prop[i,2]+prop[i,4]
    pB=prop[i,1]+prop[i,4]
    exp=n.cells*c(pA*pB,pA*(1-pB),(1-pA)*pB,(1-pA)*(1-pB))
    names(exp)=c('Off','B','A','AB')
    pval.ind[i]=1-pchisq(sum((obs-exp)^2/exp),df=1)
    if(exp['AB']<obs['AB']){non.ind.type[i]='C'}
    if((exp['A']+exp['B'])<(obs['A']+obs['B'])){non.ind.type[i]='R'}
  }
  non.ind.obj=list(pval.ind=pval.ind,non.ind.type=non.ind.type)
  return(non.ind.obj)
}