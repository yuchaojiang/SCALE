gene_classify=function(alleleA,alleleB,epsilon=NULL){
  if(is.null(epsilon)){epsilon=0.001}
  gene.category=rep(NA,nrow(alleleA))
  A.prop=B.prop=rep(NA,nrow(alleleA))
  results.list=vector(mode = "list", length = nrow(alleleA))
  n=ncol(alleleA)
  genename=rownames(alleleA)
  for(i in 1:nrow(alleleA)){
    nA=alleleA[i,]
    nB=alleleB[i,]
    k=4
    a=b=cura=curb=1
    phi=rep(1,k)/k
    diff=1
    numiters=1
    while(diff>0.0001 || numiters <= 30){
      numiters=numiters+1
      Zhat=Estep(a, b, phi, epsilon, n, nA, nB)
      curM=Mstep(Zhat, cura, n, nA, nB)
      curphi=curM$phihat
      cura=curb=curM$a.b
      diff=max(max(abs(phi-curphi)),max(abs(b-curb)))
      phi=curphi
      a=b=cura
      if(numiters>=3000) break
    }
    
    cell.post.filter=apply(Zhat,1,max)>0.8
    category=c('Off','A','B','AB')
    results=table(category[apply(Zhat[cell.post.filter,],1,which.max)])
    results.list[[i]]=results
    
    if(setequal(names(results),'Off')){
      gene.cate='Silent'
    } else if (setequal(names(results),'A') | setequal(names(results),c('A','Off'))){
      gene.cate='MonoA'
    } else if (setequal(names(results),'B') | setequal(names(results),c('B','Off'))){
      gene.cate='MonoB'
    } else {
      gene.cate='Biallelic'
    }
    
    if(!('Off' %in% names(results))){results=c(results,Off=0)}
    if(!('AB' %in% names(results))){results=c(results,AB=0)}
    if(!('A' %in% names(results))){results=c(results,A=0)}
    if(!('B' %in% names(results))){results=c(results,B=0)}
    
    A.prop[i]=round((results['A']+results['AB'])/sum(results),3)
    B.prop[i]=round((results['B']+results['AB'])/sum(results),3)
    
    if(A.prop[i]>=0.05 & B.prop[i]>=0.05 & A.prop[i]<=0.95 & B.prop[i]<=0.95){
      gene.cate='Biallelic.bursty'
    }
    
    gene.category[i]=gene.cate
    cat('Gene',i,':',genename[i],',',gene.cate,'\t')
    cat('A prop',A.prop[i],'B prop',B.prop[i],'\n') 
    
  }
  gene.category[gene.category=='Biallelic']='Biallelic.nonbursty'
  gene.class.obj=list(gene.category=gene.category,A.prop=A.prop,
                      B.prop=B.prop,results.list=results.list)
  return(gene.class.obj)
}
