allelic_plot=function(alleleA,alleleB,gene.class.obj, allelic.kinetics.obj,
                      diff.allelic.obj,non.ind.obj,i){
  
  genename=rownames(alleleA)
  gene.category=gene.class.obj$gene.category
  A.prop=gene.class.obj$A.prop
  B.prop=gene.class.obj$B.prop
  results.list=gene.class.obj$results.list
    
  konA=allelic.kinetics.obj$konA
  koffA=allelic.kinetics.obj$koffA
  sA=allelic.kinetics.obj$sA
  konB=allelic.kinetics.obj$konB
  koffB=allelic.kinetics.obj$koffB
  sB=allelic.kinetics.obj$sB
  sizeA=sA/koffA
  sizeB=sB/koffB
  
  pval.kon=diff.allelic.obj$pval.kon
  pval.size=diff.allelic.obj$pval.size
  pval.ind=non.ind.obj$pval.ind
  non.ind.type=non.ind.obj$non.ind.type
    
  N=apply(alleleA+alleleB,2,sum)
  lib.size=N/mean(N)
  Ai=(alleleA[i,])/(lib.size)
  Bi=(alleleB[i,])/(lib.size)
  pdf(file=paste('SCALE_results_',genename[i],'.pdf',sep=''),width=5,height=8) 
  par(mfrow=c(2,1))
  par(mar=c(3,4,5,4))
  par(mgp=c(2,0.6,.4))
  barplot(Ai,ylim=c(-max(Bi),max(Ai)),col='red',names.arg='',xlab='Cell',ylab='Adjusted reads')
  barplot(-Bi,ylim=c(-max(Bi),max(Ai)),col='blue',add=TRUE,names.arg='')
  title(paste('Gene',genename[i]))
  
  par(mar=c(1,4,2,4))
  plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n",xaxt = "n", yaxt = "n")
  legend('topleft',col=c('red'),legend=c('A allele'),pch=15,bty='n')
  legend('topright',y=1,col=c('blue'),legend=c('B allele'),pch=15,bty='n')
  text(x=0.5,y=0.8,paste('Gene category :',gene.category[i]))
  text(x=0.5,y=0.65,paste('Number of cells:',paste(names(results.list[[1]]),results.list[[1]],collapse = '; ')))
  if(gene.category[i]=='Biallelic.bursty'){
    text(x=0.5,y=0.5,paste('konA =',round(konA[i],3),'sizeA =',round(sizeA[i],3)))
    text(x=0.5,y=0.4,paste('konB =',round(konB[i],3),'sizeB =',round(sizeB[i],3)))
    text(x=0.5,y=0.2,paste('Test of shared burst freq: pval =',round(pval.kon[i],4)))
    text(x=0.5,y=0.1,paste('Test of shared burst size: pval =',round(pval.size[i],4)))
    text(x=0.5,y=0,paste('Test of independent bursting: pval =',round(pval.ind[i],4)))
  }
  dev.off()
}