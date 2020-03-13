allelic_kinetics=function(alleleA,alleleB,abkt,gene.category,cellsize,pdf=NULL){
  if(is.null(pdf)){pdf=TRUE}
  bursty.gene= which(gene.category=='Biallelic.bursty')
  genename=rownames(alleleA)
  alpha=abkt[1]
  beta=abkt[2]
  kappa=abkt[3]
  tau=abkt[4]
  N=apply(alleleA+alleleB,2,sum)
  lib.size=N/mean(N)

  bandwidth_output=matrix(ncol=4,nrow=10)
  colnames(bandwidth_output)=c('Bandwidth',
                               'Percent_non_neg','Cor_freq','Cor_size')
  num.cells=ncol(alleleA)
  
  # firstly determine the optimal bandwidth for the histogram repiling
  for(bandwidth in 1:10){
    cat('Bandwidth ',bandwidth,':\t')
    konA=koffA=sA=rep(NA,length(genename))
    konB=koffB=sB=rep(NA,length(genename))
    
    for(i in bursty.gene){
      # allele A
      Ai=alleleA[i,]/lib.size 
      temp = as.data.frame(table(pmax(0,ceiling(Ai/bandwidth) * bandwidth-bandwidth/2)))
      q.temp=as.numeric(as.matrix(temp[,1])) # binned expression
      c.temp=as.numeric(as.matrix(temp[,2])) # number of cells within each bin
      y.temp=round(exp((log(q.temp)-alpha)/beta)) # adjusted expression for each bin
      p.temp=expit(kappa+tau*log(pmax(0.001,y.temp))) # probability of dropout for each bin
      nc.temp=round(c.temp/p.temp) # adjusted number of cells within each adjusted expression bin
      if(sum(nc.temp[-1])>num.cells){
        loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
        nc.temp[loc.temp]=0
        nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
      } else{
        nc.temp[1]=num.cells-sum(nc.temp[-1])
      }
      pb.temp=pb.moment(rep(round(y.temp),nc.temp),cellsize)
      konA[i]=pb.temp[1]
      koffA[i]=pb.temp[2]
      sA[i]=pb.temp[3]
      
      # allele B
      Bi=alleleB[i,]/lib.size
      temp = as.data.frame(table(pmax(0,ceiling(Bi/bandwidth) * bandwidth-bandwidth/2)))
      c.temp=as.numeric(as.matrix(temp[,2]))
      q.temp=as.numeric(as.matrix(temp[,1]))
      y.temp=round(exp((log(q.temp)-alpha)/beta))
      p.temp=expit(kappa+tau*log(pmax(0.001,y.temp))) 
      nc.temp=round(c.temp/p.temp)
      if(sum(nc.temp[-1])>num.cells){
        loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
        nc.temp[loc.temp]=0
        nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
      } else{
        nc.temp[1]=num.cells-sum(nc.temp[-1])
      }
      pb.temp=pb.moment(rep(round(y.temp),nc.temp),cellsize)
      konB[i]=pb.temp[1]
      koffB[i]=pb.temp[2]
      sB[i]=pb.temp[3]
    }
    
    sel1=(sA>0)&(sB>0)&(koffA>0)&(koffB>0)
    sel2=konA>0 & konB>0
    sel=sel1&sel2
    
    percent1=sum(sel[bursty.gene])/length(bursty.gene)  # percentage of non-negative estimated kon, koff, s
    cor1.1=cor(log(konA[intersect(which(sel),bursty.gene)]),log(konB[intersect(which(sel),bursty.gene)]))
    cor1.2=cor(log((sA/koffA)[intersect(which(sel),bursty.gene)]),log((sB/koffB)[intersect(which(sel),bursty.gene)]))
    cat('% non-neg estimates',round(percent1,3),
        '\t corr. freq',round(cor1.1,3),
        '\t corr. size',round(cor1.2,3),'\n')
    bandwidth_output[bandwidth,]=c(bandwidth,round(c(percent1,cor1.1,cor1.2),3))
  }

  # now choose the optimal bandwidth and run the analysis again.
  bandwidth=which.max(rank(bandwidth_output[,2])+rank(bandwidth_output[,3])+rank(bandwidth_output[,4]))
  konA=koffA=sA=rep(NA,length(genename))
  konB=koffB=sB=rep(NA,length(genename))
  
  for(i in bursty.gene){
    # allele A
    Ai=alleleA[i,]/lib.size
    temp = as.data.frame(table(pmax(0,ceiling(Ai/bandwidth) * bandwidth-bandwidth/2)))
    q.temp=as.numeric(as.matrix(temp[,1])) # binned expression
    c.temp=as.numeric(as.matrix(temp[,2])) # number of cells within each bin
    y.temp=round(exp((log(q.temp)-alpha)/beta)) # adjusted expression for each bin
    p.temp=expit(kappa+tau*log(pmax(0.001,y.temp)))  # probability of dropout for each bin
    nc.temp=round(c.temp/p.temp) # adjusted number of cells within each adjusted expression bin
    if(sum(nc.temp[-1])>num.cells){
      loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
      nc.temp[loc.temp]=0
      nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
    } else{
      nc.temp[1]=num.cells-sum(nc.temp[-1])
    }
    pb.temp=pb.moment(rep(round(y.temp),nc.temp),cellsize)
    konA[i]=pb.temp[1]
    koffA[i]=pb.temp[2]
    sA[i]=pb.temp[3]
    
    # allele B
    Bi=alleleB[i,]/lib.size
    temp = as.data.frame(table(pmax(0,ceiling(Bi/bandwidth) * bandwidth-bandwidth/2)))
    c.temp=as.numeric(as.matrix(temp[,2]))
    q.temp=as.numeric(as.matrix(temp[,1]))
    y.temp=round(exp((log(q.temp)-alpha)/beta))
    p.temp=expit(kappa+tau*log(pmax(0.001,y.temp)))
    nc.temp=round(c.temp/p.temp)
    if(sum(nc.temp[-1])>num.cells){
      loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
      nc.temp[loc.temp]=0
      nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
    } else{
      nc.temp[1]=num.cells-sum(nc.temp[-1])
    }
    pb.temp=pb.moment(rep(round(y.temp),nc.temp),cellsize)
    konB[i]=pb.temp[1]
    koffB[i]=pb.temp[2]
    sB[i]=pb.temp[3]
  }
  
  sel1=(sA>0)&(sB>0)&(koffA>0)&(koffB>0)
  sel2=konA>0 & konB>0
  sel=sel1&sel2
  
  percent1=sum(sel[bursty.gene])/length(bursty.gene)  # percentage of non-negative estimated kon, koff, s
  cor1.1=cor(log(konA[intersect(which(sel),bursty.gene)]),log(konB[intersect(which(sel),bursty.gene)]))
  cor1.2=cor(log((sA/koffA)[intersect(which(sel),bursty.gene)]),log((sB/koffB)[intersect(which(sel),bursty.gene)]))
  
  
  if(pdf == TRUE){
      pdf(file='corr_s_koff.pdf',width=8,height=4.5)
      par(mfrow=c(1,2))
      konA.temp=log((konA)[intersect(which(sel),bursty.gene)])
      konB.temp=log((konB)[intersect(which(sel),bursty.gene)])
      plot(konA.temp,konB.temp,pch=16,cex=0.4,
           xlim=c(min(konA.temp,konB.temp),max(konA.temp,konB.temp)),
           ylim=c(min(konA.temp,konB.temp),max(konA.temp,konB.temp)),
           xlab='log(konA)',ylab='log(konB)')
      grid()
      abline(a=0,b=1,lty=2,lwd=1.5,col='blue')
      legend('topleft',paste('r =',round(cor1.1,3)),bty='n')
      
      sizeA.temp=log((sA/koffA)[intersect(which(sel),bursty.gene)])
      sizeB.temp=log((sB/koffB)[intersect(which(sel),bursty.gene)])
      plot(sizeA.temp,sizeB.temp,pch=16,cex=0.4,
           xlim=c(min(sizeA.temp,sizeB.temp),max(sizeA.temp,sizeB.temp)),
           ylim=c(min(sizeA.temp,sizeB.temp),max(sizeA.temp,sizeB.temp)),
           xlab='log(sA/koffA)',ylab='log(sB/koffB)')
      grid()
      abline(a=0,b=1,lty=2,lwd=1.5,col='blue')
      legend('topleft',paste('r =',round(cor1.2,3)),bty='n')
      dev.off()
  }
 
  allelic.kinetics.obj=list(konA=konA,koffA=koffA,sA=sA,
                        konB=konB,koffB=koffB,sB=sB,
                        bandwidth=bandwidth)
  return(allelic.kinetics.obj)
}
