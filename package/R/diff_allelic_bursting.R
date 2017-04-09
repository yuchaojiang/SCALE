diff_allelic_bursting=function(alleleA, alleleB, cellsize, gene.category,
                               abkt, allelic.kinetics.obj, 
                               nboot=NULL, mode=NULL){
  if(is.null(mode)){mode='corrected'}
  if(mode!='corrected' & mode !='raw') stop("Mode has to be either 'raw' or 'corrected'!")
  if(is.null(nboot)){nboot=10000}
    bandwidth=allelic.kinetics.obj$bandwidth
    bursty.gene= which(gene.category=='Biallelic.bursty')
    N=apply(alleleA+alleleB,2,sum)
    lib.size=N/mean(N)
    alpha=abkt[1]
    beta=abkt[2]
    kappa=abkt[3]
    tau=abkt[4]
    konA=allelic.kinetics.obj$konA
    koffA=allelic.kinetics.obj$koffA
    sA=allelic.kinetics.obj$sA
    konB=allelic.kinetics.obj$konB
    koffB=allelic.kinetics.obj$koffB
    sB=allelic.kinetics.obj$sB
    sizeA=sA/koffA
    sizeB=sB/koffB
    num.cells=ncol(alleleA)
    pval.kon=pval.size=rep(NA,nrow(alleleA))
    
    sel1=(sA>0)&(sB>0)&(koffA>0)&(koffB>0)
    sel2=konA>0 & konB>0
    sel=sel1&sel2
    sum(sel[bursty.gene])
    
    for(i in which(sel)){
        cat('Testing gene',i,': ')
        konA.samp=rep(NA,nboot)
        koffA.samp=rep(NA,nboot)
        sA.samp=rep(NA,nboot)
        konB.samp=rep(NA,nboot)
        koffB.samp=rep(NA,nboot)
        sB.samp=rep(NA,nboot)
        
        if(mode=='raw'){
          Ai.obs=(alleleA[i,])/(lib.size)
          Bi.obs=(alleleB[i,])/(lib.size)
          for(k in 1:nboot){
            ABi=c(Ai.obs,Bi.obs)
            samp.index=sample.int(length(ABi),size=length(ABi),replace=TRUE)
            ABi=ABi[samp.index]
            cellsize.temp=(c(cellsize,cellsize))[samp.index]
            Ai=ABi[1:num.cells]
            Bi=ABi[(num.cells+1):length(ABi)]
            Ai.cellsize=cellsize.temp[1:num.cells]
            Bi.cellsize=cellsize.temp[(num.cells+1):length(ABi)]
            
            temp = as.data.frame(table(pmax(0,ceiling(Ai/bandwidth) * bandwidth-bandwidth/2)))
            q.temp=as.numeric(as.matrix(temp[,1])) 
            c.temp=as.numeric(as.matrix(temp[,2])) 
            y.temp=round(exp((log(q.temp)-alpha)/beta)) 
            p.temp=expit(kappa+tau*log(y.temp))  
            nc.temp=round(c.temp/p.temp) 
            if(sum(nc.temp[-1])>num.cells){
              loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
              nc.temp[loc.temp]=0
              nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
            } else{
              nc.temp[1]=num.cells-sum(nc.temp[-1])
            }
            Ai.corrected=rep(round(y.temp),nc.temp)
            
            temp = as.data.frame(table(pmax(0,ceiling(Bi/bandwidth) * bandwidth-bandwidth/2)))
            q.temp=as.numeric(as.matrix(temp[,1]))
            c.temp=as.numeric(as.matrix(temp[,2])) 
            y.temp=round(exp((log(q.temp)-alpha)/beta))
            p.temp=expit(kappa+tau*log(y.temp))  
            nc.temp=round(c.temp/p.temp)
            if(sum(nc.temp[-1])>num.cells){
              loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
              nc.temp[loc.temp]=0
              nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
            } else{
              nc.temp[1]=num.cells-sum(nc.temp[-1])
            }
            Bi.corrected=rep(round(y.temp),nc.temp)
            
            pb.temp=pb.moment(Ai.corrected,Ai.cellsize)
            konA.samp[k]=pb.temp[1]
            koffA.samp[k]=pb.temp[2]
            sA.samp[k]=pb.temp[3]
            
            pb.temp=pb.moment(Bi.corrected,Bi.cellsize)
            konB.samp[k]=pb.temp[1]
            koffB.samp[k]=pb.temp[2]
            sB.samp[k]=pb.temp[3]
          }
        } else if (mode == 'corrected'){
          Ai=(alleleA[i,])/(lib.size)
          temp = as.data.frame(table(pmax(0,ceiling(Ai/bandwidth) * bandwidth-bandwidth/2)))
          q.temp=as.numeric(as.matrix(temp[,1])) 
          c.temp=as.numeric(as.matrix(temp[,2])) 
          y.temp=round(exp((log(q.temp)-alpha)/beta))
          p.temp=expit(kappa+tau*log(y.temp))  
          nc.temp=round(c.temp/p.temp) 
          if(sum(nc.temp[-1])>num.cells){
            loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
            nc.temp[loc.temp]=0
            nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
          } else{
            nc.temp[1]=num.cells-sum(nc.temp[-1])
          }
          Ai.corrected=rep(round(y.temp),nc.temp)
          
          Bi=(alleleB[i,])/(lib.size)
          temp = as.data.frame(table(pmax(0,ceiling(Bi/bandwidth) * bandwidth-bandwidth/2)))
          q.temp=as.numeric(as.matrix(temp[,1]))
          c.temp=as.numeric(as.matrix(temp[,2])) 
          y.temp=round(exp((log(q.temp)-alpha)/beta)) 
          p.temp=expit(kappa+tau*log(y.temp)) 
          nc.temp=round(c.temp/p.temp) 
          if(sum(nc.temp[-1])>num.cells){
            loc.temp=which(rev(cumsum(rev(nc.temp)))>num.cells)
            nc.temp[loc.temp]=0
            nc.temp[loc.temp[length(loc.temp)]]=num.cells-sum(nc.temp[-loc.temp])
          } else{
            nc.temp[1]=num.cells-sum(nc.temp[-1])
          }
          Bi.corrected=rep(round(y.temp),nc.temp)
          
          for(k in 1:nboot){
            ABi=c(Ai.corrected,Bi.corrected)
            samp.index=sample.int(length(ABi),size=length(ABi),replace=TRUE)
            ABi=ABi[samp.index]
            cellsize.temp=(c(cellsize,cellsize))[samp.index]
            Ai=ABi[1:length(Ai)]
            Bi=ABi[(length(Ai)+1):length(ABi)]
            Ai.cellsize=cellsize.temp[1:num.cells]
            Bi.cellsize=cellsize.temp[(num.cells+1):length(ABi)]
            
            pb.temp=pb.moment(Ai,Ai.cellsize)
            konA.samp[k]=pb.temp[1]
            koffA.samp[k]=pb.temp[2]
            sA.samp[k]=pb.temp[3]
            
            pb.temp=pb.moment(Bi,Bi.cellsize)
            konB.samp[k]=pb.temp[1]
            koffB.samp[k]=pb.temp[2]
            sB.samp[k]=pb.temp[3]
          }
        }
       
        
        samp.filter=konA.samp>0 & koffA.samp >0 & sA.samp >0 & konB.samp>0 & koffB.samp >0 & sB.samp >0
        samp.filter[is.na(samp.filter)]=FALSE
        konA.samp=konA.samp[samp.filter]
        koffA.samp=koffA.samp[samp.filter]
        sA.samp=sA.samp[samp.filter]
        
        konB.samp=konB.samp[samp.filter]
        koffB.samp=koffB.samp[samp.filter]
        sB.samp=sB.samp[samp.filter]
        sizeB.samp=sB.samp/koffB.samp
        sizeA.samp=sA.samp/koffA.samp
        
        pval.kon[i]=sum(abs(konA.samp-konB.samp)>=abs(konA[i]-konB[i]))/length(konA.samp)
        pval.size[i]=sum(abs(sizeA.samp-sizeB.samp)>=abs(sA[i]/koffA[i]-sB[i]/koffB[i]))/length(konA.samp)
        cat('pval_kon =',round(pval.kon[i],3),';  pval_size =',round(pval.size[i],3),'.\n')
    }
    
    diff.allelic.obj=list(pval.kon=pval.kon,pval.size=pval.size)
    return(diff.allelic.obj)
}
