output_table=function(alleleA,alleleB,gene.class.obj, allelic.kinetics.obj,
                      diff.allelic.obj,non.ind.obj){
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
    
    cellcount=matrix(nrow=nrow(alleleA),ncol=4,data=0)
    colnames(cellcount)=c('A','B','AB','Off')
    for(i in 1:nrow(alleleA)){
        cellcount[i,names(results.list[[i]])]=results.list[[i]] 
    }
    colnames(cellcount)=c('A_cell','B_cell','AB_cell','Off_cell')
    
    pval.kon=round(pval.kon,5)
    sizeA=round(sizeA,2)
    sizeB=round(sizeB,2)
    pval.size=round(pval.size,5)
    pval.ind=round(pval.ind,5)
    
    SCALE.output=cbind(genename,gene.category,konA,konB,pval.kon,
                       sizeA,sizeB,pval.size,cellcount,A.prop,B.prop,
                       pval.ind,non.ind.type)
    SCALE.output[is.na(SCALE.output)]='-'
    return(SCALE.output)
}