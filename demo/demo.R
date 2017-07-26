library(SCALE)
####################################################
####################################################
##########                                ##########
##########           data input           ##########
##########                                ##########
####################################################
####################################################

data(mouse.blastocyst)
alleleA=mouse.blastocyst$alleleA
alleleB=mouse.blastocyst$alleleB
spikein_input=mouse.blastocyst$spikein_input
genename=rownames(alleleA)
sampname=colnames(alleleA)
head(colnames(alleleA))
head(rownames(alleleA))
rownames(spikein_input)
head(colnames(spikein_input))


####################################################
####################################################
##########                                ##########
##########        quality control         ##########
##########                                ##########
####################################################
####################################################

N = apply(alleleA + alleleB, 2, sum)
lib.size = N/mean(N)
spikein_read=spikein_input[,-(1:2)]
# remove cells with extreme library size factors
lib.size.filter=(lib.size>0.5 & lib.size < 2.5)
alleleA=alleleA[,lib.size.filter]; alleleB=alleleB[,lib.size.filter]
N = apply(alleleA + alleleB, 2, sum)
lib.size = N/mean(N)
spikein_read=spikein_read[,lib.size.filter[colnames(spikein_read)]]

# remove cells with extreme ratios of spike-in versus endogenous reads
endo_readA=alleleA[,colnames(spikein_read)] # endogenous reads from allele A
endo_readB=alleleB[,colnames(spikein_read)] # endogenous reads from allele B
ratio=apply(spikein_read,2,sum)/apply(endo_readA+endo_readA,2,sum)
ratio.filter=(ratio >0.05 & ratio < 3)
alleleA=alleleA[,is.na(match(colnames(alleleA),colnames(spikein_read)[!ratio.filter]))]
alleleB=alleleB[,is.na(match(colnames(alleleB),colnames(spikein_read)[!ratio.filter]))]
N = apply(alleleA + alleleB, 2, sum)
lib.size = N/mean(N)
spikein_read=spikein_read[,ratio.filter]
spikein_input=cbind(spikein_input[,1:2],spikein_read)


####################################################
####################################################
##########                                ##########
##########          heterogeneity         ##########
##########                                ##########
####################################################
####################################################

# SCALE needs to be applied to a homogeneous population of cells.
# Below is sanity check on data heterogeneity using PCA and t-SNE.

read=alleleA+alleleB
lib.size.factor=matrix(ncol=ncol(read),nrow=nrow(read),apply(read,2,sum)/mean(apply(read,2,sum)),byrow=T)
read=read/lib.size.factor # adjust for sequencing depth
AR=alleleA/(alleleA+alleleB)

AR=AR[apply(read,1,sum)>10,] # remove silent, lowly expressed genes
read=read[apply(read,1,sum)>10,]

AR=as.matrix(AR[apply(read,1,sd)>quantile(apply(read,1,sd),0.9),]) # select for highly variable genes
read=as.matrix(read[apply(read,1,sd)>quantile(apply(read,1,sd),0.9),])

# SVD and t-SNE on total coverage
svd.read=svd(read)
plot(svd.read$v[,1], svd.read$v[,2], xlab='PC1',ylab='PC2')
scatterplot3d(svd.read$v[,1], svd.read$v[,2], svd.read$v[,3],
              xlab='PC1',ylab='PC2',zlab='PC3')
tsne_read=tsne(t(read),k=3,initial_dims=30) # If this takes too long, select highly variable genes/SNPs.
plot(tsne_read[,1],tsne_read[,2],xlab='tSNE1',ylab='tSNE2')
scatterplot3d(tsne_read[,1],tsne_read[,2],tsne_read[,3],
              xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')

# SVD and t-SNE on allelic ratio (AR)
AR=AR[!apply(is.nan(AR),1,any),]
svd.AR=svd(AR)
plot(svd.AR$v[,1], svd.AR$v[,2], xlab='PC1',ylab='PC2')
scatterplot3d(svd.AR$v[,1], svd.AR$v[,2], svd.AR$v[,3],
              xlab='PC1',ylab='PC2',zlab='PC3')
tsne_AR=tsne(t(AR),k=3,initial_dims=30) # If this takes too long, select highly variable genes/SNPs.
plot(tsne_AR[,1],tsne_AR[,2],xlab='tSNE1',ylab='tSNE2')
scatterplot3d(tsne_AR[,1],tsne_AR[,2],tsne_AR[,3],
              xlab='tSNE1',ylab='tSNE2',zlab='tSNE3')


####################################################
####################################################
##########                                ##########
##########     technical variability      ##########
##########                                ##########
####################################################
####################################################

# Make sure the column names of spikein_input correspond to (a subset of)
# those of alleleA and alleleB.
abkt=tech_bias(spikein_input=spikein_input, alleleA = alleleA, 
               alleleB = alleleB, readlength = 50, pdf=FALSE)


####################################################
####################################################
##########                                ##########
##########      gene classification       ##########
##########                                ##########
####################################################
####################################################

# Only first 200 genes are computed for demonstration.
gene.class.obj=gene_classify(alleleA = as.matrix(alleleA[1:200,]),
                             alleleB = as.matrix(alleleB[1:200,]))
A.prop=gene.class.obj$A.prop
B.prop=gene.class.obj$B.prop
gene.category=gene.class.obj$gene.category
results.list=gene.class.obj$results.list

# Below is loading pre-computed gene classification results.
data(gene.class.obj)
A.prop=gene.class.obj$A.prop
B.prop=gene.class.obj$B.prop
gene.category=gene.class.obj$gene.category
results.list=gene.class.obj$results.list


####################################################
####################################################
##########                                ##########
##########   allelic bursting kinetics    ##########
##########                                ##########
####################################################
####################################################

cellsize=rep(1,ncol(alleleA))  # cell size input
allelic.kinetics.obj=allelic_kinetics(alleleA = alleleA[1:1000,], 
                                      alleleB = alleleB[1:1000,], 
                                      abkt = abkt, 
                                      gene.category = gene.category[1:1000], 
                                      cellsize = cellsize, pdf=TRUE)

# Below is loading pre-computed allelic kinetics.
data(allelic.kinetics.obj)
bandwidth=allelic.kinetics.obj$bandwidth
konA=allelic.kinetics.obj$konA
koffA=allelic.kinetics.obj$koffA
sA=allelic.kinetics.obj$sA
konB=allelic.kinetics.obj$konB
koffB=allelic.kinetics.obj$koffB
sB=allelic.kinetics.obj$sB
sizeA=sA/koffA
sizeB=sB/koffB


####################################################
####################################################
##########                                ##########
########## differential bursting kinetics ##########
##########                                ##########
####################################################
####################################################
diff.allelic.obj=diff_allelic_bursting(alleleA = alleleA,
                                        alleleB = alleleB,
                                        cellsize = cellsize,
                                        gene.category = gene.category,
                                        abkt = abkt,
                                        allelic.kinetics.obj = allelic.kinetics.obj,
                                        mode = 'corrected')
pval.kon = diff.allelic.obj$pval.kon; pval.size = diff.allelic.obj$pval.size

# Below is loading pre-computed testing results.
data(diff.allelic.obj)
pval.kon = diff.allelic.obj$pval.kon; pval.size = diff.allelic.obj$pval.size

####################################################
####################################################
##########                                ##########
##########    non-independent bursting    ##########
##########                                ##########
####################################################
####################################################

non.ind.obj=non_ind_bursting(alleleA = alleleA, alleleB = alleleB,
                             gene.category = gene.category,
                             results.list = results.list)
pval.ind=non.ind.obj$pval.ind
non.ind.type=non.ind.obj$non.ind.type


####################################################
####################################################
##########                                ##########
##########            plot genes          ##########
##########                                ##########
####################################################
####################################################

i=which(genename=='Btf3l4')
allelic_plot(alleleA = alleleA, alleleB = alleleB,
             gene.class.obj = gene.class.obj,
             allelic.kinetics.obj = allelic.kinetics.obj,
             diff.allelic.obj = diff.allelic.obj,
             non.ind.obj = non.ind.obj, i= i)

####################################################
####################################################
##########                                ##########
##########     generate output tables     ##########
##########                                ##########
####################################################
####################################################

SCALE.output=output_table(alleleA=alleleA, alleleB=alleleB,
                          gene.class.obj = gene.class.obj,
                          allelic.kinetics.obj = allelic.kinetics.obj,
                          diff.allelic.obj = diff.allelic.obj,
                          non.ind.obj = non.ind.obj)
head(SCALE.output)
write.table(SCALE.output,file='SCALE.output.txt',col.names = T,row.names = F,quote = F,sep='\t')
