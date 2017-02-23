# SCALE
Allele-Specific Expression by Single-Cell RNA Sequencing


## Author
Yuchao Jiang, Nancy R. Zhang, Mingyao Li


## Maintainer
Yuchao Jiang <yuchaoj@upenn.edu>


## Description
SCALE is a statistical framework for **S**ingle **C**ell **AL**lelic **E**xpression analysis. SCALE estimates kinetic parameters that characterize the transcriptional bursting process at the allelic level, while accounting for technical bias and other complicating factors such as cell size. SCALE detects genes with significantly different bursting kinetics between the two alleles, as well as genes where the two alleles exhibit dependence in their bursting processes.


## Installation
```r
install.packages("devtools")
library(devtools)
install_github("yuchaojiang/SCALE/package")
```


## Demo code & Vignettes
* Bioinformatic pipeline
  * [Calling germline heterozygous SNPs by bulk-tissue sequencing](https://github.com/yuchaojiang/SCALE/blob/master/bioinfo/bulk_SNP.sh)
  * [Force-calling allele-specific read counts at germline heterozygous loci by single-cell RNA sequencing](https://github.com/yuchaojiang/SCALE/blob/master/bioinfo/scRNAseq_SNP.sh)
  * [Generating input for exogenous spike-ins for adjustment of technical variability by single-cell RNA sequencing](https://github.com/yuchaojiang/SCALE/blob/master/bioinfo/scRNAseq_spikein.sh)
* SCALE
  * [Demo code](https://github.com/yuchaojiang/SCALE/blob/master/demo/demo.R)
  * [Vignettes](https://github.com/yuchaojiang/SCALE/blob/master/demo/SCALE_vignettes.pdf)


## Google user group (Q&A)
If you have any questions with the package, please feel free to post in our Google user group https://groups.google.com/d/forum/SCALE_scRNAseq or email us at SCALE_scRNAseq@googlegroups.com. We will try our best to reply as soon as possible.


## Citation
Jiang, Y., Zhang, N.R., Li M., 2016. Modeling allele-specific gene expression by single-cell RNA sequencing. *Revision submitted*.

[Main Manuscript (bioRxiv)](http://biorxiv.org/content/early/2017/02/17/109629.full.pdf)

[Supplement (bioRxiv)](http://biorxiv.org/content/biorxiv/suppl/2017/02/17/109629.DC1/109629-1.pdf)

