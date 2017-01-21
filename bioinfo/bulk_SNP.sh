# 1. Index the genome template
bwa index ~/structure/hg19/ucsc.hg19.fasta

# 2. Align reads in .fastq file to the template (sequenced at two different lanes)
bwa mem -M -t 16 ~/structure/hg19/ucsc.hg19.fasta _EGAR00001248897_UCF_1014_NoIndex_L008_R1_001.fastq _EGAR00001248897_UCF_1014_NoIndex_L008_R2_001.fastq > _EGAR00001248897_UCF_1014_NoIndex_L008.sam

bwa mem -M -t 16 ~/structure/hg19/ucsc.hg19.fasta _EGAR00001248899_UCF1014_NoIndex_L002_R1_001.fastq _EGAR00001248899_UCF1014_NoIndex_L002_R2_001.fastq > _EGAR00001248899_UCF1014_NoIndex_L002.sam

# 3. Convert .sam to .bam and sort
samtools view -bS _EGAR00001248897_UCF_1014_NoIndex_L008.sam > _EGAR00001248897_UCF_1014_NoIndex_L008.bam
line=_EGAR00001248897_UCF_1014_NoIndex_L008.bam
java -jar ~/bin/SortSam.jar INPUT=$line OUTPUT=$line.sorted.bam SORT_ORDER=coordinate 

samtools view -bS _EGAR00001248899_UCF1014_NoIndex_L002.sam > _EGAR00001248899_UCF1014_NoIndex_L002.bam
line=_EGAR00001248899_UCF1014_NoIndex_L002.bam
java -jar ~/bin/SortSam.jar INPUT=$line OUTPUT=$line.sorted.bam SORT_ORDER=coordinate 

# 4. Add read group
line=_EGAR00001248897_UCF_1014_NoIndex_L008.bam
#line=_EGAR00001248899_UCF1014_NoIndex_L002.bam

java -jar ~/bin/AddOrReplaceReadGroups.jar INPUT=$line.sorted.bam OUTPUT=$line.sorted.rg.bam RGID=$line RGLB=WGS_UCF1014 RGPL=ILLUMINA RGPU=machine RGSM=UCF1014
samtools index $line.sorted.rg.bam

# 5. Merge two different lanes
samtools merge UCF1014.merge.bam _EGAR00001248899_UCF1014_NoIndex_L002.bam.sorted.rg.bam _EGAR00001248897_UCF_1014_NoIndex_L008.bam.sorted.rg.bam
java -jar ~/bin/SortSam.jar INPUT=UCF1014.merge.bam OUTPUT=UCF1014.merge.sorted.bam SORT_ORDER=coordinate 
samtools index UCF1014.merge.sorted.bam

# 6. Dedup
java -jar ~/bin/MarkDuplicates.jar INPUT=UCF1014.merge.sorted.bam OUTPUT=UCF1014.merge.sorted.dedup.bam  METRICS_FILE=UCF1014.merge.sorted.dedup.metrics.txt PROGRAM_RECORD_ID= MarkDuplicates PROGRAM_GROUP_VERSION=null PROGRAM_GROUP_NAME=MarkDuplicates
java -jar ~/bin/BuildBamIndex.jar INPUT=UCF1014.merge.sorted.dedup.bam

# 7. Realign
java -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /home/stat/yuchaoj/structure/hg19/ucsc.hg19.fasta -I UCF1014.merge.sorted.dedup.bam -known /home/stat/yuchaoj/structure/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /home/stat/yuchaoj/structure/hg19/1000G_phase1.indels.hg19.sites.vcf -o UCF1014.merge.sorted.dedup.target_intervals.list 

java -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner -R /home/stat/yuchaoj/structure/hg19/ucsc.hg19.fasta -I UCF1014.merge.sorted.dedup.bam -targetIntervals UCF1014.merge.sorted.dedup.target_intervals.list -known /home/stat/yuchaoj/structure/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /home/stat/yuchaoj/structure/hg19/1000G_phase1.indels.hg19.sites.vcf -o UCF1014.merge.sorted.dedup.realigned.bam

# 8. Recalibrate
java -jar ~/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/stat/yuchaoj/structure/hg19/ucsc.hg19.fasta -I UCF1014.merge.sorted.dedup.realigned.bam -knownSites /home/stat/yuchaoj/structure/hg19/dbsnp_138.hg19.vcf -knownSites /home/stat/yuchaoj/structure/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /home/stat/yuchaoj/structure/hg19/1000G_phase1.indels.hg19.sites.vcf -o UCF1014.merge.sorted.dedup.realigned.recal_data.table

java -jar ~/bin/GenomeAnalysisTK.jar -T PrintReads -R /home/stat/yuchaoj/structure/hg19/ucsc.hg19.fasta -I UCF1014.merge.sorted.dedup.realigned.bam -BQSR UCF1014.merge.sorted.dedup.realigned.recal_data.table -o UCF1014.merge.sorted.dedup.realigned.recal.bam
samtools index UCF1014.merge.sorted.dedup.realigned.recal.bam

# 9. GATK HaplotypeCaller
java -jar ~/bin/GenomeAnalysisTK.jar -R ~/structure/hg19/ucsc.hg19.fasta -T HaplotypeCaller -I UCF1014.merge.sorted.dedup.realigned.recal.bam -o UCF1014.merge.sorted.dedup.realigned.recal.raw.snps.indels.g.vcf

