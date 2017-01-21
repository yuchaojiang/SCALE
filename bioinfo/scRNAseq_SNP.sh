# 0. Get splice junction database:
wget http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/GENCODE/Old/gencode.v14.annotation.gtf.sjdb

# 1. Generate the genome using STAR, 100bp PE
genomeDir=/home/yuchaoj/hg19
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles ~/hg19/hg19.fa --sjdbFileChrStartEnd ~/hg19/gencode.v14.annotation.gtf.sjdb --sjdbOverhang 99 --runThreadN 4

# 2. Alignment jobs were executed as follows
genomeDir=/home/yuchaoj/hg19
while read samp
do
echo 'STAR --genomeDir '$genomeDir' --readFilesIn '''$samp'''_1.fastq '''$samp'''_2.fastq --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix '$samp'_ --runThreadN 4' | bsub -M 60000
done < samp.list

# 3. Convert .sam to .bam
while read samp; do
echo 'samtools view -bS '''$samp'''_Aligned.out.sam > '''$samp'''_Aligned.out.bam ' | bsub -M 60000
done < samp.list

# 4. Filter & sort
while read samp; do
echo 'perl filter_sam_v2.pl '''$samp'''_Aligned.out.bam '''$samp'''_Aligned.out.filtered.sam; samtools view -bS '''$samp'''_Aligned.out.filtered.sam > '''$samp'''_Aligned.out.filtered.bam' | bsub -M 60000
done < samp.list

while read samp; do
echo 'java -Xmx30G -jar ~/bin/SortSam.jar INPUT='''$samp'''_Aligned.out.filtered.bam OUTPUT='''$samp'''_Aligned.out.filtered.sorted.bam SORT_ORDER=coordinate ' | bsub -M 40000
done < samp.list

# 5. Add read group and index
while read samp; do
echo 'java -Xmx30G -jar ~/bin/AddOrReplaceReadGroups.jar INPUT='''$samp'''_Aligned.out.filtered.sorted.bam OUTPUT='''$samp'''_Aligned.out.filtered.sorted.rg.bam RGID=UCF1014 RGLB=scRNA_seq_UCF1014 RGPL=ILLUMINA RGPU=machine RGSM='$samp'; samtools index '''$samp'''_Aligned.out.filtered.sorted.rg.bam'  | bsub -M 40000
done < samp.list

# 6. parse file: position.txt contains all the heterogyzous loci (chr + coordinate) returned by GATK HaplotypeCaller using WGS.
while read bam; do
echo ' samtools mpileup -E -f /home/yuchaoj/hg19/hg19.fa -d 1000000 --position position.txt '$bam'  > '$bam'.mpileup;
perl pileup2base_no_strand.pl '$bam'.mpileup 30 '$bam'.parse30.txt
' | bsub -M 20000
done < rg.bam.list

