# 0. concatenate ERCC with hg19 and index the fasta file
cat ERCC92.fa hg19.fa > hg19_ERCC.fa
java -jar ~/bin/CreateSequenceDictionary.jar R= hg19_ERCC.fa O= hg19_ERCC.dict
samtools faidx hg19_ERCC.fa

# 1. Generate the genome using STAR, 50bp PE
genomeDir=/home/yuchaoj/hg19_ERCC
STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles ~/hg19_ERCC/hg19_ERCC.fa --sjdbFileChrStartEnd ~/hg19/gencode.v14.annotation.gtf.sjdb --sjdbOverhang 49 --runThreadN 4

# 2. Alignment jobs were executed as follows
genomeDir=/home/yuchaoj/hg19_ERCC
while read samp; do
echo 'STAR --genomeDir '$genomeDir' --readFilesIn '''$samp'''_R1_001.fastq '''$samp'''_R2_001.fastq --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix '$samp'_ --runThreadN 4' | bsub -M 60000
done < samp_ERCC.list

# 3. Convert .sam to .bam
while read samp; do
echo 'samtools view -bS '''$samp'''_Aligned.out.sam > '''$samp'''_Aligned.out.bam ' | bsub -M 60000
done < samp_ERCC.list

# 4. Filter & Sort
while read samp; do
echo 'perl filter_sam_v2.pl '''$samp'''_Aligned.out.bam '''$samp'''_Aligned.out.filtered.sam; samtools view -bS '''$samp'''_Aligned.out.filtered.sam > '''$samp'''_Aligned.out.filtered.bam' | bsub -M 60000
done < samp_ERCC.list

while read samp; do
echo 'java -Xmx30G -jar ~/bin/SortSam.jar INPUT='''$samp'''_Aligned.out.filtered.bam OUTPUT='''$samp'''_Aligned.out.filtered.sorted.bam SORT_ORDER=coordinate ' | bsub -M 40000
done < samp_ERCC.list

# 5. Add read group and index
while read samp; do
echo 'java -Xmx30G -jar ~/bin/AddOrReplaceReadGroups.jar INPUT='''$samp'''_Aligned.out.filtered.sorted.bam OUTPUT='''$samp'''_Aligned.out.filtered.sorted.rg.bam RGID=T2N_ERCC RGLB=scRNA_seq_T2N_ERCC RGPL=ILLUMINA RGPU=machine RGSM='$samp'; samtools index '''$samp'''_Aligned.out.filtered.sorted.rg.bam'  | bsub -M 30000
done < samp_ERCC.list

# 6. Get total read counts as well as read counts for ERCC
while read bam; do
echo $bam
samtools view -c $bam
done < rg.bam.list

while read ercc; do
echo $ercc
while read bam; do samtools view -c $bam $ercc;done < rg.bam.list | cat > $ercc.txt
done < ercc.id 

