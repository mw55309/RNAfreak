## Data and scripts in support of:

Christelle Robert and Mick Watson (2015) Errors in RNA-Seq quantification affect genes of relevance to human disease, submitted

Scripts should be considered "research quality" and are being provided for reproducibility purposes only.  Whilst we expect that these scripts may work on other datasets, this has not been tested.  We have not made any attempt to optimise the code.

## Steps to reproduce the analysis of mouse lung cancer

* 1 Download data from SRA accessions in file SRR_Acc_List.txt

* 2 Alignment to the mouse genome used STAR 2.4.0i and the Mus_musculus.GRCm38.dna.primary_assembly.fa file available from Ensembl.  STAR index built using:
```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index_mouse --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa
```
* 3 An example of the STAR command line for alignment is:
```	
STAR --runThreadN 4 --genomeDir STAR_index_mouse --readFilesIn SRR1528720_1.fastq.gz --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 50 --alignSJoverhangMin 1 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix SRR1528720_1.STAR --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif
```
* 4 BAM files were sorted by name and indexed using samtools

* 5 Read counts against a known annotation (Ensembl file Mus_musculus.GRCm38.79.gtf) were calculated using htseq-count from HTSeq version 0.6.1.  An example command is:
```
samtools view SRR1528735_1.STARAligned.sortedByName.bam | htseq-count -r name -o SRR1528735_1.STARAligned.sortedByName.sam -s no -t exon -m union -i gene_id - Mus_musculus.GRCm38.79.gtf > SRR1528735_1.STARAligned.sortedByName.htseq.u.txt
```
* 6 The above command produces two major outputs for each dataset
  * A (large) SAM file detailing all read alignments and how htseq-count has categorised the read
  * A text file containing the counts per transcript (from the GTF)

* 7 A matrix of the htseq-count results are available from: http://www.ark-genomics.org/tmp/Mouse_Lung_Cancer_counts.txt

* 8 We created an exon bed file using script gtf2bed.pl and the BED file is available from http://www.ark-genomics.org/tmp/Mus_musculus.GRCm38.79.exon.bed.  Note that this script specifically created feature names of the format exon_id.gene_id which are used further downstream in the analysis.
```
perl gtf2bed.pl Mus_musculus.GRCm38.79.gtf > Mus_musculus.GRCm38.79.exon.bed
```

* 9 If the genome file is not already indexed, then index it using samtools
```
samtools faidx Mus_musculus.GRCm38.dna.primary_assembly.fa
```

* 10 We create counts of MMG (multi-map groups) per sample using script count_from_bed.pl.  (Note this script expects feature names of the format exon_id.gene_id, and also expects you to have assigned reads to Ensembl IDs.  Specifically, the script expects uniquely mapped reads to have a SAM flag beginning XF:Z:ENS)
```
perl count_from_bed.pl SRR1528735_1.STARAligned.sortedByName.sam Mus_musculus.GRCm38.79.exon.bed Mus_musculus.GRCm38.dna.primary_assembly.fa > SRR1528735_1.STARAligned.sortedByName.sam.cts
```

* 11 We then collate all of the MMG counts across all samples using script collect_group_counts.pl.  The command below applies filters: an MMG must have a count greater than 100 in more than 50% of the .cts files to make it into the final results 
```
perl collect_group_counts.pl 100 0.5 *.gene.cts > groupcounts.txt
``` 

* 12 We then merge groups such that any group completely contained within a larger group is merged
```
perl collapse_groups.pl groupcounts.txt > mergedgroups.txt
```

The merged group counts are available from http://www.ark-genomics.org/tmp/mergedcounts.txt
