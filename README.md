## Data and scripts in support of:

Christelle Robert and Mick Watson (2015) Errors in RNA-Seq quantification affect genes of relevance to human disease, submitted

Scripts should be considered "research quality" and are being provided for reproducibility purposes only.  Whilst we expect that these scripts may work on other datasets, this has not been tested.  We have not made any attempt to optimise the code.

## Steps to reproduce the analysis of mouse lung cancer

1. Download data from SRA accessions in file SRR_Acc_List.txt

2. Alignment to the mouse genome used STAR 2.4.0i and the Mus_musculus.GRCm38.dna.primary_assembly.fa file available from Ensembl.  STAR index built using:

	STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index_mouse --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa

3. An example of the STAR command line for alignment is:
	
	STAR --runThreadN 4 --genomeDir STAR_index_mouse --readFilesIn SRR1528720_1.fastq.gz --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 50 --alignSJoverhangMin 1 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix SRR1528720_1.STAR --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif

4. BAM files were sorted by name and indexed using samtools

5. Read counts against a known annotation (Ensembl file Mus_musculus.GRCm38.79.gtf) were calculated using htseq-count from HTSeq version 0.6.1.  An example command is:

	samtools view SRR1528735_1.STARAligned.sortedByName.bam | htseq-count -r name -o SRR1528735_1.STARAligned.sortedByName.sam -s no -t exon -m union -i gene_id - Mus_musculus.GRCm38.79.gtf > SRR1528735_1.STARAligned.sortedByName.htseq.u.txt

6. 
