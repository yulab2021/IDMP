![image](https://github.com/yulab2021/IDMP/assets/133012957/004cc112-7f18-446b-83d6-99610b3f7f32)# Degradome analysis

## 0x00 Environment configuration

####

The following modules are needed to run

| module          | version |
| --------------- | ------- |
| python          | 3.9.3   |
| pandas          | 1.5.1   |
| plotnine        | 0.10.1  |
| numpy           | 1.23.4  |
| matplotlib_venn | 0.11.9  |
| seaborn         | 0.12.2  |
| matplotlib      | 3.7.0   |



## ##0x01 Data preprocessing

Before using our package, you may need to preprocess the file. Finally we need the bed file to continue the following analysis. If you already have the bed files, you can skip these steps.

We take Arabidopsis degradome data (SRP117737) as an example below.

To conduct the analysis, it's necessary to download both the Arabidopsis genome file in FASTA format and its corresponding annotation file in GFF3 format.

Please review the GFF3 file using the following criteria:
(1) Confirm that the chromosome ID matches the ID in the genome FASTA file.
(2) Ensure that the annotation in column 3 includes the "five_prime_UTR" feature, which is essential for uORF analysis.

#### Transfer sra file to fastq file

```
for file in */*sra
do
fastq-dump -O ./ $file 
done
```

#### Remove 3-prime connector
Note: verify the 3' adapter sequence to be used as input for the '-a' parameter.

```
for file in *fastq
do
cutadapt -a TGGAATTCTCGGGTGC -o $(echo $file | sed 's/.fastq/.fq/') $file ## check the adaptor sequence!
done
```

#### Filter the sequence whose length <18 nt

```
python ./Degradome_analysis/filter.py
```

#### Download Arabidopsis genome and build the star index
```
STAR --runMode genomeGenerate --genomeDir Ath.STAR.index --genomeFastaFiles TAIR10_genome.fa --runThreadN 20
```
#### Map the reads to the genome

```
mkdir map2genome

for file in filtered_data/*_trimed_filter.fq
do
STAR --runThreadN 6 --genomeDir /home/database/genome/Arabidopsis/Ath.STAR.index --readFilesIn $file --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.10 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix map2genome/$(echo $file | sed 's/filtered_data\///'| sed 's/_trimed_filter.fq//') &
done
```

#### Convert the bam file to bed files

```
for file in *bam
do
convert2bed --input=bam  <$file >$(echo $file | sed 's/.bam/.bed/') &
done
```

#### Turn bed files to a list file

```
for file in *.bed
do
echo $file >> input.list
done
```

Then you can use 'input.list' as the value of argument --input

## ##0x02 Co_translation RNA Decay

```
python ./Degradome_analysis/Degradome_analysis.py --mode 1 --input input.list --gff gff_file --genome genome_file --gene_name gene_name_annotation_file &
```

##### Output:

A lineplot shows the distribution of reads number in the range -48~2 of stop codon of all samples

A set of heatmaps

A set of tables shows the reads number in the range -48~2 of stop codon of each sample

A set of tables and barplots show the pattern_3 of each samples

A set of barplots show the reads/counts of each codon for each sample

## ##0x03 uORF identification

```
python ./Degradome_analysis/Degradome_analysis.py --mode 2 --input input.list --gff gff_file --gene_name gene_name_annotation_file --genome genome_file &
```

##### Output:

A set of heatmaps and tables show the distribution of reads

## ##0x04 Exon junction complex pausing 

```
python ./Degradome_analysis/Degradome_analysis.py --mode 3 --input input.list --gff gff_file --genome genome_file --gene_name gene_name_annotation_file --sample_info the_treatment_and_condition_of_samples &
```

##### Output:

A lineplot shows the the reads number in the range -50~0 of stop codon of each sample

A set of heatmaps

## ##0x05 5’P reads distribution along pre-miRNAs

```
python ./Degradome_analysis/Degradome_analysis.py --mode 4 --gff gff_file --input input.list --gff3_file pre-miRNA_and_mature_miRNA_annotation_file --mirna_name name_of_mirna  &
```

##### Output:

A lineplot shows the distribution of reads for a given precursor is plotted, and the relative position of the mature miRNA is given

A file shows the miRNAs have enrichment signals

## ##0x06 Identification of miRNA cleavage sites

First of all, you can use the following code to obtain the documents that need to be submitted to the website

```
python ./Degradome_analysis/Degradome_analysis.py --mode 5-1 --gff gff_file --genome genome_file --mirna_fa mature_miRNA_files --rna_id three_letter_abbreviations &
```

Then, you can download mirna.fa and cds.fa from the current working directory and submit them to the website https://www.zhaolab.org/psRNATarget/, then download the analysis result to do the next step.

```
python ./Degradome_analysis/Degradome_analysis.py --mode 5-2 --gff gff_file --input input.list --ps_file output_file --mirna_name name_of_mirna &
```

##### Output:

A file shows whether the cleavage site is significantly enriched over flanking 50nt.

A set of lineplots show the distribution of reads for the targets of the given miRNA

## ##0x07 Venn plot

Under the premise that you have run mode1, 2, 3, calling this script produces a Venn diagram showing the number of intersections that occur for genes with significant results in the three analysis modes

```
python ./Degradome_analysis/venn_plotting.py
```

