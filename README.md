# Sequencing depth bias on microdiversity estimates

Estimation of sequencing depth influence on microdiversity estimates


## Required tools: Enveomics, BBTools, inStrain, Bowtie2, bedtools and samtools programs

1) Enveomics scripts collection can be downloaded from http://enve-omics.ce.gatech.edu/enveomics/download
```
git clone git://github.com/lmrodriguezr/enveomics.git enveomics 
```
2) BBTools scripts collection can be downloaded from https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/

3) inStrain can be downloaded following the Installation instructions from https://instrain.readthedocs.io/en/latest/installation.html

```
git clone https://github.com/MrOlm/instrain.git

cd instrain

pip install .
```
4) Bowtie2 can be downloaded from https://github.com/BenLangmead/bowtie2

5) Beedtools can be downloaded from https://bedtools.readthedocs.io/en/latest/content/installation.html

6) Samtools can be downloaded from https://github.com/samtools/samtools


## Illumina metagenomic trimming (not needed for Long-reads)

```
bbduk.sh in1=Short.MG1.R1.fastq in2=Short.MG1.R2.fastq out1=Short.MG1.trimmed.R1.fastq out2=Short.MG1.trimmed.R2.fastq ref=adapters.fa ktrim=r k=28 mink=12 hdist=1 tbo=t tpe=t qtrim=rl trimq=20 minlength=100
```

## Convert metagenomes from .fastq to .fasta format

Short-reads:
```
cat Short.MG1.R1.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > Short.MG1.R1.fasta
cat Short.MG1.R2.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > Short.MG1.R2.fasta
```

Long-reads:
```
cat Long.MG1.fastq | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}' > Long.MG1.fasta
```

## Short- and long-read random subsamplings of 1,5,10,20,30,40,50,60,70,80,90% of reads

Short-reads:
```
FastA.subsample.pl -f 1,5,10,20,30,40,50,60,70,80,90 Short.MG1.R1.fasta
```
  a) Selection of same reads for both paired-end files ==> This is an example using the 10% of reads:

```
grep ">" Short.MG1.R1.fasta.10.0000-1.fa > IDS.txt
sed -i 's/-1/-2/g' IDS.txt
FastA.filter.pl IDS.txt Short.MG1.R1.fasta.10.0000-1.fa > Short.MG1.R2.fasta.10.0000-1.fa
```

Long-reads:
```
FastA.subsample.pl -f 1,5,10,20,30,40,50,60,70,80,90 Long.MG1.fasta
```

## Fragment long-reads in 200 bps length 

```
shred.sh in=Long.MG1.fasta.10.0000-1.fa out=Long.Fragmented.MG1.fasta.10.0000-1.fa length=200 minlength=0
```

## Mapping with Bowtie2 (example with short-reads and 10% of reads)

```
cat MAGs* > CAT.ALL.MAGs.fasta

bowtie2-build CAT.ALL.MAGs.fasta CAT.ALL.MAGs.db

bowtie2 --reorder --no-unal -f -p 36 -x CAT.ALL.MAGs.db -1 Short.MG1.R1.fasta.10.0000-1.fa -2 Short.MG1.R2.fasta.10.0000-1.fa > Short.MG1.10.sam

sam.filter.rb -m Short.MG1.10.sam -o filter95_Short.MG1.10.sam -t 16
```

## Mapping with Bowtie2 (example with long-reads and 10% of reads)

```
cat MAGs* > CAT.ALL.MAGs.fasta

bowtie2-build CAT.ALL.MAGs.fasta CAT.ALL.MAGs.db

bowtie2 --reorder --no-unal -f -p 36 -x CAT.ALL.MAGs.db -U Long.Fragmented.MG1.fasta.10.0000-1.fa > Long.Fragmented.MG1.10.sam

sam.filter.rb -m Long.Fragmented.MG1.10.sam -o filter95_Long.Fragmented.MG1.10.sam -t 16
```

## Nucleotide diversity (pi), rarefied nucleotide diversity at 200X, and Average Nucleotide Identity of reads (ANIr) quantification using inStrain

```
inStrain profile filter95.Long.Fragmented.MG1.sam MAG1.fasta -o IS_MAG1_Long.Fragmented.MG1 --pairing_filter non_discordant --skip_mm_profiling --rarefied_coverage 200 --skip_plot_generation
```

Desired outputs are in: IS_MAG1_Long.Fragmented.MG1/output/IS_MAG1_Long.Fragmented.MG1_genome_info.tsv

genome/MAG name: Column 1
coverage: Column 2
nucl_diversity (pi): Column 4
nucl_diversity_rarefied (pi rarefied at 200X): Column 13
reads_mean_PID (ANIr): Column 28


## Rarefaction curves and interpolation of nucleotide diversity to a standardized sequencing depth (e.g. 200X)

To generate rarefaction curves and interpolations use the Rarefaction.Interpolation.R script with coverage.table.csv and nucleotide.diversity.table.csv/anir.table.csv


## Example of rarefaction curve graph:

Each line represent the nucleotide diversity and sequencing depth obtained for one MAG across all fractions of a metagenome. The diversity ratio represents the nucleotide diversity obtained for each fraction of a metagenome (subsample) divided by the nucleotide diversity of the total (before subsampling) metagenome. Specifically:
 
a) Green lines shows increasing diversity with increasing coverage

b) Red lines shows decreasing diversity with increasing coverage

![figure](/Example.Rarefaction.svg)


## Estimation of average estimated error (%)

To generate average estimated error (%) plots use the Average.Estimated.Error.R script with coverage.table.csv and nucleotide.diversity.table.csv/anir.table.csv
