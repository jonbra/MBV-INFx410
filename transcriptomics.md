# Get data

Download raw transcriptome sequences: `wget `
Uncompress (untar): `tar -xvf `

Download the genome and gene annotation: `wget `
Uncompress: `tar -xvf `

## Quality check
module load fastqc

fastqc * # ca. 3 min. Gjøre dette på Desktopen - ser html-filene direkte

Most sequences are of really high quality. Some needs a little trimming at the end. 

module load trim-galore/0.3.3

# Trim galore
[Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
Default phred 20. Uses standard illumina sequence adapters. No singletons are kept. Runs fastqc on the trimmed reads.

Running the file like: `bash trim_galore.sh`

``
#!/bin/bash
module load trim-galore/0.3.3 # load the program

# $f1 will be R1 files and $f2 till be R2 files
for f1 in *R1.fastq.gz
do
	f2=${f1%%R1.fastq.gz}"R2.fastq.gz"

    trim_galore --fastqc -o ~/MBV-INFX410-testing/transcriptomics/RawData/TrimmedSeqs --paired $f1 $f2

done
```

19:42 - 20:00: 18 min totalt

## Mapping (4:30 per prøve => litt over halvtime)

```
#!/bin/bash
module load tophat/2.1.1
module load bowtie2/2.2.9
module load samtools/1.3.1

# Must be inside the TrimmedSeqs directory
# $f1 will be R1 files and $f2 till be R2 files
for f1 in *R1_val_1.fq.gz
do
	f2=${f1%%R1_val_1.fq.gz}"R2_val_2.fq.gz"

    tophat -G ~/MBV-INFX410-testing/transcriptomics/RawData/genome_transcriptome/ML2.2.nogene.gff3 -p 8 --library-type fr-firststrand -o ~/MBV-INFX410-testing/transcriptomics/RawData/Mapping/$f1/ ~/MBV-INFX410-testing/transcriptomics/RawData/genome_transcriptome/Ml_genome ~/MBV-INFX410-testing/transcriptomics/RawData/TrimmedSeqs/$f1 ~/MBV-INFX410-testing/transcriptomics/RawData/TrimmedSeqs/$f2

done
```

## Counting


