
Log onto Abel and create a directory called `transcriptomics`

## Get data

Download _one_ of the fastq-files (they come in pairs) from [here](http://folk.uio.no/jonbra/)

`wget http://folk.uio.no/jonbra/aboral-1_R1.fastq.gz`  
`wget http://folk.uio.no/jonbra/aboral-1_R2.fastq.gz`

Download also this file:

`wget http://folk.uio.no/jonbra/genome_transcriptome.tar`

Unpack it
`tar -xvf genome_transcriptome.tar

## Quality check

`module load fastqc`
`fastqc *.gz`

Download the `.html` files to your local computer

## Trimming

```
module load trim-galore/0.3.3

trim_galore --fastqc -o outfolder --paired R1-file R2-file
```

Run fastqc again on the trimmed files and see the changes

## Mapping

```
module load tophat/2.1.1
module load bowtie2/2.2.9
module load samtools/1.3.1

tophat -G /genome_transcriptome/ML2.2.nogene.gff3 -p 8 --library-type fr-firststrand -o mapping /genome_transcriptome/Ml_genome trimmed_R1_file trimmed_R2_file
```

## Expression count 

Sort the mapping file

```
module load samtools
samtools sort -O bam -T tmp -n mapping/accepted_hits.bam -o mapping/accepted_hits_sorted.bam
```

Install and run HTSeq

```
module load python2
pip install --user HTSeq

python -m HTSeq.scripts.count -f bam -r name -s reverse -t mRNA -i ID mapping/accepted_hits_sorted.bam mapping/genome_transcriptome/ML2.2.nogene.gff3 > sample-name.txt
```

Try this if you get an error

```
unset LC_CTYPE
unset LANG
```
