# Contents
[Exercise 1](#1)  
[Exercise 2](#2)    
[Exercise 3](#3)  


## Exercise 1 - Quality assessment of sequence data <a name="1"></a>
First log onto Freebee and create a directory called `transcriptomics` and go inside it.

### Get data

Download _one_ of the fastq-samples (they come in pairs, forward and reverse) from [here](http://folk.uio.no/jonbra/MBV-INF4410_2017/Transcriptomics/)  
For example:

```
wget http://folk.uio.no/jonbra/MBV-INF4410_2017/Transcriptomics/aboral-1_R1.fastq.gz  
wget http://folk.uio.no/jonbra/MBV-INF4410_2017/Transcriptomics/aboral-1_R2.fastq.gz
```
These are paired-end sequences (usually named R1 and R2). Most programs require that these files contains the same number of reads, and that they are sorted in the same order. Any idea how we can check this?  

Let's first count the number of sequences:

```
# We use zcat to print the contents of a gzipped file to screen, and then we pipe it directly into the wc command
zcat aboral-1_R1.fastq.gz | wc -l
zcat aboral-1_R2.fastq.gz | wc -l
```
How many sequences are there in the two files? (remember that a sequence in .fastq-format covers four lines). Do both files contain the same number?

Another way to do this is to count the number of `@` characters in each file, because each sequence in a `.fastq` begins with an `@`:
```
zcat aboral-1_R1.fastq.gz | grep -c "^@"
```

`grep` is a very useful command used to extract lines of files matching different patterns and do stuff with these lines. Such as in this case count the nuber of lines beginning with @.

But we should also check that the sequences come in the same order in the two files, i. e. that sequences from the same pair are at the same place: 

```
zcat aboral-1_R1.fastq.gz | grep "^@" | head
```
```
# @NS500336:69:H5KLLAFXX:1:11101:19827:1048 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:10714:1050 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:8501:1050 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:20889:1050 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:3254:1050 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:23010:1051 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:7060:1053 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:4582:1055 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:8956:1056 1:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:20748:1057 1:N:0:CGCTCATT+NTTCGCCT
```
```
zcat aboral-1_R2.fastq.gz | grep "^@" | head
```
```
# @NS500336:69:H5KLLAFXX:1:11101:19827:1048 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:10714:1050 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:8501:1050 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:20889:1050 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:3254:1050 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:23010:1051 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:7060:1053 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:4582:1055 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:8956:1056 2:N:0:CGCTCATT+NTTCGCCT
# @NS500336:69:H5KLLAFXX:1:11101:20748:1057 2:N:0:CGCTCATT+NTTCGCCT
```
What does these commands do?

Notice that the names of each read are the same in both files, except the number "1" or "2" which indicates the first and second read of each pair. We didn't check all the sequences, but I'm pretty confident that these files are ok!  


### Quality check  
To inspect the reads and visualize the quality, run FastQC:  

```
module load fastqc
fastqc *.gz
```  

Download the `.html` files to your local computer and look at them in a web browser.  
- How long are the reads?
- How does the quality look like? Are there any differences in the quality score between the R1 and the R2 files?
- Are there any leftover sequencing adapters?

### Trimming  
We trim the reads using [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/). Default settings is to trim nucleotides lower than phred score 20 and looks for standard Illumina sequencing adapters. The option `--paired` tells the program to expect two files of paired reads, and `--fastqc` tells it to run FastQC on the trimmed reads. The job takes a couple of minutes.

```
# First we check which versions of the program are available on Abel
module avail trim-galore

# We use the latest version
module load trim-galore/0.4.4

trim_galore --fastqc --paired R1-file R2-file
```
Inspect the trimming reports in the terminal (remember how to view files in Unix?).  
- Were any reads trimmed? 
- Were there any reads with sequencing adapters?  

Download the fastqc reports of the trimmed reads and look at them in the browser. Notice any differences compared to the untrimmed reads?

[Top](#contents)
# Exercise 2 - Mapping  <a name="2"></a>
Now we are ready to map the reads to the genome and count the gene expression. Download this file which contains the _Mnemiopsis leidyi_ genome and transcriptome (in addition to a few index files which we'll need later):

`wget http://folk.uio.no/jonbra/MBV-INF4410_2017/Transcriptomics/genome_transcriptome.tar`

and unpack it with:
`tar -xvf genome_transcriptome.tar`  

We use TopHat2 to map the trimmed reads to the genome. Tophat first tries to map the entire read to the genome. Then, for the reads that do not map it splits them into shorter pieces and tries to map those. In this way, reads which spans introns can be mapped. This will also generate information about how where exons and introns are and if any transcripts are alternatively spliced. The mapping takes about 10 min so you can take a break (but check the first minute that there are no errors!). Run the following commands:  

```
module load tophat/2.1.1
module load bowtie2/2.2.9 # bowtie2 is the actual mapper
module load samtools/1.3.1 # needed to process files

tophat -G genome_transcriptome/ML2.2.nogene.gff3 -p 8 --library-type fr-firststrand genome_transcriptome/Ml_genome trimmed_R1_file trimmed_R2_file &
```
  
### Inspecting the .bam file in IGV
Download the accepted_hits.bam, genome and gff3 file
Sort the bam-file in IGV, then make an index


[Top](#contents)
# Exercise 3 - Counting gene expression <a name="3"></a>

First, sort the mapping file (remember to load samtools if you haven't):

`samtools sort -O bam -T tmp -n tophat_out/accepted_hits.bam -o tophat_out/accepted_hits_sorted.bam` (takes 1-2 min).

Install and run HTSeq

```
module load python2
pip install --user HTSeq

python -m HTSeq.scripts.count -f bam -r name -s reverse -t mRNA -i ID tophat_out/accepted_hits_sorted.bam genome_transcriptome/ML2.2.nogene.gff3 > sample-name.txt
```
Takes about 5 min.  
Try this if you get an error

```
unset LC_CTYPE
unset LANG
```
