### The basic steps are:

    `bwa index` to index the reference genome (only needs to be done once per genome).
    'bwa aln' to align the reads
    'bwa sampe' to report the alignments.

After aligning, run SAMTools to find variations. The basic steps are

    'samtools faidx' to index the reference genome
    'samtools view' to load the alignments
    'samtools index' to index the BAM file
    'samtools mpileup' to find variations.
    'samtools view' is a useful command to inspect how the reads align to the genome at a given position

### Get the data

```
wget http://folk.uio.no/timothyh/course_vc_2016_imbv/inputData/human_g1k_v37_chr5/gatkBundle/human_g1k_v37_chr5.fasta .
wget http://folk.uio.no/timothyh/course_vc_2016_imbv/inputData/reads_agilentV1_chr5/real/real_agilentV1_chr5.R2.fastq .
wget http://folk.uio.no/timothyh/course_vc_2016_imbv/inputData/reads_agilentV1_chr5/real/real_agilentV1_chr5.R1.fastq .
```

Index the genome (takes a couple of minutes)

```
[jonbra@freebee mapping]$ bwa index human_g1k_v37_chr5.fasta
```

Align the reads
```
bwa aln human_g1k_v37_chr5.fasta real_agilentV1_chr5.R1.fastq > read1.sai
bwa aln human_g1k_v37_chr5.fasta real_agilentV1_chr5.R2.fastq > read2.sai
bwa sampe human_g1k_v37_chr5.fasta read?.sai real_agilentV1_chr5.R?.fastq > aln.sam
```


```
# We convert the SAM file to a BAM file to save space
samtools view -bT human_g1k_v37_chr5.fasta aln.sam > aln.bam

# We sort the file
samtools sort -O BAM -o aln_sorted.bam aln.bam

# And then we index it
samtools index aln_sorted.bam aln_sorted.bai
```

Download the chromosome 5 sequence and the sorted bam file and bam index file.


## Inspect the mapping in IGV browser
Start the IGV browser using Java Web Start (like JalView) from [here](http://software.broadinstitute.org/software/igv/download) (Launch with 750 MB).

Download the tiles used for the exome capture: `wget http://folk.uio.no/timothyh/course_vc_2016_imbv/inputData/human_g1k_v37_chr5/agilentV1/agilent37M.chr5.b37.bed .` and open this in IGV using Load from File.
