# Tutorial

This tutorial assumes a GNU/Linux system and that tatajuba is already installed. For installation instructions please
refer to the [README file](../README.md) in the entry page of the github repository. 

We also assume that some tools are installed or that you can handle installing missing ones.

## Downloading the data

If you download the [example_tutorial.txz](example_tutorial.txz) file and expand it with
```bash
tar Jxvf example_tutorial.txz
```
It will generate a structure similar to below:
```
example_tutorial
├── data
│   ├── borde
│   │   ├── genes.gff
│   │   ├── sequences.fa
│   │   └── snpEffectPredictor.bin
│   └── campy
│       ├── genes.gff
│       ├── sequences.fa
│       └── snpEffectPredictor.bin
├── GCF_000148705.1_ASM14870v1_genomic.fna -> data/campy/sequences.fa
├── GCF_000148705.1_ASM14870v1_genomic.fna.amb
├── GCF_000148705.1_ASM14870v1_genomic.fna.ann
├── GCF_000148705.1_ASM14870v1_genomic.fna.bwt
├── GCF_000148705.1_ASM14870v1_genomic.fna.pac
├── GCF_000148705.1_ASM14870v1_genomic.fna.sa
├── GCF_000148705.1_ASM14870v1_genomic.gff -> data/campy/genes.gff
├── outdir
│   ├── 1.vcf.gz
│   ├── 2.vcf.gz
│   ├── 4.vcf.gz
│   ├── 5.vcf.gz
│   ├── 7.vcf.gz
│   ├── 8.vcf.gz
│   ├── 9.vcf.gz
│   ├── per_sample_average_length.tsv
│   ├── per_sample_modal_frequency.tsv
│   ├── per_sample_proportional_coverage.tsv
│   ├── selected_tracts_annotated.tsv
│   ├── selected_tracts_unknown.tsv
│   ├── tract_list.tsv
│   └── variable_tracts.bed
├── reads
│   ├── ERR1701019.fastq
│   ├── ERR1701029.fastq
│   ├── ERR1701049.fastq
│   ├── ERR1701059.fastq
│   └── ERR1701079.fastq
├── snpEff.config
└── tutorial.md
```

Some files may be missing, as those generated by BWA (`GCF_000148705.1_ASM14870v1_genomic.fna.*`) and tatajuba output,
but they can be reconstructed.

The `data` directory and the `snpEff.config` file are for _snpEff_. The `reads` folder contains (severely subsampled)
read files to be used in this tutorial. The `GCF_000148705.1_ASM14870v1_genomic.*` files are the GFF and FASTA files
used as the reference genome.

## Running tatajuba

### The reference genome files
Tatajuba implements the [BWA](https://github.com/lh3/bwa) library for reference mapping, which relies on index
files it generates from the FASTA file *if they are missing*.
That is, tatajuba does not overwrite the `.amb`, `.ann`, `.bwt` etc. files every time it runs, therefore if you modify
the fasta file please delete these index files so that tatajuba can reconstruct them.

By the way, these index files are generated in the same directory as the fasta file, thus please copy the reference
files to a directory where you have write permissions. 
The FASTA file is not needed if and only if the sequences are embedded within the GFF file, as is the case for
[prokka](https://github.com/tseemann/prokka)'s output.
Tatajuba will then generate a fasta file from it (together with its index files).

For this tutorial, I have downloaded the references from the [Assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000493495.1) 
database on NCBI, since I need both the genomic GFF and the FASTA files with corresponding names.
The [RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/NC_022529.1) database also provides this, as well as obviously prokka
output.
The GFF format talks about contigs or chromosomes, which are the genomic FASTA sequences (genome or plasmid). Thus
tatajuba (and its documentation) use these words interchangeably. 

# Downstream analyses 

## Annotating the VCF files with variant effect information

As we mention in the [README file](../README.md), there are a few tools available for exploring the functional effect of
mutations: 
[bcftools consequence calling](https://samtools.github.io/bcftools/howtos/csq-calling.html), 
[Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html), and 
[snpEff](http://pcingola.github.io/SnpEff/).
They rely on the VCF file for the samples, together with the GFF3 and FASTA files for the reference genome.
We recommend `snpEff` or `VEP`, since `bcftools csq`, supports only ENSEMBL GFF3 files.

Ideally the VCF files should be generated from the BAM/SAM files, i.e. using reference-based assembly programs 
like [minimap2](https://github.com/lh3/minimap2) or [bwa](https://github.com/lh3/bwa).
[Snippy](https://github.com/tseemann/snippy) generates not only the BAM alignment and its corresponding `snps.filt.vcf` files, 
as it annotates the predicted effects with `snpEff` into the `snps.vcf` files. 
Howver keep in mind that the "core SNPs" and [further analyses from snippy-multi](https://github.com/tseemann/snippy#core-snp-phylogeny) 
by design will exclude most information from the HTs (which are indels). 

Tatajuba generates a BED file which can be used to filter the regions harbouring homopolymeric tracts. 
But it also experimentally outputs individual VCF files with minimal information (i.e. no quality information etc.).
These VCF files work well with `snpEff` and `bcftools` but remreber that the gzip format output by tatajuba is **not**
compatible with `bcftools`. 
Conversion is easy though, and if in doubt simply unconpress everything with `gunzip *.gz`. The VCF files are not big.

In this example we will describe how to use `snpEff`. 

### creating a database
For instructions on how to use your own GFF3 file, which is almost always the case, see [this discussion](https://www.biostars.org/p/50963/) and 
also [the official documentation](https://pcingola.github.io/SnpEff/se_buildingreg/#option-1-using-a-gff-file).
In a nutshell, you'll need a configuration file, which we call `snpEff.config`, with locations of the databases.
In our case the relevant lines look like

```
data.dir = ./data/
campy.genome: campy
borde.genome: borde
```

This means that there are two genome databases called "campy" and "borde", which are in fact subdirectories below
`./data/`.
As usual, we learned this [by being one with the Torstenverse](https://github.com/tseemann/snippy). 


Within each of these directories, you'll need one GFF file, named `genes.gff`, with its corresponding FASTA file,
called `sequences.fa`.
If you have been playing with tatajuba you should already have these files :wink:.
Then you tell `snpEff` to generate the database from these files:

```
zcat ../myfiles/GCF_000148705.1_ASM14870v1_genomic.gff.gz > data/campy/genes.gff
zcat ../myfiles/GCF_000148705.1_ASM14870v1_genomic.fna.gz > data/campy/sequences.fa
snpEff build -c snpEff.config -gff3 campy
```

In the example above I copy them from two zipped files from another directory, but if you downloaded the 
[example data set](example_tutorial.txz) then the files should already be there, and with links in the current directory.

### running snpEff

To actually run `snpEff` you need to tell it where the configuration file is, and which database you want to use (in our
case, the freshly created "campy"):

```bash
snpEff ann -c snpEff.config campy concatenated.vcf > concatenated.annotated.vcf
```

In this example we are annotating all variants, across all samples since we are using the file `hand_merged.vcf` [from
above](#merging-all-variants-from-the-vcf-file).


## Concatenating all variants from the VCF file

In this example we assume you want to have a summary of all variants across all samples (to get all possible mutational
effects, for instance). 
One alternative is to [merge all samples into a multi-sample VCF file](#merging-samples-with-bcftools). 
However a simpler alternative is to concatenate all VCF samples by hand, removing the duplicate events:

```
zcat somefile.vcf.gz | head -n 5 > concatenated.vcf
for i in *.vcf.gz; do zcat $i | tail -n +6; done |  sort -nk 2 | uniq >> concatenated.vcf
```
The numbers "5" and "6" above relate to the number of header lines (you nay need to check one file by eye, 
or replace the `head` by a `grep` command). 
You can also [normalise the VCF files](#merging-samples-with-bcftools) before concatenating them by hand.

In a VCF file, after the special keyword headers (starting with `##`), VCF needs a file with the column description. 
If you want to rename the sample to hightlight the fact that this is a dummy sample, a concatenation of all samples,
you thus can modify this line just before the body of the file:
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  all_samples
```

## Merging samples with bcftools

Tatajuba will output one VCF file per sample, but sometimes we are interested in having one file with all variants from
all samples together. 
This can be done with `bcftools merge`, which generated one multi-sample file from several individual ones. 

In the example below we assume that there are several `vcf.gz` files and the reference genome
`GCF_000148705.1_ASM14870v1_genomic.fna` in FASTA format in the same directory:

```bash
# uncompress the files since the standard zlib is not compatible with the expected BGZF format
gunzip *.vcf.gz
# normalise the files to make sure modifications are comparable
for i in *.vcf; do bcftools norm -f GCF_000148705.1_ASM14870v1_genomic.fna ${i} > norm.$i; done
# compress the files using the required BGZF format
for i in norm.*.vcf; do bgzip $i; done
# create an index for each file
for i in norm.*; do bcftools index $i; done
# merge the normalised files into one
bcftools merge -o bcf_merged.vcf norm.*.gz
```

In this example we also use `bcftools norm` to ["fix" each VCF file](https://samtools.github.io/bcftools/bcftools.html#norm)
(by left-aligning indels etc.).