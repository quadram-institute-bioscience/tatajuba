## Annotating the VCF files with variant effect information

As we mention in[README file](../README.md), there are a few tools available for exploring the functional effect of
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
