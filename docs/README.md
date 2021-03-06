# Tutorial

The file [210427.campy_bordetella.ipynb](210427.campy_bordetella.ipynb) contains the `jupyter` notebook for the
*Campylobacter* and *Bordetella* homopolymer tract analyses for the manuscript
https://doi.org/10.1101/2021.06.02.446710. 

# Usage

the basic help can be invoked with `tatajuba -h`, which will look something like

```
tatajuba 1.0.2
Compare histograms of homopolymeric tract lengths, within context.
The complete syntax is:

 tatajuba  [-h|--help] [-v|--version] [-p|--paired] [-k|--kmer={2,...,32}] [-m|--minsize={1,...,32}] [-i|--minreads=<int>] [-d|--maxdist=<int>] [-l|--leven=<int>] [-t|--nthreads=<int>] -g|--gff=<genome.gff3|genome.gff3.gz> [--fasta=<genome.fna|genome.fna.gz>] <fastq files> [<fastq files>]... [-o|--outdir=<file>]

  -h, --help                       print a longer help and exit
  -v, --version                    print version and exit
  -p, --paired                     paired end (pairs of) files
  -k, --kmer={2,...,32}            kmer size flanking each side of homopolymer (default=25)
  -m, --minsize={1,...,32}         minimum homopolymer tract length to be compared (default=4)
  -i, --minreads=<int>             minimum number of reads for tract+context to be considered (default=5)
  -d, --maxdist=<int>              maximum distance between kmers of a flanking region to merge them into one context (default=1)
  -l, --leven=<int>                levenshtein distance between flanking regions to merge them into one context (after ref genome mapping)
  -t, --nthreads=<int>             suggested number of threads (default is to let system decide; I may not honour your suggestion btw)
  -g, --gff=<genome.gff3|genome.gff3.gz> reference genome file in GFF3, preferencially with sequence
  --fasta=<genome.fna|genome.fna.gz> reference genome file in fasta format, if absent from GFF3
  <fastq files>                    fastq file with reads (weirdly, fasta also possible as long as contains all reads and not only contigs)
  -o, --outdir=<file>              output directory, or 'random' for generating random dir name (default=current dir '.')
'Context' is the pair of flanking k-mers. If you have paired end files, then their order should be strictly f1_R1, f1_R2, f2_R1 etc.
  that is, R1 and R2 should be consecutive. A GFF3 reference genome must also be supplied, and bwa will create a series of index files
  if tatajuba can find the DNA sequences at end of GFF3 file, after pragma '##FASTA'

You can also supply a fasta file, e.g. when you download GFF3 from NCBI it may lack the DNA genome at the end --- then downloading the
  associated fna may work (untested). The internal library updates the GFF3 structure with provided fasta sequences (so a fasta file may
  overwrite DNA sequences present in the GFF3 file. If you modify the fasta file please delete the index files generated by the BWA library
  so that they can be regenerated with the updated information.

Notice that tatajuba creates files with same prefix and in same location as the GFF3 file, which may overwrite existing ones.
  On the other hand, as suggested above, you can recreate all the generated files by deleting them and running tatajuba again.

The default values are as following:
kmer= 25         minsize=  4     minread=  5     maxdist=  1
```

The minimum required information are the GFF3 file (assuming it has fasta sequences, as e.g. generated by `prokka`) and the reads. Usually the GFF3
file will not contain the sequence in fasta format, and thus you'll have to provide the fasta file with `--fasta` for
the exact same GFF3, with matching sequence names. 
Also, if you have paired reads you need to add option `-p` and make sure the read files are given in order, with reads
from same sample together: 
always `sampleA_1.fq sampleA_2.fq sampleB_1.fq sampleB_2.fq` but never
`sampleA_1.fq sampleB_1.fq sampleA_2.fq sampleB_2.fq`.
