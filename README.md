<img src="recipe/tatajuba-text.png" height="100" alt="Tatajuba">

__Leonardo de Oliveira Martins<sup>1</sup>__, 
__Samuel Bloomfield<sup>1</sup>__,
__Emily Stoakes<sup>2</sup>__,
__Andrew Grant<sup>2</sup>__,
__Andrew Page<sup>1</sup>__,
__Alison Mather<sup>1</sup>__
<br>
<sub> 1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK; </sub>  
<sub> 2. Department of Veterinary Medicine, University of Cambridge, Madingley Road, Cambridge, CB3 0ES</sub>

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/tatajuba/blob/master/LICENSE)

# Distribution of homopolymeric tracts 

Instead of assuming a fixed length for a given homopolymer tract, tatajubá allows for the whole distribution of tract sizes to
be analysed. 
The rationale is that 1. our sequence might represent a population of non-identical individuals, with diversity of tract
lengths, and 2. sequencing errors might be more frequent near or within homopolymers (so we should not remove
uncertainty prematurely).

Tatajubá also assumes that what we call a "tract" is a homopolymeric base flanked by a specific sequence (allowing for
variability), and it can discard homopolymers absent in reverse or forward reads to minimise "strand bias". 

#### Paper and citation
"Tatajuba ― Exploring the distribution of homopolymer tracts", 
Leonardo de Oliveira Martins, Samuel Bloomfield, Emily Stoakes, Andrew Grant, Andrew J Page, Alison E Mather, (2021),
bioRxiv 2021.06.02.446710; doi: https://doi.org/10.1101/2021.06.02.446710

#### Name
Tatajuba (_Bagassa guianensis_) is a South American tree, also known as Tatajubá, Tatajuva, Garrote, Totajuba.
It means "yellow fire" in [Tupi](https://en.wikipedia.org/wiki/Tupi_language).

## Installation

Currently the software has been tested exclusively on linux systems, but hopefully you can run it on other systems through the [singularity](#singularity)
and [docker](#docker) containers. 
If you have any tips for successfull usage in other systems, [do let us know](https://github.com/quadram-institute-bioscience/tatajuba/issues).

### Conda
[![Anaconda-Server Badge](https://anaconda.org/bioconda/tatajuba/badges/platforms.svg)](https://anaconda.org/bioconda/tatajuba)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/tatajuba/badges/latest_release_date.svg)](https://anaconda.org/bioconda/tatajuba)

After you install [miniconda](https://conda.io/en/latest/miniconda.html), simply run
```[bash]
conda install -c bioconda tatajuba
```
The software tatajuba is still under development, thus the conda version may be outdated.

### Singularity

After installing [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html), you can download an executable
container with:
```bash
# check https://cloud.sylabs.io/library/leomrtns/default/tatajuba for most recent tag)
singularity pull --arch amd64 library://leomrtns/default/tatajuba:1.0.3 
```
As with conda above, the container might not have the latest improvements. In case you want the most recent version, you can
use the singularity definition file [recipe/tatajuba.def](recipe/tatajuba.def) to generate a container as in 

```bash
sudo singularity build tatajuba.sif recipe/tatajuba.def
```
If you build the container as above, the software will be up-to-date since it will download from github and compile.

### Docker
```bash
# check https://quay.io/repository/biocontainers/tatajuba?tab=tags for most recent tag 
docker pull quay.io/biocontainers/tatajuba:1.0.3--h82fb495_0
```
This is generated from bioconda, so same caveats apply.

### Compiling from source
If installing through conda/singularity is not an option, or if you want the latest version of the 
software, you can download it and compile it yourself. 
Tatajuba relies on GCC6 or newer due to assuming *OpenMP 4.5*.
This repository must be cloned with `git clone --recursive` to ensure it also downloads
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) and [our modified version of BWA](https://github.com/leomrtns/bwa).
You will need a recent version of GCC in your system. 

This sofware uses `autotools`, so you can install it with `configure` + `make`. 
You may need to define where you want it installed with `configure --prefix=DIR` which is where are your unix-like
`include/`, `lib/`, and `bin/` directories. My favourite is `~/local`. 

It will compile from the directories `biomcmc-lib`, `kalign`, and `bwa` before finally compiling `tatajuba`.
Notice that this does **not** generate the usual executables for `kalign` or `bwa`: only their libraries are used by
tatajubá.

Here is an example of its installation, please modify to better suit your needs:

```[bash]
/home/simpson/$ git clone --recursive https://github.com/quadram-institute-bioscience/tatajuba.git
/home/simpson/$ cd tatajuba && ./autogen.sh 
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../configure --prefix=${HOME}/local ## prefix is the location of your local libraries etc.
/home/simpson/$ make; make install
```
If it works, you should have `tatajuba` installed in the `${HOME}/local/bin` directory, in the example above (or
whatever you set as the `-prefix`).
You may want to [add this path to your `$PATH` variable](https://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path).

If you want, you can optionally check the installation by running a battery of unit and integration tests for both tatajuba and biomcmc-lib (the low-level C library `tatajuba` relies on):

```[bash]
/home/simpson/$ sudo apt-get install check  # preferred method, assuming you have admin priviledges on the ubuntu/debian machine
/home/simpson/$ # conda install -c conda-forge check  # alternative to apt-get get above, using conda
/home/simpson/$ make check
```

If `configure` complains about a missing library (usually `libcheck` or `zlib`), you'll need to install them before 
running `configure` again.
You will also need the `autotools` environment before running the configuration (`autogen.sh` depends on it):

```[bash]
## 'bootstrap' the configuration files (needed when cloning from github):
/home/simpson/$ apt-get install pkg-config autotools-dev autoconf automake libtool
/home/simpson/$ (cd tatajuba && autoreconf)  ## the parentheses avoid entering the directory afterwards

## install libraries possibly missing (zlib and omp are strongly suggested)
/home/simpson/$ apt-get install zlib1g-dev libomp-dev libbz2-dev check liblzma-dev
```
The libraries rely on `pkg-config` to find their location: if your `pkg-config` was installed through conda then you'd
better install the above libs via conda as well (or, you know, updating environmental variables etc).
The `zlib` library is mandatory, while `liblzma-dev` and `libbz2-dev` are called, respectively, `xz` and `bzip2` on a strict conda environment.

## Documentation 
You can find the documentation on its usage in the [`docs` folder](docs/README.md).

The program gives summary help with `tatajuba` (without arguments) and detailed help with `tatajuba -h`.
Please feel free [to report any issues](https://github.com/quadram-institute-bioscience/tatajuba/issues) or to request
clarification if anything is not clear. 

## Model
At the lowest level (C `struct`), the homopolymeric tracts are stored as the two flanking k-mers (called "context" here) and the base
comprising the homopolymer in th middle, as seen  in the figure below. 

<img src="recipe/200322_001.png" height="160" alt="context_struct" align="middle">

We define the canonical form based on the homopolymer &mdash; in the figure above the same flanking regions `CCG` and
`GAT` are stored as a completely different context b/c they flank a distinct homopolymer base. The three contextualised
tracts above are displayed by tatajubá as
```
CCG-A-GAT
ATC-A-CCG
CCG-C-GAT
```
due to the canon, we always observe the side of the homopolymers with `A` or with `C`.

Scanning through the fastq files, we now can, for each sample, generate the histograms of contextualised homopolymeric tract lengths as
depicted in the figure below.

<img src="recipe/200322_002.png" height="160" alt="context-histogram" align="middle">

Once this histogram is complete we search for this tract (i.e. homopolymer plus flanking regions) on the reference
genome, by using a typical length (a typical length would be _3_ for the figure above).
Tatajubá also tries to merge histograms if they may represent the same tract both before and after the reference genome
mapping.
*Before* mapping it tries to find contexts that are quite similar (and thus could represent the same tract). 
The parameter `maxdist` will control up to how may mismatches (per flanking region) are considered the same context.
*After* mapping we may notice very close tracts, which may in fact be the same tract but with indels in the flanking
regions. 
The parameter `leven` decides the maximum Levenshtein distance between contexts for such neighbouring tracts to be
considered the same. 
If in your results you see tracts just a few bases apart, try increasing the `leven` value. My suggestion is that both
parameters are kept as low as possible (less than two). 
As a side note, the Levenshtein distance is slower to calculate, while the pre-mapping mismatch has to be done between
all pairs and not only neighbours.

By the way, histograms with very low frequency (representing contexts+tracts observed very rarely in the fastq file) are
excluded, assuming they represent sequencing errors. This is controlled by the parameter `minreads`. The default is
currenlty _3_ (any tract observed in less than _3_ reads is discarded).

Currently our measure of dispersion (used to find tracts most variable across genomes) is the *relative difference of
ranges* (similar to the coefficient of range), defined here as (MAX-MIN)/MAX.

## Troubleshooting

As mentioned, this software is still under development; in particular the conda/singularity/docker versions might be
outdated. Here is a list of common pitfalls.

* tatajuba relies on OpenMP 4.5, which is supported [on GCC6 or newer](https://www.openmp.org/resources/openmp-compilers-tools/). 
  This means that [even the conda version might fail if your system library is older than that
](https://stackoverflow.com/questions/62098781/is-it-possible-to-use-a-different-gcc-version-inside-a-conda-environment#comment109969669_62098781).
* We have tested it exclusively on linux systems, and making it more portable to Mac or Windows is not a priority.
* the program will refuse to run on only one sample (since it compares differences between samples). However, it will
  run if more samples are given but only one has mapped HTs &mdash; in practice, equivalent to running on one sample since
  it discards samples without any HTs mapped to the reference. In this case the program will fail at the very end (you
  may have some output). I am treating this as a bug.
* The program should produce error messages; however I've seen it failing without notice. One particular case is when it
  runs out of memory (it is killed by the system). It needs at least 8GB of memory. 
* The files `selected_tracts_{annotated/unknown}.tsv` are being used for debug purposes, but will be soon replaced by a
  more useful `per_sample`-like file. In particular the locations do not correspond to the `tract_id` locations &mdash; if you
  currently want to use these files, please use the `tract_id` for mapping to the correct locations (available in files
  `per_sample*` or `tract_llist.tsv`). 
* As of 2020.08.01, the conda/singularity versions (1.0.3) use a lot of memory. This has been fixed if you use the source
  code.

Please use github [to report any issues](https://github.com/quadram-institute-bioscience/tatajuba/issues), and to see
which issues are being addressed.

## License
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2020-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

Tatajubá contains code from [bwa](https://github.com/lh3/bwa) by Heng Li and [kalign](https://github.com/TimoLassmann/kalign.git) by Timo Lassmann,
both released under a GPL-3.0 license.
(We do not currently compile or use this modified kalign, btw)

![Anurag's github stats](https://github-readme-stats.vercel.app/api?username=leomrtns&count_private=true&show_icons=true&theme=calm)
