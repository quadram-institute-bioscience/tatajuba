<img src="recipe/tatajuba-text.png" height="100" alt="Tatajuba">

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

# Distribution of homopolymeric tracts

Instead of assuming a fixed length for a given homopolymer tract, tatajubá allows for its whole distribution of sizes to
be analysed. 
The rationale is that 1. our sequence might represent a population of non-identical individuals, with diversity of tract
lengths, and 2. sequencing errors might be more frequent near or within homopolymers (so we should not remove
uncertainty prematurely).

Tatajubá also assumes that what we call a "tract" is a homopolymeric base flanked by a specific sequence (allowing for
variability).

#### Name
Tatajuba (_Bagassa guianensis_) is a South American tree, also known as Tatajubá, Tatajuva, Garrote, Totajuba.
It means "yellow fire" in [Tupi](https://en.wikipedia.org/wiki/Tupi_language).

## Installation
You should download this repository with `git clone --recursive` to ensure it also downloads
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib).

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
/home/simpson/$ ../tatajuba/configure --prefix=${HOME}/local
/home/simpson/$ make; make install
/home/simpson/$ make check  # battery of unit and integration tests for both tatajuba and biomcmc-lib
```

If `configure` complains about a missing library (usually `libcheck` or `zlib`), you'll need to install them before 
running `configure` again.
If there is no `configure` file at all in the distribution, or it
You will most likely need to install the `autotools` before running the configuration (`autogen.sh` depends on it). 
Both cases are shown below, if you can install them system-wide:

```[bash]
## 'bootstrap' the configuration files (needed when cloning from github):
/home/simpson/$ apt-get install pkg-config autotools-dev autoconf automake libtool
/home/simpson/$ (cd tatajuba && autoreconf)  ## the parentheses avoid entering the directory afterwards

## install libraries possibly missing (only check is mandatory, but zlib and omp are strongly suggested)
/home/simpson/$ apt-get install zlib1g-dev libomp-dev libbz2-dev check liblzma-dev
```
The libraries rely on `pkg-config` to find their location: if your `pkg-config` was installed through conda then you'd
better install the above libs via conda as well (or, you know, updating environmental variables etc)

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

## License
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

Tatajubá contains code from [bwa](https://github.com/lh3/bwa) by Heng Li and [kalign](https://github.com/TimoLassmann/kalign.git) by Timo Lassmann,
both released under a GPL-3.0 license.
