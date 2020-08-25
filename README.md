<img src="recipe/tatajuba-text.png" height="100" alt="Tatajuba">
### Distribution of homopolymeric tracts


__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Name
Tatajuba (_Bagassa guianensis_) is a South American tree, also known as Tatajubá, Tatajuva, Garrote, Totajuba.
If I'm not mistaken it means "yellow fire" in [Tupi](https://en.wikipedia.org/wiki/Tupi_language).

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
/home/simpson/$ git clone --recursive git@github.com:quadram-institute-bioscience/tatajuba.git
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../tatajuba/configure --prefix=${HOME}/local
/home/simpson/$ make; make install
/home/simpson/$ make check  # battery of unit and integration tests for both tatajuba and biomcmc-lib
```

If `configure` complains about a missing library (usually `libcheck` or `zlib`), you'll need to install them before 
running `configure` again.
If there is no `configure` file at all in the distribution, or it doesn't run, then you you need to install the
`autotools` and rerun the configuration. 
Both cases are shown below:

```[bash]
## install libraries possibly missing:
/home/simpson/$ apt-get install zlib1g-dev check

## 'bootstrap' the configuration files:
/home/simpson/$ apt-get install automake autoconf libtool 
/home/simpson/$ (cd tatajuba && autoreconf)  ## the parentheses avoid entering the directory afterwards
```

## Model
At the lowest level, the homopolymeric tracts are stored as the two flanking k-mers (called "context" here) and the base
comprising the homopolymer in th middle, as seen  in the figure below. 

<img src="recipe/200322_001.png" height="160" alt="Tatajuba">

We define the canonical form based on the homopolymer &mdash; in the figure above the same flanking regions `CCG` and
`GAT` are stored as a completely different context b/c they flank a distinct homopolymer base. The three contextualised
tracts above are displayed by tatajubá as
```
CCG-A-GAT
ATC-A-CCG
CCG-C-GAT
```
due to the canon, we only observe `A` or `C` as the homopolymers. 


Scanning through the fastq files, we can now for each sample generate a histogram of homopolymeric tract lengths,
depicted in the figure below.
<img src="recipe/200322_002.png" height="100" alt="Tatajuba">

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
