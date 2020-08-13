<img src="recipe/tatajuba-text.png" height="100" alt="Tatajuba">
Distribution of homopolymeric tracts


__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Name
Tatajuba (_Bagassa guianensis_) is a South American tree, also known as Tatajubá, Tatajuva, Garrote, Totajuba.
If I'm not mistaken it means "yellow fire" in [Tupi](https://en.wikipedia.org/wiki/Tupi_language).

## Installation
You should download this repository with `git clone --recursive` to ensure it also downloads
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib).

This sofware uses `autotools`, so you can install it with `configure` + `make`. I'll add more instructions later, but
remember to define where you want it installed with `configure --prefix=DIR`. 
It will compile from the directories `biomcmc-lib`, `kalign`, and `bwa` before finally compiling `tatajuba`.
Notice that this does **not** generate the usual executables for `kalign` or `bwa` (only their libraries are used here).

## License
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2019-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

Tatajubá contains code from [bwa](https://github.com/lh3/bwa) by Heng Li and [kalign](https://github.com/TimoLassmann/kalign.git) by Timo Lassmann,
both released under a GPL-3.0 license.
