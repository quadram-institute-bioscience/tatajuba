/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#ifndef BISECTINGKMEANS_H
#define BISECTINGKMEANS_H

#include "global.h"
#include "rng.h"
#include "alignment_parameters.h"
#include "sequence_distance.h"

void build_tree_kmeans(struct msa* msa, struct aln_param* ap);
#endif
