/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

/*! \file
 *  \brief library based on [kalign3 multiple sequence alignment program](https://github.com/TimoLassmann/kalign.git)
 */

#ifndef _cumaru_kalign_h
#define _cumaru_kalign_h

#include "alignment.h"
#include "alignment_parameters.h"
#include "bisectingKmeans.h"
#include "bpm.h"
#include "global.h"
#include "rng.h"
#include "sequence_distance.h"
//#include "tldevel.h"

char_vector kalign3_from_char_vector (char_vector dna);
#endif
