/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#ifndef RNG_H
#define RNG_H

#include "tldevel.h"
#include "stdint.h"

struct rng_state{
        uint64_t s[4];
};

extern struct rng_state* init_rng(uint64_t seed);
extern void free_rng(struct rng_state* rng);

extern double tl_random_double(struct rng_state* rng);
extern int tl_random_int(struct rng_state* rng,int a);

#endif
