/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#ifndef BPM_H
#define BPM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

/* Must be called before bpm_256!!!!  */
extern void set_broadcast_mask(void);
uint8_t bpm_256 (const uint8_t* t,const uint8_t* p, uint32_t n, uint32_t m);
uint8_t bpm (const uint8_t* t,const uint8_t* p, uint32_t n, uint32_t m);

#endif
