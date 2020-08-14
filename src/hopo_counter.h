/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _hopo_counter_h_
#define _hopo_counter_h_

#include <kalign.h>

typedef struct hopo_counter_struct* hopo_counter;

typedef struct
{ 
  uint64_t context[2]; // flanking kmers (bitstring, not hashed)
  /* bit fields below are signed to faciliate arithm comparisons, thus we lose one bit for signal */
  int64_t base:2,    // base: 0=AT 1=CG (forward or reverse, we use canonical which is A side or C side)  
          length:8,  // (former base_size) length of homopolymeric tract (in bases)
          count:24,  // frequency of homopolymer in this context (due to coverage)
          extra:30;  // unused |  7bit=128; 23bit=8mi 
} hopo_element;

struct hopo_counter_struct
{
  hopo_element *elem;
  char *name;
  int n_elem, n_alloc;
  double coverage[2], variance[2];
  int *idx, n_idx;
  int ref_counter;
};

hopo_counter new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, int kmer_size, int min_hopo_size);
void del_hopo_counter (hopo_counter hc);
void finalise_hopo_counter (hopo_counter hc);
void compare_hopo_counters (hopo_counter hc1, hopo_counter hc2, double *result);
void print_debug_hopo_counter (hopo_counter hc);

#endif
