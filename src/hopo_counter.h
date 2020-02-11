/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _hopo_counter_h_
#define _hopo_counter_h_

#include <kalign.h>

# define MAXFILENAMELENGTH 64

typedef struct hopo_counter_struct* hopo_counter;

typedef struct
{ // in future I may have vect[3] to help memory alignment
  uint64_t context; // flanking kmers (bitstring or hashed)
  uint8_t base;     // 0=AT 1=CG (forward or reverse, we use canonical)  
  int base_size;    // length of homopolymeric tract (in bases)
  int count;          // frequency of homopolymer in this context (due to coverage)
} hopo_element;

struct hopo_counter_struct
{
  hopo_element *elem;
  char name[MAXFILENAMELENGTH+1];
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
