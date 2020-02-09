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
{ // in future I may have vect[3] to help memory alignment
  uint64_t context; // flanking kmers (bitstring or hashed)
  uint8_t base;     // 0=AT 1=CG (forward or reverse, we use canonical)  
  int base_size;    // length of homopolymeric tract (in bases)
  int count;          // frequency of homopolymer in this context (due to coverage)
} hopo_element;

struct hopo_counter_struct
{
  hopo_element *elem;
  int n_elem, n_alloc, coverage[2];
  int ref_counter;
};

hopo_counter new_hopo_counter_from_file (const char *filename);
void del_hopo_counter (hopo_counter hc);
void print_debug_hopo_counter (hopo_counter hc);


#endif
