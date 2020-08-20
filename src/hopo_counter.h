/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _hopo_counter_h_
#define _hopo_counter_h_

#include <kalign.h> // biomcmc is included here
#include <wrapper_bwa.h>

typedef struct hopo_counter_struct* hopo_counter;
typedef struct context_histogram_struct* context_histogram_t;
typedef struct genomic_context_list_struct* genomic_context_list_t;

typedef struct
{ 
  uint64_t context[2]; // flanking kmers (bitstring, not hashed)
  /* bit fields below are signed to faciliate arithm comparisons, thus we lose one bit for signal */
  int64_t base:2,    // base: 0=AT 1=CG (forward or reverse, we use canonical which is A side or C side)  
          length:16, // (former base_size) length of homopolymeric tract (in bases)
          count:32;  // frequency of homopolymer in this context (due to coverage)
} hopo_element;

struct hopo_counter_struct
{
  hopo_element *elem;
  char *name;
  int n_elem, n_alloc, kmer_size;
  double coverage[2], variance[2];
  int *idx, n_idx;
  int ref_counter;
};

struct context_histogram_struct
{
  uint64_t *context; /*! \brief now a vector since we store all within distance */
  int8_t base;       /*! \brief homopolymer base (AT or CG) */
  char *name;        /*! \brief context name is flanking kmers with tract base in the middle */
  int n_context,     /*! \brief vector size (of neighbourhood) */
      integral,      /*! \brief sum of frequencies */
      location,      /*! \brief genomic location(s) of context */
      mode_context_count,   /*! \brief frequency of reads, defining "best homopolymer+context" */
      mode_context_length,  /*! \brief tract length of best homopolymer+context */
      mode_context_id;      /*! \brief which context (from neighboUrhood) has best homopolymer+context */
  int *tmp_count, *tmp_length, n_tmp;
  empfreq h;
};

struct genomic_context_list_struct
{
  context_histogram_t *hist;
  int n_hist, kmer_size;
};

hopo_counter new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, int kmer_size, int min_hopo_size);
void del_hopo_counter (hopo_counter hc);
void finalise_hopo_counter (hopo_counter hc, const char *reference_genome_filename);
void compare_hopo_counters (hopo_counter hc1, hopo_counter hc2, double *result);
void print_debug_hopo_counter (hopo_counter hc);

genomic_context_list_t new_genomic_context_list (hopo_counter hc, int max_distance_per_flank, int min_coverage);
void del_genomic_context_list (genomic_context_list_t genome);
void finalise_genomic_context_hist (genomic_context_list_t genome,  const char *reference_genome_filename);

#endif
