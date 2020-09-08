/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _context_histogram_h_
#define _context_histogram_h_

#include <kalign.h> // biomcmc is included here
#include <wrapper_bwa.h>

extern uint8_t dna_in_2_bits[256][2];
extern char bit_2_dna[];

typedef struct hopo_counter_struct* hopo_counter;
typedef struct context_histogram_struct* context_histogram_t;
typedef struct genomic_context_list_struct* genomic_context_list_t;

typedef struct 
{
  char *reference_genome_filename; 
  bool paired_end;
  gff3_t gff;
  int max_distance_per_flank, 
      kmer_size,
      min_tract_size,
      levenshtein_distance,
      min_coverage,
      n_threads;
} tatajuba_options_t;

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
  int n_elem, n_alloc, kmer_size, coverage;
  int *idx, n_idx;
  tatajuba_options_t opt;
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
      coverage, n_tracts,   /*! values related to genome (not histogram) */
      mode_context_count,   /*! \brief frequency of reads, defining "best homopolymer+context" */
      mode_context_length,  /*! \brief tract length of best homopolymer+context */
      mode_context_id;      /*! \brief which context (from neighbourhood) has best homopolymer+context */
  int *tmp_count, *tmp_length, index; // index in genome_set, first used temporarily as counter
  empfreq h; // h.idx = tract length; h.freq = count (number of reads supporting this length)
  int ref_counter;
};

struct genomic_context_list_struct
{
  context_histogram_t *hist;
  char *name;
  tatajuba_options_t opt;
  int n_hist, coverage, ref_start;  /*! \brief ref_start is index of first hist found on ref genome (all before were not found) */
};

void print_tatajuba_options (tatajuba_options_t opt);
hopo_counter new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, tatajuba_options_t opt);
void del_hopo_counter (hopo_counter hc);
void del_context_histogram (context_histogram_t ch);

void print_debug_genomic_context_hist (genomic_context_list_t genome);
genomic_context_list_t  new_genomic_context_list (hopo_counter hc);
void del_genomic_context_list (genomic_context_list_t genome);
void finalise_genomic_context_hist (genomic_context_list_t genome);

int distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result); // return is not distance
bool context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt);

#endif
