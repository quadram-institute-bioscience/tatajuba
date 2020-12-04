/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _hopo_counter_h
#define _hopo_counter_h

#include <biomcmc.h>
#include <wrapper_bwa.h>

extern uint8_t dna_in_2_bits[256][2];
extern char bit_2_dna[];

typedef struct hopo_counter_struct* hopo_counter;

typedef struct 
{
  char *reference_fasta_filename, *outdir; 
  bool paired_end;
  gff3_t gff;
  int max_distance_per_flank, 
      kmer_size,
      min_tract_size,
      levenshtein_distance,
      min_coverage,
      n_samples,
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

int compare_hopo_element_decreasing (const void *a, const void *b);
int compare_hopo_context (hopo_element a, hopo_element b);
int distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist);

void print_tatajuba_options (tatajuba_options_t opt);
hopo_counter new_hopo_counter (int kmer_size);
hopo_counter new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, tatajuba_options_t opt);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_tract_size);
void del_hopo_counter (hopo_counter hc);
char* leftmost_hopo_name_and_length_from_string (char *seq, size_t len, int kmer_size, int min_tract_size, int *tract_length);

int hopo_counter_histogram_integral (hopo_counter hc, int start);
void finalise_hopo_counter (hopo_counter hc);
char* generate_name_from_flanking_contexts (uint64_t *context, int8_t base, int kmer_size);

#endif
