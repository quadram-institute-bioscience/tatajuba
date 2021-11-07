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
  bool paired_end, remove_biased;
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
  uint64_t context[2]; /*! \brief flanking kmers (bitstring, not hashed) */
  /* bit fields below are signed to faciliate arithm comparisons, thus we lose one bit for signal */
  int64_t base:2,    /*! \brief base: 0=AT 1=CG (forward or reverse, we use canonical which is A side or C side) */
          length:10, /*! \brief  (former base_size) length of homopolymeric tract (in bases) */
          count:20,  /*! \brief frequency of homopolymer in this context (due to coverage) */
          mismatches:12,  /*! \brief edit distance NM (could also be mismatches n_mm plus indels) from bwa */
          multi:3,        /*! \brief more than one match */
          neg_strand:2,   /*! \brief if maps to neg strand of the reference genome (_not_ read strand) */
          canon_flag:3; /*! \brief original read in canonical form (1) or not (2); if tract was seen in both, then both forw and rev are present: flag=3 */
  int32_t read_offset, /*! \brief at beginning, used as start of context+tract in read (when searching in reference fasta); later, becomes 1D flattened location from bwa */
          loc_ref_id,
          loc_pos,     /*! \brief 2D BWA location [ref_id,position] which are ref sequence ID and site position within this refseq */
          loc_last;    /*! \brief position of last base in reference matching HT (so that we can have the whole mathing segment in reference) */
} hopo_element;

struct hopo_counter_struct
{
  hopo_element *elem;
  char *name;
  int ref_start, n_elem, n_alloc, kmer_size, coverage;
  int *idx_initial, *idx_final, n_idx;
  tatajuba_options_t opt;
  int ref_counter;
};

int compare_hopo_element_decreasing (const void *a, const void *b);
int compare_hopo_context (hopo_element a, hopo_element b);
int distance_between_single_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist);
int distance_between_context_kmer_pair (uint64_t *c1, uint64_t *c2);
/*! \brief short edit distance between contexts by bitwise shift of contexts; stores best shift operation per context */
int distance_between_context_kmer_pair_with_edit_shift (uint64_t *c1, uint64_t *c2, int *best_shift); 

void print_tatajuba_options (tatajuba_options_t opt);
hopo_counter new_hopo_counter (int kmer_size);
hopo_counter new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, tatajuba_options_t opt);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_tract_size);
/*! \brief creates a dummy counter, with all monomers (i.e. all contexts and single bases in middle), for cases no HT is found within seq[] */
void update_hopo_counter_from_seq_all_monomers (hopo_counter hc, char *seq, int seq_length);
void del_hopo_counter (hopo_counter hc);
char* leftmost_hopo_name_and_length_from_string (char *seq, size_t len, int kmer_size, int min_tract_size, int *tract_length);

int hopo_counter_histogram_integral (hopo_counter hc, int start);
void finalise_hopo_counter (hopo_counter hc);
char* generate_tract_as_string (uint64_t *context, int8_t base, int kmer_size, int tract_length, bool neg_strand);
char* generate_name_from_flanking_contexts (uint64_t *context, int8_t base, int kmer_size, bool neg_strand);
char* protein_from_dna_string (char *dna, size_t n_dna, bool reverse);

#endif
