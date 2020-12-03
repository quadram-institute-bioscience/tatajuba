/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _context_histogram_h_
#define _context_histogram_h_

#include "hopo_counter.h" 

typedef struct context_histogram_struct* context_histogram_t;
typedef struct genomic_context_list_struct* genomic_context_list_t;

struct context_histogram_struct
{
  uint64_t *context; /*! \brief now a vector since we store all within distance */
  int8_t base;       /*! \brief homopolymer base (AT or CG) */
  char *name;        /*! \brief context name is flanking kmers with tract base in the middle */
  int n_context,     /*! \brief vector size (of neighbourhood) */
      integral,      /*! \brief sum of frequencies */
      location,      /*! \brief genomic location(s) of context in 1D flattened contigs/genomes */
      loc2d[2],      /*! \brief 2D genomic location [ref genome, location within genome] */
      coverage, n_tracts,   /*! values related to genome (not histogram) */
      mode_context_count,   /*! \brief frequency of reads, defining "best homopolymer+context" */
      mode_context_length,  /*! \brief tract length of best homopolymer+context */
      mode_context_id;      /*! \brief which context (from neighbourhood) has best homopolymer+context */
  int *tmp_count, *tmp_length, index; // index in genome_set (i.e. fastq pair), first used temporarily as counter
  empfreq h; // h.idx = tract length; h.freq = count (number of reads supporting this length)
  gff3_fields gffeature; /*! \brief could be a list (even for a single position on a single genome) but here we store first belonging */
  int tract_id; /*! \brief tract_id used in genome->concat, can be same for several contexts if similar enough */
  int ref_counter;
};

struct genomic_context_list_struct
{
  context_histogram_t *hist;
  char *name;
  tatajuba_options_t opt;
  int n_hist, coverage, ref_start;  /*! \brief ref_start is index of first hist found on ref genome (all before were not found) */
};

int distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int *idx_match);
int distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result);
int compare_context_histogram_for_qsort (const void *a, const void *b);
bool context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt);

void del_context_histogram (context_histogram_t ch);
void print_debug_genomic_context_hist (genomic_context_list_t genome);
genomic_context_list_t  new_genomic_context_list (hopo_counter hc);
void del_genomic_context_list (genomic_context_list_t genome);
void finalise_genomic_context_hist (genomic_context_list_t genome);

int distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result); // return is not distance
bool context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt);

#endif
