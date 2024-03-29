/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _context_histogram_h_
#define _context_histogram_h_

#include "hopo_counter.h" 
#define CH_MAX_DIST 0xffff

typedef struct context_histogram_struct* context_histogram_t;
typedef struct genomic_context_list_struct* genomic_context_list_t;

struct context_histogram_struct
{
  uint64_t *context;    /*! \brief now a vector since we store all within distance */
  int32_t base:2,       /*! \brief homopolymer base (AT or CG) */
         multi:3,       /*! \brief  0,1,2 (from hopo_counter) if all (1) or some (2) contexts have paralogs */
         indel:2,       /*! \brief if indel is closer than existing distances for at least one context pair */
         neg_strand:1,  /*! \brief strand in the _reference_ genome */
         mismatches:12; /*! \brief edit (Levenstein) distance between modal HT for this sample and reference */
  char *name;        /*! \brief context name is flanking kmers with tract base in the middle */
  int n_context,     /*! \brief vector size (of neighbourhood) */
      integral,      /*! \brief sum of frequencies */
      location,      /*! \brief genomic location(s) of context in 1D flattened contigs/genomes */
      loc2d[3],      /*! \brief 2D genomic location [ref genome, location within genome, last location within genome] */
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
  char *name; /*! \brief simplified sample name, excluding prefix (directory) and suffix (e.g. ".fastq.gz") */
  tatajuba_options_t opt;
  int n_hist, coverage, ref_start;  /*! \brief ref_start is index of first hist found on ref genome (all before were not found) */
};

int indel_distance_between_context_histogram_and_hopo_context (context_histogram_t ch, char *name);
int distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int location_difference, int *idx_match);
int distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result);
int compare_context_histogram_for_qsort (const void *a, const void *b);
bool context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt);
/*! \brief store lengths of prefix and suffix in common into l[]. Returns true if seqs are distinct, false o.w. */
bool common_prefix_suffix_lengths_from_strings (char *s1, size_t n1, char *s2, size_t n2, int *l);

void del_context_histogram (context_histogram_t ch);
void print_debug_genomic_context_hist (genomic_context_list_t genome);
genomic_context_list_t  new_genomic_context_list (hopo_counter hc);
void del_genomic_context_list (genomic_context_list_t genome);
void finalise_genomic_context_hist (genomic_context_list_t genome);

int distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result); // return is not distance
bool context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt);

#endif
