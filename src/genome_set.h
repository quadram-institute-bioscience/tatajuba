/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief set of homopolymer counters (set of single end or paired end fastq files)
 */

#ifndef _tata_genome_set_h_
#define _tata_genome_set_h_

#define N_SUMMARY_TABLES 5 

#include "context_histogram.h" 

typedef struct genome_set_struct* genome_set_t;
typedef struct g_tract_vector_struct* g_tract_vector_t;

typedef struct
{
  int location, n_dist, lev_distance, id_in_concat; // levenstein distance across genomes in neighbour locations (while still belonging to same tract)
  double *tab0, *d1, *d2; // d2 is pointer to middle of d1; only d1 is alloced/feed || tab0 is new, for tables
  double *gentab[N_SUMMARY_TABLES], reldiff[N_SUMMARY_TABLES]; // reldiff are relative differences of extreme values from genome table
  int n_genome_total, n_genome_id, *genome_id;
  context_histogram_t example;
} g_tract_s;

typedef struct
{
  uint64_t tract_length:10, /*! \brief length of tract in reference */
           max_length:10,   /*! \brief max length of tract amongst samples (used only to alloc memory search) */
           tract_location:27, /*! \brief location (from left/start of contig in fasta) of first tract site */
           contig_location:27, /*! \brief location in contig (from fasta) of beginning of context (2^24=1e8) */ 
           concat_idx:27;   /*! \brief index in concat[] of longest length tract (used as example) */   
  char *contig_name; // contig_name is just a pointer to char_vector
  char *seq;
  char *tract_name; /*! \brief string representation of homopolymer plus context (equiv to context_histogram[]->name) */
} tract_in_reference_s;

struct g_tract_vector_struct
{
  g_tract_s *summary; // summary of tracts over genomes
  context_histogram_t *concat; // tracts pooled over genomes, which can be found in reference
  int n_summary, n_concat;
};

struct genome_set_struct 
{
  genomic_context_list_t *genome;
  g_tract_vector_t tract;
  tract_in_reference_s *tract_ref; /*! \brief all tracts found in reference (idx mapped to concat[]->tract_id) with extra reference info */
  char_vector ref_names; 
  int n_genome, n_tract_ref;
  double secs[4];
  int ref_counter;
};

genome_set_t new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt);
void del_genome_set (genome_set_t g);
void print_selected_g_tract_vector (genome_set_t g);
void print_debug_g_tract_vector (genome_set_t g);
void print_tract_list (genome_set_t g);

#endif
