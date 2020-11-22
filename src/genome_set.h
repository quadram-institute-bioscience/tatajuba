/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief set of homopolymer counters (set of single end or paired end fastq files)
 */

#ifndef _genome_set_h_
#define _genome_set_h_

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
  int tract_length;
} tract_in_reference_s;

struct g_tract_vector_struct
{
  g_tract_s *summary; // summary of tracts over genomes
  context_histogram_t *concat; // tracts pooled over genomes, which can be found in reference
  tract_in_reference_t *ref;
  int n_summary, n_concat, n_trait;
};

struct genome_set_struct 
{
  genomic_context_list_t *genome;
  g_tract_vector_t tract;
  int n_genome;
  double secs[3];
  int ref_counter;
};

genome_set_t new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt);
void del_genome_set (genome_set_t g);
void print_selected_g_tract_vector (genome_set_t g);
void print_debug_g_tract_vector (genome_set_t g);

#endif
