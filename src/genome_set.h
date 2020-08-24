/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief set of homopolymer counters (set of single end or paired end fastq files)
 */

#ifndef _genome_set_h_
#define _genome_set_h_

#include "context_histogram.h" 

typedef struct genome_set_struct* genome_set_t;
typedef struct g_tract_vector_struct* g_tract_vector_t;

typedef struct
{
  int location, n_dist, lev_distance; // levenstein distance across genomes in neighbour locations (while still belonging to same tract)
  double *d1, *d2; // d2 is pointer to middle of d1; only d1 is alloced/feed
  empfreq mode;    // idx = index of genome; freq = model tract length
  context_histogram_t example;
} g_tract_t;

struct g_tract_vector_struct
{
  g_tract_t *distinct;
  context_histogram_t *concat; // tracts pooled over genomes that can be found in reference
  int n_distinct, n_concat;
};

struct genome_set_struct 
{
  genomic_context_list_t *genome;
  g_tract_vector_t tract;
  int n_genome;
  double secs_read, secs_finalise, secs_comparison;
  int ref_counter;
};

genome_set_t new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt);
void del_genome_set (genome_set_t g);
void print_interesting_tracts (genome_set_t g);
void print_debug_g_tract_vector (g_tract_vector_t tract);

#endif
