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

struct genome_set_struct 
{
  genomic_context_list_t *genome;
  context_histogram_t *tract; // tracts pooled over genomes that can be found in reference
  int n_genome, n_tract;
  double secs_read, secs_finalise, secs_comparison;
  distance_generator generator;
  int ref_counter;
};

genome_set_t new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt);
void del_genome_set (genome_set_t g);
// distance_generator new_distance_generator_from_hopo_set (hopo_set hs);

#endif
