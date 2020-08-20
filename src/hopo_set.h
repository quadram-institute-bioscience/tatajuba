/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief set of homopolymer counters (set of single end or paired end fastq files)
 */

#ifndef _hopo_set_h_
#define _hopo_set_h_

#include "hopo_counter.h" 

typedef struct genome_set_struct* genome_set_t;

struct genome_set_struct 
{
  genomic_context_list_t *genome;
  int n_genome;
  double secs_read, secs_finalise, secs_comparison;
  distance_generator generator;
  int ref_counter;
};

genome_set new_genome_set_from_files (const char **filenames, int n_filenames, bool paired_end, int kmer_size, int min_hopo_size, const char *reference_genome_filename);
void del_hopo_set (hopo_set hs);
distance_generator new_distance_generator_from_hopo_set (hopo_set hs);

#endif
