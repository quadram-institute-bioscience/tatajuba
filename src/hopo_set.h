/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief set of homopolymer counters (set of single end or paired end fastq files)
 */

#ifndef _hopo_set_h_
#define _hopo_set_h_

#include "hopo_counter.h" 


typedef struct hopo_set_struct* hopo_set;

struct hopo_set_struct
{
  hopo_counter *hc;
  int n_hc;
  double secs_read, secs_finalise, secs_comparison;
  distance_generator generator;
  int ref_counter;
};

hopo_set new_hopo_set_from_files (const char **filenames, int n_filenames, bool paired_end, int kmer_size, int min_hopo_size);
void del_hopo_set (hopo_set hs);
distance_generator new_distance_generator_from_hopo_set (hopo_set hs);


#endif
