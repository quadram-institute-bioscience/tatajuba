/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "hopo_set.h"
  
void hopo_set_distance_wrapper (void *data, int s1, int s2, double *result);

hopo_set
new_hopo_set_from_files (const char **filenames, int n_filenames, bool paired_end, int kmer_size, int min_hopo_size)
{
  clock_t time0, time1;
  int i;
  hopo_set hs = (hopo_set) biomcmc_malloc (sizeof (struct hopo_set_struct));
  hs->n_hc = n_filenames;

  if (paired_end) hs->n_hc /= 2;
  if (kmer_size < 2)  kmer_size = 2;
  if (kmer_size > 15) kmer_size = 15; 
  if (min_hopo_size < 1)  min_hopo_size = 1;
  if (min_hopo_size > 32) min_hopo_size = 32; 
  
  time0 = clock ();
  hs->hc = (hopo_counter*) biomcmc_malloc (hs->n_hc * sizeof (hopo_counter));
  if (paired_end) for (i = 0; i < hs->n_hc; i++) {
    hs->hc[i] = new_or_append_hopo_counter_from_file (NULL,      filenames[2*i],   kmer_size, min_hopo_size);
    hs->hc[i] = new_or_append_hopo_counter_from_file (hs->hc[i], filenames[2*i+1], kmer_size, min_hopo_size);
  }
  else for (i = 0; i < hs->n_hc; i++) {
    hs->hc[i] = new_or_append_hopo_counter_from_file (NULL, filenames[i], kmer_size, min_hopo_size);
  }
  time1 = clock (); hs->secs_read = (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 
  for (i = 0; i < hs->n_hc; i++) finalise_hopo_counter (hs->hc[i]); 
  time1 = clock (); hs->secs_finalise = (double)(time1-time0)/(double)(CLOCKS_PER_SEC);

  hs->secs_comparison = 0.;
  hs->generator = NULL;
  hs->ref_counter = 1;

  return hs;
}

void
del_hopo_set (hopo_set hs)
{
  int i;
  if (!hs) return;
  if (--hs->ref_counter) return;
  for (i = hs->n_hc - 1; i >= 0; i--) del_hopo_counter (hs->hc[i]);
  if (hs->hc) free (hs->hc);
  del_distance_generator (hs->generator);
  free (hs);
}

distance_generator
new_distance_generator_from_hopo_set (hopo_set hs)
{
  distance_generator dg = new_distance_generator (hs->n_hc, 2); // 2 distances available (2 rescalings)
  distance_generator_set_function_data (dg, &hopo_set_distance_wrapper, (void*) hs); // hs will have any extra data needed to calculate distances
  hs->generator = dg; dg->ref_counter++;
  return dg;
}

void
hopo_set_distance_wrapper (void *data, int s1, int s2, double *result)
{
  clock_t time0, time1;
  time0 = clock ();
  compare_hopo_counters ( ((hopo_set)data)->hc[s1], ((hopo_set)data)->hc[s2], result);
  time1 = clock (); ((hopo_set)data)->secs_comparison += (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  return;
}

