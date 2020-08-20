/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "hopo_set.h"
  
void hopo_set_distance_wrapper (void *data, int s1, int s2, double *result);

genome_set_t
new_genome_set_from_files (const char **filenames, int n_filenames, bool paired_end, int kmer_size, int min_tract_size, const char *reference_genome_filename)
{
  clock_t time0, time1;
  int i;
  hopo_counter hc;
  genome_set_t g = (genome_set_t) biomcmc_malloc (sizeof (struct genome_set_struct));
  g->n_genomes = n_filenames;

  if (paired_end) g->n_genome /= 2;
  if (kmer_size < 2)  kmer_size = 2;
  if (kmer_size > 32) kmer_size = 32; 
  if (min_tract_size < 2)  min_tract_size = 2;
  if (min_tract_size > 32) min_tract_size = 32; 
  
  time0 = clock ();
  g->genome = (genomic_context_list_t*) biomcmc_malloc (g>n_genome * sizeof (genomic_context_list));
  if (paired_end) for (i = 0; i < g->n_genome; i++) {
    hc = new_or_append_hopo_counter_from_file (NULL, filenames[2*i],   kmer_size, min_tract_size);
    hc = new_or_append_hopo_counter_from_file (hc,   filenames[2*i+1], kmer_size, min_tract_size);
    time1 = clock (); g->secs_read = (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

    g->genome[i] = new_genomic_context_list (hc, reference_genome_filename, max_distance_per_flank, min_coverage);
    del_hopo_counter (hc); hc = NULL;
    time1 = clock (); g->secs_finalise = (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  }
  else for (i = 0; i < g->n_genome; i++) {
    hc = new_or_append_hopo_counter_from_file (NULL, filenames[i], kmer_size, min_tract_size);
    time1 = clock (); g->secs_read = (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

    g->genome[i] = new_genomic_context_list (hc, reference_genome_filename, max_distance_per_flank, min_coverage);
    del_hopo_counter (hc); hc = NULL;
    time1 = clock (); g->secs_finalise = (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  }

  g->secs_comparison = 0.;
  g->generator = NULL;
  g->ref_counter = 1;

  return g;
}

void
del_genome_set (genome_set_t g)
{
  int i;
  if (!g) return;
  if (--g->ref_counter) return;
  for (i = g->n_genome - 1; i >= 0; i--) del_genomic_context_list (g->genome[i]);
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

