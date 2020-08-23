/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "genome_set.h"
  
//void hopo_set_distance_wrapper (void *data, int s1, int s2, double *result);
void  genome_set_concatenate_tracts (genome_set_t g);

genome_set_t
new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt) 
{
  clock_t time0, time1;
  int i;
  hopo_counter hc;
  genome_set_t g = (genome_set_t) biomcmc_malloc (sizeof (struct genome_set_struct));
  g->n_genome = n_filenames;
  g->tract = NULL;
  g->n_tract = 0;

  if (opt.paired_end) g->n_genome /= 2;
  
  time0 = clock ();
  g->genome = (genomic_context_list_t*) biomcmc_malloc (g->n_genome * sizeof (genomic_context_list_t));
  if (opt.paired_end) for (i = 0; i < g->n_genome; i++) {
    hc = new_or_append_hopo_counter_from_file (NULL, filenames[2*i],   opt);
    hc = new_or_append_hopo_counter_from_file (hc,   filenames[2*i+1], opt);
    time1 = clock (); g->secs_read = (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

    g->genome[i] = new_genomic_context_list (hc);
    del_hopo_counter (hc); hc = NULL;
    time1 = clock (); g->secs_finalise = (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  }
  else for (i = 0; i < g->n_genome; i++) {
    hc = new_or_append_hopo_counter_from_file (NULL, filenames[i], opt);
    time1 = clock (); g->secs_read = (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

    g->genome[i] = new_genomic_context_list (hc);
    del_hopo_counter (hc); hc = NULL;
    time1 = clock (); g->secs_finalise = (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  }
  g->secs_comparison = 0.;
  g->generator = NULL;
  g->ref_counter = 1;

  genome_set_concatenate_tracts (g);
  return g;
}

void
del_genome_set (genome_set_t g)
{
  int i;
  if (!g) return;
  if (--g->ref_counter) return;
  for (i = g->n_genome - 1; i >= 0; i--) del_genomic_context_list (g->genome[i]);
  if (g->genome) free (g->genome);
  del_distance_generator (g->generator);
  free (g);
}

void 
genome_set_concatenate_tracts (genome_set_t g)
{
  int i, j, k1, k2, dif, n_new_h;
  context_histogram_t *new_h = NULL;

  /* 1. label each histogram with index of genome they belong to */
  for (i = 0; i < g->n_genome; i++) for (j = 0; j < g->genome[i]->n_hist; j++) g->genome[i]->hist[j]->index = i;

  /* 2. skip genomes without (genome-location) annotated histograms and start new histogram */ 
  for (i = 0; g->genome[i]->ref_start == g->genome[i]->n_hist; i++); // do nothing besides loop
  g->n_tract = g->genome[i]->n_hist - g->genome[i]->ref_start;
  g->tract   = (context_histogram_t*) biomcmc_malloc (sizeof (context_histogram_t) * g->n_tract);
  for (j = g->genome[i]->ref_start; j < g->genome[i]->n_hist; j++) g->tract[j] = g->genome[i]->hist[j];
  new_h = g->tract;
  n_new_h = g->n_tract;

  /* 3. tract[] will be new_h and genome[i]->hist merged, in order (both are ordered) */
  for (i++; i < g->n_genome; i++) if (g->genome[i]->n_hist > g->genome[i]->ref_start) {
    printf ("DBG::%4d::%5d\n",i, g->genome[i]->n_hist - g->genome[i]->ref_start);
    g->n_tract = n_new_h + g->genome[i]->n_hist - g->genome[i]->ref_start; // sum of lengths
    g->tract = (context_histogram_t*) biomcmc_malloc (sizeof (context_histogram_t) * g->n_tract);
    j = 0;
    for (k1 = 0, k2 = g->genome[i]->ref_start; (k1 < n_new_h) && (k2 < g->genome[i]->n_hist); ) { 
      dif = new_h[k1]->location - g->genome[i]->hist[k2]->location;
      if (dif > 0) g->tract[j++] = g->genome[i]->hist[k2++]; // k2 steps forward 
      if (dif < 0) g->tract[j++] = new_h[k1++];              // k1 steps forward
      else { // both are the same, both step forward and tract increases by two
        g->tract[j++] = g->genome[i]->hist[k2++]; 
        g->tract[j++] = new_h[k1++];
      }
    } // loop ends when one of them finishes, now we need to complete with the other
    for (; k2 < g->genome[i]->n_hist; k2++, j++) g->tract[j] = g->genome[i]->hist[k2]; 
    for (; k1 < n_new_h; k1++, j++) g->tract[j] = new_h[k1];
    free (new_h); // only place where we free is here (o.w. we only have pointers)
    new_h = g->tract;
    n_new_h = g->n_tract;
  }
  printf ("DBG::total::%5d\n", g->n_tract); 
}

/*
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
*/
