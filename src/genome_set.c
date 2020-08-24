/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "genome_set.h"
  
g_tract_vector_t new_g_tract_vector_from_genomic_context_list (genomic_context_list_t *genome, int n_genome);
void del_g_tract_vector (g_tract_vector_t tract);
void g_tract_vector_concatenate_tracts  (g_tract_vector_t tract, genomic_context_list_t *genome, int n_genome);
void generate_g_tract_from_context_histogram (g_tract_vector_t tract, int prev, int curr);

genome_set_t
new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt) 
{
  clock_t time0, time1;
  int i, j;
  hopo_counter hc;
  genome_set_t g = (genome_set_t) biomcmc_malloc (sizeof (struct genome_set_struct));
  g->n_genome = n_filenames;
  g->ref_counter = 1;
  g->secs_read = g->secs_finalise = g->secs_comparison = 0.;

  if (opt.paired_end) g->n_genome /= 2;
  
  time0 = clock ();
  g->genome = (genomic_context_list_t*) biomcmc_malloc (g->n_genome * sizeof (genomic_context_list_t));
  if (opt.paired_end) for (i = 0; i < g->n_genome; i++) {
    hc = new_or_append_hopo_counter_from_file (NULL, filenames[2*i],   opt);
    hc = new_or_append_hopo_counter_from_file (hc,   filenames[2*i+1], opt);
    time1 = clock (); g->secs_read += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

    g->genome[i] = new_genomic_context_list (hc);
    del_hopo_counter (hc); hc = NULL;
    time1 = clock (); g->secs_finalise += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 
  }
  else for (i = 0; i < g->n_genome; i++) {
    hc = new_or_append_hopo_counter_from_file (NULL, filenames[i], opt);
    time1 = clock (); g->secs_read += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

    g->genome[i] = new_genomic_context_list (hc);
    del_hopo_counter (hc); hc = NULL;
    time1 = clock (); g->secs_finalise += (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 
  }

  /* label each histogram with index of genome they belong to */
  for (i = 0; i < g->n_genome; i++) for (j = 0; j < g->genome[i]->n_hist; j++) g->genome[i]->hist[j]->index = i;
  /* generate comparisons between genomes */
  g->tract = new_g_tract_vector_from_genomic_context_list (g->genome, g->n_genome);
  time1 = clock (); g->secs_comparison = (double)(time1-time0)/(double)(CLOCKS_PER_SEC); time0 = time1; 

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
  del_g_tract_vector (g->tract);
  free (g);
}

g_tract_vector_t
new_g_tract_vector_from_genomic_context_list (genomic_context_list_t *genome, int n_genome)
{
  int i, prev;
  g_tract_vector_t tract = (g_tract_vector_t) biomcmc_malloc (sizeof (struct g_tract_vector_struct));
  tract->distinct = NULL;
  tract->concat = NULL;
  tract->n_distinct = tract->n_concat = 0;
  g_tract_vector_concatenate_tracts  (tract, genome, n_genome); // updates tract->concat

  for (prev=0, i = 1; i < tract->n_concat; i++) if (tract->concat[i]->location != tract->concat[prev]->location) {
    generate_g_tract_from_context_histogram (tract, prev, i); // updates tract->distinct
    prev = i;
  }
  generate_g_tract_from_context_histogram (tract, prev, i-1); //  last block i == tract->n_concat

  return tract;
}

void
del_g_tract_vector (g_tract_vector_t tract)
{
  if (!tract) return;
  if (tract->concat) free (tract->concat);
  if (tract->distinct) {
    for (int i = tract->n_distinct - 1; i >= 0; i--) {
      if (tract->distinct[i].dist) free (tract->distinct[i].dist);
      del_empfreq (tract->distinct[i].mode);
    }
    free (tract->distinct);
  }
  free (tract);
}

void 
g_tract_vector_concatenate_tracts (g_tract_vector_t tract, genomic_context_list_t *genome, int n_genome)
{
  int i, j, k1, k2, dif, n_new_h;
  context_histogram_t *new_h = NULL;
  /** genomic_context_list is always sorted by position (with position=-1 at beginning, before ref_start) */

  /* 1. skip genomes without (genome-location) annotated histograms and start new histogram */ 
  for (i = 0; genome[i]->ref_start == genome[i]->n_hist; i++); // do nothing besides loop
  tract->n_concat = genome[i]->n_hist - genome[i]->ref_start;
  tract->concat   = (context_histogram_t*) biomcmc_malloc (sizeof (context_histogram_t) * tract->n_concat);
  for (j=0, k1 = genome[i]->ref_start; k1 < genome[i]->n_hist; j++, k1++) tract->concat[j] = genome[i]->hist[k1];
  new_h = tract->concat;
  n_new_h = tract->n_concat;

  /* 2. concat[] will be new_h and genome[i]->hist merged, in order (both are ordered) */
  for (i++; i < n_genome; i++) if (genome[i]->n_hist > genome[i]->ref_start) {
    tract->n_concat = n_new_h + genome[i]->n_hist - genome[i]->ref_start; // sum of lengths
    printf ("DBG::adding genome::%4d::%5d\n",i, genome[i]->n_hist - genome[i]->ref_start);
    tract->concat = (context_histogram_t*) biomcmc_malloc (tract->n_concat * sizeof (context_histogram_t));
    j = 0;
    for (k1 = 0, k2 = genome[i]->ref_start; (k1 < n_new_h) && (k2 < genome[i]->n_hist); ) { 
      dif = new_h[k1]->location - genome[i]->hist[k2]->location;
      if (dif > 0) tract->concat[j++] = genome[i]->hist[k2++]; // k2 steps forward 
      if (dif < 0) tract->concat[j++] = new_h[k1++];           // k1 steps forward
      else { // both are the same, both step forward and tract increases by two
        tract->concat[j++] = genome[i]->hist[k2++]; 
        tract->concat[j++] = new_h[k1++];
      }
    } // loop ends when one of them finishes, now we need to complete with the other
    for (; k2 < genome[i]->n_hist; k2++, j++) tract->concat[j] = genome[i]->hist[k2]; 
    for (; k1 < n_new_h; k1++, j++) tract->concat[j] = new_h[k1];
    free (new_h); // only place where we free is here (o.w. we only have pointers)
    new_h = tract->concat;
    n_new_h = tract->n_concat;
  }
  printf ("DBG::total::%5d\n", tract->n_concat); 
}

void
generate_g_tract_from_context_histogram (g_tract_vector_t tract, int prev, int curr)
{
  int *modlen=NULL, *index=NULL;
  int i, j, k, this = tract->n_distinct;
  tract->n_distinct++;
  tract->distinct = (g_tract_t*) biomcmc_realloc ((g_tract_t*) tract->distinct, tract->n_distinct * sizeof (g_tract_t));
  tract->distinct[this].location = tract->concat[prev]->location;
  tract->distinct[this].dist = NULL;

  /* 1. empirical frequency of modal lengths */
  modlen = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  index  = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  for (i = 0; i < (curr - prev); i++) {
    modlen[i]  = tract->concat[i]->h->i[0].idx; // modal tract length for this genome id
    index[i] = tract->concat[i]->index; // genome id (index)
  }
  tract->distinct[this].mode = new_empfreq_from_int_weighted (index, (curr - prev), modlen);
  if (modlen) free (modlen);
  if (index) free (index);

  /* 2. matrix of pairwise distances (1D vector sorted from highest to smallest distance) */
  tract->distinct[this].n_dist = ((curr-prev) * (curr-prev-1))/2;
  if (!tract->distinct[this].n_dist) return;
  tract->distinct[this].dist = (double*) biomcmc_malloc (tract->distinct[this].n_dist * sizeof (double));
  for (k = i = 0; i < (curr - prev); i++) for (j = 0; j < i; j++) distance_between_context_histograms (tract->concat[i], tract->concat[j], tract->distinct[this].dist + (k++));
  qsort (tract->distinct[this].dist, tract->distinct[this].n_dist, sizeof (double), compare_double_decreasing);
  return;
}

void
print_debug_g_tract_vector (g_tract_vector_t tract)
{
  int i;
  empfreq e;
  printf ("location n_genomes min_length max_length max_dist min_dist\n"); 
  for (i = 0; i < tract->n_distinct; i++) { 
    e = tract->distinct[i].mode;
    printf ("%-5d %-5d %-5d %-5d ", tract->distinct[i].location, e->n, e->i[0].freq, e->i[e->n-1].freq);
    if (tract->distinct[i].n_dist>0) printf ("%9e %9e\n", tract->distinct[i].dist[0], tract->distinct[i].dist[tract->distinct[i].n_dist-1]);
    else printf ("0.0 0.0\n");
  }
}
