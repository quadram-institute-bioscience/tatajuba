/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "genome_set.h"
  
g_tract_vector_t new_g_tract_vector_from_genomic_context_list (genomic_context_list_t *genome, int n_genome);
void del_g_tract_vector (g_tract_vector_t tract);
void g_tract_vector_concatenate_tracts  (g_tract_vector_t tract, genomic_context_list_t *genome, int n_genome);
void generate_g_tract_from_context_histogram (g_tract_vector_t tract, int prev, int curr, int lev_distance);

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
  int i, prev, lev_distance = 0, max_lev_distance = 0;
  bool is_same;
  g_tract_vector_t tract = (g_tract_vector_t) biomcmc_malloc (sizeof (struct g_tract_vector_struct));
  tract->distinct = NULL;
  tract->concat = NULL; // big vector with all tracts from all genomes in linear order
  tract->n_distinct = tract->n_concat = 0;
  g_tract_vector_concatenate_tracts  (tract, genome, n_genome); // updates tract->concat

  for (prev=0, i = 1; i < tract->n_concat; i++) { // overlap() is true if tracts are the same, false if distinct 
    is_same = context_histograms_overlap (tract->concat[i], tract->concat[prev], &lev_distance, genome[0]->opt);
    if (!is_same) {
      generate_g_tract_from_context_histogram (tract, prev, i, max_lev_distance); // updates tract->distinct
      prev = i;
      max_lev_distance = 0;
    }
    else if (max_lev_distance < lev_distance) max_lev_distance = lev_distance; // tracts are similar enough 
  }
  // call context() just to calculate levenshtein distance 
  context_histograms_overlap (tract->concat[i-1], tract->concat[prev], &lev_distance, genome[0]->opt);
  if (max_lev_distance < lev_distance) max_lev_distance = lev_distance; // tracts are similar enough 
  generate_g_tract_from_context_histogram (tract, prev, i-1, max_lev_distance); //  last block i == tract->n_concat

  return tract;
}

void
del_g_tract_vector (g_tract_vector_t tract)
{
  if (!tract) return;
  if (tract->concat) free (tract->concat);
  if (tract->distinct) {
    for (int i = tract->n_distinct - 1; i >= 0; i--) {
      if (tract->distinct[i].d1) free (tract->distinct[i].d1);
      del_empfreq (tract->distinct[i].mode);
      del_context_histogram (tract->distinct[i].example);
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
  for (i+=1; i < n_genome; i++) if (genome[i]->n_hist > genome[i]->ref_start) {
    tract->n_concat = n_new_h + genome[i]->n_hist - genome[i]->ref_start; // sum of lengths
    tract->concat = (context_histogram_t*) biomcmc_malloc (tract->n_concat * sizeof (context_histogram_t));
    for (j = 0, k1 = 0, k2 = genome[i]->ref_start; (k1 < n_new_h) && (k2 < genome[i]->n_hist); ) { 
      dif = new_h[k1]->location - genome[i]->hist[k2]->location;
      if (dif > 0) tract->concat[j++] = genome[i]->hist[k2++]; // k2 steps forward (if equal both step forward)
      else if (dif < 0) tract->concat[j++] = new_h[k1++];           // k1 steps forward
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
}

void
generate_g_tract_from_context_histogram (g_tract_vector_t tract, int prev, int curr, int lev_distance)
{
  int *modlen=NULL, *index=NULL;
  int i, j, k, this = tract->n_distinct;
  double result[2];
  tract->n_distinct++;
  tract->distinct = (g_tract_t*) biomcmc_realloc ((g_tract_t*) tract->distinct, tract->n_distinct * sizeof (g_tract_t));
  tract->distinct[this].location = tract->concat[prev]->location;
  tract->distinct[this].d1 = NULL;
  tract->distinct[this].lev_distance = lev_distance;

  tract->distinct[this].example = tract->concat[prev]; // example of a context_histogram_t
  tract->distinct[this].example->ref_counter++;

  /* 1. empirical frequency of modal lengths */
  modlen = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  index  = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  for (i = 0, j = prev; j < curr; j++, i++) {
    modlen[i]  = tract->concat[j]->h->i[0].idx; // modal tract length for this genome id
//    for(k=0;k<tract->concat[j]->h->n;k++) { printf (" %3d (%4d)", tract->concat[j]->h->i[k].idx, tract->concat[j]->h->i[k].freq);} printf ("::%5d::DBG\n", prev);
    index[i] = tract->concat[i]->index; // genome id (index), not very useful now :)
  }
  tract->distinct[this].mode = new_empfreq_from_int_weighted (index, (curr - prev), modlen);
  if (modlen) free (modlen);
  if (index) free (index);

  /* 2. matrix of pairwise distances (1D vector sorted from highest to smallest distance) */
  tract->distinct[this].n_dist = ((curr-prev) * (curr-prev-1))/2;
  if (!tract->distinct[this].n_dist) return;
  tract->distinct[this].d1 = (double*) biomcmc_malloc (2 * tract->distinct[this].n_dist * sizeof (double)); // d1 and d2
  tract->distinct[this].d2 = tract->distinct[this].d1 + tract->distinct[this].n_dist;

  for (k = 0, i = prev+1; i < curr; i++) for (j = prev; j < i; j++) {
    distance_between_context_histograms (tract->concat[i], tract->concat[j], result);
    tract->distinct[this].d1[k]   = result[0]; 
    tract->distinct[this].d2[k++] = result[1]; 
  }
  qsort (tract->distinct[this].d1, tract->distinct[this].n_dist, sizeof (double), compare_double_decreasing);
  qsort (tract->distinct[this].d2, tract->distinct[this].n_dist, sizeof (double), compare_double_decreasing);
  return;
}

void
print_interesting_tracts (genome_set_t g)
{
  int i;
  g_tract_t *t;
  bool to_print = false;
  printf ("tract location n_genomes | max_tract_length min_tract_length | lev_distance | max_chi2_dist max_avge_dist\n"); 
  for (i = 0; i < g->tract->n_distinct; i++) { 
    to_print = false;
    t = g->tract->distinct + i;
    if (t->mode->n < g->n_genome) to_print = true;
    if ((!to_print) && (t->mode->i[0].freq != t->mode->i[t->mode->n-1].freq)) to_print = true;
    if ((!to_print) && (t->lev_distance > 0)) to_print = true;
    if ((!to_print) && (t->d1[0] > 0.)) to_print = true;
    if ((!to_print) && (t->d2[0] > 0.)) to_print = true;
    if (to_print) {
      printf ("%-60s %-5d %-5d |%-5d %-5d | %-5d | ", t->example->name, t->location, t->mode->n, t->mode->i[0].freq, t->mode->i[t->mode->n-1].freq, t->lev_distance);
      if (t->n_dist > 0) printf ("%12e %12e\n", t->d1[0], t->d2[0]);
      printf ("%12e %12e\n", 0., 0.);
    }
  }
}

void
print_debug_g_tract_vector (g_tract_vector_t tract)
{
  int i;
  empfreq e;
  printf ("tract location n_genomes | min_length max_length | lev_distance | max_chi2_dist max_avge_dist\n"); 
  for (i = 0; i < tract->n_distinct; i++) { 
    e = tract->distinct[i].mode;
    printf ("%-60s %-5d %-5d |%-5d %-5d | %-5d | ", tract->distinct[i].example->name, tract->distinct[i].location, e->n, e->i[0].freq, e->i[e->n-1].freq, tract->distinct[i].lev_distance);
    if (tract->distinct[i].n_dist>0) printf ("%12e %12e\n", tract->distinct[i].d1[0], tract->distinct[i].d2[0]);
    else printf ("0.0 0.0\n");
  }
}
