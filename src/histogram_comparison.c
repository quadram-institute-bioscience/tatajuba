/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#ifndef _histogram_comparison_headers_
#define _histogram_comparison_headers_

#include "context_histogram.h"

int compare_hopo_element_decreasing (const void *a, const void *b);
int compare_hopo_context (hopo_element a, hopo_element b);
int distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist);
int distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int *idx_match);
int compare_context_histogram_for_qsort (const void *a, const void *b);

#endif // as with regular header files, functions must be declared once  

int
compare_hopo_element_decreasing (const void *a, const void *b)
{
  int result = ((hopo_element *)b)->base - ((hopo_element *)a)->base;
  if (result) return result; // sort first by homopolymer (A/T or G/C)
  if (((hopo_element *)b)->context[0] > ((hopo_element *)a)->context[0]) return 1;// sort by left kmers
  if (((hopo_element *)b)->context[0] < ((hopo_element *)a)->context[0]) return -1; // unsigned is never negative!
  if (((hopo_element *)b)->context[1] > ((hopo_element *)a)->context[1]) return 1;// sort by right kmers
  if (((hopo_element *)b)->context[1] < ((hopo_element *)a)->context[1]) return -1;
  return ((hopo_element *)b)->length - ((hopo_element *)a)->length; // same context and homopolymer base, thus sort by tract length
}

int
compare_hopo_context (hopo_element a, hopo_element b)
{
  int result = b.base - a.base;
  if (result) return result;
  if (b.context[0] > a.context[0]) return 1;
  if (b.context[0] < a.context[0]) return -1;
  if (b.context[1] > a.context[1]) return 1;
  if (b.context[1] < a.context[1]) return -1;
  return 0; 
}

/*! \brief max_dist must be positive, and is the max allowed distance _per_ flanking region */
int
distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist)
{
  uint64_t d = *c1 ^ *c2; // XOR is one if bits are different, zero o.w. 
  int dist = 0;
  while (d && (dist < max_dist)) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference  
  return dist;
}

int
distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int *idx_match)
{ 
  int distance = 0, this_max = 0, i;
  *idx_match = -1;
  if (ch->base != he.base) return 2 * max_distance + 1; // homopolymer tracts are different
  for (i = 0; i < ch->n_context; i++) {
    distance = distance_between_context_kmer (&(ch->context[2*i]), &(he.context[0]), 2 * max_distance);
    if (distance >= 2 * max_distance) return distance;
    distance += distance_between_context_kmer (&(ch->context[2*i + 1]), &(he.context[1]), 2 * max_distance - distance);
    if (distance >= 2 * max_distance) return distance;
    if (distance > this_max) this_max = distance;
    if (distance == 0) {
      *idx_match = i;
      return 0;
    }
  }
  return this_max;
}

// TODO: distance between contexts (for each context in c1, smallest dist to all contexts in c2)
int
distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result)
{ 
  int i, min_n;
  double x,y;
  min_n = BIOMCMC_MIN (c1->h->n, c2->h->n);
  result[0] = 0.; // chi-square based distance, but 'bin-bin' with modes only 
  for (i = 0; i < min_n; i++) {
    x = (double)(1 + abs (c1->h->i[i].idx - c2->h->i[i].idx)); // absolute distance between nth-mode plus one
    y = (double)(c1->h->i[i].freq)/(double)(c1->integral) - (double)(c2->h->i[i].freq)/(double)(c2->integral); // absolute distance between nth-mode plus one
    result[0] += x * fabs (y);
  }
  result[1] = 0.; x = y = 0.; // weighted average
  for (i = 0; i < c1->h->n; i++) x += (double)(c1->h->i[i].freq * c1->h->i[i].idx)/(double)(c1->integral);
  for (i = 0; i < c2->h->n; i++) y += (double)(c2->h->i[i].freq * c2->h->i[i].idx)/(double)(c2->integral);
  x = x - y;
  result[1] = x * x; 

  return 2; // number of results
}

int
compare_context_histogram_for_qsort (const void *a, const void *b) // increasing, resolve ties with integral (e.g. location==-1)
{
  int result = (*(context_histogram_t*)a)->location - (*(context_histogram_t*)b)->location; // increasing
  if (result) return result;
  result = (*(context_histogram_t*)b)->integral - (*(context_histogram_t*)a)->integral;  // decreasing
  if (result) return result;
  /* solve ties by lexicographic order of contexts (like hopo_element) */
  if ((*(context_histogram_t*)b)->context[0] > (*(context_histogram_t*)a)->context[0]) return 1; 
  if ((*(context_histogram_t*)b)->context[0] < (*(context_histogram_t*)a)->context[0]) return -1; 
  if ((*(context_histogram_t*)b)->context[1] > (*(context_histogram_t*)a)->context[1]) return 1; 
  if ((*(context_histogram_t*)b)->context[1] < (*(context_histogram_t*)a)->context[1]) return -1; 
  return 0;
}

bool 
context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt)
{
  int tract_length, n_bases_apart;

  n_bases_apart = abs (c1->location - c2->location);
  if (!n_bases_apart) {
    if (distance) *distance = 0;
    return true; // same location
  }

  tract_length = BIOMCMC_MIN (c1->mode_context_length, c2->mode_context_length);
//  printf ("DBG::LEN::%6d :: %6d %6d\n", c1->location, n_bases_apart, tract_length);
  if (n_bases_apart > (tract_length - 1) || (c1->base != c2->base)) { // bases SHOULD be the same, but you never know...
    if (distance) *distance = -1; // should NOT be used
    return false; // locations too different
  }

  uint32_t i1, i2;
  i1 = strlen (c1->name); // name = "ATTGC-A-CCCAG"
  i2 = biomcmc_levenshtein_distance (c1->name, i1, c2->name, i1, 1, 1, true); // allows for indels
  if (distance) *distance = (int) i2;
//  printf ("DBG::DIST::%6d :: %8u\t", c1->location, i2);
  if (i2 <= (uint32_t)(opt.levenshtein_distance)) return true;
  return false;
}
