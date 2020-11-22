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
int distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result);
int compare_context_histogram_for_qsort (const void *a, const void *b);
bool context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt);

#endif
