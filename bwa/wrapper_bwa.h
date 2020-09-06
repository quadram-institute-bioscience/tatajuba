/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief wrapper around BWA library developed by Heng Li
 *  https://github.com/lh3/bwa commit 13b5637 (dowloaded 2020.08.05)
 */

#ifndef _wrapper_bwa_h_
#define _wrapper_bwa_h_

#include <stddef.h>
#include <biomcmc.h>

typedef struct bwase_match_struct* bwase_match_t;
typedef struct
{
  int query_id, ref_id, position, score;
  uint32_t gapo:8, gape:8, mm:8, neg_strand:1, is_best_hit:1; 
  uint64_t top1_hits:28, top2_hits:28, mapQ:8; // c1 c2 seQ from bwtaln
} bwase_elem_t;

struct bwase_match_struct
{
  bwase_elem_t *m;
  char_vector ref; /*! \brief reference genomes from index files */
  int n_m;
};

char *save_bwa_index (const char *genome_filename, const char *suffix, char overwrite);
/*! \brief "index + aln + bwase"; creates indices if absent. do not forget to match_list=NULL the first time (if you are not appending) 
 *  \param print_to_stdout zero if you want to fill match_list[] (most common use), otherwise just prints SAM format to stdout
 *  \result number n of successful matches; most important is however one-dimensional match_list[] with five columns (consecutive values) per match.
 * */
int bwa_aln_bwase (const char *index_filename, char **seqname, char **dnaseq, char **qual, size_t *seq_len, int n_dnaseq, int n_occurrences, int **match_list, char print_to_stdout);

#endif
