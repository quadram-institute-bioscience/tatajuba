/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief wrapper around BWA library developed by Heng Li
 *  https://github.com/lh3/bwa commit 13b5637 (dowloaded 2020.08.05)
 */

#ifndef _wrapper_bwa_h_
#define _wrapper_bwa_h_

#include "bwa.h"
#include "bntseq.h"
#include "bwtaln.h"

typedef struct bwase_match_struct* bwase_match_t;
typedef struct
{
  int query_id, ref_id, position;
  uint32_t gapo:8, gape:8, mm:8, neg_strand:1, is_best_hit:1; 
  uint64_t top1_hits:28, top2_hits:28, mapQ:8; // c1 c2 seQ from bwtaln
  char *cigar;
} bwase_elem_t;

struct bwase_match_struct
{
  bwase_elem_t *m;
  bntseq_t *bns; /*! \brief index files with ref genome info */
  char *prefix;
  int n_m, ref_counter;
};

typedef struct
{
  float fnr;                 //-n NUM  missing prob under 0.02 err rate (bwa-aln allows integer for max_diff) 
  int64_t max_gapo:12,       //-o maximum number or fraction of gap opens [%d]\n", opt->max_gapo
          max_gape: 12,      //-e maximum number of gap extensions, -1 for disabling long gaps [-1]
          s_mm:12,           //-M mismatch penalty [%d]\n", opt->s_mm
          s_gapo:12,         //-O gap open penalty [%d]\n", opt->s_gapo
          s_gape:12,         //-E gap extension penalty [%d]\n", opt->s_gape
          max_seed_diff:4;   //-k maximum differences in the seed [%d]\n", opt->max_seed_diff
  int64_t max_top2:16,       //-R stop searching when there are >INT equally best hits [%d]\n", opt->max_top2
          max_del_occ:16,    //-d maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ
          max_entries:32;    //-m maximum entries in the queue [%d]\n", opt->max_entries
  int64_t trim_qual:15,      //-q quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual
          seed_len:15,       //-l seed length [%d]\n", opt->seed_len
          indel_end_skip:15, //-i do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip
          logscaled:2,       //-L log-scaled gap penalty for long deletions
          all_hits:2,        //-N non-iterative mode: search for all n-difference hits (slooow)
          n_threads:15;
} bwase_options_t;

char *save_bwa_index (const char *genome_filename, const char *suffix, char overwrite);
/*! \brief "index + aln + bwase"; creates indices if absent. do not forget to match_list=NULL the first time (if you are not appending) 
 *  \param print_to_stdout zero if you want to fill match_list[] (most common use), otherwise just prints SAM format to stdout
 *  \result number n of successful matches; most important is however one-dimensional match_list[] with five columns (consecutive values) per match.
 * */
int bwa_aln_bwase (const char *index_filename, char **seqname, char **dnaseq, char **qual, size_t *seq_len, int n_dnaseq, int n_occurrences, int **match_list, char print_to_stdout);

/*! \brief new param options with default values, mandatory for bwase_to_match() */ 
bwase_options_t new_bwase_options_t (int level); 
/*! \brief new version, storing results to struct (from bwase.c) [this is low level component, assumes a bwase-match_t exists] */
bwa_seq_t *bwase_to_match_t (bwase_match_t match, bwa_seq_t *seqs, int n_dnaseq, int n_occ, gap_opt_t *opt);
char *bwase_match_ref_genome_name (bwase_match_t match, int i);

void del_bwase_match_t (bwase_match_t match);
/*! \brief high-level function that creates the bwase_match_t struct and finds matches to index files (defined by FASTA name */
bwase_match_t new_bwase_match_from_bwa_and_char_vector (const char *index_filename, char_vector seqname, char_vector dnaseq, int n_occurrences, bwase_options_t bopt);

#include "bwamem.h"
typedef struct bwmem_match_struct* bwmem_match_t;

typedef struct
{
  int query_id, ref_id, score;
	int64_t position, rb, re; // [rb,re): reference sequence in the alignment
	int qb, qe;     // [qb,qe): query sequence in the alignment
  uint32_t mapQ:8, edit_distance:22, neg_strand:1, is_primary:1; 
  char *cigar;
} bwmem_elem_t;

struct bwmem_match_struct
{
  bwmem_elem_t *m;
  bwaidx_t *idx; /*! \brief index files with ref genome info */ 
  char *prefix;
  int n_m, ref_counter;
};

bwmem_match_t new_bwmem_match_t (const char *index_filename);
bwmem_match_t new_bwmem_match_from_bwa_and_char_vector (const char *index_filename, char_vector dnaseq);
void del_bwmem_match_t (bwmem_match_t match);
void update_bwmem_match (bwmem_match_t match, int query_id, char *seq, size_t len, mem_opt_t *opt);
char * bwmem_match_ref_genome_name (bwmem_match_t match, int i);
#endif
