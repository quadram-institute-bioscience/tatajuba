/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file 
 *  \brief wrapper around BWA library developed by Heng Li
 *  https://github.com/lh3/bwa commit 13b5637 (dowloaded 2020.08.05)
 */

#include <unistd.h>  // access() returns zero in sucess, -1 if failure;  R_OK | W_OK | X_OK or exists: F_OK
#include "bwtindex.c"
#include "bwase.c"
#include "bwtaln.c"
#include "wrapper_bwa.h"

char*
save_bwa_index (const char *genome_filename, const char *suffix, char overwrite)
{
	int algo_type = BWTALGO_AUTO, block_size = 10000000;
  char *prefix;
  if (suffix) {
    prefix = (char*) biomcmc_malloc ((strlen(genome_filename) + strlen(suffix) + 2) * sizeof (char));
    strcpy (prefix, genome_filename); strcat (prefix, "."); strcat (prefix, suffix);
  }
  else {
    prefix = (char*) biomcmc_malloc ((strlen(genome_filename) + 1) * sizeof (char));
    strcpy (prefix, genome_filename);
  }
  if (overwrite == 0) { // verify if all index files exist
    char *idxname, *term[] = {".amb",".ann",".bwt",".pac",".sa"};
    char idx_file_exists = 1;
    int i;
    idxname = (char*) biomcmc_malloc ((strlen(prefix) + 5) * sizeof (char));
    for (i = 0; (i < 5) && idx_file_exists; i++) {
      strcpy(idxname, prefix); strcat(idxname, term[i]);
      if (access(idxname, F_OK) == -1) idx_file_exists = 0;
    }
    if (idxname) free (idxname);
    if (idx_file_exists) return prefix; // all files are present
  }
  bwa_idx_build (genome_filename, prefix, algo_type, block_size);
  return prefix;
}   

int
bwa_aln_bwase (const char *index_filename, char **seqname, char **dnaseq, char **qual, size_t *seq_len, int n_dnaseq, int n_occurrences, int **match_list, char sam_to_stdout)
{
  bwa_seq_t *seqs;
  int n_matches = 0;
  gap_opt_t *opt = gap_init_opt();
  char *prefix = save_bwa_index (index_filename, NULL, false);
  if (n_occurrences < 1) n_occurrences = 1;
  if (sam_to_stdout != 0) n_matches = -1; // bwa_sai2sam() will then print SAM to stdout instead of creating match_list[]

  // TODO: bwa works with chunks of 0x40000 (262k) query seqs 
  seqs = bwa_read_seq_from_vector (seqname, dnaseq, qual, seq_len, n_dnaseq, opt->trim_qual);
  seqs = bwa_aln_from_vector (prefix, seqs, n_dnaseq, opt);
  seqs = bwa_sai2sam_se_from_vector (prefix, seqs, n_dnaseq, n_occurrences, opt, match_list, &n_matches);
  bwa_free_read_seq (n_dnaseq, seqs);
  if (opt) free (opt);
  if (prefix) free (prefix);
  return n_matches; 
}

bwase_match_t
new_bwase_match_t (const char *index_filename)
{
  bwase_match_t match = (bwase_match_t) biomcmc_malloc (sizeof (struct bwase_match_struct));
  match->prefix = save_bwa_index (index_filename, NULL, false); // generates index files
  match->bns = bns_restore (match->prefix);
  match->m = NULL;
  match->n_m = 0;
  match->ref_counter = 1;
  return match;
}

void
del_bwase_match_t (bwase_match_t match)
{
  if (!match) return;
  if (--match->ref_counter) return;
  for (int i = 0; i < match->n_m; i++) if (match->m[i].cigar)  free (match->m[i].cigar); 
  if (match->m) free (match->m);
  if (match->bns) bns_destroy (match->bns);
  free (match);
  return;
}

bwase_match_t
new_bwase_match_from_bwa_and_char_vector (const char *index_filename, char_vector seqname, char_vector dnaseq, int n_occurrences)
//char *index_filename, char **seqname, char **dnaseq, char **qual, size_t *seq_len, int n_dnaseq, int n_occurrences
{
  int i, c, n_reads, chunks = 0x40000;
  bwa_seq_t *seqs;
  gap_opt_t *opt = gap_init_opt();
  if (seqname->nstrings != dnaseq->nstrings) biomcmc_error ("Vector length of read names do not match the one with reads in bwase_match()");
  bwase_match_t match = new_bwase_match_t (index_filename);

  for (c = 0; c < seqname->nstrings; c += chunks) { 
    n_reads = chunks;
    if ((c + n_reads) > seqname->nstrings) n_reads = seqname->nstrings - c; // last loop
    seqs = bwa_read_seq_from_vector (seqname->string + c, dnaseq->string + c, NULL, dnaseq->nchars + c, n_reads, opt->trim_qual);
    seqs = bwa_aln_from_vector (match->prefix, seqs, n_reads, opt);
    seqs = bwase_to_match_t (match, seqs, n_reads, n_occurrences, opt); // updates bwase_match_t
    bwa_free_read_seq (n_reads, seqs);
  }
  if (opt) free (opt);
  return match; 
}

