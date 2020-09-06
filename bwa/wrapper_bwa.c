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
