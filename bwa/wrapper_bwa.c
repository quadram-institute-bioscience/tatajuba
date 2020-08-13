/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file 
 *  \brief wrapper around BWA library developed by Heng Li
 *  https://github.com/lh3/bwa commit 13b5637 (dowloaded 2020.08.05)
 */

#include <biomcmc.h>
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

void
bwa_aln_bwase (const char *index_filename, char **seqname, char **dnaseq, char **qual, size_t *seq_len, int n_dnaseq, int n_occurrences)
{
  bwa_seq_t *seqs;
  gap_opt_t *opt = gap_init_opt();
  char *prefix = save_bwa_index (index_filename, NULL, true);
  if (n_occurrences < 1) n_occurrences = 1;

  seqs = bwa_read_seq_from_vector (seqname, dnaseq, qual, seq_len, n_dnaseq, opt->trim_qual);
  seqs = bwa_aln_from_vector (prefix, seqs, n_dnaseq, opt);
  seqs = bwa_sai2sam_se_from_vector (prefix, seqs, n_dnaseq, n_occurrences, opt);
  bwa_free_read_seq (n_dnaseq, seqs);
  if (opt) free (opt);
  if (prefix) free (prefix);
  return; // FIXME: should return hits etc
}
