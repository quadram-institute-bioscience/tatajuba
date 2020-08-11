/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file 
 *  \brief wrapper around BWA library developed by Heng Li
 *  https://github.com/lh3/bwa commit 13b5637 (dowloaded 2020.08.05)
 */

#include <biomcmc.h>
#include "bwtindex.c"

void
save_bwa_index (const char *genome_filename, const char *suffix)
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
  bwa_idx_build(genome_filename, prefix, algo_type, block_size);
  if (prefix) free (prefix);
  return;
}   
  
