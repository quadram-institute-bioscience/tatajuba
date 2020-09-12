/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

/*! \file
 *  \brief 
 *
 *  Algorithms in this file based on [kalign3](https://github.com/TimoLassmann/kalign.git)(GPL-3.0-or-later  &copy;2006, 2019 Timo Lassmann
 */

#include <biomcmc.h>
#include "global.h"
#include "alignment_parameters.h"
#include "bisectingKmeans.h"
#include "alignment.h"

bool charvector_is_protein (char_vector seq);

char_vector kalign3_from_char_vector (char_vector seq)
{
  struct msa* msa = NULL;
  struct aln_param* ap = NULL;
  int** map = NULL; /* holds all alignment paths  */
  uint32_t i;
  bool is_protein;
  char_vector aligned_charvector = NULL;

  is_protein = charvector_is_protein (seq); // once protein alignment is finished
  if (is_protein) biomcmc_warning ("protein detected\n");
  else  biomcmc_warning ("DNA detected\n");
  /* Step 1: read all input sequences & figure out output  */
  msa = read_char_vector_to_msa (seq, is_protein);
  ap = init_ap (msa->numseq, is_protein);/* allocate aln parameters  */
  /* Start bi-secting K-means sequence clustering */
  build_tree_kmeans (msa, ap);
  if (is_protein) convert_msa_to_internal (msa, 2); // is_protein = 2 means to make standard (full) protein 
  /* Start alignment stuff */ 
  map = hirschberg_alignment(msa, ap);
  weave (msa, map, ap->tree); /* it's aligned already */
  /* clean up map */
  for (i = 0; i < msa->num_profiles ;i++) if(map[i]) free (map[i]); 
  if (map) free (map);
  /* We are done. */
  aligned_charvector = aligned_msa_to_charvector (msa);
  free_msa(msa);
  free_ap(ap);
  return aligned_charvector;
}

bool
charvector_is_protein (char_vector seq)
{
  char  dna_letters[] = "acgtuACGTUnN"; // MRWSYK are ambiguous
  char prot_letters[] = "defhiklmpqrsvwyDEFHIKLMPQRSVWY"; // excluding acgtn
  size_t i, len_dna, len_prot;
  uint8_t dna[256], prot[256];
  int j, count_dna = 0, count_prot = 0, evidence_dna = 0;

  len_dna  = strlen (dna_letters);
  len_prot = strlen (prot_letters);
  for (i = 0; i < 256; i++) dna[i] = prot[i] = 0;
  for (i = 0; i < len_dna;  i++)  dna[ (int)dna_letters[i] ] = 1;
  for (i = 0; i < len_prot; i++) prot[ (int)prot_letters[i] ] = 1;

  for (j = 0; j < seq->nstrings; j++) {
    for (i = 0; i < seq->nchars[i]; i++) {
      count_dna  +=  dna[ (int)seq->string[j][i] ];
      count_prot += prot[ (int)seq->string[j][i] ]; // indels etc are not counted by either
    }
    if ((count_dna > 0xfffffff) || (count_prot < 0xfffffff)) { // zeroes counter to avoid overflow 
      evidence_dna += (count_dna > count_prot? 1:-1);
      count_dna = count_prot = 0;
    }
  } // for strings
  evidence_dna += (count_dna > count_prot? 1:-1);
  if (evidence_dna > 0) return false;
  else return true; // if in doubt, treat as protein
}
