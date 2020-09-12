/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#ifndef ALIGNMENT_PARAMETERS_H
#define ALIGNMENT_PARAMETERS_H

#include <biomcmc.h>
#include <float.h>
#include "global.h"
#include "tldevel.h"
#include "rng.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


struct msa_seq {
  char* seq;   // original (string) representation of sequence
  uint32_t id; // leo: we don't store seq names here, external char_vector should have it
  uint8_t* s;  // compact representation of sequence
  uint32_t* gaps;
  uint32_t len;
};

struct msa {
  struct msa_seq** sequences;
  int** sip;
  int* nsip;
  uint32_t* plen;
  uint32_t numseq;
  uint32_t num_profiles;
  bool is_protein; 
};

struct aln_param {
  struct rng_state* rng;
  float** subm;
  float gpo;
  float gpe;
  float tgpe;
  uint32_t* tree;
};

struct alphabet {
  int8_t to_internal[128];
  int8_t to_external[32];
  int type;
};

/**< main input function **/
struct msa* read_char_vector_to_msa (char_vector dna, bool is_protein);
/**< main output function **/
char_vector aligned_msa_to_charvector (struct msa* msa);

struct aln_param* init_ap (int numseq, bool is_protein);
void free_ap (struct aln_param* ap);
uint32_t* pick_anchor (struct msa* msa, uint32_t* n);
int make_aliged_seq (uint8_t* aligned, uint8_t* unaligned, int* gaps,int len);

void convert_msa_to_internal(struct msa* msa, bool is_protein);/* convert */
int write_msa(struct msa* msa, char* outfile, int type);
void free_msa(struct msa* msa);

void weave (struct msa* msa, int** map, uint32_t* tree);
int clean_aln(struct msa* msa);
// euclidean_distance.h
extern  int edist_256(const float* a,const float* b, const int len, float* ret);
extern int edist_serial(const float* a,const float* b,const int len, float* ret);
extern int edist_serial_d(const double* a,const double* b,const int len, double* ret);

extern int shuffle_arr_r(int* arr,int n, struct rng_state* rng);
/* Steinegger, Martin, and Johannes Söding. "Clustering huge protein sequence sets in linear time." Nature communications 9.1 (2018): 2542. */
// (c) 2017 Johannes Soeding & Martin Steinegger, Gnu Public License version 3
uint16_t circ_hash(const uint8_t* x, const uint8_t length);
uint16_t circ_hash_next(const uint8_t * x,const uint8_t length,const uint8_t x_first, uint16_t h);
// kmeans.h
double** kmeans(double** data,int* cluster_assignment, int len_a,int len_b, int k);

#endif
