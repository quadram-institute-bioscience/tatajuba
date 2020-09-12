/* SPDX-License-Identifier: GPL-3.0-or-later
* Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
* 
* This file is based on the Kalign3  program, commit 
* [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
* TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
*/


#include <xmmintrin.h>
#include "sequence_distance.h"
#include "alignment_parameters.h"
#include "alignment.h"
#include "bpm.h"

#define NODESIZE 16

/* small hash implementation */
struct bignode{
  struct bignode *next;
  unsigned int pos[NODESIZE];
  unsigned int num;
};

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos);
void big_remove_nodes(struct bignode *n);
float calc_distance (uint8_t* seq_a, uint8_t* seq_b, uint32_t len_a, uint32_t len_b, bool is_protein);
float protein_wu_distance_calculation (struct bignode* hash[], const uint8_t* seq, const int seqlen, const int diagonals, const float mode);
float dna_distance_calculation(struct bignode* hash[],const uint8_t * p,const uint32_t seqlen, uint32_t diagonals, float mode);

float** d_estimation (struct msa* msa, uint32_t* samples, uint32_t num_samples,int pair)
{
  float** dm = NULL;
  uint8_t* seq_a;
  uint8_t* seq_b;
  float dist;
  uint32_t len_a, len_b, i, j;
#if HAVE_AVX2
  set_broadcast_mask();
#endif

  if(pair){
    // dm = galloc(dm, num_samples, num_samples, 0.0f); // DEBUG (-g -Wall) doesn't like _Generic()...
    dm = alloc_2D_array_size_float (dm, num_samples, num_samples, 0.0f);
    for(i = 0; i < num_samples;i++){
      seq_a = msa->sequences[samples[i]]->s;// aln->s[samples[i]];
      len_a = msa->sequences[samples[i]]->len;//aln->sl[samples[i]];
      for(j = 0;j < num_samples;j++){
        seq_b = msa->sequences[samples[j]]->s; //aln->s[ samples[j]];
        len_b = msa->sequences[samples[j]]->len;//aln->sl[selection[j]];
        dist = calc_distance (seq_a, seq_b, len_a, len_b, msa->is_protein);
        dm[i][j] = dist;//*dist;
        dm[j][i] = dm[i][j];
      }
    }
  } else {
    uint32_t a, numseq = msa->numseq;
    dm = (float**) biomcmc_malloc (sizeof(float*) * numseq);
    a = num_samples / 8;
    if (num_samples%8) a++;
    a = a << 3;

    for(i = 0; i < numseq; i++) {
      dm[i] = NULL;
      dm[i] = _mm_malloc(sizeof(float) * a,32);
      for(j = 0; j < a; j++) dm[i][j] = 0.0f;
    }

    for(i = 0; i < numseq;i++){
      seq_a = msa->sequences[i]->s;// aln->s[i];
      len_a = msa->sequences[i]->len;//  aln->sl[i];
      for(j = 0;j < num_samples;j++){
        seq_b = msa->sequences[samples[j]]->s;// aln->s[ seeds[j]];
        len_b = msa->sequences[samples[j]]->len;// aln->sl[seeds[j]];
        dist = calc_distance (seq_a, seq_b, len_a, len_b, msa->is_protein);
        dm[i][j] = dist;
      }
    }
  }
  return dm;
}

float calc_distance (uint8_t* seq_a, uint8_t* seq_b, uint32_t len_a, uint32_t len_b, bool is_protein __attribute__((unused)))
{
#ifdef HAVE_AVX2
  uint8_t dist;
  if (len_a > len_b) dist = bpm_256 (seq_a, seq_b, len_a, len_b);
  else               dist = bpm_256 (seq_b, seq_a, len_b, len_a);
  return (float) dist;
#else
  struct bignode* hash[1024];
  uint32_t i;
  float dist;
  unsigned int hv;
  for (i = 0;i < 1024;i++) hash[i] = 0;

  if (is_protein) {
    for (i = len_a-2;i--;){
      hv = (seq_a[i] << 5) + seq_a[i+1];
      hash[hv] = big_insert_hash(hash[hv],i);
      hv = (seq_a[i] << 5) + seq_a[i+2];
      hash[hv] = big_insert_hash(hash[hv],i);
    }
    dist = protein_wu_distance_calculation(hash,seq_b,len_b,len_a+len_b, 58.9);
  } else {
    for (i = len_a-5; i--;){
      hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+2]&3)<<4) + ((seq_a[i+3]&3)<<2) + (seq_a[i+4]&3);//ABCDE
      hash[hv] = big_insert_hash(hash[hv],i);
      hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+2]&3)<<4) + ((seq_a[i+3]&3)<<2) + (seq_a[i+5]&3);//ABCDF
      hash[hv] = big_insert_hash(hash[hv],i);
      hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+2]&3)<<4) + ((seq_a[i+4]&3)<<2) + (seq_a[i+5]&3);//ABCEF
      hash[hv] = big_insert_hash(hash[hv],i);
      hv = ((seq_a[i]&3)<<8) + ((seq_a[i+1]&3)<<6) + ((seq_a[i+3]&3)<<4) + ((seq_a[i+4]&3)<<2) + (seq_a[i+5]&3);//ABDEF
      hash[hv] = big_insert_hash(hash[hv],i);
      hv = ((seq_a[i]&3)<<8) + ((seq_a[i+2]&3)<<6) + ((seq_a[i+3]&3)<<4) + ((seq_a[i+4]&3)<<2) + (seq_a[i+5]&3);//ACDEF
      hash[hv] = big_insert_hash(hash[hv],i);
    }
    dist = dna_distance_calculation (hash, seq_b, len_b, len_a+len_b, 61.08);
  } // is not protein

  for (i = 1024;i--;) if (hash[i]) {
    big_remove_nodes(hash[i]);
    hash[i] = 0;
  }
  return dist;
#endif
}

float protein_wu_distance_calculation (struct bignode* hash[], const uint8_t* seq, const int seqlen, const int diagonals, const float mode)
{
  struct bignode* node_p;
  unsigned int* d = NULL;
  unsigned int* tmp = NULL;
  float out = 0.0;
  register int i,j;
  register int c;
  register int num;
  register unsigned int hv;

  d = malloc(sizeof(unsigned int)*diagonals);
  for (i = 0;i < diagonals;i++) d[i] = 0;
  for (i = seqlen-2;i--;) {
    hv = (seq[i] << 5) + seq[i+1];
    node_p = hash[hv];
    while(node_p) {
      tmp = node_p->pos;
      num = node_p->num;
      for(j = 0;j < num;j++) {
        c = tmp[j];
        d[c]++;
        c++;
        d[c]++;
      }
      node_p = node_p->next;
    }
    hv = (seq[i] << 5) + seq[i+2];

    node_p = hash[hv];

    while(node_p) {
      tmp = node_p->pos;
      num = node_p->num;
      for(j = 0;j < num;j++){
        c = tmp[j];
        d[c]++;
      }
      node_p = node_p->next;
    }
    d++;
  } // for i in seqlen-2
  d -= (seqlen-2);
  for (i = diagonals;i--;) if(d[i] > mode) out += d[i];
  free(d);
  return out;
}

float dna_distance_calculation (struct bignode* hash[], const uint8_t * p, const uint32_t seqlen, uint32_t diagonals, float mode)
{
  struct bignode* node_p;
  float out = 0.0;
  unsigned int* tmp = NULL;
  unsigned int* d = NULL;
  uint32_t i,j;
  unsigned int hv;

  d = (unsigned int*) biomcmc_malloc (sizeof(unsigned int) * diagonals);
  for (i = 0; i < diagonals; i++) d[i] = 0;
  for (i = seqlen-5;i--;){
    hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+4]&3);//ABCDE
    if (hash[hv]){
      node_p = hash[hv];
      while(node_p){
        tmp = node_p->pos;
        for(j = 0;j < node_p->num;j++) d[tmp[j]]++;
        node_p = node_p->next;
      }
    }

    hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+3]&3)<<2) + (p[i+5]&3);//ABCDF
    if (hash[hv]){
      node_p = hash[hv];
      while(node_p){
        tmp = node_p->pos;
        for(j = 0;j < node_p->num;j++) d[tmp[j]]++;
        node_p = node_p->next;
      }
    }
    hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+2]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABCEF
    if (hash[hv]){
      node_p = hash[hv];
      while(node_p){
        tmp = node_p->pos;
        for(j = 0;j < node_p->num;j++) d[tmp[j]]++;
        node_p = node_p->next;
      }
    }
    hv = ((p[i]&3)<<8) + ((p[i+1]&3)<<6) + ((p[i+3]&3)<<4)  + ((p[i+4]&3)<<2) + (p[i+5]&3);//ABDEF
    if (hash[hv]){
      node_p = hash[hv];
      while(node_p){
        tmp = node_p->pos;
        for(j = 0;j < node_p->num;j++) d[tmp[j]]++;
        node_p = node_p->next;
      }
    }
    hv = ((p[i]&3)<<8) + ((p[i+2]&3)<<6) + ((p[i+3]&3)<<4) + ((p[i+4]&3)<<2) + (p[i+5]&3);//ACDEF
    if (hash[hv]){
      node_p = hash[hv];
      while(node_p){
        tmp = node_p->pos;
        for(j = 0;j < node_p->num;j++) d[tmp[j]]++;
        node_p = node_p->next;
      }
    }

    d++;
  }
  d -= (seqlen-5);

  for (i = diagonals;i--;) if(d[i] > mode) out += d[i];

  if (d) free(d);
  return out;
}

struct bignode* big_insert_hash (struct bignode *n,const unsigned int pos)
{
  struct bignode* p = NULL;
  if (n) {
    if (n->num < NODESIZE) {
      n->pos[n->num] = pos;
      n->num++;
      return n;
    } else {
      p = (struct bignode*) biomcmc_malloc (sizeof(struct bignode));
      p->pos[0] = pos;
      p->num = 1;
      p->next = n;
    }
  } else {
    p = (struct bignode*) biomcmc_malloc (sizeof(struct bignode));
    p->pos[0] = pos;
    p->num = 1;
    p->next = n;
  }
  return p;
}

void big_remove_nodes(struct bignode *n)
{
  struct bignode* p = NULL;
  while (n) {
    p = n;
    n = n->next;
    if (p) free(p);
  }
}

