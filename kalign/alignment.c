/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#include <xmmintrin.h>
#include "alignment.h"
#include "alignment_parameters.h"

#define MAX3(a,b,c) MAX(MAX(a,b),c)

struct states{
  float a;
  float ga;
  float gb;
};

struct hirsch_mem{
  struct states* f;
  struct states* b;
  int size;
  int starta, startb, enda, endb;
  uint32_t len_a, len_b;
};

struct dp_matrix{
  struct states* s;
  void* tb_mem;
  char** tb;
  int x;
  int y;
};

/* Memory allocation for forward and backward slices  */
struct hirsch_mem* hirsch_mem_alloc(int x);
int hirsch_mem_realloc(struct hirsch_mem* hm,int x);
void hirsch_mem_free(struct hirsch_mem* hm);

/* setting up fast data structures for alignment */
float* make_profile(struct aln_param* ap,const uint8_t* seq,const uint32_t len);
int set_gap_penalties(float* prof, uint32_t len,int nsip);

/* Main dyn programming functions */
int hirsch_ss_dyn_score(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, float* score);
void hirsch_align_two_ss_vector_score(const struct aln_param* ap, struct hirsch_mem* hm, int32_t old_cor[],float* score);
/* Align 2 sequences  */
int hirsch_ss_dyn(const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2, struct hirsch_mem* hm, int* hirsch_path);
void hirsch_align_two_ss_vector(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path, float input_states[], int32_t old_cor[]);
int foward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm);
int backward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm);

/* Align one sequence to a profile */
int hirsch_ps_dyn(const struct aln_param* ap, const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip);
int hirsch_align_two_ps_vector(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[], int32_t old_cor[],int sip);
int foward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip);
int backward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip);

/* Align 2 profiles  */
int hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
int hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[], int32_t old_cor[]);
int foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
int backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);

/* auxiliary dyn prog functions  */
int* mirror_hirsch_path(int* hirsch_path, uint32_t len_a, uint32_t len_b);
int* add_gap_info_to_hirsch_path(int* hirsch_path, uint32_t len_a, uint32_t len_b);
int update(const float* profa, const float* profb,float* newp, struct aln_param*ap, int* path,int sipa,int sipb);
int calc_seq_id(const int* path,int* a, int*b,float* dist);

int** hirschberg_alignment(struct msa* msa, struct aln_param* ap)
{
  struct hirsch_mem *hm = hirsch_mem_alloc (2048);
  uint32_t i,j,g,a,b,c, len_a, len_b, numseq, *tree = NULL;
  float** profile = NULL;
  int** map = NULL;  // not unsigned since initial value is -1

  g = msa->num_profiles;
  numseq = msa->numseq;
  tree = ap->tree;
  profile = (float**) biomcmc_malloc (sizeof(float*) * g);
  map = (int**) biomcmc_malloc (sizeof(int*) * g);
  for ( i = 0;i< g;i++) { profile[i] = NULL; map[i] = NULL;  }

  for (i = 0; i < (numseq-1); i++) {
    a = tree[i*3];
    b = tree[i*3+1];
    c = tree[i*3+2];
    if(a < numseq) len_a = msa->sequences[a]->len;//  aln->sl[a];
    else           len_a = msa->plen[a];
    if(b < numseq)  len_b = msa->sequences[b]->len;// aln->sl[b];
    else            len_b = msa->plen[b];

    g = (len_a > len_b)? len_a: len_b;
    map[c] = (int*) biomcmc_malloc (sizeof(int) * (g+2));
    hirsch_mem_realloc (hm, g);

    for (j = 0; j < (g+2); j++) map[c][j] = -1;

    if (a < numseq) profile[a] = make_profile(ap,msa->sequences[a]->s,len_a);
    else set_gap_penalties(profile[a],len_a,msa->nsip[b]);

    if (b < numseq) profile[b] = make_profile(ap,msa->sequences[b]->s,len_b);
    else set_gap_penalties(profile[b],len_b,msa->nsip[a]);

    hm->starta = hm->startb = 0;
    hm->enda = (int) len_a;
    hm->endb = (int) len_b;
    hm->len_a = len_a;
    hm->len_b = len_b;
    hm->f[0].a = hm->b[0].a = 0.0;
    hm->f[0].ga = hm->b[0].ga = -FLT_MAX;
    hm->f[0].gb = hm->b[0].gb = -FLT_MAX;
    if(a < numseq) {
      if(b < numseq) {
        hirsch_ss_dyn (ap, msa->sequences[a]->s, msa->sequences[b]->s, hm, map[c]);
      } else {
        hm->enda = (int) len_b;
        hm->endb = (int) len_a;
        hm->len_a = len_b;
        hm->len_b = len_a;
        hirsch_ps_dyn (ap,profile[b], msa->sequences[a]->s, hm, map[c], msa->nsip[b]);
        map[c] = mirror_hirsch_path(map[c],len_a,len_b);
      } 
    } else {
      if (b < numseq) {
        hirsch_ps_dyn (ap, profile[a], msa->sequences[b]->s , hm, map[c], msa->nsip[a]);
      } else {
        if (len_a < len_b) {
          hirsch_pp_dyn (profile[a], profile[b], hm, map[c]);
        } else {
          hm->enda = (int) len_b;
          hm->endb = (int) len_a;
          hm->len_a = len_b;
          hm->len_b = len_a;
          hirsch_pp_dyn (profile[b], profile[a], hm, map[c]);
          map[c] = mirror_hirsch_path (map[c],len_a,len_b);
        }
      }
    }

    map[c] = add_gap_info_to_hirsch_path (map[c],len_a,len_b);

    if(i != numseq-2){
      profile[c] = (float*) biomcmc_malloc (sizeof(float) * 64 * (map[c][0]+2));
      update (profile[a],profile[b],profile[c],ap,map[c],msa->nsip[a],msa->nsip[b]);
    }

    msa->plen[c] = map[c][0];
    msa->nsip[c] = msa->nsip[a] + msa->nsip[b];
    msa->sip[c] = (int*) biomcmc_malloc (sizeof(int) * (msa->nsip[a] + msa->nsip[b]));
    g =0;
    for (j = msa->nsip[a];j--;) msa->sip[c][g++] = msa->sip[a][j];
    for (j = msa->nsip[b];j--;) msa->sip[c][g++] = msa->sip[b][j];
    if (profile[a]) free (profile[a]);
    if (profile[b]) free (profile[b]);
  }
  if (profile) free (profile);
  hirsch_mem_free(hm);
  return map;
}

int hirsch_ss_dyn (const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path)
{
  int32_t mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
  float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
  int32_t old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

  if (hm->starta  >= hm->enda) return OK;
  if (hm->startb  >= hm->endb) return OK;

  hm->enda = mid;
  foward_hirsch_ss_dyn(ap,seq1,seq2,hm);
  hm->starta = mid;
  hm->enda = old_cor[1];
  backward_hirsch_ss_dyn(ap,seq1,seq2,hm);
  hirsch_align_two_ss_vector (ap,seq1,seq2,hm,hirsch_path,input_states,old_cor);
  return  OK;
}

int hirsch_ss_dyn_score (const struct aln_param* ap, const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm, float* score)
{
  int32_t mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
  int32_t old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

  if(hm->starta  >= hm->enda) return OK;
  if(hm->startb  >= hm->endb) return OK;

  hm->enda = mid;
  foward_hirsch_ss_dyn(ap,seq1,seq2,hm);
  hm->starta = mid;
  hm->enda = old_cor[1];
  backward_hirsch_ss_dyn(ap,seq1,seq2,hm);
  hirsch_align_two_ss_vector_score (ap,hm,old_cor,score);
  return  OK;
}

void hirsch_align_two_ss_vector_score (const struct aln_param* ap, struct hirsch_mem* hm, int32_t old_cor[], float* score)
{
  struct states* f = hm->f;
  struct states* b = hm->b;
  int32_t i;
  float gpo,gpe,tgpe;
  float max = -FLT_MAX;
  float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
  float sub = 0.0;
  gpo = ap->gpo;
  gpe = ap->gpe;
  tgpe = ap->tgpe;

  i = old_cor[2];
  for(i = old_cor[2]; i < old_cor[3];i++) {
    sub = fabsf(middle - (float)(i));
    sub /= 1000.;
    if(f[i].a+b[i].a-sub > max)        max = f[i].a+b[i].a  - sub;
    if(f[i].a+b[i].ga-gpo-sub > max)   max = f[i].a+b[i].ga - gpo-sub;
    if(f[i].a+b[i].gb -gpo-sub > max)  max = f[i].a+b[i].gb - gpo-sub;
    if(f[i].ga+b[i].a - gpo-sub > max) max = f[i].ga+b[i].a - gpo-sub;

    if(hm->startb == 0){
      if(f[i].gb+b[i].gb - tgpe-sub > max) max = f[i].gb+b[i].gb - tgpe-sub;
    }else{
      if(f[i].gb+b[i].gb - gpe -sub> max)  max = f[i].gb+b[i].gb - gpe-sub;
    }
    if(f[i].gb+b[i].a - gpo-sub > max)     max = f[i].gb+b[i].a  - gpo-sub;
  }
  i = old_cor[3];
  sub = fabsf(middle -(float)(i));
  sub /= 1000.;

  if (f[i].a+b[i].gb-gpo-sub > max) max = f[i].a+b[i].gb - gpo - sub;
  if (hm->endb == (int32_t)(hm->len_b)) {
    if(f[i].gb+b[i].gb -tgpe-sub > max) max = f[i].gb+b[i].gb - tgpe-sub;
  } else {
    if(f[i].gb+b[i].gb - gpe-sub > max) max = f[i].gb+b[i].gb - gpe-sub;
  }
  *score  = max;
}

void hirsch_align_two_ss_vector (const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[], int32_t old_cor[])
{
  struct states* f = hm->f;
  struct states* b = hm->b;
  int32_t i;
  int c;
  int transition = -1;
  float gpo,gpe,tgpe;

  gpo = ap->gpo;
  gpe = ap->gpe;
  tgpe = ap->tgpe;

  //code: // a -> a = 1 // a -> ga = 2 // a -> gb = 3 // ga ->ga = 4 // ga -> a = 5 //gb->gb = 6; //gb->a = 7;
  float max = -FLT_MAX;
  float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
  float sub = 0.0;

  i = old_cor[2];
  c = -1;
  for (i = old_cor[2]; i < old_cor[3]; i++) {
    sub = fabsf(middle - (float)(i));
    sub /= 1000.;
    if(f[i].a+b[i].a-sub > max)        { max = f[i].a+b[i].a  - sub;     transition = 1; c = i; }
    if(f[i].a+b[i].ga-gpo-sub > max)   { max = f[i].a+b[i].ga - gpo-sub; transition = 2; c = i; }
    if(f[i].a+b[i].gb -gpo-sub > max)  { max = f[i].a+b[i].gb - gpo-sub; transition = 3; c = i; }
    if(f[i].ga+b[i].a - gpo-sub > max) { max = f[i].ga+b[i].a - gpo-sub; transition = 5; c = i; }

    if(hm->startb == 0){
      if(f[i].gb+b[i].gb - tgpe-sub > max) { max = f[i].gb+b[i].gb - tgpe - sub; transition = 6; c = i; }
    } else {
      if(f[i].gb+b[i].gb - gpe -sub> max)  { max = f[i].gb+b[i].gb - gpe  - sub; transition = 6; c = i; }
    }
    if(f[i].gb+b[i].a - gpo-sub > max)     { max = f[i].gb+b[i].a  - gpo  - sub; transition = 7; c = i; }
  } // for (old_cor[2] ... old_cor[3])

  i = old_cor[3];
  sub = fabsf(middle - (float)(i));
  sub /= 1000.;

  if (f[i].a+b[i].gb-gpo-sub > max)      { max = f[i].a  + b[i].gb - gpo  - sub; transition = 3; c = i; }
  if (hm->endb == (int32_t)(hm->len_b)) {
    if (f[i].gb+b[i].gb -tgpe-sub > max) { max = f[i].gb + b[i].gb - tgpe - sub; transition = 6; c = i; }
  } else {
    if (f[i].gb+b[i].gb - gpe-sub > max) { max = f[i].gb + b[i].gb - gpe  - sub; transition = 6; c = i; }
  }

  switch(transition){
    case 1: //a -> a = 1
      hirsch_path[old_cor[4]] = c;
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      break;
    case 2:// a -> ga = 2
      hirsch_path[old_cor[4]] = c;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4];
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = hm->f[0].gb = -FLT_MAX;
      hm->f[0].ga = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      break;
    case 3:// a -> gb = 3
      hirsch_path[old_cor[4]] = c;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c;
      hm->endb = old_cor[3];
      hm->f[0].a = hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      break;
    case 5://ga -> a = 5
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = hm->b[0].gb = -FLT_MAX;
      hm->b[0].ga = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4];
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      break;
    case 6://gb->gb = 6;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c;
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c;
      hm->endb = old_cor[3];
      hm->f[0].a = hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      break;
    case 7://gb->a = 7;
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c;
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ss_dyn(ap,seq1,seq2,hm,hirsch_path);
      break;
  }
}

int foward_hirsch_ss_dyn (const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm)
{
  struct states* s = hm->f;
  float** subm =  ap->subm;
  float *subp = 0;
  const int32_t starta = hm->starta;
  const int32_t enda = hm->enda;
  const int32_t startb =hm->startb;
  const int32_t endb = hm->endb;
  register int32_t i = 0;
  register int32_t j = 0;
  register float pa = 0;
  register float pga = 0;
  register float pgb = 0;
  register float ca = 0;
  register float xa = 0;
  register float xga = 0;
  float gpo,gpe,tgpe;
  gpo = ap->gpo;
  gpe = ap->gpe;
  tgpe = ap->tgpe;

  s[startb].a = s[0].a;
  s[startb].ga = s[0].ga;
  s[startb].gb = s[0].gb;
  if(startb){
    for (j = startb+1; j < endb;j++){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);
      s[j].gb = -FLT_MAX;
    }
  }else {
    for (j = startb+1; j < endb;j++){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;
      s[j].gb = -FLT_MAX;
    }
  }
  s[endb].a = -FLT_MAX;
  s[endb].ga = -FLT_MAX;
  s[endb].gb = -FLT_MAX;

  seq2--;
  for (i = starta;i < enda;i++){
    subp = subm[seq1[i]];

    pa = s[startb].a;
    pga = s[startb].ga;
    pgb = s[startb].gb;
    s[startb].a = -FLT_MAX;
    s[startb].ga = -FLT_MAX;

    xa = s[startb].a;
    xga = s[startb].ga;

    if(startb) s[startb].gb = MAX(pgb - gpe,pa - gpo);
    else       s[startb].gb = MAX(pgb,pa) - tgpe;

    for (j = startb+1; j < endb;j++){
      ca = s[j].a;
      pa = MAX3(pa,pga-gpo,pgb-gpo);
      pa += subp[seq2[j]];
      s[j].a = pa;
      pga = s[j].ga;
      s[j].ga = MAX(xga-gpe,xa-gpo);
      pgb = s[j].gb;
      s[j].gb = MAX(pgb-gpe ,ca-gpo);
      pa = ca;
      xa = s[j].a;
      xga = s[j].ga;
    }
    ca = s[j].a;
    pa = MAX3(pa,pga-gpo,pgb-gpo);
    pa += subp[seq2[j]];
    s[j].a = pa;
    s[j].ga = -FLT_MAX;
    if (endb != (int32_t) (hm->len_b)) s[j].gb = MAX(s[j].gb-gpe ,ca-gpo);
    else                               s[j].gb = MAX(s[j].gb,ca) - tgpe;
  }
  return OK;
}

int backward_hirsch_ss_dyn(const struct aln_param* ap,const uint8_t* seq1,const uint8_t* seq2,struct hirsch_mem* hm)
{
  struct states* s = hm->b;
  float** subm = ap->subm;
  float *subp = NULL;
  const int32_t starta = hm->starta;
  const int32_t enda = hm->enda;
  const int32_t startb =hm->startb;
  const int32_t endb = hm->endb;
  register int32_t i = 0;
  register int32_t j = 0;
  register float pa = 0;
  register float pga = 0;
  register float pgb = 0;
  register float ca = 0;
  register float xa = 0;
  register float xga = 0;
  float gpo,gpe,tgpe;

  gpo = ap->gpo;
  gpe = ap->gpe;
  tgpe = ap->tgpe;
  s[endb].a = s[0].a ;
  s[endb].ga = s[0].ga;
  s[endb].gb = s[0].gb;
  //init of first row;
  if(endb != (int32_t)(hm->len_b)) {
    for(j = endb-1;j > startb;j--) {
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
      s[j].gb = -FLT_MAX;
    }
  }else{
    for(j = endb-1;j > startb;j--){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;
      s[j].gb = -FLT_MAX;
    }
  }

  s[startb].a = -FLT_MAX;
  s[startb].ga = -FLT_MAX;
  s[startb].gb = -FLT_MAX;

  i = enda-starta;
  seq1 += starta; // FIXME: does it work for more than 64k seqs?
  while(i--){
    subp = subm[seq1[i]];
    pa = s[endb].a;
    pga = s[endb].ga;
    pgb = s[endb].gb;
    s[endb].a = -FLT_MAX;
    s[endb].ga = -FLT_MAX;

    xa = s[endb].a;
    xga = s[endb].ga;

    if (endb != (int32_t)(hm->len_b)) s[endb].gb = MAX(pgb-gpe,pa-gpo);
    else                              s[endb].gb = MAX(pgb,pa)-tgpe;

    for(j = endb-1;j > startb;j--){
      ca = s[j].a;
      pa = MAX3(pa,pga - gpo,pgb-gpo);
      pa += subp[seq2[j]];
      s[j].a = pa;
      pga = s[j].ga;
      s[j].ga = MAX(xga-gpe,xa-gpo);
      pgb = s[j].gb;
      s[j].gb = MAX(pgb-gpe,ca-gpo);
      pa = ca;
      xa = s[j].a;
      xga = s[j].ga;
    }
    ca = s[j].a;
    pa = MAX3(pa,pga - gpo,pgb-gpo);
    pa += subp[seq2[j]];
    s[j].a = pa;
    s[j].ga = -FLT_MAX;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
    if(startb) s[j].gb = MAX(s[j].gb-gpe,ca-gpo);
    else       s[j].gb = MAX(s[j].gb,ca) - tgpe;
  }
  return OK;
}

int hirsch_ps_dyn (const struct aln_param* ap, const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip)
{
  int32_t mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
  float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
  int32_t old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

  if(hm->starta  >= hm->enda) return OK;
  if(hm->startb  >= hm->endb) return OK;

  hm->enda = mid;
  foward_hirsch_ps_dyn(ap,prof1,seq2,hm,sip);
  hm->starta = mid;
  hm->enda = old_cor[1];
  backward_hirsch_ps_dyn(ap,prof1,seq2,hm,sip);
  hirsch_align_two_ps_vector(ap,prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
  return OK;
}

int hirsch_align_two_ps_vector(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[], int32_t old_cor[],int sip)
{
  struct states* f = hm->f;
  struct states* b = hm->b;
  int32_t i, c; 
  int transition = -1;

  const float open = ap->gpo * sip;
  //code: // a-> a = 1 // a-> ga = 2 // a-> gb = 3 // ga->ga = 4 // ga-> a = 5 // gb->gb = 6 // gb->a = 7
  float max = -FLT_MAX;
  float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
  float sub = 0.0;
  prof1+= ((old_cor[4]+1)<<6);
  i = old_cor[2];
  c = -1;
  for(i = old_cor[2]; i < old_cor[3]; i++){
    sub = fabsf(middle - (float)(i));
    sub /= 1000.;
    if(f[i].a+b[i].a-sub> max)             { max = f[i].a+b[i].a -sub;           transition = 1; c = i; }
    if(f[i].a+b[i].ga-open-sub > max)      { max = f[i].a+b[i].ga-open-sub;      transition = 2; c = i; }
    if(f[i].a+b[i].gb+prof1[27]-sub > max) { max = f[i].a+b[i].gb+prof1[27]-sub; transition = 3; c = i; }
    if(f[i].ga+b[i].a-open-sub > max)      { max = f[i].ga+b[i].a-open-sub;      transition = 5; c = i; }

    if(hm->startb == 0){
      if(f[i].gb+b[i].gb+prof1[29]-sub > max) { max = f[i].gb+b[i].gb+prof1[29]-sub; transition = 6; c = i; }
    }else{
      if(f[i].gb+b[i].gb+prof1[28]-sub > max) { max = f[i].gb+b[i].gb+prof1[28]-sub; transition = 6; c = i;}
    }
    if(f[i].gb+b[i].a+prof1[-37]-sub > max)   { max = f[i].gb+b[i].a+prof1[-37]-sub; transition = 7;c = i; }
  }
  i = old_cor[3];
  sub = fabsf(middle - (float)(i));
  sub /= 1000.;
  if(f[i].a+b[i].gb+prof1[27]-sub > max)    { max = f[i].a+b[i].gb+prof1[27]-sub;  transition = 3; c = i; }
  if(hm->endb == (int32_t) hm->len_b){
    if(f[i].gb+b[i].gb+prof1[29]-sub > max) { max = f[i].gb+b[i].gb+prof1[29]-sub; transition = 6;c = i; }
  } else {
    if(f[i].gb+b[i].gb+prof1[28]-sub > max) { max = f[i].gb+b[i].gb+prof1[28]-sub; transition = 6; c = i; }
  }
  prof1-= ((old_cor[4]+1)<<6);

  switch(transition){
    case 1: //a -> a = 1
      hirsch_path[old_cor[4]] = c;
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      break;
    case 2:// a -> ga = 2
      hirsch_path[old_cor[4]] = c;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      //backward:
      hm->starta = old_cor[4];
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = -FLT_MAX;
      hm->f[0].ga = 0.0;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      break;
    case 3:// a -> gb = 3
      hirsch_path[old_cor[4]] = c;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c;
      hm->endb = old_cor[3];
      hm->f[0].a = -FLT_MAX;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      break;
    case 5://ga -> a = 5
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = -FLT_MAX;
      hm->b[0].ga = 0.0;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4];
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      break;
    case 6://gb->gb = 6;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = -FLT_MAX;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c;
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c;
      hm->endb = old_cor[3];
      hm->f[0].a = -FLT_MAX;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      break;
    case 7://gb->a = 7;
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = -FLT_MAX;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c;
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_ps_dyn(ap,prof1,seq2,hm,hirsch_path,sip);
      break;
  }
  return OK;
}

int foward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip)
{
  struct states* s = hm->f;
  register float pa = 0;
  register float pga = 0;
  register float pgb = 0;
  register float ca = 0;
  register float xa = 0;
  register float xga = 0;
  register int32_t i = 0;
  register int32_t j = 0;
  const float open = ap->gpo * sip;
  const float ext = ap->gpe *sip;
  const float text = ap->tgpe * sip;

  prof1 += (hm->starta)<< 6;
  s[hm->startb].a = s[0].a;
  s[hm->startb].ga = s[0].ga;
  s[hm->startb].gb = s[0].gb;
  if(hm->startb){
    for (j = hm->startb+1; j < hm->endb;j++){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
      s[j].gb = -FLT_MAX;
    }
  }else{
    for (j = hm->startb+1; j < hm->endb;j++){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
      s[j].gb = -FLT_MAX;
    }
  }
  s[hm->endb].a = -FLT_MAX;
  s[hm->endb].ga = -FLT_MAX;
  s[hm->endb].gb = -FLT_MAX;
  seq2--;
  for (i = hm->starta;i < hm->enda;i++){
    prof1 += 64;
    pa = s[hm->startb].a;
    pga = s[hm->startb].ga;
    pgb = s[hm->startb].gb;
    s[hm->startb].a = -FLT_MAX;
    s[hm->startb].ga = -FLT_MAX;
    xa = s[hm->startb].a;
    xga = s[hm->startb].ga;

    if(hm->startb) s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
    else           s[hm->startb].gb = MAX(pgb,pa)+prof1[29];
    for (j = hm->startb+1; j < hm->endb;j++){
      ca = s[j].a;
      pa = MAX3(pa,pga -open,pgb + prof1[-37]);
      pa += prof1[32 + seq2[j]];
      s[j].a = pa;
      pga = s[j].ga;
      s[j].ga = MAX(xga-ext,xa-open);
      pgb = s[j].gb;
      s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
      pa = ca;
      xa = s[j].a;
      xga = s[j].ga;
    }
    ca = s[j].a;
    pa = MAX3(pa,pga -open,pgb + prof1[-37]);
    pa += prof1[32 + seq2[j]];
    s[j].a = pa;
    s[j].ga = -FLT_MAX;//MAX(s[j-1].ga-ext,s[j-1].a-open);
    if (hm->endb != (int32_t)(hm->len_b)) s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
    else                                  s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
  }
  prof1 -= hm->enda << 6;
  return OK;
}

int backward_hirsch_ps_dyn(const struct aln_param* ap,const float* prof1,const uint8_t* seq2,struct hirsch_mem* hm,int sip)
{
  struct states* s = hm->b;
  register float pa = 0;
  register float pga = 0;
  register float pgb = 0;
  register float ca = 0;
  register float xa = 0;
  register float xga = 0;
  register int32_t i = 0;
  register int32_t j = 0;
  const float open = ap->gpo * sip;
  const float ext = ap->gpe *sip;
  const float text = ap->tgpe * sip;

  prof1 += (hm->enda+1) << 6;
  s[hm->endb].a = s[0].a;
  s[hm->endb].ga = s[0].ga;
  s[hm->endb].gb = s[0].gb;
  if(hm->endb != (int32_t)(hm->len_b)){
    for(j = hm->endb-1;j > hm->startb;j--){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
      s[j].gb = -FLT_MAX;
    }
  }else{
    for(j = hm->endb-1;j > hm->startb;j--){
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;
      s[j].gb = -FLT_MAX;
    }
  }
  s[hm->startb].a = -FLT_MAX;
  s[hm->startb].ga = -FLT_MAX;
  s[hm->startb].gb = -FLT_MAX;
  i = hm->enda-hm->starta;
  while(i--){
    prof1 -= 64;
    pa = s[hm->endb].a;
    pga = s[hm->endb].ga;
    pgb = s[hm->endb].gb;
    s[hm->endb].a = -FLT_MAX;
    s[hm->endb].ga = -FLT_MAX;
    xa = s[hm->endb].a;
    xga = s[hm->endb].ga;

    if(hm->endb != (int32_t)(hm->len_b)) s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
    else                                 s[hm->endb].gb = MAX(pgb,pa) +prof1[29];

    for(j = hm->endb-1;j > hm->startb;j--){
      ca = s[j].a;
      pa = MAX3(pa,pga - open,pgb +prof1[91]);
      pa += prof1[32 + seq2[j]];
      s[j].a = pa;
      pga = s[j].ga;
      s[j].ga = MAX(xga-ext,xa-open);
      pgb = s[j].gb;
      s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
      pa = ca;
      xa = s[j].a;
      xga = s[j].ga;
    }
    ca = s[j].a;
    pa = MAX3(pa,pga - open,pgb +prof1[91]);
    pa += prof1[32 + seq2[j]];
    s[j].a = pa;
    s[j].ga = -FLT_MAX;//MAX(s[j+1].ga-ext,s[j+1].a-open);
    if(hm->startb) s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
    else           s[j].gb = MAX(s[j].gb,ca)+prof1[29];
  }
  return OK;
}

int hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)
{
  int32_t mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
  float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
  int32_t old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

  if(hm->starta >= hm->enda) return OK;
  if(hm->startb >= hm->endb) return OK;
  hm->enda = mid;
  foward_hirsch_pp_dyn(prof1,prof2,hm);
  hm->starta = mid;
  hm->enda = old_cor[1];
  backward_hirsch_pp_dyn(prof1,prof2,hm);
  hirsch_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
  return OK;
}

int hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[], int32_t old_cor[])
{
  struct states* f = hm->f;
  struct states* b = hm->b;
  int32_t i;
  int c;
  int transition = -1;

  //code: // a-> a = 1 // a-> ga = 2 // a-> gb = 3 // ga->ga = 4 // ga-> a = 5 // gb->gb = 6; // gb->a = 7;
  float max = -FLT_MAX;
  float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
  float sub = 0.0;
  prof1+= ((old_cor[4]+1) << 6);
  prof2 += old_cor[2] << 6;
  i = old_cor[2];
  c = -1;
  for(i = old_cor[2]; i < old_cor[3]; i++) {
    sub = fabsf(middle - (float)(i));
    sub /= 1000.;
    prof2 += 64;
    if(f[i].a+b[i].a-sub > max){ max = f[i].a+b[i].a-sub; transition = 1; c = i; }
    if(f[i].a+b[i].ga+prof2[27]-sub > max){ max = f[i].a+b[i].ga+prof2[27]-sub; transition = 2; c = i; }
    if(f[i].a+b[i].gb+prof1[27] -sub> max){ max = f[i].a+b[i].gb+prof1[27]-sub; transition = 3; c = i; }
    if(f[i].ga+b[i].a+prof2[-37]-sub > max){ max = f[i].ga+b[i].a+prof2[-37]-sub; transition = 5; c = i; }
    if(hm->startb == 0){
      if(f[i].gb+b[i].gb+prof1[29]-sub > max){ max = f[i].gb+b[i].gb+prof1[29]-sub; transition = 6; c = i; }
    }else{
      if(f[i].gb+b[i].gb+prof1[28]-sub > max){ max = f[i].gb+b[i].gb+prof1[28]-sub; transition = 6; c = i; }
    }
    if(f[i].gb+b[i].a+prof1[-37]-sub > max){ max = f[i].gb+b[i].a+prof1[-37]-sub; transition = 7; c = i; }
  }
  i = old_cor[3];
  sub = fabsf(middle - (float)(i));
  sub /= 1000.;
  if(f[i].a+b[i].gb+prof1[27]-sub > max){ max = f[i].a+b[i].gb+prof1[27]-sub; transition = 3; c = i; }
  if(hm->endb == (int32_t)(hm->len_b)) {
    if(f[i].gb+b[i].gb+prof1[29]-sub > max){ max = f[i].gb+b[i].gb+prof1[29]-sub; transition = 6;  c = i; }
  }else{
    if(f[i].gb+b[i].gb+prof1[28]-sub > max){ max = f[i].gb+b[i].gb+prof1[28]-sub; transition = 6; c = i; }
  }

  prof1-= (old_cor[4]+1)<<6;
  prof2 -= old_cor[3] << 6;
  switch(transition){
    case 1: //a -> a = 1
      hirsch_path[old_cor[4]] = c;
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      break;
    case 2:// a -> ga = 2
      hirsch_path[old_cor[4]] = c;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4];
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = -FLT_MAX;
      hm->f[0].ga = 0.0;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      break;
    case 3:// a -> gb = 3
      hirsch_path[old_cor[4]] = c;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = 0.0;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c;
      hm->endb = old_cor[3];
      hm->f[0].a = -FLT_MAX;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      break;
    case 5://ga -> a = 5
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = -FLT_MAX;
      hm->b[0].ga = 0.0;
      hm->b[0].gb = -FLT_MAX;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4];
      hm->startb = old_cor[2];
      hm->endb = c-1;
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      break;
    case 6://gb->gb = 6;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = -FLT_MAX;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c;
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c;
      hm->endb = old_cor[3];
      hm->f[0].a = -FLT_MAX;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = 0.0;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      break;
    case 7://gb->a = 7;
      hirsch_path[old_cor[4]+1] = c+1;
      //foward:
      hm->f[0].a = input_states[0];
      hm->f[0].ga = input_states[1];
      hm->f[0].gb = input_states[2];
      hm->b[0].a = -FLT_MAX;
      hm->b[0].ga = -FLT_MAX;
      hm->b[0].gb = 0.0;
      hm->starta = old_cor[0];
      hm->enda = old_cor[4]-1;
      hm->startb = old_cor[2];
      hm->endb = c;
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      //backward:
      hm->starta = old_cor[4]+1;
      hm->enda = old_cor[1];
      hm->startb = c+1;
      hm->endb = old_cor[3];
      hm->f[0].a = 0.0;
      hm->f[0].ga = -FLT_MAX;
      hm->f[0].gb = -FLT_MAX;
      hm->b[0].a = input_states[3];
      hm->b[0].ga = input_states[4];
      hm->b[0].gb = input_states[5];
      hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
      break;
  }
  return OK;
}

int foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
  uint32_t freq[21];
  struct states* s = hm->f;
  register float pa = 0;
  register float pga = 0;
  register float pgb = 0;
  register float ca = 0;
  register float xa = 0;
  register float xga = 0;
  register int32_t i = 0;
  register int32_t j = 0;
  register int32_t c = 0;
  register int32_t f = 0;
  prof1 += (hm->starta) << 6;
  prof2 += (hm->startb) << 6;
  s[hm->startb].a = s[0].a;
  s[hm->startb].ga = s[0].ga;
  s[hm->startb].gb = s[0].gb;
  if(hm->startb){
    for (j = hm->startb+1; j < hm->endb;j++){
      prof2+=64;
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
      s[j].gb = -FLT_MAX;
    }
    prof2+=64;
  }else{
    for (j = hm->startb+1; j < hm->endb;j++){
      prof2+=64;
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];
      s[j].gb = -FLT_MAX;
    }
    prof2+=64;
  }
  prof2 -= (hm->endb-hm->startb) << 6;
  s[hm->endb].a = -FLT_MAX;
  s[hm->endb].ga = -FLT_MAX;
  s[hm->endb].gb = -FLT_MAX;

  for (i = hm->starta;i < hm->enda;i++){
    prof1 += 64;
    f = 0;
    for (j = 0;j < 21; j++) if(prof1[j]){ freq[f] = (uint32_t) j; f++; }
    f--;
    pa = s[hm->startb].a;
    pga = s[hm->startb].ga;
    pgb = s[hm->startb].gb;
    s[hm->startb].a = -FLT_MAX;
    s[hm->startb].ga = -FLT_MAX;
    xa = s[hm->startb].a;
    xga = s[hm->startb].ga;

    if(hm->startb) s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
    else s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];
    for (j = hm->startb+1; j < hm->endb;j++) {
      prof2 += 64;
      ca = s[j].a;
      pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);
      prof2 += 32;
      for (c = f+1;c--;) pa += prof1[freq[c]]*prof2[freq[c]]; // f+1 since first loop iteration already substracts one
      prof2 -= 32;
      s[j].a = pa;
      pga = s[j].ga;
      s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);
      pgb = s[j].gb;
      s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);
      pa = ca;
      xa = s[j].a;
      xga = s[j].ga;
    }
    prof2 += 64;
    ca = s[j].a;
    pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);
    prof2 += 32;
    for (c = f+1;c--;) pa += prof1[freq[c]]*prof2[freq[c]];
    prof2 -= 32;
    s[j].a = pa;
    s[j].ga = -FLT_MAX;
    if (hm->endb != (int32_t)(hm->len_b)) s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
    else s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
    prof2 -= (hm->endb-hm->startb) << 6;
  }
  prof1 -=  (hm->enda) << 6;
  return OK;
}

int backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
  uint32_t freq[21];
  struct states* s = hm->b;
  register float pa = 0;
  register float pga = 0;
  register float pgb = 0;
  register float ca = 0;
  register float xa = 0;
  register float xga = 0;
  register int32_t i = 0;
  register int32_t j = 0;
  register int32_t c = 0;
  register int32_t f = 0;

  prof1 += (hm->enda+1) << 6;
  prof2 += (hm->endb+1) << 6;
  s[hm->endb].a = s[0].a;
  s[hm->endb].ga = s[0].ga;
  s[hm->endb].gb = s[0].gb;
  if(hm->endb != (int32_t)(hm->len_b)) {
    for(j = hm->endb-1;j > hm->startb;j--){
      prof2 -= 64;
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);
      s[j].gb = -FLT_MAX;
    }
    prof2 -= 64;
  }else{
    for(j = hm->endb-1;j > hm->startb;j--){
      prof2 -= 64;
      s[j].a = -FLT_MAX;
      s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];
      s[j].gb = -FLT_MAX;
    }
    prof2 -= 64;
  }
  s[hm->startb].a = -FLT_MAX;
  s[hm->startb].ga = -FLT_MAX;
  s[hm->startb].gb = -FLT_MAX;
  i = hm->enda-hm->starta;
  while(i--){
    prof1 -= 64;
    f = 0;
    for (j = 0;j < 21; j++) if(prof1[j]) { freq[f] = (uint32_t) j; f++; }
    f--;

    pa = s[hm->endb].a;
    pga = s[hm->endb].ga;
    pgb = s[hm->endb].gb;
    s[hm->endb].a = -FLT_MAX;
    s[hm->endb].ga = -FLT_MAX;
    xa = s[hm->endb].a;
    xga = s[hm->endb].ga;
    if(hm->endb != (int32_t)(hm->len_b)) s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);
    else s[hm->endb].gb = MAX(pgb,pa)+prof1[29];

    prof2 += (hm->endb-hm->startb) << 6;
    for(j = hm->endb-1;j > hm->startb;j--){
      prof2 -= 64;
      ca = s[j].a;
      pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
      prof2 += 32;
      for (c = f+1;c--;) pa += prof1[freq[c]]*prof2[freq[c]];
      prof2 -= 32;
      s[j].a = pa;
      pga = s[j].ga;
      s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);
      pgb = s[j].gb;
      s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);
      pa = ca;
      xa = s[j].a;
      xga = s[j].ga;
    }
    prof2 -= 64;
    ca = s[j].a;

    pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
    prof2 += 32;
    for (c = f+1;c--;) pa += prof1[freq[c]]*prof2[freq[c]];
    prof2 -= 32;
    s[j].a = pa;
    s[j].ga = -FLT_MAX;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
    if(hm->startb) s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
    else s[j].gb = MAX(s[j].gb,ca)+prof1[29];
  }
  return OK;
}

int* mirror_hirsch_path(int* hirsch_path, uint32_t len_a, uint32_t len_b)
{
  uint32_t i;
  int *np = (int*) biomcmc_malloc (sizeof(int) * (len_a+2));
  for(i =0; i < len_a+2; i++) np[i] = -1;
  for(i = 1; i <= len_b; i++) if(hirsch_path[i] != -1) np[hirsch_path[i]] = (int) i;
  if (hirsch_path) free (hirsch_path);
  return np;
}

int* add_gap_info_to_hirsch_path(int* hirsch_path, uint32_t len_a, uint32_t len_b)
{
  int b = 0, a = 0;
  uint32_t i, j;
  int *np = (int*) biomcmc_malloc (sizeof(int) * (len_b + len_a + 2));
  for(i =0; i < len_a + len_b +2; i++) np[i] = 0;
  j = 1;
  b = -1;
  if(hirsch_path[1] == -1){
    np[j] = 2;
    j++;
  }else{
    if(hirsch_path[1] != 1){
      for ( a = 0; a < hirsch_path[1] -1; a++){
        np[j] = 1;
        j++;
      }
      np[j] = 0;
      j++;
    }else{
      np[j] = 0;
      j++;
    }
  }
  b = hirsch_path[1];

  for(i = 2; i <= len_a; i++){
    if(hirsch_path[i] == -1){
      np[j] = 2;
      j++;
    }else{
      if(hirsch_path[i]-1 != b && b != -1){
        for ( a = 0; a < hirsch_path[i] - b - 1; a++){
          np[j] = 1;
          j++;
        }
        np[j] = 0;
        j++;
      }else{
        np[j] = 0;
        j++;
      }
    }
    b = hirsch_path[i];
  }

  if(hirsch_path[len_a] < (int) (len_b) && hirsch_path[len_a] != -1) for ( a = 0; a < (int)(len_b) - hirsch_path[len_a]; a++) {
    np[j] = 1;
    j++;
  }
  np[0] = j-1;
  np[j] = 3;

  np = (int*) biomcmc_realloc ((int*)np, sizeof(int) * (np[0]+2));
  if (hirsch_path) free (hirsch_path);

  //add gap info..
  i = 2;
  while(np[i] != 3){
    if ((np[i-1] &3) && !(np[i] & 3)){
      if(np[i-1] & 8) np[i-1] += 8;
      else            np[i-1] |= 16;
    }else if (!(np[i-1] & 3) &&(np[i] &3)){
      np[i] |= 4;
    }else if ((np[i-1] & 1) && (np[i] & 1)){
      np[i] |= 8;
    }else if ((np[i-1] & 2) && (np[i] & 2)){
      np[i] |= 8;
    }
    i++;
  }
  //add terminal gap...
  a = 1;
  while(np[a] != 0) { np[a] |= 32; a++; }
  a = np[0];
  while(np[a] != 0) { np[a] |= 32; a--; }
  return np;
}

float* make_profile (struct aln_param* ap,const uint8_t* seq,const uint32_t len)
{
  uint32_t i,j,c;
  float** subm =  ap->subm;
  float* prof = NULL;
  float gpo,gpe,tgpe;
  gpo = ap->gpo;
  gpe = ap->gpe;
  tgpe = ap->tgpe;
  prof = (float*) biomcmc_malloc (sizeof(float)*(len+2) * 64);
  prof += (64 *(len+1));
  for (i = 0;i < 64;i++) prof[i] = 0;
  prof[23+32] = -gpo;
  prof[24+32] = -gpe;
  prof[25+32] = -tgpe;

  i = len;
  while(i--) {
    prof -= 64;
    for (j = 0;j < 64;j++) prof[j] = 0;
    c = seq[i];
    prof[c] += 1;
    prof += 32;
    for(j = 21;j--;) prof[j] = subm[c][j];
    prof[23] = -gpo;
    prof[24] = -gpe;
    prof[25] = -tgpe;
    prof -= 32;
  }
  prof -= 64;
  for (i = 0;i < 64;i++) prof[i] = 0;
  prof[23+32] = -gpo;
  prof[24+32] = -gpe;
  prof[25+32] = -tgpe;
  return prof;
}

int set_gap_penalties(float* prof, uint32_t len,int nsip)
{
  uint32_t i;
  prof +=  (64 *(len+1));
  prof[27] = prof[55] * (float)(nsip);//gap open or close  23
  prof[28] = prof[56] * (float)(nsip);//gap extention 24
  prof[29] = prof[57] * (float)(nsip);//gap open or close 25
  i = len+1;
  while(i--){
    prof -= 64;
    prof[27] = prof[55] * (float)(nsip);//gap open or close
    prof[28] = prof[56] * (float)(nsip);//gap extention
    prof[29] = prof[57] * (float)(nsip);//gap open or close
  }
  return OK;
}

int update (const float* profa, const float* profb,float* newp, struct aln_param *ap, int* path, int sipa, int sipb)
{
  int i,j,c;
  for (i = 64; i--;) newp[i] = profa[i] + profb[i];
  profa += 64;
  profb += 64;
  newp += 64;
  c = 1;
  while(path[c] != 3){
    if (!path[c]){
      for (i = 64; i--;) newp[i] = profa[i] + profb[i];
      profa += 64;
      profb += 64;
    }
    if (path[c] & 1){
      for (i = 64; i--;) newp[i] = profb[i];
      profb += 64;
      if(!(path[c] & 20)){
        if(path[c] & 32){
          newp[25] += sipa;//1;
          i = ap->tgpe*sipa;
        }else{
          newp[24] += sipa;//1;
          i = ap->gpe*sipa;
        }

        for (j = 32; j < 55;j++) newp[j] -=i;
      }else{
        if (path[c] & 16){
          if(path[c] & 32){
            newp[25] += sipa;//1;
            i = ap->tgpe*sipa;
            newp[23] += sipa;//1;
            i += ap->gpo*sipa;
          }else{
            newp[23] += sipa;//1;
            i = ap->gpo*sipa;
          }
          for (j = 32; j < 55;j++) newp[j] -=i;
        }
        if (path[c] & 4){
          if(path[c] & 32){
            newp[25] += sipa;//1;
            i = ap->tgpe*sipa;
            newp[23] += sipa;//1;
            i += ap->gpo*sipa;
          }else{
            newp[23] += sipa;//1;
            i = ap->gpo*sipa;
          }
          for (j = 32; j < 55;j++) newp[j] -=i;
        }
      }
    }
    if (path[c] & 2) {
      for (i = 64; i--;) newp[i] = profa[i];
      profa+=64;
      if(!(path[c] & 20)){
        if(path[c] & 32){
          newp[25] += sipb;//1;
          i = ap->tgpe*sipb;
        }else{
          newp[24] += sipb;//1;
          i = ap->gpe*sipb;
        }
        for (j = 32; j < 55;j++) newp[j] -=i;
      }else{
        if (path[c] & 16){
          if(path[c] & 32){
            newp[25] += sipb;//1;
            i =  ap->tgpe*sipb;
            newp[23] += sipb;//1;
            i +=  ap->gpo*sipb;
          }else{
            newp[23] += sipb;//1;
            i =  ap->gpo*sipb;
          }
          for (j = 32; j < 55;j++) newp[j] -=i;
        }
        if (path[c] & 4){
          if(path[c] & 32){
            newp[25] += sipb;//1;
            i = ap->tgpe*sipb;
            newp[23] += sipb;//1;
            i += ap->gpo*sipb;
          }else{
            newp[23] += sipb;//1;
            i = ap->gpo*sipb;
          }
          for (j = 32; j < 55;j++) newp[j] -=i;
        }
      }
    }
    newp += 64;
    c++;
  }
  for (i = 64; i--;) newp[i] =  profa[i] + profb[i];
  newp -= (path[0]+1) *64;
  return OK;
}

int calc_seq_id(const int *path,int *a, int *b, float* dist)
{
  int i,j,c;
  float id = 0.0f;
  float len = 0.0f;
  i = 0;
  j = 0;
  c = 1;

  while(path[c] != 3){
    if (!path[c]){
      if(a[i] == b[j]) id++;
      i++;
      j++;
    }
    if (path[c] & 1) j++;
    if (path[c] & 2) i++;
    len += 1.0f;
    c++;
  }
  *dist = id / len;
  return OK;
}

struct hirsch_mem* hirsch_mem_alloc(int x)
{
  struct hirsch_mem* hm = NULL;
  hm = (struct hirsch_mem*) biomcmc_malloc (sizeof(struct hirsch_mem));
  hm->starta = 0;
  hm->startb = 0;
  hm->enda = 0;
  hm->endb = 0;
  hm->size = x+1;
  hm->len_a = 0;
  hm->len_b = 0;
  hm->f = (struct states*) biomcmc_malloc (sizeof(struct states) * (x+1));
  hm->b = (struct states*) biomcmc_malloc (sizeof(struct states) * (x+1));
  return hm;
}

int hirsch_mem_realloc(struct hirsch_mem* hm, int x)
{
  if((x+1) > hm->size) {
    hm->size = x+1;
    hm->f = (struct states*) biomcmc_realloc ((struct states*) (hm->f), sizeof(struct states)* (x+1));
    hm->b = (struct states*) biomcmc_realloc ((struct states*) (hm->b), sizeof(struct states)* (x+1));
  }
  return OK;
}

void hirsch_mem_free(struct hirsch_mem* hm)
{
  if(!hm) return;
  if (hm->f) free (hm->f);
  if (hm->b) free (hm->b);
  free (hm);
}
