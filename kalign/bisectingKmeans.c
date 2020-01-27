/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#include <xmmintrin.h>
#include "bisectingKmeans.h"
#include "alignment.h"
#include "alignment_parameters.h"

struct node{
  struct node* left;
  struct node* right;
  int id;
};

struct kmeans_result{
  uint32_t *sl, *sr, nl, nr;
  float score;
};

static struct kmeans_result* alloc_kmeans_result(uint32_t num_samples);
static void free_kmeans_results(struct kmeans_result* k);
struct node* upgma (float **dm, uint32_t* samples, uint32_t numseq);
struct node* alloc_node(void);

int label_internal (struct node*n, int label);
uint32_t* readbitree (struct node* p, uint32_t* tree);
struct node* bisecting_kmeans (struct msa* msa, struct node* n, float** dm, uint32_t* samples, uint32_t numseq, uint32_t num_anchors, uint32_t num_samples, struct rng_state* rng);

void build_tree_kmeans (struct msa* msa, struct aln_param* ap)
{
  struct node* root = NULL;
  float** dm = NULL;
  uint32_t *tree = NULL;
  uint32_t i, *samples = NULL, *anchors = NULL, num_anchors, numseq;

  ASSERT (msa != NULL, "No alignment.");
  ASSERT (ap != NULL, "No alignment parameters.");
  tree = ap->tree;
  numseq = msa->numseq;
  anchors = pick_anchor (msa, &num_anchors);
  dm = d_estimation (msa, anchors, num_anchors, 0);
  if (anchors) free (anchors);
  samples = (uint32_t*) biomcmc_malloc (sizeof(uint32_t)* numseq);
  for(i = 0; i < numseq; i++) samples[i] = i;
  root = bisecting_kmeans (msa, root, dm, samples, numseq, num_anchors, numseq, ap->rng);
  label_internal (root, numseq);
  ap->tree[0] = 1;
  ap->tree = readbitree(root, ap->tree);
  for (i = 0; i < (numseq*3);i++) tree[i] = tree[i+1];
  if (root) free (root);
  for(i =0 ; i < msa->numseq;i++)  _mm_free(dm[i]);
  if (dm) free (dm);
}

struct node* bisecting_kmeans (struct msa* msa, struct node* n, float** dm, uint32_t* samples, uint32_t numseq, uint32_t num_anchors, uint32_t num_samples, struct rng_state* rng)
{
  struct kmeans_result* res_tmp = NULL;
  struct kmeans_result* best = NULL;
  struct kmeans_result* res_ptr = NULL;

  uint32_t tries = 50, t_iter, r, *sl = NULL, *sr = NULL, num_l, num_r;
  uint32_t num_var, i,j,s, stop = 0;
  float* w = NULL;
  float* wl = NULL;
  float* wr = NULL;
  float* cl = NULL;
  float* cr = NULL;
  float dl = 0.0f;
  float dr = 0.0f;
  float score;

  if (num_samples < 100) {
    float** dm = NULL;
    dm = d_estimation (msa, samples, num_samples, 1);// anchors, num_anchors,1));
    n = upgma (dm, samples, num_samples);
    gfree (dm);
    if (samples) free (samples);
    return n;
  }

  num_var = num_anchors / 8;
  if( num_anchors%8) num_var++;
  num_var = num_var << 3;

  wr = _mm_malloc (sizeof(float) * num_var, 32);
  wl = _mm_malloc (sizeof(float) * num_var, 32);
  cr = _mm_malloc (sizeof(float) * num_var, 32);
  cl = _mm_malloc (sizeof(float) * num_var, 32);

  best = alloc_kmeans_result (num_samples);
  res_tmp = alloc_kmeans_result (num_samples);
  best->score = FLT_MAX;

  for (t_iter = 0; t_iter < tries; t_iter++) {
    res_tmp->score = FLT_MAX;
    sl = res_tmp->sl;
    sr = res_tmp->sr;
    w = _mm_malloc(sizeof(float) * num_var,32);
    for (i = 0; i < num_var; i++) w[i] = wr[i] = wl[i] = cr[i] = cl[i] = 0.0f;
    for(i = 0; i < num_samples;i++){
      s = samples[i];
      for(j = 0; j < num_anchors;j++) w[j] += dm[s][j];
    }
    for(j = 0; j < num_anchors;j++) w[j] /= num_samples;
    r = (uint32_t) tl_random_int (rng, num_samples);
    s = samples[r];
    for(j = 0; j < num_anchors;j++) cl[j] = dm[s][j];
    for(j = 0; j < num_anchors;j++) cr[j] = w[j] - (cl[j] - w[j]);
    _mm_free(w);

    /* check if cr == cl - we have identical sequences  */
    s = 0;
    for(j = 0; j < num_anchors; j++) if(fabsf(cl[j]-cr[j]) >  1.0E-6) {
      s = 1;
      break;
    }

    if(!s) {
      score = num_l = num_r = 0;
      sl[num_l] = samples[0];
      num_l++;
      for(i = 1 ; i <num_samples; i++) {
        sr[num_r] = samples[i];
        num_r++;
      }
    } else {
      w = NULL;
      while(1){
        stop++;
        if(stop == 10000) biomcmc_error ("kalign3 bisectingKmeans algorithm failed\n");
        num_l = num_r = 0;
        for(i = 0; i < num_anchors;i++) wr[i] = wl[i] = 0.0f;
        score = 0.0f;
        for(i = 0; i < num_samples;i++) {
          s = samples[i];
#ifdef HAVE_AVX2
          edist_256(dm[s], cl, num_anchors, &dl);
          edist_256(dm[s], cr, num_anchors, &dr);
#else
          edist_serial(dm[s], cl, num_anchors, &dl);
          edist_serial(dm[s], cr, num_anchors, &dr);
#endif
          score += MIN(dl,dr);

          if(dr < dl) {
            w = wr;
            sr[num_r] = s;
            num_r++;
          } else {
            w = wl;
            sl[num_l] = s;
            num_l++;
          }
          for(j = 0; j < num_anchors;j++) w[j] += dm[s][j];
        }

        for(j = 0; j < num_anchors;j++){
          wl[j] /= num_l;
          wr[j] /= num_r;
        }
        s = 0;

        for(j = 0; j < num_anchors;j++) {
          if(wl[j] != cl[j]) {
            s = 1;
            break;
          }
          if(wr[j] != cr[j]) {
            s = 1;
            break;

          }
        }
        if(s) {
          w = cl;
          cl = wl;
          wl = w;
          w = cr;
          cr = wr;
          wr = w;
        } else {
          break;
        }
      }
    }
    res_tmp->nl =  num_l;
    res_tmp->nr =  num_r;
    res_tmp->score = score;

    if(res_tmp->score < best->score){
      res_ptr = res_tmp;
      res_tmp = best;
      best = res_ptr;
    }
  }
  free_kmeans_results (res_tmp);

  sl = best->sl;
  sr = best->sr;
  num_l = best->nl;
  num_r = best->nr;
  if (best) free (best);

  _mm_free(wr);
  _mm_free(wl);
  _mm_free(cr);
  _mm_free(cl);
  if (samples) free (samples);
  n = alloc_node();
  n->left =  bisecting_kmeans(msa,n->left, dm, sl, numseq, num_anchors, num_l,rng);
  n->right = bisecting_kmeans(msa,n->right, dm, sr, numseq, num_anchors, num_r,rng);
  return n;
}

struct node* upgma (float **dm, uint32_t* samples, uint32_t numseq)
{
  struct node** tree = NULL;
  struct node* tmp = NULL;
  uint32_t i, j, node_a = 0, node_b = 0, cnode = numseq, numprofiles, *as = NULL;
  float max;

  numprofiles = (numseq << 1) - 1;
  as = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * numseq);
  for (i = numseq; i--;) as[i] = i+1;

  tree = (struct node**) biomcmc_malloc (sizeof(struct node*) * numseq);
  for (i = 0;i < numseq;i++){
    tree[i] = NULL;
    tree[i] = alloc_node();
    tree[i]->id = samples[i];
  }

  while (cnode != numprofiles) {
    max = FLT_MAX;
    for (i = 0;i < numseq-1; i++) if (as[i]) for ( j = i + 1;j < numseq;j++) if (as[j]) if (dm[i][j] < max) {
      max = dm[i][j];
      node_a = i;
      node_b = j;
    }
    tmp = NULL;
    tmp = alloc_node();
    tmp->left = tree[node_a];
    tmp->right = tree[node_b];

    tree[node_a] = tmp;
    tree[node_b] = NULL;

    /*deactivate  sequences to be joined*/
    as[node_a] = cnode+1;
    as[node_b] = 0;
    cnode++;

    /*calculate new distances*/
    for (j = numseq;j--;) if (j != node_b) dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])*0.5f;
    dm[node_a][node_a] = 0.0f;
    for (j = numseq;j--;) dm[j][node_a] = dm[node_a][j];
  }
  tmp = tree[node_a];
  if (tree) free (tree);
  if (as) free (as);
  return tmp;
}

struct node* alloc_node (void)
{
  struct node *n = (struct node*) biomcmc_malloc (sizeof(struct node));
  n->left = n->right = NULL;
  n->id = -1;
  return n;
}

int label_internal (struct node*n, int label)
{
  if(n->left)  label = label_internal(n->left, label);
  if(n->right) label = label_internal(n->right, label);
  if(n->id == -1){
    n->id = label;
    label++;
  }
  return label;
}

uint32_t* readbitree (struct node* p, uint32_t* tree)
{
  if(p->left)  tree = readbitree(p->left,tree);
  if(p->right) tree = readbitree(p->right,tree);

  if(p->left) if(p->right){
    tree[tree[0]] = p->left->id;
    tree[tree[0]+1] = p->right->id;
    tree[tree[0]+2] = p->id;
    tree[0] +=3;
    if (p->left)  free (p->left);
    if (p->right) free (p->right);
  }
  return tree;
}

struct kmeans_result* alloc_kmeans_result (uint32_t num_samples)
{
  struct kmeans_result* k = NULL;
  ASSERT(num_samples != 0, "No samples???");
  k =  (struct kmeans_result*) biomcmc_malloc (sizeof(struct kmeans_result));
  k->nl = k->nr = 0;
  k->sl = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * num_samples);
  k->sr = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * num_samples);
  k->score = FLT_MAX;
  return k;
}

void free_kmeans_results(struct kmeans_result* k)
{
  if(!k) return;
  if(k->sl) free (k->sl);
  if(k->sr) free (k->sr);
  free (k);
}
