/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * 
 * This file is based on the Kalign3  program, commit 
 * [7b3395a30d](https://github.com/TimoLassmann/kalign/tree/7b3395a30d60e994c9f2101bd2055cc3a426b7f7).
 * TimoLassmann/kalign is licensed under the GNU General Public License v3.0 or later.
 */

#include <xmmintrin.h>
#include <immintrin.h>
#include <stdalign.h>
#include "alignment_parameters.h"
#include "rng.h"
#include "float.h"

struct sort_struct{
  int len;
  int id;
};

/* only local; */
struct line_buffer {
  struct out_line** lines;
  int max_line_len;
  int alloc_num_lines;
  int num_line;
};

struct out_line {
  char* line;
  int block;
  int seq_id;
};

struct msa_seq* msa_seq_from_char_vector_string (char *string, size_t nchars, uint32_t id);
struct msa* read_char_vector_to_msa (char_vector dna);

void set_subm_gaps_DNA(struct aln_param* ap);
int clean_and_set_to_extern(struct alphabet* a);
void merge_codes (struct alphabet*a,const int X, const int Y);
int sort_by_len (const void *a, const void *b);
uint32_t* select_seqs (struct msa* msa, uint32_t num_anchor);

/* from rwalign.c */
void free_msa_seq(struct msa_seq* seq);
struct line_buffer* alloc_line_buffer(int max_line_len);
int resize_line_buffer(struct line_buffer* lb);
void free_line_buffer(struct line_buffer* lb);
static void set_sip_nsip(struct msa* msa);/* local helper functions  rwalign */
/* end rwalign.c */

struct aln_param* init_ap (int numseq)
{
  struct aln_param* ap =  (struct aln_param*) biomcmc_malloc (sizeof(struct aln_param));
  int i,j;

  ap->tree = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * (numseq*3+1));

  for(i = 0;i < (numseq*3+1);i++) ap->tree[i] = 0;

  ap->rng = init_rng(42);
  ap->subm = (float**) biomcmc_malloc (sizeof(float*) * 21);

  for (i = 21;i--;){
    ap->subm[i] = NULL;
    ap->subm[i] = (float*) biomcmc_malloc (sizeof(float) * 21);
    for (j = 21;j--;) ap->subm[i][j] = 0.0f;
  }
  set_subm_gaps_DNA (ap);
  return ap;
}

void free_ap (struct aln_param* ap)
{
  int i;
  if(!ap) return;
  if(ap->subm){
    for (i = 21;i--;) MFREE(ap->subm[i]);
    MFREE(ap->subm);
  }
  if(ap->rng) free_rng(ap->rng);
  if(ap->tree) free (ap->tree);
  MFREE(ap);
}

/* Timo: These are old parameters from kalign 2 */
void set_subm_gaps_DNA (struct aln_param* ap)
{
  int i,j;
  for(i = 0; i < 5; i++) for(j =0; j < 5;j++) ap->subm[i][j] = 283;
  //	A   91 -114  -31 -123    0  -43
  ap->subm[0][0] += 91;
  ap->subm[0][1] += -114;
  ap->subm[0][2] += -31;
  ap->subm[0][3] += -123;
  //	C -114  100 -125  -31    0  -43
  ap->subm[1][0] += -114;
  ap->subm[1][1] += 100;
  ap->subm[1][2] += -125;
  ap->subm[1][3] += -31;
  //	G  -31 -125  100 -114    0  -43
  ap->subm[2][0] += -31;
  ap->subm[2][1] += -125;
  ap->subm[2][2] += 100;
  ap->subm[2][3] += -114;
  //	T -123  -31 -114   91    0  -43
  ap->subm[3][0] += -123;
  ap->subm[3][1] += -31;
  ap->subm[3][2] += -114;
  ap->subm[3][3] += 91;

  ap->gpo = 217;
  ap->gpe = 39.4;
  ap->tgpe =  292.6;
  //param->secret = 28.3;
}

// alphabet.c

struct alphabet* create_dna_alphabet (void)
{
  struct alphabet* a = (struct alphabet*) biomcmc_malloc (sizeof(struct alphabet));
  int i;
  char dnacode[16] = "ACGTUNRYSWKMBDHV";
  int code = 0;

  for(i = 0; i < 128;i++) a->to_internal[i] = -1;
  for(i = 0; i < 32;i++)  a->to_external[i] = -1;
  for(i = 0; i < 16;i++) {
    a->to_internal[(int) dnacode[i]] = code;
    code++;
  }

  merge_codes(a,'U','T');
  merge_codes(a,'N','R');  /* R.................A or G */
  merge_codes(a,'N','Y');  /* Y.................C or T */
  merge_codes(a,'N','S');  /* S.................G or C */
  merge_codes(a,'N','W');  /* W.................A or T */
  merge_codes(a,'N','K');  /* K.................G or T */
  merge_codes(a,'N','M');  /* M.................A or C */
  merge_codes(a,'N','B');  /* B.................C or G or T */
  merge_codes(a,'N','D');  /* D.................A or G or T */
  merge_codes(a,'N','H');  /* H.................A or C or T */
  merge_codes(a,'N','V');  /* V.................A or C or G */

  clean_and_set_to_extern(a);
  return a;
}

void merge_codes (struct alphabet*a, const int X, const int Y)
{
  int min;
  min = MACRO_MIN(a->to_internal[X],a->to_internal[Y]);
  ASSERT(min != -1, "code not set!");
  a->to_internal[X] = min;
  a->to_internal[Y] = min;
}

int clean_and_set_to_extern (struct alphabet* a)
{
  int i;
  int code = 0;
  int8_t trans[32];
  for(i = 0; i < 32;i++) trans[i] = -1;

  for(i = 64; i < 96;i++) if(a->to_internal[i] != -1) trans[a->to_internal[i]] = 1;
  code = 0;
  for(i = 0; i < 32;i++) if(trans[i] == 1) {
    trans[i] = code;
    code++;
  }
  for(i = 64; i < 96;i++) if(a->to_internal[i] != -1) {
    a->to_internal[i] = trans[a->to_internal[i]];//a->to_internal[i]];
    a->to_internal[i+32] = a->to_internal[i];
  }
  for(i = 64;i < 96;i++) if(a->to_internal[i] != -1) a->to_external[a->to_internal[i]] = i;
  return OK;
}

// pick_anchor.c

uint32_t* pick_anchor(struct msa* msa, uint32_t* n)
{
  uint32_t* anchors = NULL, num_anchor = 0, powlog2;
  ASSERT(msa != NULL, "No alignment.");
  powlog2 = (uint32_t) pow(log2((double) msa->numseq), 2.0);
  num_anchor = MAX(MIN(32, msa->numseq), powlog2);
  anchors = select_seqs (msa, num_anchor);
  *n = num_anchor;
  return anchors;
}

uint32_t* select_seqs(struct msa* msa, uint32_t num_anchor)
{
  struct sort_struct** seq_sort = (struct sort_struct**) biomcmc_malloc (sizeof(struct sort_struct*) * msa->numseq);
  uint32_t* anchors = NULL;
  uint32_t i,stride;
  for(i = 0; i < msa->numseq;i++) {
    seq_sort[i] = (struct sort_struct*) biomcmc_malloc (sizeof(struct sort_struct));
    seq_sort[i]->id = i;
    seq_sort[i]->len = msa->sequences[i]->len;//  aln->sl[i];
  }

  qsort (seq_sort, msa->numseq, sizeof(struct sort_struct*),sort_by_len);
  anchors = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * num_anchor);
  stride = msa->numseq / num_anchor;
  for(i = 0; i < num_anchor;i++) anchors[i] = seq_sort[i*stride]->id;
  ASSERT(i == num_anchor,"Cound not select all anchors\tnum_anchor:%d\t numseq:%d",num_anchor,msa->numseq);
  for(i = 0; i < msa->numseq;i++) MFREE(seq_sort[i]);
  MFREE(seq_sort);
  return anchors;
}

int sort_by_len(const void *a, const void *b)
{
  struct sort_struct* const *one = a;
  struct sort_struct* const *two = b;
  if((*one)->len > (*two)->len) return -1;
  else  return 1;
}

// rwalign.c 

/**< only DNA sequences and integer IDs, no names */
struct msa* read_char_vector_to_msa (char_vector dna) 
{
  struct msa* msal = (struct msa*) biomcmc_malloc (sizeof (struct msa));
  uint32_t i;
  msal->sequences = NULL;
  msal->numseq = dna->nstrings;
  msal->plen = NULL;
  msal->sip = NULL;
  msal->nsip = NULL;
  msal->sequences = (struct msa_seq**) biomcmc_malloc (sizeof(struct msa_seq*) * msal->numseq);
  for (i = 0; i < msal->numseq; i++) msal->sequences[i] = msa_seq_from_char_vector_string (dna->string[i], dna->nchars[i], (uint32_t) i);
  convert_msa_to_internal (msal);
  set_sip_nsip (msal);
  return msal;
}

void set_sip_nsip (struct msa* msa)
{
  uint32_t i;
  ASSERT(msa!= NULL, "No msa");
  if(msa->plen){
    for (i = msa->num_profiles;i--;) if(msa->sip[i]) MFREE(msa->sip[i]);
    if(msa->plen) MFREE(msa->plen);
    if(msa->sip)  MFREE(msa->sip);
    if(msa->nsip) MFREE(msa->nsip);
    msa->plen = NULL;
    msa->sip = NULL;
    msa->nsip = NULL;
  }

  msa->num_profiles = (msa->numseq << 1 ) - 1;
  msa->sip = (int**) biomcmc_malloc (sizeof(int*) * msa->num_profiles);
  msa->nsip = (int*) biomcmc_malloc (sizeof(int) * msa->num_profiles);
  msa->plen = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * msa->num_profiles);
  for (i =0; i < msa->num_profiles; i++){
    msa->sip[i] = NULL;
    msa->nsip[i] = 0;
  }
  for(i = 0; i < msa->numseq; i++) {
    msa->sip[i] = (int*) biomcmc_malloc (sizeof(int) * msa->num_profiles);
    msa->nsip[i] = 1;
    msa->sip[i][0] = i;
    msa->plen[i] = 0;
  }
}

void convert_msa_to_internal (struct msa* msa)
{
  struct alphabet* a = NULL;
  struct msa_seq* seq = NULL;
  int8_t* t = NULL;
  uint32_t i, j;

  a = create_dna_alphabet();
  t = a->to_internal;
  for(i = 0; i <  msa->numseq;i++){
    seq = msa->sequences[i];
    for(j =0 ; j < seq->len;j++){
      if(t[(int) seq->seq[j]] == -1){
        WARNING_MSG("there should be no character not matching the alphabet");
        WARNING_MSG("offending character: >>>%c<<<", seq->seq[j]);
      }else{
        seq->s[j] = t[(int) seq->seq[j]];
      }
    }
  }
  if (a) free (a);
}

char_vector aligned_msa_to_charvector (struct msa* msa)
{
  uint32_t c, f, i, j, id;
  size_t aln_len = msa->sequences[0]->len;
  char_vector aligned;

  for (j = 0; j <= msa->sequences[0]->len; j++) aln_len += msa->sequences[0]->gaps[j];

  aligned = new_char_vector_fixed_length (msa->numseq, aln_len);

  for(i = 0; i < msa->numseq;i++) {
    id = msa->sequences[i]->id; /* I don't know if the order can change so better safe than sorry */
    f = 0;
    for (j = 0; j < msa->sequences[i]->len; j++) {
      for (c = 0; c < msa->sequences[i]->gaps[j]; c++)  aligned->string[id][f++] = '-';
      aligned->string[id][f++] = msa->sequences[i]->seq[j];
    }
    for (c = 0; c < msa->sequences[i]->gaps[j]; c++) aligned->string[id][f++] = '-';
  }
  return aligned;
}

void free_msa (struct msa* msa)
{
  int i;
  if(!msa) return;
  for(i = msa->numseq - 1; i>=0; --i) free_msa_seq (msa->sequences[i]);
  for (i = msa->num_profiles;i--;) if(msa->sip[i]) free(msa->sip[i]);
  if (msa->plen) free (msa->plen);
  if (msa->sip)  free (msa->sip);
  if (msa->nsip) free (msa->nsip);
  if (msa->sequences) free (msa->sequences);
}

struct msa_seq* msa_seq_from_char_vector_string (char *string, size_t nchars, uint32_t id)
{
  struct msa_seq* seq = (struct msa_seq*) biomcmc_malloc (sizeof (struct msa_seq));
  size_t i;
  seq->id = id;
  seq->len = 0;
  nchars++; // we incorporate null_terminate_sequences() whereby we add final 0

  seq->seq  = (char*) biomcmc_malloc (sizeof(char) * nchars);
  seq->gaps = (uint32_t*) biomcmc_malloc (sizeof(uint32_t) * (nchars + 1)); // in kalign3 it's one more than the dna (to allow trailing indels?)
  for(i = 0;i < nchars + 1; i++) seq->gaps[i] = 0;
  for(i = 0;i < nchars - 1; i++) { // minus one since we added one 
    if( isalpha((int) string[i])) seq->seq[seq->len++] = string[i];
    if( ispunct((int) string[i])) seq->gaps[seq->len]++; // no increment len++ here
  }
  seq->seq[seq->len] = 0; // null_terminate_sequences() in kalign3
  if (seq->len < nchars - 1) {
    seq->seq  = (char*) biomcmc_realloc ((char*)seq->seq, sizeof(char) * seq->len + 1);
    seq->gaps = (uint32_t*)  biomcmc_realloc ((uint32_t*)seq->gaps, sizeof(uint32_t) * (nchars + 2)); // in kalign3 it's one more than the dna (to allow trailing indels?)
  }
  seq->s  = (uint8_t*) biomcmc_malloc (sizeof(uint8_t) * seq->len + 1);
  return seq;
}

void free_msa_seq(struct msa_seq* seq)
{
  if (!seq) return;
  if (seq->seq)  free (seq->seq);
  if (seq->s)    free (seq->s);
  if (seq->gaps) free (seq->gaps);
  if (seq) free (seq);
}

struct line_buffer* alloc_line_buffer(int max_line_len)
{
  struct line_buffer* lb = NULL;
  int i;
  ASSERT(max_line_len > 60, "max_line_len:%d too small", max_line_len);
  MMALLOC(lb, sizeof(struct line_buffer));
  lb->alloc_num_lines = 1024;

  lb->num_line = 0;
  lb->lines = NULL;
  lb->max_line_len = max_line_len;
  MMALLOC(lb->lines, sizeof(struct out_line*) * lb->alloc_num_lines);

  for(i = 0; i < lb->alloc_num_lines;i++){
    lb->lines[i] = NULL;
    MMALLOC(lb->lines[i], sizeof(struct out_line));
    lb->lines[i]->block = 0;
    lb->lines[i]->seq_id =0;
    lb->lines[i]->line = NULL;
    MMALLOC(lb->lines[i]->line,sizeof(char) * lb->max_line_len);
  }

  return lb;
ERROR:
  return NULL;
}


int resize_line_buffer(struct line_buffer* lb)
{
  int old_len = 0;
  int i;
  old_len = lb->alloc_num_lines;
  lb->alloc_num_lines = lb->alloc_num_lines + 1024;

  MREALLOC(lb->lines, sizeof(struct out_line*) * lb->alloc_num_lines);

  for(i = old_len; i < lb->alloc_num_lines;i++){
    lb->lines[i] = NULL;
    MMALLOC(lb->lines[i], sizeof(struct out_line));
    lb->lines[i]->block = 0;
    lb->lines[i]->seq_id =0;
    lb->lines[i]->line = NULL;
    MMALLOC(lb->lines[i]->line,sizeof(char) * lb->max_line_len);
  }
  return OK;
ERROR:
  return FAIL;
}


void free_line_buffer(struct line_buffer* lb)
{
  int i;

  if(lb){
    for(i = 0; i < lb->alloc_num_lines;i++){
      MFREE(lb->lines[i]->line);
      MFREE(lb->lines[i]);
    }
    MFREE(lb->lines);
    MFREE(lb);
  }
}

// weave_alignment.c

void make_seq (struct msa* msa, uint32_t a, uint32_t b, int* path);
int update_gaps (uint32_t old_len, uint32_t *gis, uint32_t *newgaps);

void weave (struct msa* msa, int** map, uint32_t* tree)
{
  uint32_t i, a,b; 
  for (i = 0; i < (msa->numseq-1)*3;i +=3){
    a = tree[i];
    b = tree[i+1];
    make_seq (msa, a, b, map[tree[i+2]]);
  }
}

void make_seq (struct msa* msa, uint32_t a, uint32_t b, int* path)
{
  uint32_t i, *gap_a, *gap_b, c, posa = 0, posb = 0;

  gap_a = (uint32_t*) biomcmc_malloc ((path[0]+1)*sizeof(uint32_t));
  gap_b = (uint32_t*) biomcmc_malloc ((path[0]+1)*sizeof(uint32_t));

  for (c = (uint32_t) (path[0] + 1); c--;) gap_a[c] = gap_b[c] = 0;
  c = 1;
  while(path[c] != 3) {
    if (!path[c]) {
      posa++;
      posb++;
    } else {
      if (path[c] & 1){
        gap_a[posa] += 1;
        posb++;
      } else if (path[c] & 2){
        gap_b[posb] += 1;
        posa++;
      }
    }
    c++;
  }
  for (i = msa->nsip[a]; i--;) update_gaps (msa->sequences[msa->sip[a][i]]->len, msa->sequences[msa->sip[a][i]]->gaps, gap_a);
  for (i = msa->nsip[b]; i--;) update_gaps (msa->sequences[msa->sip[b][i]]->len, msa->sequences[msa->sip[b][i]]->gaps, gap_b);
  if (gap_a) free (gap_a);
  if (gap_b) free (gap_b);
}

int update_gaps (uint32_t old_len, uint32_t *gis, uint32_t *newgaps)
{
  uint32_t i, j, add = 0, rel_pos = 0;
  for (i = 0; i <= old_len; i++){
    add = 0;
    for (j = rel_pos; j <= rel_pos + gis[i]; j++) if (newgaps[j] != 0) add += newgaps[j];
    rel_pos += gis[i]+1;
    gis[i] += add;
  }
  return OK;
}

/* euclidean_distances.c -- functions taken from https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-float-vector-sum-on-x86 */
float hsum256_ps_avx(__m256 v);
float hsum_ps_sse3(__m128 v);

int edist_serial(const float* a,const float* b,const int len, float* ret)
{
  int i;
  float d = 0.0f;
  float t;

  for(i = 0; i < len;i++){
    t = (a[i] - b[i]);
    d += t *t;
  }
  *ret = sqrtf(d);
  return OK;
}

int edist_serial_d(const double* a,const double* b,const int len, double* ret)
{
  int i;
  double d = 0.0f;
  double t;

  for(i = 0; i < len;i++){
    t = (a[i] - b[i]);
    d += t *t;
  }
  *ret = sqrt(d);
  return OK;
}

#ifdef HAVE_AVX2
int edist_256(const float* a,const float* b, const int len, float* ret)
{
  float d = 0.0f;
  register int i;
  __m256 xmm1;// = _mm256_load_ps(a);
  __m256 xmm2;// = _mm256_load_ps(b);
  __m256 r = _mm256_set1_ps(0.0f);
  for(i = 0;i < len;i+=8){
    xmm1 = _mm256_load_ps(a);
    xmm2 = _mm256_load_ps(b);
    xmm1 =  _mm256_sub_ps(xmm1, xmm2);
    xmm1 = _mm256_mul_ps(xmm1, xmm1);
    r = _mm256_add_ps(r, xmm1);
    a+=8;
    b+=8;
  }
  d = hsum256_ps_avx(r);
  *ret = sqrtf(d);
  return OK;
}

float hsum256_ps_avx(__m256 v)
{
  __m128 vlow  = _mm256_castps256_ps128(v);
  __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
  vlow  = _mm_add_ps(vlow, vhigh);     // add the low 128
  return hsum_ps_sse3(vlow);         // and inline the sse3 version, which is optimal for AVX
  // (no wasted instructions, and all of them are the 4B minimum)
}

float hsum_ps_sse3(__m128 v)
{
  __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
  __m128 sums = _mm_add_ps(v, shuf);
  shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
  sums        = _mm_add_ss(sums, shuf);
  return        _mm_cvtss_f32(sums);
}

#endif

// (c) 2017 Johannes Soeding & Martin Steinegger, Gnu Public License version 3
// Rotate left macro: left circular shift by numbits within 16 bits
#define RoL(val, numbits) (val << numbits) ^ (val >> (16 - numbits))
// Transform each letter x[i] to a fixed random number RAND[x[i]]
// to ensure instantaneous mixing into the 16 bits  Do XOR with RAND[x[i]] and 5-bit rotate left for each i from 1 to k
uint16_t circ_hash(const uint8_t* x, const uint8_t length)
{
  const uint16_t RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec,
    0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146, 0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};// 16 bit random numbers
  register uint16_t h = 0x0;
  h = h^ RAND[x[0]];// XOR h and ki
  for (register uint8_t i = 1; i < length; ++i){
    h = RoL(h, 5);
    h ^= RAND[x[i]];// XOR h and ki
  }
  return h;
}

// Rolling hash variant for previous hash function:
// Computes hash value for next key x[0:length-1] from previous hash value hash( x[-1:length-2] ) and x_first = x[-1]
uint16_t circ_hash_next(const uint8_t * x,const uint8_t length,const uint8_t x_first, uint16_t h)
{
  const uint16_t RAND[21] = {0x4567, 0x23c6, 0x9869, 0x4873, 0xdc51, 0x5cff, 0x944a, 0x58ec,
    0x1f29, 0x7ccd, 0x58ba, 0xd7ab, 0x41f2, 0x1efb, 0xa9e3, 0xe146, 0x007c, 0x62c2, 0x0854, 0x27f8, 0x231b};// 16 bit random numbers
  h ^= RoL(RAND[x_first], (5*(length-1)) % 16);// undo INITIAL_VALUE and first letter x[0] of old key
  // circularly permute all letters x[1:length-1] to 5 positions to left
  h =  RoL(h, 5);// add new, last letter of new key x[1:length]
  h ^= RAND[x[length-1]];
  return h;
}

int shuffle_arr_r(int* arr,int n, struct rng_state* rng)
{
  int r;
  int i,j;
  int tmp;
  for (i = 0; i < n - 1; i++) {
    r = tl_random_int(rng,n);
    j = i +  r % (n-i);
    tmp = arr[j];
    arr[j] = arr[i];
    arr[i] = tmp;
  }
  return OK;
}

// kmeans.c

double** kmeans(double** data,int* cluster_assignment, int len_a, int len_b, int k)
{
  double** means = NULL;
  double** tmp = NULL;
  double** tmp2 = NULL;
  double** tmp_ptr = NULL;
  double best_solution;
  double score;
  double old_score;
  double min;
  double d;
  int min_index;
  int i,j;
  int num_attempts = 10000;
  int a_item;
  struct rng_state* rng;
  int* sel = NULL;
  double* n_item = NULL;
  ASSERT(k > 0,"K needs to be greater than one");
  ASSERT(k < len_a,"K larger than number of items");
  rng = init_rng(0);
  sel = galloc(sel,len_a);
  for(i = 0; i < len_a;i++) sel[i] = i;
  n_item = galloc(n_item,k);
  for(i = 0; i < k;i++) n_item[i] = 0.0;
  tmp = galloc(tmp,k,len_b,0.0f);
  tmp2= galloc(tmp2,k,len_b,0.0f);
  means= galloc(means,k,len_b,0.0f);
  best_solution = DBL_MAX;
  for(a_item = 0; a_item < num_attempts;a_item++){
    for(i = 0; i < k;i++) for(j = 0; j < len_b;j++) tmp[i][j] = 0.0;
    /* shuffle to select first k means */
    shuffle_arr_r(sel, len_a,rng);
    /* initial selection  */
    for(i = 0; i < k;i++) for(j = 0; j < len_b;j++){
      tmp[i][j] = data[sel[i]][j];
      tmp2[i][j] = 0.0;
    }
    score =1.0;
    old_score = 0.0;
    while(score != old_score) {
      for(i = 0; i < k;i++) n_item[i] = 0.0;
      for(i = 0; i < k;i++) for(j = 0; j < len_b;j++) tmp2[i][j] = 0.0;
      score = 0.0;
      for(i= 0; i < len_a;i++){
        min = DBL_MAX;
        min_index = -1;
        for(j = 0; j < k;j++) {
          edist_serial_d(data[i], tmp[j], len_b, &d);
          if(d < min) {
            min = d;
            min_index = j;
          }
        }
        n_item[min_index]++;
        score += min;
        for(j = 0; j < len_b;j++) tmp2[min_index][j] += data[i][j];
        if(cluster_assignment) cluster_assignment[i] = min_index;
      }
      /* calculate new means  */
      for(i = 0; i < k;i++) for(j = 0; j < len_b;j++) tmp2[i][j] /= n_item[i];
      /* switch */
      tmp_ptr = tmp;
      tmp= tmp2;
      tmp2 = tmp_ptr;
      if(old_score == score) break;
      old_score = score;
    }

    if(score < best_solution){
      fprintf(stdout,"%f better than %f\n", score,best_solution);
      best_solution = score;
      for(i = 0; i < k;i++) for(j = 0; j < len_b;j++) means[i][j] = tmp[i][j];
    }
  }
  for(i = 0; i < k;i++) n_item[i] = 0.0;
  for(i= 0; i < len_a;i++) if(cluster_assignment) n_item[cluster_assignment[i]]++;
  for(i = 0; i < k;i++) fprintf(stdout,"%d %f\n",i,n_item[i]);

  gfree(tmp);
  gfree(tmp2);
  gfree(n_item);
  gfree(sel);
  return means;
}
