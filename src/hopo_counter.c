/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "hopo_counter.h"
#include "kseq.h"

#define MIN_TRACT_LENGTH 2
#define MAX_TRACT_LENGTH 34
#define TRACT_LENGTH_RANGE 32  // (MAX_TRACT_LENGTH - MIN_TRACT_LENGTH)

BMC2_KSEQ_INIT(gzFile, gzread);

static uint8_t dna_in_2_bits[256][2] = {{0xff}};
static char bit_2_dna[] = {'A', 'C', 'G', 'T'};

static void initialize_dna_to_bit_tables (void);
int compare_hopo_element_decreasing (const void *a, const void *b);
int compare_hopo_context (hopo_element a, hopo_element b);
int compare_hopo_context_within_distance (hopo_element a, hopo_element b, int max_distance);
int distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist);
hopo_counter new_hopo_counter (void);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int kmer_size, int min_hopo_size);
void add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, int context_size, uint8_t hopo_base_int, int hopo_size);
void copy_hopo_element_start_count (hopo_element to, hopo_element from);
void bin_similarity_one (hopo_counter hc, int start, double *result);
void bin_similarity_two (hopo_counter hc1, int start1, hopo_counter hc2, int start2, double *result);
void estimate_coverage_hopo_counter (hopo_counter hc);
void estimate_variance_hopo_counter (hopo_counter hc);

int
compare_hopo_element_decreasing (const void *a, const void *b)
{
  int result = ((hopo_element *)b)->base - ((hopo_element *)a)->base;
  if (result) return result; // sort first by homopolymer (A/T or G/C)
  result = ((hopo_element *)b)->context[0] - ((hopo_element *)a)->context[0];
  if (result) return result; // sort second by kmers 
  result = ((hopo_element *)b)->context[1] - ((hopo_element *)a)->context[1];
  if (result) return result; // sort second by kmers 
  return ((hopo_element *)b)->length - ((hopo_element *)a)->length; // same context and homopolymer base, thus sort by tract length
}

int
compare_hopo_context (hopo_element a, hopo_element b)
{
  int result = b.base - a.base;
  if (result) return result;
  result = b.context[0] - a.context[0];
  if (result) return result;
  return b.context[1] - a.context[1];
}

int
compare_hopo_context_within_distance (hopo_element a, hopo_element b, int max_distance)
{ 
  if ((b.base - a.base) != 0) return 2 * max_distance; // homopolymer tracts are different
  int distance = distance_between_context_kmer (&(b.context[0]), &(a.context[0]), 2 * max_distance);
  if (distance >= 2 * max_distance) return distance; 
  return distance + distance_between_context_kmer (&(b.context[1]), &(a.context[1]), 2 * max_distance - distance);
}

/*! \brief max_dist must be positive, and is the max allowed distance _per_ flanking region */
int
distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist)
{
  uint64_t d = *c1 ^ *c2; // XOR is one if bits are different, zero o.w. 
  int dist = 0;
  while (d && (dist < max_dist)) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference  
  return dist;
}

hopo_counter
new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, int kmer_size, int min_hopo_size)
{
  int i;
  hopo_counter hc_local = hc;

  gzFile fp = gzopen (filename, "r");
  bmc2_kseq_t *seq = bmc2_kseq_init (fp);
  if (!hc) {
    hc_local = new_hopo_counter ();
    name_length = strlen(filename);
    hc_local->name = (char*) biomcmc_malloc (sizeof (char) * (name_length + 1));
    strncpy (hc_local->name, filename, name_length);
    hc_local->name[name_length] = '\0';
  }
  if (hc_local->idx) biomcmc_error ("This counter has been compared to another; cannot add more reads to it");
  while ((i = bmc2_kseq_read (seq)) >= 0) update_hopo_counter_from_seq (hc_local, seq->seq.s, seq->seq.l, kmer_size, min_hopo_size); 
  bmc2_kseq_destroy(seq); // other kseq_t parameters: seq->name.s, seq->seq.l, seq->qual.l
  gzclose(fp);
  return hc_local;
}

hopo_counter
new_hopo_counter (void)
{
  hopo_counter hc = (hopo_counter) biomcmc_malloc (sizeof (struct hopo_counter_struct));
  hc->n_alloc = 32;
  hc->n_idx = hc->n_elem = 0;
  hc->variance[0] = hc->variance[1] = 1.;
  hc->coverage[0] = hc->coverage[1] = 0.;
  hc->elem = (hopo_element*) biomcmc_malloc (hc->n_alloc * sizeof (hopo_element));
  hc->ref_counter = 1;
  hc->name = hc->idx = NULL;
  if (dna_in_2_bits[0][0] == 0xff) initialize_dna_to_bit_tables (); 

  return hc;
}

void
del_hopo_counter (hopo_counter hc)
{
  if (!hc) return;
  if (--hc->ref_counter) return;
  if (hc->elem) free (hc->elem);
  if (hc->name) free (hc->name);
  if (hc->idx) free (hc->idx);
  free (hc);
  return;
}

static void
initialize_dna_to_bit_tables (void)
{
  int i;
  /* The ACGT is PAUP convention (and maybe DNAml, fastDNAml); PAML uses TCAG ordering */
  for (i = 0; i < 256; i++) dna_in_2_bits[i][0] = dna_in_2_bits[i][1] = 4;
  dna_in_2_bits['A'][0] = dna_in_2_bits['a'][0] = 0; dna_in_2_bits['A'][1] = dna_in_2_bits['a'][1] = 3;  /*  A  <-> T  */
  dna_in_2_bits['C'][0] = dna_in_2_bits['c'][0] = 1; dna_in_2_bits['C'][1] = dna_in_2_bits['c'][1] = 2;  /*  C  <-> G  */
  dna_in_2_bits['G'][0] = dna_in_2_bits['g'][0] = 2; dna_in_2_bits['G'][1] = dna_in_2_bits['g'][1] = 1;  /*  G  <-> C  */
  dna_in_2_bits['T'][0] = dna_in_2_bits['t'][0] = 3; dna_in_2_bits['T'][1] = dna_in_2_bits['t'][1] = 0;  /*  T  <-> A  */
  dna_in_2_bits['U'][0] = dna_in_2_bits['u'][0] = 3; dna_in_2_bits['U'][1] = dna_in_2_bits['u'][1] = 0;  /*  U  <-> A  */
}

// check if doesnt handle overlapping hopos (e.g. TGACCCAAAATGC -> CCC and AAAA)
void
update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int kmer_size, int min_hopo_size)
{
  int i, j, k, count_same = 0, start_mono = -1;
  uint8_t context[2 * kmer_size]; 
  char prev_char = '$';
  printf ("\nDBG:: %s (%d)\n", seq, seq_length); // DEBUG
  count_same = 0; // zero because previous is "$" 
  for (i = 0; i < (seq_length - kmer_size); i++) {
    if (seq[i] == prev_char) {
      count_same++;
      if ((count_same > min_hopo_size) && (start_mono >= kmer_size)) {
        while (i < (seq_length-kmer_size-1) && (seq[i+1] == prev_char)) { i++; count_same++; }
        for (j=start_mono; j <= i; j++) printf ("%c", seq[j]); //DEBUG
        if (dna_in_2_bits[(int)prev_char][0] < dna_in_2_bits[(int)prev_char][1]) { // A or C : forward strand
          k = 0;
          for (j = start_mono - kmer_size; j < start_mono; j++) context[k++] = dna_in_2_bits[ (int)seq[j] ][0]; 
          for (j = i + 1; j <= i + kmer_size; j++) context[k++] = dna_in_2_bits[ (int)seq[j] ][0]; 
          add_kmer_to_hopo_counter (hc, context, kmer_size, dna_in_2_bits[(int)prev_char][0], count_same);
        }
        else if (dna_in_2_bits[(int)prev_char][0] > dna_in_2_bits[(int)prev_char][1]) { // T/U or G : reverse strand
          k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
          for (j = i + kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          for (j = start_mono - 1; j >= start_mono - kmer_size; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          add_kmer_to_hopo_counter (hc, context, kmer_size, dna_in_2_bits[(int)prev_char][1], count_same);
        } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
        printf ("\t : ");
        for (j=0; j < 2*kmer_size; j++) printf ("%c", bit_2_dna[context[j]]); //DEBUG
        printf ("\n");
      } // if (count_same>2) [i.e. we found valid homopolymer]
    } else { 
      count_same = 1;
      prev_char = seq[i];
      start_mono = i;
    } // else (seq[i] == prev_char)
  } // for (i in seq[])
}

void
add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, int context_size, uint8_t hopo_base_int, int hopo_size)
{
  int i;
  if (hc->n_elem == hc->n_alloc) {
    hc->n_alloc *= 2;
    hc->elem = (hopo_element*) biomcmc_realloc ((hopo_element*)hc->elem, hc->n_alloc * sizeof (hopo_element));
  }
  hc->elem[hc->n_elem].base   = hopo_base_int; // zero (AT) or one (CG) but always store direction with A or C
  hc->elem[hc->n_elem].length = hopo_size;     // homopolymer track length, in bases
  hc->elem[hc->n_elem].count  = 1;
  hc->elem[hc->n_elem].context[0] = context[0]; // left context 
  for (i = 1; i < context_size; i++) hc->elem[hc->n_elem].context[0] |= (context[i] << (2 * i));
  hc->elem[hc->n_elem].context[1] = context[context_size]; // right context
  for (i = 1; i < context_size; i++) hc->elem[hc->n_elem].context[1] |= (context[context_size + i] << (2 * i));
  hc->n_elem++;
}

void
copy_hopo_element_start_count (hopo_element to, hopo_element from)
{
  to->base    =  from->base; 
  to->length  =  from->length;
  to->context[0] = from->context[0];
  to->context[1] = from->context[1];
  to->count = 1;
}

void
finalise_hopo_counter (hopo_counter hc)
{
  hopo_element *pivot, *efreq;
  int i, n1 = 0;
  qsort (hc->elem, hc->n_elem, sizeof (hopo_element), compare_hopo_element_decreasing);
  efreq = (hopo_element*) biomcmc_malloc (hc->n_elem * sizeof (hopo_element));
  copy_hopo_element_start_count (&(efreq[0]), &(hc->elem[0]));
  n1 = 0; hc->n_idx = 0;

  /* frequency of each polymer, in context, into efreq[] vector of elements  */
  for (i=1; i < hc->n_elem; i++) {
    if (compare_hopo_element_decreasing ((const void*) &(hc->elem[i-1]), (const void*) &(hc->elem[i]))) 
      copy_hopo_element_start_count (&(efreq[++n1]), &(hc->elem[i]));
    else efreq[n1].count++; // same context _and_tract length
  }
  /* new efreq[] becomes hc->elem[] */
  pivot = hc->elem;
  hc->n_alloc = hc->n_elem = n1 + 1;
  hc->elem = (hopo_element*) biomcmc_realloc ((hopo_element*) efreq, hc->n_alloc * sizeof (hopo_element));
  free (pivot);

  /* idx[] will have indices of distinct contexts (i.e. histogram of polymer lengths) */
  hc->idx = (int*) biomcmc_malloc ((hc->n_elem + 1) * sizeof (int));
  hc->idx[hc->n_idx++] = 0; 
  for (i=1; i < hc->n_elem; i++) if (compare_hopo_context (hc->elem[i-1], hc->elem[i])) hc->idx[hc->n_idx++] = i;
  hc->idx[hc->n_idx++] = i; // last element (plus one) is also on index (to avoid conditional checking */
  hc->idx = (int*) biomcmc_realloc ((int*) hc->idx, hc->n_idx* sizeof (int));

  // TODO: remove singletons
  estimate_coverage_hopo_counter (hc);
  estimate_variance_hopo_counter (hc); 
}

/* coverage is estimated through kmer frequency (kmer=each flanking region) */
void
estimate_coverage_hopo_counter (hopo_counter hc)
{
  int sum_count = 0, i, *kmer, *count;
  empfreq ef;
  kmer  = (int*) biomcmc_malloc (2 * hc->n_elem * sizeof (int));
  count = (int*) biomcmc_malloc (2 * hc->n_elem * sizeof (int));
  for (i=0; i < hc->n_elem; i++) { // several elements may share (integer version of) context
    kmer[ i]  = (int) (hc->elem[i].context[0] & ((1ULL << 31) -1)); // use 30 bits
    count[i] = hc->elem[i].count;
    kmer[ i + hc->n_elem] = (int) (hc->elem[i].context[1] & ((1ULL << 31) -1)); // use 30 bits
    count[i + hc->n_elem] = hc->elem[i].count;
    sum_count += count[i];
  }
  ef = new_empfreq_from_int_weighted (kmer, 2 * hc->n_elem, count);
  hc->coverage[0] = (double)(ef->i[0].freq);
  //hc->coverage[1] = (int)((double)(sum_count)/(double)(ef->n));
  hc->coverage[1] = (double)(sum_count);
  if (hc->coverage[0] < 1e-15) hc->coverage[0] = 1e-15;
  if (hc->coverage[1] < 1e-15) hc->coverage[1] = 1e-15;
  if (kmer) free (kmer);
  if (count) free (count);
  del_empfreq (ef);
}

void
estimate_variance_hopo_counter (hopo_counter hc)
{
  int i;
  hc->variance[0] = hc->variance[1] = 0.;
  for (i = 0; i < hc->n_idx -1; i++) bin_similarity_two (hc, i, hc, i, hc->variance);
  if (hc->variance[0] < 1e-15) hc->variance[0] = 1e-15;
  if (hc->variance[1] < 1e-15) hc->variance[1] = 1e-15;
  hc->variance[0] = sqrt (hc->variance[0]);
  hc->variance[1] = sqrt (hc->variance[1]);
}

void 
compare_hopo_counters (hopo_counter hc1, hopo_counter hc2, double *result)
{
  int order, j1, j2;

  if (!hc1->idx) finalise_hopo_counter (hc1); 
  if (!hc2->idx) finalise_hopo_counter (hc2); 

  result[0] = result[1] = 0.; /* TODO: result should have distances for all contexts (i.e. hc1->n_idx + hc2->n_idx) */
  for (j1 = j2 = 0; (j1 < hc1->n_idx-1) && (j2 < hc2->n_idx-1);) { // both are in decreasing order
    order = compare_hopo_context (hc1->elem[ hc1->idx[j1] ], hc2->elem[ hc2->idx[j2] ]); // 0 = they're the same; o.w. which comes first
    if      (order < 0) { bin_similarity_one (hc1, j1, result); j1++; }
    else if (order > 0) { bin_similarity_one (hc2, j2, result); j2++; }
    else {  bin_similarity_two (hc1, j1, hc2, j2, result); j1++; j2++;}
  }
  /* loop above ends when one of the two vectors finish; the other still has unmatched contexts */
  for (;j1 < hc1->n_idx-1; j1++) bin_similarity_one (hc1, j1, result);
  for (;j2 < hc2->n_idx-1; j2++) bin_similarity_one (hc2, j2, result);
  result[0] = fabs(result[0] / (hc1->variance[0] * hc2->variance[0]) - 1.);
  result[1] = fabs(result[1] / (hc1->variance[1] * hc2->variance[1]) - 1.);
}

void
bin_similarity_one (hopo_counter hc, int start, double *result)
{ // alternatives: find closest context (hard); assume unseen bin < min_tract_length (easier) (TRUE bin_simi_matrices dont have this issue) 
  int i;
  for (i = hc->idx[start]; i < hc->idx[start+1]; i++) {
    result[0] += (1 + hc->elem[i].length) * (double)(hc->elem[i].count)/hc->coverage[0]; 
    result[1] += (1 + hc->elem[i].length) * (double)(hc->elem[i].count)/hc->coverage[1]; 
  }
}

void
bin_similarity_two (hopo_counter hc1, int start1, hopo_counter hc2, int start2, double *result)
{
  int i1, i2;
  double x, y;
  /* bin_hist = srqt[ t(P-Q) A (P-Q) ] where A=1/cov_mat but we just calculate a 'weighted' Euclidean distance */
  for (i1 = hc1->idx[start1]; i1 < hc1->idx[start1+1]; i1++)  for (i2 = hc2->idx[start2]; i2 < hc2->idx[start2+1]; i2++) {
    x = (double)(1 + abs (hc1->elem[i1].length - hc2->elem[i2].length)); // absolute distance between hist bins plus one 
    y = (double)(hc1->elem[i1].count)/hc1->coverage[0] - (double)(hc2->elem[i2].count)/hc2->coverage[0]; 
    result[0] += x * fabs(y);
    y = (double)(hc1->elem[i1].count)/hc1->coverage[1] - (double)(hc2->elem[i2].count)/hc2->coverage[1]; 
    result[1] += x * fabs(y);
  }
}

void
print_debug_hopo_counter (hopo_counter hc)
{
  int i;
  for (i = 0; i < 4; i++) printf ("%3d %5d | ", hc->elem[i].length, hc->elem[i].count);
  printf ("%9.4lf %9.4lf %9.4lf %9.4lf\n", hc->coverage[0], hc->coverage[1], hc->variance[0], hc->variance[1]);
}

context_histogram_t
new_context_histogram_from_hopo_elem (hopo_element he)
{
  context_histogram_t ch = (context_histogram) biomcmc_malloc (sizeof (struct context_histogram_struct));
  ch->context = (uint64_t*) biomcmc_malloc (2 * sizeof (uint64_t));
  ch->n_context = 1; // how many context pairs 
  ch->count = (uint32_t*) biomcmc_malloc (TRACT_LENGTH_RANGE * sizeof (uint32_t));
  for (int i = 0; i < TRACT_LENGTH_RANGE; i++) ch->count = 0;
  ch->mode_context_id = 0; // location, in context[], of best context (the one with length of highest frequency)
  ch->genome_location = -1; // bwa location
  /* first context is best context: */
  ch->base = he.base;
  ch->l_1 = ch->l_2 = he.length - MIN_TRACT_LENGTH; 
  ch->mode_context_count = he.count; // frequency of best context
  ch->mode_context_length = he.length; // tract length of best context
  ch->context[0] = he.context[0];
  ch->context[1] = he.context[1];
  ch->count[he.length - MIN_TRACT_LENGTH] = he.count;
  return ch;
}

void
del_context_histogram (context_histogram_t ch)
{
  if (!ch) return;
  if (ch->context) free (ch->context);
  if (ch->count) free (ch->count);
  free (ch);
}

void
context_histogram_update_if_close (context_histogram_t ch, hopo_element he)
{
  if (/* he is close to all */) context_histogram_add_hopo_elem (ch, he);
}


void
context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he)
{
  int i, this_id = -1;
  /* Check if context is already on ch */
  for (i = 0; i < ch->n_context; i++) if ((he.context[0] == ch->context[2*i]) && (he.context[1] == ch->context[2*i+1])) {this_id = i; break; } 
  if (this_id < 0) { // new context
    ch->context = (uint64_t*) biomcmc_realloc ((uint64_t*) ch->context, (2 * (ch->n_context+1)) * sizeof (uint64_t));
    this_id = ch->n_context++;
    ch->context[2 * this_id]     = he.context[0];
    ch->context[2 * this_id + 1] = he.context[1];
  }
  /* if context+tract more frequent than observed so far, then this context is best */
  if (ch->mode_context_count < he.count) {
    ch->mode_context_count = he.count;
    ch->mode_context_length = he.length; // may be same length, diff context
    ch->mode_context_id = this_id; 
  }
  if (ch->l_1 > (he.length - MIN_TRACT_LENGTH)) ch->l_1 = he.length - MIN_TRACT_LENGTH; 
  if (ch->l_2 < (he.length - MIN_TRACT_LENGTH)) ch->l_2 = he.length - MIN_TRACT_LENGTH; 
  ch->count[he.length - MIN_TRACT_LENGTH] += he.count;
}
