/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "hopo_counter.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);

static uint8_t dna_in_2_bits[256][2] = {{0xff}};
//static char bit_2_dna[] = {'A', 'C', 'G', 'T'};

static void initialize_dna_to_bit_tables (void);
int compare_hopo_element_decreasing (const void *a, const void *b);
int compare_hopo_context (hopo_element a, hopo_element b);
hopo_counter new_hopo_counter (void);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int kmer_size, int min_hopo_size);
void add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, int context_size, uint8_t hopo_base_int, int hopo_size);
void finalise_hopo_counter (hopo_counter hc);
void estimate_coverage_hopo_counter (hopo_counter hc);

int
compare_hopo_element_decreasing (const void *a, const void *b)
{
  int result = ((hopo_element *)b)->base - ((hopo_element *)a)->base;
  if (result) return result; // sort first by homopolymer (A/T or G/C)
  result = ((hopo_element *)b)->context - ((hopo_element *)a)->context;
  if (result) return result; // sort second by kmers 
  return ((hopo_element *)b)->base_size - ((hopo_element *)a)->base_size; // same context and homopolymer base, thus sort by tract length
} 

int
compare_hopo_context (hopo_element a, hopo_element b)
{
  int result = b.base - a.base;
  if (result) return result;
  return b.context - a.context;
}

hopo_counter
new_hopo_counter_from_file (const char *filename)
{
  int i;
  gzFile fp = gzopen (filename, "r");
  kseq_t *seq = kseq_init (fp);
  hopo_counter hc = new_hopo_counter ();
  while ((i = kseq_read (seq)) >= 0) update_hopo_counter_from_seq (hc, seq->seq.s, seq->seq.l, 5, 4); // seq_l, kmer_size, min_homopol_size 
  finalise_hopo_counter (hc);
  kseq_destroy(seq); // other kseq_t parameters: seq->name.s, seq->seq.l, seq->qual.l
  gzclose(fp);
  return hc;
}

hopo_counter
new_hopo_counter (void)
{
  hopo_counter hc = (hopo_counter) biomcmc_malloc (sizeof (struct hopo_counter_struct));
  hc->n_alloc = 32;
  hc->n_idx = hc->n_elem = hc->coverage[0] = hc->coverage[1] = 0;
  hc->elem = (hopo_element*) biomcmc_malloc (hc->n_alloc * sizeof (hopo_element));
  hc->ref_counter = 1;
  hc->idx = NULL;
  if (dna_in_2_bits[0][0] == 0xff) initialize_dna_to_bit_tables (); 

  return hc;
}

void
del_hopo_counter (hopo_counter hc)
{
  if (!hc) return;
  if (--hc->ref_counter) return;
  if (hc->elem) free (hc->elem);
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

void
update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int kmer_size, int min_hopo_size)
{
  int i, j, k, count_same = 0, start_mono = -1;
  uint8_t context[2 * kmer_size]; 
  char prev_char = '$';
  // printf ("\nDBG:: %s (%d)\n", seq, seq_length);
  count_same = 0;
  for (i = 0; i < (seq_length - kmer_size); i++) {
    if (seq[i] == prev_char) {
      count_same++;
      if ((count_same > min_hopo_size) && (start_mono >= kmer_size)) {
        while (i < (seq_length-kmer_size-1) && (seq[i+1] == prev_char)) { i++; count_same++; }
        // for (j=start_mono; j <= i; j++) printf ("%c", seq[j]); //DEBUG
        if (dna_in_2_bits[(int)prev_char][0] < dna_in_2_bits[(int)prev_char][1]) { // A or C : forward strand
          k = 0;
          for (j = start_mono - kmer_size; j < start_mono; j++) context[k++] = dna_in_2_bits[ (int)seq[j] ][0]; 
          for (j = i + 1; j <= i + kmer_size; j++) context[k++] = dna_in_2_bits[ (int)seq[j] ][0]; 
          add_kmer_to_hopo_counter (hc, context, 2 * kmer_size, dna_in_2_bits[(int)prev_char][0], count_same);
        }
        else if (dna_in_2_bits[(int)prev_char][0] > dna_in_2_bits[(int)prev_char][1]) { // T/U or G : reverse strand
          k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
          for (j = i + kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          for (j = start_mono - 1; j >= start_mono - kmer_size; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          add_kmer_to_hopo_counter (hc, context, 2 * kmer_size, dna_in_2_bits[(int)prev_char][1], count_same);
        } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
        // for (j=0; j < 2*kmer_size; j++) printf ("%c", bit_2_dna[ context[j] ]); //DEBUG
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
  hc->elem[hc->n_elem].base      = hopo_base_int; // zero (AT) or one (CG)
  hc->elem[hc->n_elem].base_size = hopo_size;     // homopolymer track length, in bases
  hc->elem[hc->n_elem].count     = 1;
  // if (context_size < 32)  // assuming kmer < 16 bases, otherwise we need to create hash function
  hc->elem[hc->n_elem].context   = context[0];
  for (i = 1; i < context_size; i++) hc->elem[hc->n_elem].context |= (context[i] << (2 * i));
  hc->n_elem++;
}

void
finalise_hopo_counter (hopo_counter hc)
{
  hopo_element *pivot, *efreq;
  int i, n1 = 0;
  qsort (hc->elem, hc->n_elem, sizeof (hopo_element), compare_hopo_element_decreasing);
  efreq = (hopo_element*) biomcmc_malloc (hc->n_elem * sizeof (hopo_element));
  efreq[0].base      = hc->elem[0].base;
  efreq[0].base_size = hc->elem[0].base_size;
  efreq[0].context   = hc->elem[0].context;
  efreq[0].count = 1; n1 = 0; hc->n_idx = 0;

  /* frequency of each polymer, in context, into efreq[] vector of elements  */
  for (i=1; i < hc->n_elem; i++) {
    if (compare_hopo_element_decreasing ((const void*) &(hc->elem[i-1]), (const void*) &(hc->elem[i]))) {
      efreq[++n1].base   = hc->elem[i].base;
      efreq[n1].base_size = hc->elem[i].base_size;
      efreq[n1].context   = hc->elem[i].context;
      efreq[n1].count = 1;
    }
    else efreq[n1].count++;
  }
  /* new efreq[] becomes hc->elem[] */
  pivot = hc->elem;
  hc->n_alloc = hc->n_elem = n1 + 1;
  hc->elem = (hopo_element*) biomcmc_realloc ((hopo_element*) efreq, hc->n_alloc * sizeof (hopo_element));
  free (pivot);

  /* idx[] will have indices of distinct contexts (i.e. histogram of polymer lengths */
  hc->idx = (int*) biomcmc_malloc ((hc->n_elem + 1) * sizeof (int));
  hc->idx[hc->n_idx++] = 0; 
  for (i=1; i < hc->n_elem; i++) 
    if ((hc->elem[i-1].base != hc->elem[i].base) || (hc->elem[i-1].context != hc->elem[i].context)) hc->idx[hc->n_idx++] = i;
  hc->idx[hc->n_idx++] = i; // last element (plus one) is also on index (to avoid conditional checking */
  hc->idx = (int*) biomcmc_realloc ((int*) hc->idx, hc->n_idx* sizeof (int));

  estimate_coverage_hopo_counter (hc);
}

void
estimate_coverage_hopo_counter (hopo_counter hc)
{
  int sum_count = 0, i, *kmer, *count;
  empfreq ef;
  kmer  = (int*) biomcmc_malloc (hc->n_elem * sizeof (int));
  count = (int*) biomcmc_malloc (hc->n_elem * sizeof (int));
  for (i=0; i < hc->n_elem; i++) { // several elements may share (integer version of) context
    kmer[i]  = (int) (hc->elem[i].context & ((1ULL << 31) -1)); // use 30 bits
    count[i] = hc->elem[i].count;
    sum_count += count[i];
  }
  ef = new_empfreq_from_int_weighted (kmer, hc->n_elem, count);
  hc->coverage[0] = ef->i[0].freq;
  //hc->coverage[1] = (int)((double)(sum_count)/(double)(ef->n));
  hc->coverage[1] = sum_count;
  assert (hc->coverage[0] > 0);
  assert (hc->coverage[1] > 0);
  if (kmer) free (kmer);
  if (count) free (count);
  del_empfreq (ef);
}

void 
compare_hopo_counters (hopo_counter hc1, hopo_counter hc2, double *result)
{
  int order, j1, j2;
  double scale1[] = {hc1->coverage[0]/hc2->coverage[0], hc1->coverage[1]/hc2->coverage[1]};
  double scale2[] = {1./scale1[0], 1./scale1[1]};
  result[0] = result[1] = 0.;
  for (j1 = j2 = 0; (j1 < hc1->n_idx-1) && (j2 < hc2->n_idx-1);) { // both are in decreasing order
    order = compare_hopo_context (hc1->elem[ hc1->idx[j1] ], hc2->elem[ hc2->idx[j2] ]);
    if (order < 0)      { bin_similarity_one (hc1, hc1->idx[j1], hc1->idx[j1+1], scale1, result); j1++; }
    else if (order > 0) { bin_similarity_one (hc2, hc2->idx[j2], hc2->idx[j2+1], scale2, result); j2++; }
    else {  bin_similarity_two (hc1, hc1->idx[j1], hc1->idx[j1+1], hc2, hc2->idx[j2], hc2->idx[j2+1], scale1, result); j1++; j2++;}
  }
  // for(j1) bin_one; for (j2) bin_one;
}



void
print_debug_hopo_counter (hopo_counter hc)
{
  int i;
  for (i = 0; i < 4; i++) printf ("%3d %5d | ", hc->elem[i].base_size, hc->elem[i].count);
  printf ("%12d %12d\n", hc->coverage[0], hc->coverage[1]);
}
