/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "hopo_counter.h"
#include "kseq.h"

BMC2_KSEQ_INIT(gzFile, gzread);

static uint8_t dna_in_2_bits[256][2] = {{0xff}};
static char bit_2_dna[] = {'A', 'C', 'G', 'T'};

static void initialize_dna_to_bit_tables (void);
int compare_hopo_element_decreasing (const void *a, const void *b);
int compare_hopo_context (hopo_element a, hopo_element b);
int compare_hopo_context_within_distance (hopo_element a, hopo_element b, int max_distance);
int distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist);

hopo_counter new_hopo_counter (int kmer_size);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_hopo_size);
void add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, uint8_t hopo_base_int, int hopo_size);
void copy_hopo_element_start_count (hopo_element *to, hopo_element *from);
void bin_similarity_one (hopo_counter hc, int start, double *result);
void bin_similarity_two (hopo_counter hc1, int start1, hopo_counter hc2, int start2, double *result);
int hopo_counter_hist_integral (hopo_counter hc, int start);
void estimate_coverage_hopo_counter (hopo_counter hc);
void estimate_variance_hopo_counter (hopo_counter hc);

context_histogram_t new_context_histogram_from_hopo_elem (hopo_element he);
void del_context_histogram (context_histogram_t ch);
int  distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int *idx_match);
void context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, int idx_match);

char* context_histogram_tract_as_string (context_histogram_t ch, int kmer_size);
char* context_histogram_generate_name (context_histogram_t ch, int kmer_size);
void genomic_context_find_reference_location (genomic_context_list_t genome, const char *reference_genome_filename);

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
  size_t name_length;

  gzFile fp = gzopen (filename, "r");
  bmc2_kseq_t *seq = bmc2_kseq_init (fp);
  if (!hc) {
    hc_local = new_hopo_counter (kmer_size);
    name_length = strlen(filename);
    hc_local->name = (char*) biomcmc_malloc (sizeof (char) * (name_length + 1));
    strncpy (hc_local->name, filename, name_length);
    hc_local->name[name_length] = '\0';
  }
  if (hc_local->idx) biomcmc_error ("This counter has been compared to another; cannot add more reads to it");
  while ((i = bmc2_kseq_read (seq)) >= 0) update_hopo_counter_from_seq (hc_local, seq->seq.s, seq->seq.l, min_hopo_size); 
  bmc2_kseq_destroy(seq); // other kseq_t parameters: seq->name.s, seq->seq.l, seq->qual.l
  gzclose(fp);
  return hc_local;
}

hopo_counter
new_hopo_counter (int kmer_size)
{
  hopo_counter hc = (hopo_counter) biomcmc_malloc (sizeof (struct hopo_counter_struct));
  hc->n_alloc = 32;
  hc->kmer_size = kmer_size;
  hc->n_idx = hc->n_elem = 0;
  hc->variance[0] = hc->variance[1] = 1.;
  hc->coverage[0] = hc->coverage[1] = 0.;
  hc->elem = (hopo_element*) biomcmc_malloc (hc->n_alloc * sizeof (hopo_element));
  hc->ref_counter = 1;
  hc->name = NULL; hc->idx = NULL;
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
update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_hopo_size)
{
  int i, j, k, count_same = 0, start_mono = -1;
  uint8_t context[2 * hc->kmer_size]; 
  char prev_char = '$';
  //printf ("\nDBG:: %s (%d)\n", seq, seq_length); // DEBUG
  count_same = 0; // zero because previous is "$" 
  for (i = 0; i < (seq_length - hc->kmer_size); i++) {
    if (seq[i] == prev_char) {
      count_same++;
      if ((count_same >= min_hopo_size) && (start_mono >= hc->kmer_size)) {
        while (i < (seq_length-hc->kmer_size-1) && (seq[i+1] == prev_char)) { i++; count_same++; }
        //for (j=start_mono; j <= i; j++) printf ("%c", seq[j]); //DEBUG
        if (dna_in_2_bits[(int)prev_char][0] < dna_in_2_bits[(int)prev_char][1]) { // A or C : forward strand
          k = 0;
          for (j = start_mono - hc->kmer_size; j < start_mono; j++) context[k++] = dna_in_2_bits[ (int)seq[j] ][0]; 
          for (j = i + 1; j <= i + hc->kmer_size; j++) context[k++] = dna_in_2_bits[ (int)seq[j] ][0]; 
          add_kmer_to_hopo_counter (hc, context, dna_in_2_bits[(int)prev_char][0], count_same);
        }
        else if (dna_in_2_bits[(int)prev_char][0] > dna_in_2_bits[(int)prev_char][1]) { // T/U or G : reverse strand
          k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
          for (j = i + hc->kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          for (j = start_mono - 1; j >= start_mono - hc->kmer_size; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          add_kmer_to_hopo_counter (hc, context, dna_in_2_bits[(int)prev_char][1], count_same);
        } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
        // printf ("\t : "); for (j=0; j < 2*hc->kmer_size; j++) { printf ("%c", bit_2_dna[context[j]]); } printf ("\n"); // DEBUG
      } // if (count_same>2) [i.e. we found valid homopolymer]
    } else { 
      count_same = 1;
      prev_char = seq[i];
      start_mono = i;
    } // else (seq[i] == prev_char)
  } // for (i in seq[])
}

void
add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, uint8_t hopo_base_int, int hopo_size)
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
  for (i = 1; i < hc->kmer_size; i++) hc->elem[hc->n_elem].context[0] |= (context[i] << (2 * i));
  hc->elem[hc->n_elem].context[1] = context[hc->kmer_size]; // right context
  for (i = 1; i < hc->kmer_size; i++) hc->elem[hc->n_elem].context[1] |= (context[hc->kmer_size + i] << (2 * i));
  hc->n_elem++;
}

void
copy_hopo_element_start_count (hopo_element *to, hopo_element *from)
{
  to->base    =  from->base; 
  to->length  =  from->length;
  to->context[0] = from->context[0];
  to->context[1] = from->context[1];
  to->count = 1;
}

void
finalise_hopo_counter (hopo_counter hc, const char *reference_genome_filename)
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

  estimate_coverage_hopo_counter (hc);
  estimate_variance_hopo_counter (hc); 

  printf ("DBG::start_genomic_context\n");
  genomic_context_list_t genome = new_genomic_context_list (hc, 1, 1);
  finalise_genomic_context_hist (genome, reference_genome_filename);
  for (i=0; i < genome->n_hist; i++) 
    printf ("DBG:: %4d minmax=(%4d, %4d) modal_tract_length= %-4d (with freq=%-3d) total_number_reads= %-6d location= %6d\t%s\n", 
            genome->hist[i]->n_context, genome->hist[i]->h->min, genome->hist[i]->h->max, 
            genome->hist[i]->h->i[0].idx,  genome->hist[i]->h->i[0].freq,  genome->hist[i]->integral, 
            genome->hist[i]->location, genome->hist[i]->name);
  del_genomic_context_list (genome);
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

//  if (!hc1->idx) finalise_hopo_counter (hc1);  // FIXME: must refactor to include reference_genome_filename
//  if (!hc2->idx) finalise_hopo_counter (hc2); 

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

int hopo_counter_hist_integral (hopo_counter hc, int start)
{
  int i, cov = 0;
  for (i = hc->idx[start]; i < hc->idx[start+1]; i++) cov += hc->elem[i].count;
  return cov;
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
  context_histogram_t ch = (context_histogram_t) biomcmc_malloc (sizeof (struct context_histogram_struct));
  ch->context = (uint64_t*) biomcmc_malloc (2 * sizeof (uint64_t));
  ch->n_context = 1; // how many context pairs 
  ch->mode_context_id = 0; // location, in context[], of best context (the one with length of highest frequency)
  ch->location = -1; // bwa location, calculate at the end
  /* since this is first context, it is best context: */
  ch->base = he.base;
  ch->mode_context_count = he.count; // frequency of best context
  ch->mode_context_length = he.length; // tract length of best context
  ch->context[0] = he.context[0];
  ch->context[1] = he.context[1];
  ch->integral = he.count;
  ch->name = NULL; // defined on finalise()
  ch->h= NULL; // empfreq created at the end (finalise_genomic_context)
  ch->n_tmp = 1;
  ch->tmp_count  = (int*) biomcmc_malloc (sizeof (int));
  ch->tmp_length = (int*) biomcmc_malloc (sizeof (int));
  ch->tmp_count[0]  = he.count;
  ch->tmp_length[0] = he.length;
  return ch;
}

void
del_context_histogram (context_histogram_t ch)
{
  if (!ch) return;
  if (ch->context) free (ch->context);
  if (ch->name)    free (ch->name);
  if (ch->tmp_count)  free (ch->tmp_count);
  if (ch->tmp_length) free (ch->tmp_length);
  del_empfreq (ch->h);
  free (ch);
}

int
distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int *idx_match)
{ 
  int distance = 0, this_max = 0, i;
  *idx_match = -1;
  if ((ch->base - he.base) != 0) return 2 * max_distance + 1; // homopolymer tracts are different
  for (i = 0; i < ch->n_context; i++) {
    distance = distance_between_context_kmer (&(ch->context[2*i]), &(he.context[0]), 2 * max_distance);
    if (distance >= 2 * max_distance) return distance;
    distance += distance_between_context_kmer (&(ch->context[2*i + 1]), &(he.context[1]), 2 * max_distance - distance);
    if (distance >= 2 * max_distance) return distance;
    if (distance > this_max) this_max = distance;
    if (distance == 0) {
      *idx_match = i;
      return 0;
    }
  }
  return this_max;
}

void
context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, int idx_match)
{
  /* Check if context is already on ch -- currently done at distance_between_context_()  */
  // for (i = 0; i < ch->n_context; i++) if ((he.context[0] == ch->context[2*i]) && (he.context[1] == ch->context[2*i+1])) {this_id = i; break; } 
  if (idx_match < 0) { // new context, but still within distance boundary
    // printf ("DBG::realloc::%6d\n", ch->n_context); // short kmer and big files cause realloc corruption
    ch->context = (uint64_t*) biomcmc_realloc ((uint64_t*) ch->context, 2 * (ch->n_context + 1) * sizeof (uint64_t));
    idx_match = ch->n_context++;
    ch->context[2 * idx_match]     = he.context[0];
    ch->context[2 * idx_match + 1] = he.context[1];
  }
  /* if context+tract more frequent than observed so far, then this context is best */
  if (ch->mode_context_count < he.count) {
    ch->mode_context_count = he.count;
    ch->mode_context_length = he.length; // may be same length, diff context
    ch->mode_context_id = idx_match; 
  }
  ch->integral += he.count;

  ch->tmp_count  = (int*) biomcmc_realloc ((int*) ch->tmp_count,  (ch->n_tmp +1) * sizeof (int));
  ch->tmp_length = (int*) biomcmc_realloc ((int*) ch->tmp_length, (ch->n_tmp +1) * sizeof (int));
  ch->tmp_count [ch->n_tmp  ] = he.count;
  ch->tmp_length[ch->n_tmp++] = he.length;
}

genomic_context_list_t
new_genomic_context_list (hopo_counter hc, int max_distance_per_flank, int min_coverage)
{
  int i1, i2, idx_match, j, distance, best_dist=0xffff, best_h=-1, coverage;
  genomic_context_list_t genome = (genomic_context_list_t) biomcmc_malloc (sizeof (struct genomic_context_list_struct));
  genome->hist = NULL;
  genome->n_hist = 0;
  genome->kmer_size = hc->kmer_size; // size on each flanking region (i.e. < 32)
  if (min_coverage < 2) min_coverage = 2;

  for (i1 = 0; i1 < hc->n_idx - 1; i1++) {
    coverage = hopo_counter_hist_integral (hc, i1);
    if (coverage < min_coverage) continue; // too few reads, skip this context+homopolymer
    idx_match = -1; // will become context element of j if contexts match exactly 
    best_dist=0xffff; distance = 0xffff; best_h=-1;
    for (j = 0; (j < genome->n_hist) && (idx_match < 0); j++) {
      distance = distance_between_context_histogram_and_hopo_context (genome->hist[j], hc->elem[hc->idx[i1]], max_distance_per_flank, &idx_match);
      if (distance < best_dist) { best_dist = distance; best_h = j; } // even if distance = 0 we need best_j 
    }
    if (distance >= 2 * max_distance_per_flank) { // no similar context found (or first context ever)
      genome->hist = (context_histogram_t*) biomcmc_realloc ((context_histogram_t*) genome->hist, (genome->n_hist+1) * sizeof (context_histogram_t));
      genome->hist[genome->n_hist] = new_context_histogram_from_hopo_elem (hc->elem[hc->idx[i1]]);
      idx_match = 0; // he->idx gives list with same context so we dont need to compare until idx[i1+1]
      best_h = genome->n_hist++; // outside if/else we need both idx_match, best_h, and i2
      i2 = hc->idx[i1] + 1; // skip first element, which was used to create context_histogram_t 
    }
    else i2 = hc->idx[i1]; // do not skip first element, add it to context_histogram_t in for() loop below 

    for (; i2 < hc->idx[i1+1]; i2++) context_histogram_add_hopo_elem (genome->hist[best_h], hc->elem[i2], idx_match);
  }
  return genome;
}

void
del_genomic_context_list (genomic_context_list_t genome)
{
  if (!genome) return;
  if (genome->hist) {
    for (int i = genome->n_hist-1; i >=0; i--) del_context_histogram (genome->hist[i]);
    free (genome->hist);
  }
  free (genome);
}

char*
context_histogram_tract_as_string (context_histogram_t ch, int kmer_size)
{
  int i, j=0, length = 2 * kmer_size + ch->mode_context_length;
  char *s = (char*) biomcmc_malloc (sizeof (char) * (length + 1));
  uint64_t ctx = ch->context[2 * ch->mode_context_id]; 

  for (i = 0; i < kmer_size; i++) s[i] = bit_2_dna[ (ctx >> (2 * i)) & 3 ]; // left context
  for (; i < kmer_size + ch->mode_context_length; i++) s[i] = bit_2_dna[ch->base]; // homopolymer tract
  ctx = ch->context[2 * ch->mode_context_id + 1]; // right context
  for (j = 0; i < 2 * kmer_size + ch->mode_context_length; j++, i++) s[i] =  bit_2_dna[ (ctx >> (2 * j)) & 3 ]; // i,j not a typo 
  s[i] = '\0';
  return s;
}

char*
context_histogram_generate_name (context_histogram_t ch, int kmer_size)
{
  int i, j=0, length = 2 * kmer_size + 4;
  char *s = (char*) biomcmc_malloc (sizeof (char) * length);
  uint64_t ctx = ch->context[2 * ch->mode_context_id]; 

  for (i = 0; i < kmer_size; i++) s[i] = bit_2_dna[ (ctx >> (2 * i)) & 3 ]; // left context
  s[i++] = '-'; s[i++] = bit_2_dna[ch->base]; s[i++] = '-';// homopolymer tract represented as "-A-" or "-T-"
  ctx = ch->context[2 * ch->mode_context_id + 1]; // right context
  for (j = 0; j < kmer_size; j++, i++) s[i] =  bit_2_dna[ (ctx >> (2 * j)) & 3 ]; 
  s[i] = '\0';
  return s;
}

void
genomic_context_find_reference_location (genomic_context_list_t genome) // new version, but OBSOLETE (now it's done with hopo_counter)
{
  char_vector readseqs;
  char *read;
  int i, j, n_features;
  uint32_t *refseq_offset;
  gff3_fields *features;
  bwase_match_t match;
  bwase_options_t bopt = new_bwase_options_t (0);

  readseqs = new_char_vector (genome->n_hist);

  for (i = 0; i < genome->n_hist; i++) {
    read = context_histogram_tract_as_string (genome->hist[i], genome->opt.kmer_size);
    char_vector_link_string_at_position (readseqs, read, readseqs->next_avail); // just link 
    genome->hist[i]->name =  context_histogram_generate_name (genome->hist[i], genome->opt.kmer_size);
  }
  match = new_bwase_match_from_bwa_and_char_vector (genome->opt.reference_fasta_filename, readseqs, readseqs, 1, bopt);
  del_char_vector(readseqs);

  /* use info from BWA's index for each reference sequence (contig/genome) */
  refseq_offset = (uint32_t*) biomcmc_malloc (sizeof (uint32_t) * match->bns->n_seqs);
  refseq_offset[0] = 0; // equiv to concatenating ref contigs (one-dim "location")
  for (i = 1; i < match->bns->n_seqs; i++) refseq_offset[i] = refseq_offset[i-1] + match->bns->anns[i-1].len;

  for (i = 0; i < match->n_m; i++) { // FIXME: we assume only one bwa hit per histogram but paralogs exist; break ties by leftmost location?
    genome->hist[match->m[i].query_id]->loc2d[0] = match->m[i].ref_id;   // genome/contig ID 
    genome->hist[match->m[i].query_id]->loc2d[1] = match->m[i].position; // location within genome/contig
    genome->hist[match->m[i].query_id]->location = refseq_offset[match->m[i].ref_id] + match->m[i].position; // one-dimensional (flat) index 
    // find_gff3_fields_within_position (gff struct, name of reference genome exactly as in gff, position within ref genome, number of features found)
    features = find_gff3_fields_within_position (genome->opt.gff, match->bns->anns[match->m[i].ref_id].name, match->m[i].position, &n_features);
    for (j = 0; j < n_features; j++) if (features[j].type.id != GFF3_TYPE_region) {
      genome->hist[match->m[i].query_id]->gffeature = features[j]; // store any feature (may be a gene, cds, ...); however ... 
      //printf ("DBG::%6d %6d | (%6d-%6d:%6d)  %s [%d]\n", i, j, features[j].start, features[j].end, match->m[i].position, features[j].attr_id.str, features[j].type.id);
      if (features[j].type.id == GFF3_TYPE_cds) {j = n_features; break; } // .. CDS have priority (overwrite previous and leave for loop)
    }
    if (features) free (features); 
  }
  del_bwase_match_t (match);
  if (refseq_offset) free (refseq_offset);
}

void
genomic_context_find_reference_location (genomic_context_list_t genome, const char *ref_genome_filename) // FIRST VERSION
{
  char_vector readname, readseqs;
  char *read;
  int i, j1, j2, n_matches, *match_list = NULL;

  readname = new_char_vector (genome->n_hist);
  readseqs = new_char_vector (genome->n_hist);

  for (i = 0; i < genome->n_hist; i++) {
    read = context_histogram_tract_as_string (genome->hist[i], genome->kmer_size);
    char_vector_link_string_at_position (readseqs, read, readseqs->next_avail); // just link 
    genome->hist[i]->name =  context_histogram_generate_name (genome->hist[i], genome->kmer_size);
    char_vector_add_string (readname, genome->hist[i]->name); // unlike _link_ above, this will alloc memory and copy name[]  
  }
  // rightmost zero means to store in int[], instead of plotting SAM to stdout
  n_matches = bwa_aln_bwase (ref_genome_filename, readname->string, readseqs->string, NULL, readseqs->nchars, genome->n_hist, 1, &match_list, 0);
  del_char_vector(readname);
  del_char_vector(readseqs);

  for (i = 0; i < n_matches; i++) {
    j1 = match_list[5 * i]; // read index
    j2 = match_list[5 * i + 2]; // location in reference
    genome->hist[j1]->location = j2;
  }
  if (match_list) free (match_list);
}

void
finalise_genomic_context_hist (genomic_context_list_t genome,  const char *reference_genome_filename)
{  
  int i;
  context_histogram_t ch;
  /* 1. empirical frequency histogram of tract lengths */
  for (i = 0; i < genome->n_hist; i++) {
    ch = genome->hist[i];
    ch->h = new_empfreq_from_int_weighted (ch->tmp_length, ch->n_tmp, ch->tmp_count); // histogram, from high to low count
    if (ch->tmp_length) free (ch->tmp_length);
    if (ch->tmp_count) free (ch->tmp_count);
    ch->tmp_length = ch->tmp_count = NULL;
  }
  /* find reference location for each context */
  genomic_context_find_reference_location (genome, reference_genome_filename);
}

typedef struct {context_histogram_t ch; int location; int integral; } cont_hist_ptr_t;
int
compare_cont_hist_ptr_by_location (const void *a, const void *b) // increasing, resolve ties with integral (e.g. location==-1)
{
  int result = (int)(((cont_hist_ptr_t *)a)->location - ((cont_hist_ptr_t *)b)->location); // increasing
  if (result) return result;
  return (int)(((ch_ptr *)b)->integral - ((ch_ptr *)a)->integral);  // decreasing
}
void
genomic_context_sort_context_histogram (genomic_context_list_t genome)
{
  int i; cont_hist_ptr_t *chp;
  chp = (cont_hist_ptr_t*) biomcmc_malloc (genome->n_hist * sizeof (cont_hist_ptr_t));
  for (i = 0; i < genome->n_hist; i++) chp[i] = {.ch = genome->hist[i], .location = genome->hist[i]->location, .integral = genome->hist[i]->integral};
}


distance_generator
new_distance_generator_from_hopo_set (hopo_set hs)
{
  distance_generator dg = new_distance_generator (hs->n_hc, 2); // 2 distances available (2 rescalings)
  distance_generator_set_function_data (dg, &hopo_set_distance_wrapper, (void*) hs); // hs will have any extra data needed to calculate distances
  hs->generator = dg; dg->ref_counter++;
  return dg;
}

void
hopo_set_distance_wrapper (void *data, int s1, int s2, double *result)
{
  clock_t time0, time1;
  time0 = clock ();
  compare_hopo_counters ( ((hopo_set)data)->hc[s1], ((hopo_set)data)->hc[s2], result);
  time1 = clock (); ((hopo_set)data)->secs_comparison += (double)(time1-time0)/(double)(CLOCKS_PER_SEC);
  return;
}

void // 2020.08.24 (with empirical freq)
update_g_tract_summary_from_context_histogram (g_tract_vector_t tract, int prev, int curr, int lev_distance, int n_genome)
{
  int *modlen=NULL, *index=NULL;
  int i, j, k, this = tract->n_summary;
  double result[2];
  tract->n_summary++;
  tract->summary = (g_tract_t*) biomcmc_realloc ((g_tract_t*) tract->summary, tract->n_summary * sizeof (g_tract_t));
  tract->summary[this].location = tract->concat[prev]->location;
  tract->summary[this].d1 = NULL;
  tract->summary[this].lev_distance = lev_distance;
  tract->summary[this].id_in_concat = prev;
  tract->summary[this].n_g = n_genome;

  tract->summary[this].example = tract->concat[prev]; // example of a context_histogram_t
  tract->summary[this].example->ref_counter++;

  tract->summary[this].tab0 = (double*) biomcmc_malloc (n_genome * N_SUMMARY_TABLES * sizeof (double));
  for (i = 0; i < n_genome * N_SUMMARY_TABLES; i++) tract->summary[this].tab0[i] = -FLT_MAX; // smallest negative value
  for (i = 0; i < N_SUMMARY_TABLES; i++) tract->summary[this].gentab[i] = tract->summary[this].tab0 + (n_genome * i); // pointers
  fill_g_tract_summary_tables (tract->summary[this].gentab, tract->concat, prev, curr, n_genome);

  /* 1. empirical frequency of modal lengths */
  modlen = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  index  = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  for (i = 0, j = prev; j < curr; j++, i++) {
    modlen[i]  = tract->concat[j]->h->i[0].idx; // modal tract length for this genome id
    index[i] = tract->concat[j]->index; // genome id (index), not very useful now :)
    //for(k=0;k<tract->concat[j]->h->n;k++) { printf (" %3d (%4d)", tract->concat[j]->h->i[k].idx, tract->concat[j]->h->i[k].freq);} printf ("::%5d::DBG\n", prev);
  }
  tract->summary[this].mode = new_empfreq_from_int_weighted (index, (curr - prev), modlen);
  if (modlen) free (modlen);
  if (index) free (index);

  /* 2. matrix of pairwise distances (1D vector sorted from highest to smallest distance) */
  tract->summary[this].n_dist = ((curr-prev) * (curr-prev-1))/2;
  if (!tract->summary[this].n_dist) return;
  tract->summary[this].d1 = (double*) biomcmc_malloc (2 * tract->summary[this].n_dist * sizeof (double)); // d1 and d2
  tract->summary[this].d2 = tract->summary[this].d1 + tract->summary[this].n_dist;

  for (k = 0, i = prev+1; i < curr; i++) for (j = prev; j < i; j++) {
    distance_between_context_histograms (tract->concat[i], tract->concat[j], result);
    tract->summary[this].d1[k]   = result[0];
    tract->summary[this].d2[k++] = result[1];
  }
  qsort (tract->summary[this].d1, tract->summary[this].n_dist, sizeof (double), compare_double_decreasing);
  qsort (tract->summary[this].d2, tract->summary[this].n_dist, sizeof (double), compare_double_decreasing);
  return;
}

void
create_tract_in_reference_structure (genome_set_t g)
{
  int i, tid;
  /* 1. create a list of reference genome names in same order as bwa's index */
  bwase_match_t match = new_bwase_match_t(g->genome[0]->opt.reference_fasta_filename); // obtain ref seq names
  g->ref_names = new_char_vector (match->bns->n_seqs);
  for (i = 0; i < match->bns->n_seqs; i++) char_vector_add_string (g->ref_names, match->bns->anns[i].name);
  del_bwase_match_t (match);

  /* 2. initialise tract_ref[] with reference genomes tract lengths etc. */
  g->n_tract_ref = g->tract->concat[g->tract->n_concat-1]->tract_id + 1; // if last tract_id is 2 then there are 3 distinct tracts
  g->tract_ref = (tract_in_reference_s*) biomcmc_malloc (g->n_tract_ref * sizeof (tract_in_reference_s));
  for (i = 0; i < g->n_tract_ref; i++) { 
    g->tract_ref[i].tract_length = -1;
    g->tract_ref[i].contig_name = g->tract_ref[i].seq = g->tract_ref[i].tract_name = NULL;
  }

  /* 3. store information from context_histogram concat[] which will be used in copying reference sequence segments */
  for (i = 0; i < g->tract->n_concat; i++) {
    tid = g->tract->concat[i]->tract_id;
    if (!g->tract_ref[tid].contig_name) { // first time tid is met
      g->tract_ref[tid].contig_name = g->ref_names->string[ g->tract->concat[i]->loc2d[0] ];
      g->tract_ref[tid].contig_location = g->tract->concat[i]->loc2d[1];
      g->tract_ref[tid].max_length = g->tract->concat[i]->mode_context_length;
      g->tract_ref[tid].first_idx = i; 
    }
    else { // tid is present more than once: sequence in reference will be longest 
      if (g->tract_ref[tid].contig_location > g->tract->concat[i]->loc2d[1]) 
        g->tract_ref[tid].contig_location = g->tract->concat[i]->loc2d[1]; // leftmost location 
      if (g->tract_ref[tid].max_length < g->tract->concat[i]->mode_context_length) 
        g->tract_ref[tid].max_length = g->tract->concat[i]->mode_context_length;
    }
  }

  /* 4. read fasta file with references and copy tract info from it */
  alignment aln = read_fasta_alignment_from_file (g->genome[0]->opt.reference_fasta_filename, false); 
 // printf ("DEBUG::%s\n", aln->character->string[0]);
  size_t len;
  int kmer_size = g->genome[0]->opt.kmer_size;
  int min_tract_size = g->genome[0]->opt.min_tract_size - 1; // less stringent for ref, and also allows for extra context size
  int start_location = 0;

  for (tid = 0; tid < g->n_tract_ref; tid++) {
    i = lookup_hashtable (aln->taxlabel_hash, g->tract_ref[tid].contig_name);
    if (i < 0) i = lookup_bruteforce (aln->taxlabel, g->tract_ref[tid].contig_name);
    if (i < 0) biomcmc_error ("Contig/genome sequence %s not found in fasta file", g->tract_ref[tid].contig_name);
    start_location = g->tract_ref[tid].contig_location - min_tract_size; // left shift (allow for mismatches) must be smaller than min tract length (to avoid finding spurious)
    if (start_location < 0) start_location = 0;
    len = g->tract_ref[tid].max_length + 2 * g->genome[0]->opt.kmer_size + 2 * min_tract_size;
    if (len > aln->character->nchars[i]) len = aln->character->nchars[i];
    //g->tract_ref[tid].seq = (char*) biomcmc_malloc (sizeof (char) * len + 1);   // NOT NEEDED
    //strncpy (g->tract_ref[tid].seq, aln->character->string[i] + start_location, len);
    //g->tract_ref[tid].seq[len] = '\0';
    /* 4.1 create context+hopo name for reference as in samples */
    //g->tract_ref[tid].tract_name = leftmost_hopo_name_and_length_from_string (g->tract_ref[tid].seq, len+1, kmer_size, min_tract_size, &g->tract_ref[tid].tract_length);
    g->tract_ref[tid].tract_name = leftmost_hopo_name_and_length_from_string (aln->character->string[i] + start_location, len, kmer_size, min_tract_size, &g->tract_ref[tid].tract_length);
    if (!g->tract_ref[tid].tract_length) { // homopolymer not found (e.g. AAAAA in sample is AATAA in reference)
      g->tract_ref[tid].tract_name = (char*) biomcmc_malloc (sizeof (char) * 4);
      strcpy (g->tract_ref[tid].tract_name, "---");
    }
  }
  del_alignment (aln);
}
