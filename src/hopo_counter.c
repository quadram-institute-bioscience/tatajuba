/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "hopo_counter.h"
#include "kseq.h"

BMC2_KSEQ_INIT(gzFile, gzread);

uint8_t dna_in_2_bits[256][2] = {{0xff}};
char bit_2_dna[] = {'A', 'C', 'G', 'T'}; // {00, 01, 10, 11}
// vector from iqtree
//char symbols_protein[] = "ARNDCQEGHILKMFPSTWYVX*"; 
// Base1:               AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
// Base2:               AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
// Base3:               ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
char genetic_code[]  = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLFX"; // X added 


static void initialize_dna_to_bit_tables (void);

void add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, uint8_t hopo_base_int, int hopo_size, int offset, uint8_t reverse_forward_flag);
void copy_hopo_element_start_count_at (hopo_element *to, hopo_element *from, int count);
void copy_hopo_element_locations (hopo_element *to, hopo_element *from);
void estimate_coverage_hopo_counter (hopo_counter hc);
void find_reference_location_and_sort_hopo_counter (hopo_counter hc);

int
compare_hopo_element_decreasing (const void *a, const void *b)
{
  int result = ((hopo_element *)b)->base - ((hopo_element *)a)->base;
  if (result) return result; // sort first by homopolymer (A/T or G/C)
  if (((hopo_element *)b)->context[0] > ((hopo_element *)a)->context[0]) return 1;// sort by left kmers
  if (((hopo_element *)b)->context[0] < ((hopo_element *)a)->context[0]) return -1; // unsigned is never negative!
  if (((hopo_element *)b)->context[1] > ((hopo_element *)a)->context[1]) return 1;// sort by right kmers
  if (((hopo_element *)b)->context[1] < ((hopo_element *)a)->context[1]) return -1;
  return ((hopo_element *)b)->length - ((hopo_element *)a)->length; // same context and homopolymer base, thus sort by tract length
}

int
compare_hopo_element_location (const void *a, const void *b)
{
  int result = ((hopo_element *)a)->read_offset - ((hopo_element *)b)->read_offset;
  if (result) return result; 
  return ((hopo_element *)b)->length - ((hopo_element *)a)->length; // same context and homopolymer base, thus sort by tract length
}

int
compare_hopo_context (hopo_element a, hopo_element b)
{
  int result = b.base - a.base;
  if (result) return result;
  if (b.context[0] > a.context[0]) return 1;
  if (b.context[0] < a.context[0]) return -1;
  if (b.context[1] > a.context[1]) return 1;
  if (b.context[1] < a.context[1]) return -1;
  return 0; 
}

/*! \brief max_dist must be positive, and is the max allowed distance _per_ flanking region */
int
distance_between_single_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist)
{
  uint64_t d = *c1 ^ *c2; // XOR is one if bits are different, zero o.w. 
  int dist = 0;
  while (d && (dist < max_dist)) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference  
  return dist;
}

int // pointer *c1 since we look at both flanking contexts
distance_between_context_kmer_pair (uint64_t *c1, uint64_t *c2)
{
  uint64_t d = c1[0] ^ c2[0]; // XOR is one if bits are different, zero o.w. 
  int dist = 0;
  while (d) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference  
  d = c1[1] ^ c2[1]; 
  while (d) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference  
  return dist;
}

int // pointer *c1 since we look at both flanking contexts
distance_between_context_kmer_pair_with_edit_shift (uint64_t *c1, uint64_t *c2, int *best_shift)
{
  uint64_t d = c1[0] ^ c2[0];  
  uint64_t shift[][3] = {{0,0,0},{0,2,2},{0,4,4},{0,6,6},{2,0,2},{4,0,4},{6,0,6}}; // shift of c1, c2, and mask to be applied (max shift); notice two bits each shift
  int i, best1 = 0xffffff, best2 = 0xffffff, dist = 0, n_shift = 7;

  for (i = 0; (i < n_shift) && (best1 > 0); i++) { // left context
    dist  = (int) shift[i][3]/2; // edit cost of shift
    d = ((c1[0] >> shift[i][0]) ^ (c2[0] >> shift[i][1])) & (~0ULL >> shift[i][3]); // XOR is one if bits are different, zero o.w. 
    while (d) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference
    if (best1 > dist) {
      best1 = dist;
      if (best_shift) {
        best_shift[0] = (int) shift[i][0]/2; // store best shifts for use by VCF for instance
        best_shift[1] = (int) shift[i][1]/2; // shift is in bits (2 bits = one base), but best_shift is in chars (bases)
      }
    }
  }
  for (i = 0; (i < n_shift) && (best2 > 0); i++) { // right context
    dist  = (int) shift[i][3]/2; // edit cost of shift
    d = ((c1[1] >> shift[i][0]) ^ (c2[1] >> shift[i][1])) & (~0ULL >> shift[i][3]); // XOR is one if bits are different, zero o.w. 
    while (d) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference
    if (best2 > dist) {
      best2 = dist;
      if (best_shift) {
        best_shift[2] = (int) shift[i][0]/2; // store best shifts for use by VCF for insance
        best_shift[3] = (int) shift[i][1]/2; // shift is in bits (2 bits = one base), but best_shift is in chars (bases)
      }
    }
  }
  return best1 + best2;
}

void
print_tatajuba_options (tatajuba_options_t opt)
{
  biomcmc_fprintf_colour (stderr, 0, 2, PACKAGE_STRING, "\n");
  fprintf (stderr, "Reference genome fasta file: %s\n", opt.reference_fasta_filename);
  fprintf (stderr, "Reference GFF3 file prefix:  %s\n", opt.gff->file_basename);
  fprintf (stderr, "Output directory:            %s\n", opt.outdir);
  fprintf (stderr, "Number of samples:           %5d (%s)\n", opt.n_samples, (opt.paired_end?"paired-end":"single-end"));
  fprintf (stderr, "Max distance per flanking k-mer:  %6d\n", opt.max_distance_per_flank);
  fprintf (stderr, "Levenshtein distance for merging: %6d\n", opt.levenshtein_distance);
  fprintf (stderr, "Flanking k-mer size (context):    %6d\n", opt.kmer_size);
  fprintf (stderr, "Min tract length to consider:     %6d\n", opt.min_tract_size);
  fprintf (stderr, "Min depth of tract lengths:       %6d\n", opt.min_coverage);
  fprintf (stderr, "Remove biased tracts:             %s\n", (opt.remove_biased?"yes":"no"));
  if (opt.n_threads) fprintf (stderr, "Number of threads (requested or optimised): %3d\n", opt.n_threads);
  else fprintf (stderr, "Software compiled without multithreaded support\n");
  if (opt.paired_end) fprintf (stderr, "Assuming paired-end samples: their file names should be consecutive (no file name check is conducted)\n");
  else fprintf (stderr, "Assuming single-end samples\n");
}

hopo_counter
new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, tatajuba_options_t opt)
{
  int i;
  hopo_counter hc_local = hc;
  size_t name_length;

  gzFile fp = gzopen (filename, "r");
  bmc2_kseq_t *seq = bmc2_kseq_init (fp);
  if (!hc) {
    hc_local = new_hopo_counter (opt.kmer_size);
    name_length = strlen(filename);
    hc_local->name = (char*) biomcmc_malloc (sizeof (char) * (name_length + 1));
    strncpy (hc_local->name, filename, name_length);
    hc_local->name[name_length] = '\0';
    hc_local->opt = opt;
  }
  if (hc_local->idx_initial) biomcmc_error ("This counter has been compared to another; cannot add more reads to it");
  while ((i = bmc2_kseq_read (seq)) >= 0) update_hopo_counter_from_seq (hc_local, seq->seq.s, seq->seq.l, opt.min_tract_size); 
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
  hc->coverage = 0.;
  hc->elem = (hopo_element*) biomcmc_malloc (hc->n_alloc * sizeof (hopo_element));
  hc->ref_counter = 1;
  hc->name = NULL; hc->idx_initial = hc->idx_final = NULL;
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
  if (hc->idx_initial) free (hc->idx_initial);
  if (hc->idx_final) free (hc->idx_final);
  free (hc);
  return;
}

char*
leftmost_hopo_name_and_length_from_string (char *seq, size_t len, int kmer_size, int min_tract_size, int *tract_length)
{
  hopo_counter hc = new_hopo_counter (kmer_size);
  update_hopo_counter_from_seq (hc, seq, (int) len, min_tract_size);
  if (!hc->n_elem) {
    //biomcmc_warning ("No homopolymer was found on sequence (from reference): %s", seq);
    del_hopo_counter (hc);
    *tract_length = 0;
    return NULL;
  }
  *tract_length = hc->elem[0].length;
  char *s = generate_name_from_flanking_contexts (hc->elem[0].context, hc->elem[0].base, kmer_size, false); // false since ref genome is positive strand by design
  del_hopo_counter (hc);
  return s;
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
update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_tract_size)
{
  int i, j, k, count_same = 0, start_mono = -1, offset = -1; // offset is leftmost in seq, not in canonic (used only when seq is from fasta ref) 
  uint8_t context[2 * hc->kmer_size], hopo_base_int, reverse_forward_flag = 0; 
  char prev_char = '$';
  count_same = 0; // zero because previous is "$" 
  for (i = 0; i < (seq_length - hc->kmer_size); i++) {
    if (seq[i] == prev_char) {
      count_same++;
      if ((count_same >= min_tract_size) && (start_mono >= hc->kmer_size)) {
//        while (i < (seq_length-hc->kmer_size-1) && (seq[i+1] == prev_char)) { i++; count_same++; }
        while (i < (seq_length - 1) && (seq[i+1] == prev_char)) { i++; count_same++; }
        if (i >= seq_length - hc->kmer_size) return; // premature end, since no space left for right kmer
        if (dna_in_2_bits[(int)prev_char][0] < dna_in_2_bits[(int)prev_char][1]) { // A or C : forward strand
          for (k =0, j = start_mono - hc->kmer_size; j < start_mono; k++, j++) context[k] = dna_in_2_bits[ (int)seq[j] ][0]; 
          for (j = i + 1; j <= i + hc->kmer_size; k++, j++) context[k] = dna_in_2_bits[ (int)seq[j] ][0]; 
          hopo_base_int = dna_in_2_bits[(int)prev_char][0];
          reverse_forward_flag = 1;
        }
        else if (dna_in_2_bits[(int)prev_char][0] > dna_in_2_bits[(int)prev_char][1]) { // T/U or G : reverse strand
          k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
          for (j = i + hc->kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          for (j = start_mono - 1; j >= start_mono - hc->kmer_size; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          hopo_base_int = dna_in_2_bits[(int)prev_char][1];
          reverse_forward_flag = 2;
        //offset = i + hc->kmer_size;
        } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
        offset = start_mono - hc->kmer_size; // always leftmost, irrespective of revcomplement
        add_kmer_to_hopo_counter (hc, context, hopo_base_int, count_same, offset, reverse_forward_flag);
        //for (j=0; j < hc->kmer_size; j++) printf ("%c", bit_2_dna[context[j]]); printf ("-%c-", bit_2_dna[hopo_base_int]); //DEBUG
        //for (; j < 2 * hc->kmer_size; j++) printf ("%c", bit_2_dna[context[j]]); //DEBUG
      } // if (count_same>2) [i.e. we found valid homopolymer]
    } else { 
      count_same = 1;
      prev_char = seq[i];
      start_mono = i;
    } // else (seq[i] == prev_char)
  } // for (i in seq[])
}

void
update_hopo_counter_from_seq_all_monomers (hopo_counter hc, char *seq, int seq_length) 
{ // does not count HTs, but stores all contexts (used in case HT is missing from reference)
  int i, j, k;
  uint8_t context[2 * hc->kmer_size], hopo_base_int, reverse_forward_flag = 0; 
  for (i = hc->kmer_size; i < (seq_length - hc->kmer_size); i++) {
    if ((seq[i] != seq[i-1]) && (seq[i] != seq[i+1])) {
      if (dna_in_2_bits[(int)seq[i]][0] < dna_in_2_bits[(int)seq[i]][1]) { // A or C : forward strand
        for (k = 0, j = i - hc->kmer_size; j < i; k++, j++) context[k] = dna_in_2_bits[ (int)seq[j] ][0]; 
        for (j = i + 1; j <= i + hc->kmer_size; k++, j++) context[k] = dna_in_2_bits[ (int)seq[j] ][0]; 
        hopo_base_int = dna_in_2_bits[(int)seq[i]][0];
        reverse_forward_flag = 1;
      }
      else if (dna_in_2_bits[(int)seq[i]][0] > dna_in_2_bits[(int)seq[i]][1]) { // T/U or G : reverse strand
        k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
        for (j = i + hc->kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
        for (j = i - 1; j >= i - hc->kmer_size; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
        hopo_base_int = dna_in_2_bits[(int)seq[i]][1];
        reverse_forward_flag = 2;
      } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
      add_kmer_to_hopo_counter (hc, context, hopo_base_int, 1, i - hc->kmer_size, reverse_forward_flag);
    } // if seq[i] not a polymer (o.w. we use standard function) 
  } // for (i in seq[])
}

void
add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, uint8_t hopo_base_int, int hopo_size, int offset, uint8_t reverse_forward_flag)
{
  int i;
  if (hc->n_elem == hc->n_alloc) {
    hc->n_alloc *= 2;
    hc->elem = (hopo_element*) biomcmc_realloc ((hopo_element*)hc->elem, hc->n_alloc * sizeof (hopo_element));
  }
  hc->elem[hc->n_elem].base   = hopo_base_int; // zero (AT) or one (CG) but always store direction with A or C
  hc->elem[hc->n_elem].length = hopo_size;     // homopolymer track length, in bases
  hc->elem[hc->n_elem].count  = 1;
  hc->elem[hc->n_elem].neg_strand = 0; // unless we know better, we use the canonical notation
  hc->elem[hc->n_elem].loc_ref_id =  hc->elem[hc->n_elem].loc_pos = hc->elem[hc->n_elem].loc_last = -1;
  hc->elem[hc->n_elem].mismatches = 0xffe; // 12 bits 
  hc->elem[hc->n_elem].multi = false; 
  
  hc->elem[hc->n_elem].read_offset  = offset;
  hc->elem[hc->n_elem].canon_flag = reverse_forward_flag;
  hc->elem[hc->n_elem].context[0] = hc->elem[hc->n_elem].context[1] = 0ULL; // left context 
  for (i = 0; i < hc->kmer_size; i++) hc->elem[hc->n_elem].context[0] |= ((context[i] & 3ULL) << (2 * i));
  for (i = 0; i < hc->kmer_size; i++) hc->elem[hc->n_elem].context[1] |= ((context[hc->kmer_size + i] & 3ULL) << (2 * i));
  hc->n_elem++;
}

void
copy_hopo_element_start_count_at (hopo_element *to, hopo_element *from, int count)
{
  to->base        = from->base; 
  to->length      = from->length;
  to->read_offset = from->read_offset;
  to->loc_ref_id  = from->loc_ref_id;
  to->loc_pos     = from->loc_pos; 
  to->loc_last    = from->loc_last; 
  to->mismatches  = from->mismatches;
  to->multi       = from->multi;
  to->neg_strand  = from->neg_strand;
  to->context[0] = from->context[0];
  to->context[1] = from->context[1];
  to->count = count;
  to->canon_flag = from->canon_flag;
}

void
copy_hopo_element_locations (hopo_element *to, hopo_element *from)
{
  to->read_offset = from->read_offset;
  to->loc_ref_id  = from->loc_ref_id;
  to->loc_pos     = from->loc_pos; 
  to->loc_last    = from->loc_last; 
  to->mismatches  = from->mismatches;
  to->multi       = from->multi;
  to->neg_strand  = from->neg_strand;
}

void
finalise_hopo_counter (hopo_counter hc)
{
  hopo_element *pivot, *efreq;
  int i, j, coverage = 0, n1 = 0;

  qsort (hc->elem, hc->n_elem, sizeof (hopo_element), compare_hopo_element_decreasing);
  efreq = (hopo_element*) biomcmc_malloc (hc->n_elem * sizeof (hopo_element));
  copy_hopo_element_start_count_at (&(efreq[0]), &(hc->elem[0]), 1);

  /* 1.  frequency of each polymer, in context, into efreq[] vector of elements  */
  n1 = 0;
  for (i=1; i < hc->n_elem; i++) {
    if (compare_hopo_element_decreasing ((const void*) &(hc->elem[i-1]), (const void*) &(hc->elem[i]))) 
      copy_hopo_element_start_count_at (&(efreq[++n1]), &(hc->elem[i]), 1);
    else { // same context _and_tract length
      efreq[n1].count++;
      efreq[n1].canon_flag |= hc->elem[i].canon_flag; 
    }
  }
  n1++; // n1 is (zero-based) index of existing hopo_element (and below we use it as _number_of_elements)

  /* 2. remove tracts seen only on one strand or seen only once */
  hc->n_elem = n1;
  if (hc->opt.remove_biased) {// exclude tract lengths appearing only in one strand (rev or forward)
    for (i = 0, n1 = 0; i < hc->n_elem; i++) if (efreq[i].canon_flag == 3) copy_hopo_element_start_count_at (&(efreq[n1++]), &(efreq[i]), efreq[i].count);
  }
  else { // even if we don't care about biased tracts, remove context_tract_lengths seen only once (but keep original depth whenever count > 1) 
    for (i = 0, n1 = 0; i < hc->n_elem; i++) if (efreq[i].count > 1) copy_hopo_element_start_count_at (&(efreq[n1++]), &(efreq[i]), efreq[i].count);
  }

  /* 3.  hc->elem[] will now point to efreq above, that is, non-identical reads with depth > 1 */
  pivot = hc->elem;
  hc->n_alloc = hc->n_elem = n1; 
  hc->elem = (hopo_element*) biomcmc_realloc ((hopo_element*) efreq, hc->n_alloc * sizeof (hopo_element));
  free (pivot); // free original elem[] 

  /* 4.  idx_initial[] final[] will have indices of distinct contexts (i.e. histogram of polymer lengths) exceeding coverage dept threshold */
  hc->n_idx = 0;
  hc->idx_initial = (int*) biomcmc_malloc (hc->n_elem * sizeof (int));
  hc->idx_final   = (int*) biomcmc_malloc (hc->n_elem * sizeof (int));
  hc->idx_initial[hc->n_idx] = 0; 
  for (i=1; i < hc->n_elem; i++) if (compare_hopo_context (hc->elem[i-1], hc->elem[i])) {
    coverage = 0;
    for (j = hc->idx_initial[hc->n_idx]; j < i; j++) coverage += hc->elem[j].count;
    if (coverage >= hc->opt.min_coverage) { 
      hc->idx_final[hc->n_idx++] = i;
      hc->idx_initial[hc->n_idx] = i; // n_idx was incremented, we are adding a new line
    }
    else hc->idx_initial[hc->n_idx] = i; 
  }
  coverage = 0; // last block c in "aaa.bbbb.cccc" since last compare_context will be "b.c" 
  for (j = hc->idx_initial[hc->n_idx]; j < i; j++) coverage += hc->elem[j].count;
  if (coverage >= hc->opt.min_coverage) hc->idx_final[hc->n_idx++] = i;

  hc->idx_initial = (int*) biomcmc_realloc ((int*) hc->idx_initial, hc->n_idx * sizeof (int));
  hc->idx_final   = (int*) biomcmc_realloc ((int*) hc->idx_final,   hc->n_idx * sizeof (int));

  estimate_coverage_hopo_counter (hc);
  find_reference_location_and_sort_hopo_counter (hc);
}

/* coverage is estimated through kmer frequency (kmer=each flanking region) */
void
estimate_coverage_hopo_counter (hopo_counter hc)
{
  int i, *kmer, *count;
  empfreq ef;
  kmer  = (int*) biomcmc_malloc (2 * hc->n_elem * sizeof (int));
  count = (int*) biomcmc_malloc (2 * hc->n_elem * sizeof (int));
  for (i=0; i < hc->n_elem; i++) { // several elements may share (integer version of) context
    kmer[ i]  = (int) (hc->elem[i].context[0] & ((1ULL << 31) -1)); // use 30 bits
    count[i] = hc->elem[i].count;
    kmer[ i + hc->n_elem] = (int) (hc->elem[i].context[1] & ((1ULL << 31) -1)); // use 30 bits
    count[i + hc->n_elem] = hc->elem[i].count;
  }
  ef = new_empfreq_from_int_weighted (kmer, 2 * hc->n_elem, count);
  hc->coverage = ef->i[0].freq;
  if (kmer) free (kmer);
  if (count) free (count);
  del_empfreq (ef);
}

int hopo_counter_histogram_integral (hopo_counter hc, int start) // OBSOLETE
{
  int i, cov = 0;
  for (i = hc->idx_initial[start]; i < hc->idx_final[start]; i++) cov += hc->elem[i].count;
  return cov;
}

char*
generate_tract_as_string (uint64_t *context, int8_t base, int kmer_size, int tract_length, bool neg_strand)
{
  int i = 0, j = 0, length = (2 * kmer_size) + tract_length + 1;
  char *s = (char*) biomcmc_malloc (sizeof (char) * length);
  uint64_t ctx = context[0]; 

  if (neg_strand) {
    s[length - 1] = '\0'; // A=00, C=01, G=10, T=11 therefore ~A=T etc. 
    for (j = length-2, i = 0; i < kmer_size; i++, j--) s[j] = bit_2_dna[ (~(ctx >> (2 * i))) & 3ULL ]; // left context at the end;
    for (i = 0; i < tract_length; j--, i++) s[j] = bit_2_dna[(~base) & 3ULL]; // homopolymer tract
    ctx = context[1];
    for (i = 0; i < kmer_size; j--, i++) s[j] =  bit_2_dna[ (~(ctx >> (2 * i))) & 3ULL ]; 
  }
  else {
    for (i = 0; i < kmer_size; i++) s[i] = bit_2_dna[ (ctx >> (2 * i)) & 3ULL ]; // left context
    for (; i < kmer_size + tract_length; i++) s[i] = bit_2_dna[base]; // homopolymer tract
    ctx = context[1]; // right context
    for (j = 0; j < kmer_size; j++, i++) s[i] =  bit_2_dna[ (ctx >> (2 * j)) & 3ULL ]; 
    s[i] = '\0';
  }
  return s;
}

char*
generate_name_from_flanking_contexts (uint64_t *context, int8_t base, int kmer_size, bool neg_strand)
{
  int i, j=0, length = 2 * kmer_size + 4;
  char *s = (char*) biomcmc_malloc (sizeof (char) * length);
  uint64_t ctx = context[0];

  if (neg_strand) {
    s[length - 1] = '\0'; // A=00, C=01, G=10, T=11 therefore ~A=T etc. 
    for (j = length-2, i = 0; i < kmer_size; i++, j--) s[j] = bit_2_dna[ (~(ctx >> (2 * i))) & 3ULL ]; // left context at the end;
    s[j--] = '.'; s[j--] = bit_2_dna[(~base) & 3ULL]; s[j--] = '.';// homopolymer tract represented as "-A-" or "-T-"
    ctx = context[1];
    for (i = 0; i < kmer_size; j--, i++) s[j] =  bit_2_dna[ (~(ctx >> (2 * i))) & 3ULL ]; 
  }
  else {
    for (i = 0; i < kmer_size; i++) s[i] = bit_2_dna[ (ctx >> (2 * i)) & 3ULL ]; // left context
    s[i++] = '.'; s[i++] = bit_2_dna[base]; s[i++] = '.';// homopolymer tract represented as "-A-" or "-T-"
    ctx = context[1];
    for (j = 0; j < kmer_size; j++, i++) s[i] =  bit_2_dna[ (ctx >> (2 * j)) & 3ULL ]; 
    s[i] = '\0';
  }
  return s;
}

void
find_reference_location_and_sort_hopo_counter (hopo_counter hc)
{
  char_vector readseqs;
  char *read;
  int i, j, qid;
  uint32_t *refseq_offset;
  uint16_t mismatch;
  int32_t read_offset; 
  bool skip_match;
  bwase_match_t match;
  bwase_options_t bopt = new_bwase_options_t (0);

  readseqs = new_char_vector (hc->n_idx);

  /* 0.  preprocessing: make sure read_offset (originally with location info within fastq read) can be used as flattened location */
  for (i = 0; i < hc->n_elem; i++) hc->elem[i].read_offset = -1;  
  /* 0.  preprocessing: generating DNA segments with context+tract fror BWA, and dummy sequence names */
  for (i = 0; i < hc->n_idx; i++) { // link_string() assumes char* was allocated outside, but assumes control of if (i.e. will free() it)
    read = generate_tract_as_string (hc->elem[hc->idx_initial[i]].context, 
                                     hc->elem[hc->idx_initial[i]].base,
                                     hc->opt.kmer_size, 
                                     hc->elem[hc->idx_initial[i]].length, // longest
                                     false); // do not flip to negative strand (BWA will tell if flipping is necessary)
    char_vector_link_string_at_position (readseqs, read, readseqs->next_avail); // just link 
  }

  /* 1.  run BWA; notice that first "readseqs" is usually the sequence name but we dont care */ 
  match = new_bwase_match_from_bwa_and_char_vector (hc->opt.reference_fasta_filename, readseqs, readseqs, 1, bopt);
  del_char_vector(readseqs);

  /* 2.  use info from BWA's index for each reference sequence (contig/genome) */
  refseq_offset = (uint32_t*) biomcmc_malloc (sizeof (uint32_t) * match->bns->n_seqs);
  refseq_offset[0] = 0; // equiv to concatenating ref contigs (one-dim "location")
  for (i = 1; i < match->bns->n_seqs; i++) refseq_offset[i] = refseq_offset[i-1] + match->bns->anns[i-1].len;

  /* 3. fill first hc->elem[] from each distinct context+tract with info from BWA (best match, and leftmost) */
  for (i = 0; i < match->n_m; i++) { 
    skip_match = true; // we assume that this match is worse than one already seen
    qid = hc->idx_initial[match->m[i].query_id]; 
    // mismatch = match->m[i].mm + match->m[i].gape + match->m[i].gapo; // mismatches plus extra indels 
    mismatch = match->m[i].nm; //  nm (edit distance)
    read_offset = refseq_offset[match->m[i].ref_id] + match->m[i].position; // one-dimensional (flat) index 

    
    if (hc->elem[qid].loc_pos < 0) skip_match = false; // first time this element is seen 
    if (skip_match && (hc->elem[qid].mismatches > mismatch)) { hc->elem[qid].multi = true; skip_match = false; } // better match than existing 
    if (skip_match && (hc->elem[qid].read_offset > read_offset)) { hc->elem[qid].multi = true; skip_match = false; } // same match but leftmost to existing

    //if (hc->elem[qid].multi) printf ("DEBUG::hopo::%d mism=%d %d offset=%d %d\n", match->m[i].ref_id,hc->elem[qid].mismatches, mismatch, hc->elem[qid].read_offset, read_offset);

    if (!skip_match) {
      hc->elem[qid].mismatches = mismatch;
      hc->elem[qid].loc_ref_id = match->m[i].ref_id;   // genome/contig ID 
      hc->elem[qid].loc_pos    = match->m[i].position; // location within genome/contig where HT starts
      hc->elem[qid].loc_last   = match->m[i].position + match->m[i].ref_length - 1; // location within genome/contig of last position of HT
      hc->elem[qid].read_offset = read_offset;         // one-dimensional (flat) index 
      hc->elem[qid].neg_strand = match->m[i].neg_strand; // if negative strand on reference genome (unrelated to for/rev read strand, BTW) 
    }
  }
  del_bwase_match_t (match);

  /* 4.  copy location info to all elements from same context+tract */
  for (i = 0; i < hc->n_idx; i++) for (j = hc->idx_initial[i] + 1; j < hc->idx_final[i]; j++) 
    copy_hopo_element_locations (&(hc->elem[j]), &(hc->elem[hc->idx_initial[i]]));
  
  /* 5.  sort according to location, with unknown elements first */
  qsort (hc->elem, hc->n_elem, sizeof (hopo_element), compare_hopo_element_location);
  for (i = 0; (i < hc->n_elem) && (hc->elem[i].read_offset < 0); i++); // just scan i
  hc->ref_start = i; // downstream analysis will use only these
  biomcmc_fprintf_colour (stderr, 0,2, hc->name, ": %6d found and %6d context+tracts were not found in reference\n", hc->n_elem - i, i);
  if (i > hc->n_elem/2) biomcmc_warning ("%6d out of %6d (more than half) context+tracts were not found in reference for sample %s", i, hc->n_elem, hc->name);

  if (refseq_offset)   free (refseq_offset);
  if (hc->idx_initial) free (hc->idx_initial);
  if (hc->idx_final)   free (hc->idx_final);
  hc->idx_initial = hc->idx_final = NULL;
}

char*
protein_from_dna_string (char *dna, size_t n_dna, bool reverse)
{
  int i, j, k, codon, n_codon = n_dna/3;
  char *prot = NULL;
  prot = (char*) biomcmc_malloc ((n_codon + 1) * sizeof (char));
  if (dna_in_2_bits[0][0] == 0xff) initialize_dna_to_bit_tables (); // should not happen since new_hopo_counter calls it, but... 

  if (reverse) for (i = codon = 0; i < n_codon; i++) {
    for (j = 0; (j < 3) && ((3*i+j) < (int) n_dna); j++) {
      k = dna_in_2_bits[ (int) dna[n_dna - 1 - (3*i + j)] ][1];
      if (k < 4) codon |= k << 2 * j; 
      else       codon = 64;
    }
    if (codon > 64) codon = 64;
    prot[i] = genetic_code[codon];
  }
  else for (i = codon = 0; i < n_codon; i++) {
    for (j = 0; (j < 3) && ((3*i+j) < (int) n_dna); j++) {
      k = dna_in_2_bits[ (int) dna[3*i + j] ][0];
      if (k < 4) codon |= k << 2 * j; 
      else       codon = 64;
    }
    if (codon > 64) codon = 64;
    prot[i] = genetic_code[codon];
  }
  return prot;
}
