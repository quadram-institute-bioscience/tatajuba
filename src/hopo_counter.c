/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "hopo_counter.h"
#include "kseq.h"

BMC2_KSEQ_INIT(gzFile, gzread);

uint8_t dna_in_2_bits[256][2] = {{0xff}};
char bit_2_dna[] = {'A', 'C', 'G', 'T'};

static void initialize_dna_to_bit_tables (void);

hopo_counter new_hopo_counter (int kmer_size);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_tract_size);
void add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, uint8_t hopo_base_int, int hopo_size);
void copy_hopo_element_start_count_at (hopo_element *to, hopo_element *from, int count);
void estimate_coverage_hopo_counter (hopo_counter hc);

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
distance_between_context_kmer (uint64_t *c1, uint64_t *c2, int max_dist)
{
  uint64_t d = *c1 ^ *c2; // XOR is one if bits are different, zero o.w. 
  int dist = 0;
  while (d && (dist < max_dist)) { if (d & 3) dist++; d >>= 2; } // every two bits, check if there is any difference  
  return dist;
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
  if (hc_local->idx) biomcmc_error ("This counter has been compared to another; cannot add more reads to it");
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
  char *s = generate_name_from_flanking_contexts (hc->elem[0].context, hc->elem[0].base, kmer_size);
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
  int i, j, k, count_same = 0, start_mono = -1;
  uint8_t context[2 * hc->kmer_size], hopo_base_int; 
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
        }
        else if (dna_in_2_bits[(int)prev_char][0] > dna_in_2_bits[(int)prev_char][1]) { // T/U or G : reverse strand
          k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
          for (j = i + hc->kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          for (j = start_mono - 1; j >= start_mono - hc->kmer_size; j--) context[k++] = dna_in_2_bits[ (int)seq[j] ][1]; 
          hopo_base_int = dna_in_2_bits[(int)prev_char][1];
        } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
        add_kmer_to_hopo_counter (hc, context, hopo_base_int, count_same);
/*
        printf ("BDG: ");
        for (j=0; j < hc->kmer_size; j++) printf ("%c", bit_2_dna[context[j]]); //DEBUG
        printf ("-%c-", bit_2_dna[hopo_base_int]); //DEBUG
        for (; j < 2 * hc->kmer_size; j++) printf ("%c", bit_2_dna[context[j]]); //DEBUG
        printf ("\n");
*/
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
  hc->elem[hc->n_elem].context[0] = hc->elem[hc->n_elem].context[1] = 0ULL; // left context 
  for (i = 0; i < hc->kmer_size; i++) hc->elem[hc->n_elem].context[0] |= ((context[i] & 3ULL) << (2 * i));
  for (i = 0; i < hc->kmer_size; i++) hc->elem[hc->n_elem].context[1] |= ((context[hc->kmer_size + i] & 3ULL) << (2 * i));
  hc->n_elem++;
}

void
copy_hopo_element_start_count_at (hopo_element *to, hopo_element *from, int count)
{
  to->base    =  from->base; 
  to->length  =  from->length;
  to->context[0] = from->context[0];
  to->context[1] = from->context[1];
  to->count = count;
}

void
finalise_hopo_counter (hopo_counter hc)
{
  hopo_element *pivot, *efreq;
  int i, n1 = 0;
  qsort (hc->elem, hc->n_elem, sizeof (hopo_element), compare_hopo_element_decreasing);
  efreq = (hopo_element*) biomcmc_malloc (hc->n_elem * sizeof (hopo_element));
  copy_hopo_element_start_count_at (&(efreq[0]), &(hc->elem[0]), 1);
  n1 = 0; hc->n_idx = 0;

  /* frequency of each polymer, in context, into efreq[] vector of elements  */
  for (i=1; i < hc->n_elem; i++) {
    if (compare_hopo_element_decreasing ((const void*) &(hc->elem[i-1]), (const void*) &(hc->elem[i]))) 
      copy_hopo_element_start_count_at (&(efreq[++n1]), &(hc->elem[i]), 1);
    else efreq[n1].count++; // same context _and_tract length
  }
  n1++; // n1 is (zero-based) index of existing hopo_element (and below we use it as _number_of_elements)
  /* remove context_tract_lengths seen only once (but keep original depth whenever count > 1) */

  hc->n_elem = n1;
  for (i = 0, n1 = 0; i < hc->n_elem; i++) if (efreq[i].count > 1) copy_hopo_element_start_count_at (&(efreq[n1++]), &(efreq[i]), efreq[i].count);

  /* new efreq[] becomes hc->elem[] */
  pivot = hc->elem;
  hc->n_alloc = hc->n_elem = n1; 
  hc->elem = (hopo_element*) biomcmc_realloc ((hopo_element*) efreq, hc->n_alloc * sizeof (hopo_element));
  free (pivot);

  /* idx[] will have indices of distinct contexts (i.e. histogram of polymer lengths) */
  hc->idx = (int*) biomcmc_malloc ((hc->n_elem + 1) * sizeof (int));
  hc->idx[hc->n_idx++] = 0; 
  for (i=1; i < hc->n_elem; i++) if (compare_hopo_context (hc->elem[i-1], hc->elem[i])) hc->idx[hc->n_idx++] = i;
  hc->idx[hc->n_idx++] = i; // last element (plus one) is also on index (to avoid conditional checking */
  hc->idx = (int*) biomcmc_realloc ((int*) hc->idx, hc->n_idx* sizeof (int));

  estimate_coverage_hopo_counter (hc);
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

int hopo_counter_histogram_integral (hopo_counter hc, int start)
{
  int i, cov = 0;
  for (i = hc->idx[start]; i < hc->idx[start+1]; i++) cov += hc->elem[i].count;
  return cov;
}

char*
generate_name_from_flanking_contexts (uint64_t *context, int8_t base, int kmer_size)
{
  int i, j=0, length = 2 * kmer_size + 4;
  char *s = (char*) biomcmc_malloc (sizeof (char) * length);
  uint64_t ctx = context[0];

  for (i = 0; i < kmer_size; i++) s[i] = bit_2_dna[ (ctx >> (2 * i)) & 3 ]; // left context
  s[i++] = '-'; s[i++] = bit_2_dna[base]; s[i++] = '-';// homopolymer tract represented as "-A-" or "-T-"
  ctx = context[1];
  for (j = 0; j < kmer_size; j++, i++) s[i] =  bit_2_dna[ (ctx >> (2 * j)) & 3 ]; 
  s[i] = '\0';
  return s;
}
