/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "histogram_comparison.h"
#include "kseq.h"

BMC2_KSEQ_INIT(gzFile, gzread);

uint8_t dna_in_2_bits[256][2] = {{0xff}};
char bit_2_dna[] = {'A', 'C', 'G', 'T'};

static void initialize_dna_to_bit_tables (void);

hopo_counter new_hopo_counter (int kmer_size);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int min_tract_size);
void add_kmer_to_hopo_counter (hopo_counter hc, uint8_t *context, uint8_t hopo_base_int, int hopo_size);
void copy_hopo_element_start_count_at (hopo_element *to, hopo_element *from, int count);
int hopo_counter_histogram_integral (hopo_counter hc, int start);
void estimate_coverage_hopo_counter (hopo_counter hc);
void finalise_hopo_counter (hopo_counter hc);

context_histogram_t new_context_histogram_from_hopo_elem (hopo_element he);
void context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, int idx_match);

char* context_histogram_tract_as_string (context_histogram_t ch, int kmer_size);
char* context_histogram_generate_name (context_histogram_t ch, int kmer_size);
char* generate_name_from_flanking_contexts (uint64_t *context, int8_t base, int kmer_size);
void genomic_context_find_reference_location (genomic_context_list_t genome);
void genomic_context_merge_histograms_at_same_location (genomic_context_list_t genome);
void accumulate_from_context_histogram (context_histogram_t to, context_histogram_t from);

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

context_histogram_t
new_context_histogram_from_hopo_elem (hopo_element he)
{
  context_histogram_t ch = (context_histogram_t) biomcmc_malloc (sizeof (struct context_histogram_struct));
  ch->context = (uint64_t*) biomcmc_malloc (2 * sizeof (uint64_t));
  ch->ref_counter = 1;
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
  ch->index = 1; // index in genome_set (here used provisorily as n_tmp 
  ch->tmp_count  = (int*) biomcmc_malloc (sizeof (int));
  ch->tmp_length = (int*) biomcmc_malloc (sizeof (int));
  ch->tmp_count[0]  = he.count;
  ch->tmp_length[0] = he.length;
  ch->gffeature = return_null_gff3_field(); // zero, but for struct
  ch->tract_id = -1; // will be set only in tract->concat, and used to summarise context histograms
  return ch;
}

void
del_context_histogram (context_histogram_t ch)
{
  if (!ch) return;
  if (--ch->ref_counter) return;
  if (ch->context) free (ch->context);
  if (ch->name)    free (ch->name);
  if (ch->tmp_count)  free (ch->tmp_count);
  if (ch->tmp_length) free (ch->tmp_length);
  del_empfreq (ch->h);
  free (ch);
}

void
context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, int idx_match)
{
  if (idx_match < 0) { // new context, but still within distance boundary (since this function was called...)
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

  ch->tmp_count  = (int*) biomcmc_realloc ((int*) ch->tmp_count,  (ch->index +1) * sizeof (int));
  ch->tmp_length = (int*) biomcmc_realloc ((int*) ch->tmp_length, (ch->index +1) * sizeof (int));
  ch->tmp_count [ch->index]   = he.count;
  ch->tmp_length[ch->index++] = he.length;
}

genomic_context_list_t
new_genomic_context_list (hopo_counter hc)
{
  int i1, i2, idx_match, j, distance, best_dist=0xffff, best_h=-1, read_coverage;
  genomic_context_list_t genome = (genomic_context_list_t) biomcmc_malloc (sizeof (struct genomic_context_list_struct));
  genome->hist = NULL;
  genome->n_hist = 0;
  
  /* 1. sort hopo_counter, calculate coverage etc. Remember that hopo_counter handles only identical contexts */
  finalise_hopo_counter (hc);
  genome->opt = hc->opt;
  genome->coverage = hc->coverage;
  genome->name = hc->name;
  hc->name = NULL;

  int dbg_count = -1;
  /* 2. accumulate histograms of 'equivalent' (almost identical) contexts */
  for (i1 = 0; i1 < hc->n_idx - 1; i1++) {
    read_coverage = hopo_counter_histogram_integral (hc, i1);
    if (read_coverage < genome->opt.min_coverage) continue; // too few reads, skip this context+homopolymer
    dbg_count++;
    idx_match = -1; // will become context element of j if contexts match exactly 
    best_dist=0xffff; distance = 0xffff; best_h=-1;
    for (j = 0; (j < genome->n_hist) && (idx_match < 0); j++) {
      distance = distance_between_context_histogram_and_hopo_context (genome->hist[j], hc->elem[hc->idx[i1]], genome->opt.max_distance_per_flank, &idx_match);
      if (distance < best_dist) { best_dist = distance; best_h = j; } // even if distance = 0 we need best_j 
    }
    if (distance >= 2 * genome->opt.max_distance_per_flank) { // no similar context found (or first context ever)
      genome->hist = (context_histogram_t*) biomcmc_realloc ((context_histogram_t*) genome->hist, (genome->n_hist+1) * sizeof (context_histogram_t));
      genome->hist[genome->n_hist] = new_context_histogram_from_hopo_elem (hc->elem[hc->idx[i1]]);
      idx_match = 0; // he->idx gives list with same context so we dont need to compare until idx[i1+1]
      best_h = genome->n_hist++; // outside if/else we need both idx_match, best_h, and i2
      i2 = hc->idx[i1] + 1; // skip first element, which was used to create context_histogram_t 
    }
    else i2 = hc->idx[i1]; // do not skip first element, add it to context_histogram_t in for() loop below 

    for (; i2 < hc->idx[i1+1]; i2++) context_histogram_add_hopo_elem (genome->hist[best_h], hc->elem[i2], idx_match);
  }
  printf ("DEBUG::hopo_counters=  %12d dbg_count =   %12d  n_hist=  %12d\n", hc->n_idx, dbg_count, genome->n_hist);
  /* 3. generate empfreq histograms, name each histogram and find location in ref genome */
  finalise_genomic_context_hist (genome);
  return genome;
}

void
finalise_genomic_context_hist (genomic_context_list_t genome)
{  
  int i;
  context_histogram_t ch;
  /* 1. empirical frequency histogram of tract lengths */
  for (i = 0; i < genome->n_hist; i++) {
    ch = genome->hist[i];
    ch->h = new_empfreq_from_int_weighted (ch->tmp_length, ch->index, ch->tmp_count); // histogram, from high to low count
    if (ch->tmp_length) free (ch->tmp_length);
    if (ch->tmp_count) free (ch->tmp_count);
    ch->tmp_length = ch->tmp_count = NULL; ch->index = -1;
  }

  /* 2. find reference location for each context */
  genomic_context_find_reference_location (genome);

  /* 3. sort context_histograms based on genomic location, ties broken with more frequent first. Ties are found when 
   *    location == -1 i.e. not found on reference, and ultimately ties are sorted by context */
  qsort (genome->hist, genome->n_hist, sizeof (context_histogram_t), compare_context_histogram_for_qsort);
  for (i = 0; (i < genome->n_hist) && (genome->hist[i]->location < 0); i++); // just scan i
  genome->ref_start = i;
  biomcmc_fprintf_colour (stderr, 0,2, genome->name, ": %6d out of %6d context+tracts were not found in reference\n", i, genome->n_hist);
  if (i > genome->n_hist/2) 
    biomcmc_warning ("%6d out of %6d (more than half) context+tracts were not found in reference for sample %s\n", i, genome->n_hist, genome->name);

  /* 4. merge context_histograms mapped to same ref genome location.BWA may detect that slightly different contexts are 
   *    actually the same, specially when max_flank_distance is too strict */
  genomic_context_merge_histograms_at_same_location (genome);
  //print_debug_genomic_context_hist (genome);
  
  /* 5. add genome-wide information to each histogram (genome coverage and number of histograms), useful to summarise tracts */
  for (i = 0; i < genome->n_hist; i++) { genome->hist[i]->coverage = genome->coverage; genome->hist[i]->n_tracts = genome->n_hist; } 

  return;
}

void
del_genomic_context_list (genomic_context_list_t genome)
{
  if (!genome) return;
  if (genome->hist) {
    for (int i = genome->n_hist-1; i >=0; i--) del_context_histogram (genome->hist[i]);
    free (genome->hist);
  }
  if (genome->name) free (genome->name);
  free (genome);
}

char*
context_histogram_tract_as_string (context_histogram_t ch, int kmer_size)
{
  int i = 0, j = 0, best_length;
  char *s; 
  uint64_t ctx = ch->context[2 * ch->mode_context_id]; 

  best_length = ch->h->i[0].idx; // alternative is ch->mode_context_length
  s = (char*) biomcmc_malloc (sizeof (char) * (2 * kmer_size + best_length + 1));

  for (i = 0; i < kmer_size; i++) s[i] = bit_2_dna[ (ctx >> (2 * i)) & 3 ]; // left context
  for (; i < kmer_size + best_length; i++) s[i] = bit_2_dna[ch->base]; // homopolymer tract
  ctx = ch->context[2 * ch->mode_context_id + 1]; // right context
  for (j = 0; j < kmer_size; j++, i++) s[i] =  bit_2_dna[ (ctx >> (2 * j)) & 3 ]; 
  s[i] = '\0';
  return s;
}

char*
context_histogram_generate_name (context_histogram_t ch, int kmer_size)
{
  return generate_name_from_flanking_contexts (ch->context + (2 * ch->mode_context_id), ch->base, kmer_size); // assuming consecutive context[]
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

char*
context_histogram_generate_name_bkp (context_histogram_t ch, int kmer_size)
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
genomic_context_find_reference_location (genomic_context_list_t genome)
{
  char_vector readname, readseqs;
  char *read;
  int i, j, n_features;
  uint32_t *refseq_offset;
  gff3_fields *features;
  bwase_match_t match;
  bwase_options_t bopt = new_bwase_options_t (0);

  readname = new_char_vector (genome->n_hist);
  readseqs = new_char_vector (genome->n_hist);

  for (i = 0; i < genome->n_hist; i++) {
    read = context_histogram_tract_as_string (genome->hist[i], genome->opt.kmer_size);
    char_vector_link_string_at_position (readseqs, read, readseqs->next_avail); // just link 
    genome->hist[i]->name =  context_histogram_generate_name (genome->hist[i], genome->opt.kmer_size);
    char_vector_add_string (readname, genome->hist[i]->name); // unlike _link_ above, this will alloc memory and copy name[]  
  }
  match = new_bwase_match_from_bwa_and_char_vector (genome->opt.reference_fasta_filename, readname, readseqs, 1, bopt);
  del_char_vector(readname);
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
print_debug_genomic_context_hist (genomic_context_list_t genome)
{
  printf ("DBG::start_genomic_context::%s\n", genome->name);
  for (int i=0; i < genome->n_hist; i++) 
    printf ("DBG::number_contexts= %4d |minmax= %5d %5d |modal_tract_length= %-5d |freq= %-6d |total_number_reads= %-7d |location= %7d\t%s\n", 
            genome->hist[i]->n_context, genome->hist[i]->h->min, genome->hist[i]->h->max, 
            genome->hist[i]->h->i[0].idx,  genome->hist[i]->h->i[0].freq,  genome->hist[i]->integral, 
            genome->hist[i]->location, genome->hist[i]->name);
}

void  // FIXME: CTCT.3xA.CCCC at position 1 and GCTC.4xA.CCCC at position zero are the same (w/ T->A subst)
genomic_context_merge_histograms_at_same_location (genomic_context_list_t genome)
{
  context_histogram_t *new_h;
  int i, j;

  if (genome->ref_start == genome->n_hist) return; // nothing to do if no contexts were found in ref genome 
  new_h = (context_histogram_t*) biomcmc_malloc (genome->n_hist * sizeof (context_histogram_t));

  for (j = 0; j < genome->ref_start; j++) new_h[j] = genome->hist[j]; // negative locations are neglected
  new_h[j] = genome->hist[j]; 
  for (i = j + 1; i < genome->n_hist; i++) {
    // NULL is pointer, it means that we don't want to store the actual levenshtein distance 
    if (context_histograms_overlap (new_h[j], genome->hist[i], NULL, genome->opt)) accumulate_from_context_histogram (new_h[j], genome->hist[i]);
    else new_h[++j] = genome->hist[i];
  }

  if (genome->n_hist > j) {
    genome->n_hist = j;
    new_h = (context_histogram_t*) biomcmc_realloc ((context_histogram_t*)new_h, genome->n_hist * sizeof (context_histogram_t));
  }
  free (genome->hist);
  genome->hist = new_h;
}

void
accumulate_from_context_histogram (context_histogram_t to, context_histogram_t from)
{
  uint64_t *c1; 
  int i, j, n1;
  char *s1;
  /* 1. merge empirical frequencies */
  to->h = new_empfreq_merge_empfreqs (to->h, from->h);
  /* 2. find context_histogram with "best" tract (exact context+tract with highest frequency/depth) */
  if (to->mode_context_count < from->mode_context_count) { // needs to update modal values
    to->location = from->location; // is the same except when one has indels (found with levenshtein distance)
    to->mode_context_count  = from->mode_context_count;
    to->mode_context_length = from->mode_context_length;
    to->mode_context_id     = from->mode_context_id;
    n1 = to->n_context; to->n_context = from->n_context; from->n_context = n1; // swap
    c1 = to->context;   to->context =   from->context;   from->context = c1;
    s1 = to->name;      to->name =      from->name;      from->name = s1;
  }
  /* 3. assume *to has best tract; update it with info from *from */
  to->integral += from->integral;
  n1 = to->n_context; // size before expanding with elements from *from
  for (i = 0; i < from->n_context; i++) for (j = 0; j < n1; j++)
    if ((to->context[2*j] != from->context[2*i]) || (to->context[2*j+1] != from->context[2*i+1])) {
      to->context = (uint64_t*) biomcmc_realloc ((uint64_t*) to->context, 2 * (to->n_context + 1) * sizeof (uint64_t));
      to->context[2 * to->n_context] = from->context[2 * i];
      to->context[2 * to->n_context + 1] = from->context[2 * i + 1];
      to->n_context++;
    } // if contexts are identical, do nothing

  del_context_histogram (from);
}
