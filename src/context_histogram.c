/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "context_histogram.h"

context_histogram_t new_context_histogram_from_hopo_elem (hopo_element he);
void context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, int idx_match);

char* context_histogram_tract_as_string (context_histogram_t ch, int kmer_size);
char* context_histogram_generate_name (context_histogram_t ch, int kmer_size);
void genomic_context_find_reference_location (genomic_context_list_t genome);
void genomic_context_merge_histograms_at_same_location (genomic_context_list_t genome);
void accumulate_from_context_histogram (context_histogram_t to, context_histogram_t from);

int
distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int *idx_match)
{ 
  int distance = 0, this_max = 0, i;
  *idx_match = -1;
  if (ch->base != he.base) return 2 * max_distance + 1; // homopolymer tracts are different
  for (i = 0; i < ch->n_context; i++) {
    distance  = distance_between_single_context_kmer (&(ch->context[2*i]), &(he.context[0]), 2 * max_distance);
    if (distance >= 2 * max_distance) return distance;
    distance += distance_between_single_context_kmer (&(ch->context[2*i + 1]), &(he.context[1]), 2 * max_distance - distance);
    if (distance >= 2 * max_distance) return distance;
    if (distance > this_max) this_max = distance;
    if (distance == 0) {
      *idx_match = i;
      return 0;
    }
  }
  return this_max;
}

// TODO: distance between contexts (for each context in c1, smallest dist to all contexts in c2)
int
distance_between_context_histograms (context_histogram_t c1, context_histogram_t c2, double *result)
{ 
  int i, min_n;
  double x,y;
  min_n = BIOMCMC_MIN (c1->h->n, c2->h->n);
  result[0] = 0.; // chi-square based distance, but 'bin-bin' with modes only 
  for (i = 0; i < min_n; i++) {
    x = (double)(1 + abs (c1->h->i[i].idx - c2->h->i[i].idx)); // absolute distance between nth-mode plus one
    y = (double)(c1->h->i[i].freq)/(double)(c1->integral) - (double)(c2->h->i[i].freq)/(double)(c2->integral); // absolute distance between nth-mode plus one
    result[0] += x * fabs (y);
  }
  result[1] = 0.; x = y = 0.; // weighted average
  for (i = 0; i < c1->h->n; i++) x += (double)(c1->h->i[i].freq * c1->h->i[i].idx)/(double)(c1->integral);
  for (i = 0; i < c2->h->n; i++) y += (double)(c2->h->i[i].freq * c2->h->i[i].idx)/(double)(c2->integral);
  x = x - y;
  result[1] = x * x; 

  return 2; // number of results
}

int
compare_context_histogram_for_qsort (const void *a, const void *b) // increasing, resolve ties with integral (e.g. location==-1)
{
  int result = (*(context_histogram_t*)a)->location - (*(context_histogram_t*)b)->location; // increasing
  if (result) return result;
  result = (*(context_histogram_t*)b)->integral - (*(context_histogram_t*)a)->integral;  // decreasing
  if (result) return result;
  /* solve ties by lexicographic order of contexts (like hopo_element) */
  if ((*(context_histogram_t*)b)->context[0] > (*(context_histogram_t*)a)->context[0]) return 1; 
  if ((*(context_histogram_t*)b)->context[0] < (*(context_histogram_t*)a)->context[0]) return -1; 
  if ((*(context_histogram_t*)b)->context[1] > (*(context_histogram_t*)a)->context[1]) return 1; 
  if ((*(context_histogram_t*)b)->context[1] < (*(context_histogram_t*)a)->context[1]) return -1; 
  return 0;
}

bool 
context_histograms_overlap (context_histogram_t c1, context_histogram_t c2, int *distance, tatajuba_options_t opt)
{
  int tract_length, n_bases_apart;

  n_bases_apart = abs (c1->location - c2->location);
  if (!n_bases_apart) {
    if (distance) *distance = 0;
    return true; // same location
  }

  tract_length = BIOMCMC_MIN (c1->mode_context_length, c2->mode_context_length);
//  printf ("DBG::LEN::%6d :: %6d %6d\n", c1->location, n_bases_apart, tract_length);
  if (n_bases_apart > (tract_length - 1) || (c1->base != c2->base)) { // bases SHOULD be the same, but you never know...
    if (distance) *distance = -1; // should NOT be used
    return false; // locations too different
  }

  uint32_t i1, i2;
  i1 = strlen (c1->name); // name = "ATTGC-A-CCCAG"
  i2 = biomcmc_levenshtein_distance (c1->name, i1, c2->name, i1, 1, 1, true); // allows for indels
  if (distance) *distance = (int) i2;
//  printf ("DBG::DIST::%6d :: %8u\t", c1->location, i2);
  if (i2 <= (uint32_t)(opt.levenshtein_distance)) return true;
  return false;
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
