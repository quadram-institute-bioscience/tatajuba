/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2020-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://leomrtns.github.io ]
 */

#include "context_histogram.h"

void add_new_context_histogram_from_hopo_elem (genomic_context_list_t genome, hopo_element he, char *name);
context_histogram_t new_context_histogram_from_hopo_elem (hopo_element he, char *name);
void context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, char *name, int idx_match);

char* context_histogram_tract_as_string (context_histogram_t ch, int kmer_size);
char* context_histogram_generate_name (context_histogram_t ch, int kmer_size);
void genomic_context_find_features (genomic_context_list_t genome);
void genomic_context_find_reference_location (genomic_context_list_t genome);
void genomic_context_merge_histograms_at_same_location (genomic_context_list_t genome);
void accumulate_from_context_histogram (context_histogram_t to, context_histogram_t from);

int
indel_distance_between_context_histogram_and_hopo_context (context_histogram_t ch, char *name)
{ 
  int len = strlen (ch->name);
  return biomcmc_levenshtein_distance (ch->name, len, name, len, 1, 1, true); // allows for indels
}

int
distance_between_context_histogram_and_hopo_context (context_histogram_t ch, hopo_element he, int max_distance, int location_difference, int *idx_match)
{ 
  int distance = 0, loc_diff = 0, this_max = 0, i;
  *idx_match = -1;
  if (ch->base != he.base) return CH_MAX_DIST; // homopolymer tracts are different

  loc_diff = he.read_offset - ch->location;
  if (loc_diff < 0) loc_diff = -loc_diff;
  if (loc_diff > location_difference) return CH_MAX_DIST;

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
  if (n_bases_apart > (tract_length - 1) || (c1->base != c2->base)) { // bases SHOULD be the same, but you never know...
    if (distance) *distance = -1; // should NOT be used
    return false; // locations too different
  }

  uint32_t i1, i2;
  i1 = strlen (c1->name); // name = "ATTGC-A-CCCAG"
  i2 = biomcmc_levenshtein_distance (c1->name, i1, c2->name, i1, 1, 1, true); // allows for indels
  if (distance) *distance = (int) i2;
  if (i2 <= (uint32_t)(opt.levenshtein_distance)) return true;
  return false;
}

void
add_new_context_histogram_from_hopo_elem (genomic_context_list_t genome, hopo_element he, char *name)
{ 
  genome->hist = (context_histogram_t*) biomcmc_realloc ((context_histogram_t*) genome->hist, (genome->n_hist+1) * sizeof (context_histogram_t));
  genome->hist[genome->n_hist++] = new_context_histogram_from_hopo_elem (he, name);
}

context_histogram_t
new_context_histogram_from_hopo_elem (hopo_element he, char *name)
{
  context_histogram_t ch = (context_histogram_t) biomcmc_malloc (sizeof (struct context_histogram_struct));

  ch->context = (uint64_t*) biomcmc_malloc (2 * sizeof (uint64_t));
  ch->ref_counter = 1;
  ch->n_context = 1; // how many context pairs 
  ch->mode_context_id = 0; // location, in context[], of best context (the one with length of highest frequency)
  /* since this is first context, it is best context: */
  ch->base = he.base;
  ch->indel = false;
  ch->multi = he.multi;
  ch->mode_context_count = he.count; // frequency of best context
  ch->mode_context_length = he.length; // tract length of best context
  ch->context[0] = he.context[0];
  ch->context[1] = he.context[1];
  ch->location = he.read_offset; 
  ch->loc2d[0] = he.loc_ref_id;
  ch->loc2d[1] = he.loc_pos;
  ch->neg_strand = he.neg_strand;
  ch->integral = he.count;
  ch->name = name;  
  ch->h = NULL; // empfreq created at the end (finalise_genomic_context)
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
context_histogram_add_hopo_elem (context_histogram_t ch, hopo_element he, char *name, int idx_match)
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
    ch->location = he.read_offset;
    ch->loc2d[0] = he.loc_ref_id;
    ch->loc2d[1] = he.loc_pos;
    if (ch->name) free (ch->name);
    ch->name = name;
  }
  else if (name) free (name);
/*
  if (ch->location > he.read_offset) {  // leftmost location
    ch->location = he.read_offset;
    ch->loc2d[0] = he.loc_ref_id;
    ch->loc2d[1] = he.loc_pos;
  }
*/
  if ((ch->multi ^ he.multi) == 1) ch->multi = 2; // if>0, then some are multi 
  ch->integral += he.count;
  // we assume ch->neg_strand is unchanged (since it's within distance boundary or just different length)

  ch->tmp_count  = (int*) biomcmc_realloc ((int*) ch->tmp_count,  (ch->index +1) * sizeof (int));
  ch->tmp_length = (int*) biomcmc_realloc ((int*) ch->tmp_length, (ch->index +1) * sizeof (int));
  ch->tmp_count [ch->index]   = he.count;
  ch->tmp_length[ch->index++] = he.length;
}

genomic_context_list_t
new_genomic_context_list (hopo_counter hc)
{
  int i, idx_match, j, distance;
  char *histname;
  
  /* 1. sort hopo_counter, calculate coverage etc. Remember that hopo_counter handles only identical contexts */
  finalise_hopo_counter (hc);
  if (hc->ref_start == hc->n_elem) {
    biomcmc_warning ("Sample %s doesn't contain any HT mapped to reference: it will be excluded from analysis", hc->name);
    return NULL;
  }

  genomic_context_list_t genome = (genomic_context_list_t) biomcmc_malloc (sizeof (struct genomic_context_list_struct));
  genome->hist = NULL;
  genome->n_hist = 0;
  genome->opt = hc->opt;
  genome->coverage = hc->coverage;
  genome->name = hc->name;
  hc->name = NULL;

  i = hc->ref_start;
  histname = generate_name_from_flanking_contexts (hc->elem[i].context, hc->elem[i].base, genome->opt.kmer_size, hc->elem[i].neg_strand);
  add_new_context_histogram_from_hopo_elem (genome, hc->elem[i], histname);

  /* 2. accumulate histograms of 'equivalent' contexts */
  for (++i; i < hc->n_elem; i++) {
    j = genome->n_hist - 1; 
    histname = generate_name_from_flanking_contexts (hc->elem[i].context, hc->elem[i].base, genome->opt.kmer_size, hc->elem[i].neg_strand);
    distance = distance_between_context_histogram_and_hopo_context (genome->hist[j], hc->elem[i], genome->opt.max_distance_per_flank, 
                                                                    genome->opt.min_tract_size, &idx_match);
    if (distance < genome->opt.max_distance_per_flank) { 
      context_histogram_add_hopo_elem (genome->hist[j], hc->elem[i], histname, idx_match);
    }
    else {
      if (distance < CH_MAX_DIST) // try again, now using indels 
        distance = indel_distance_between_context_histogram_and_hopo_context (genome->hist[j], histname);

      if (distance < genome->opt.levenshtein_distance) { 
        context_histogram_add_hopo_elem (genome->hist[j], hc->elem[i], histname, idx_match);
        genome->hist[j]->indel = true;
      }
      else  add_new_context_histogram_from_hopo_elem (genome, hc->elem[i], histname); 
    }
  }
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

  /* 2. find gff3 feature of each context */
  genomic_context_find_features (genome);

  /* 3. sort context_histograms based on genomic location, ties broken with more frequent first. Ties are found when 
   *    location == -1 i.e. not found on reference, and ultimately ties are sorted by context */
 // qsort (genome->hist, genome->n_hist, sizeof (context_histogram_t), compare_context_histogram_for_qsort);
  genome->ref_start = 0;

  /* 4. merge context_histograms mapped to same ref genome location. BWA may detect that slightly different contexts are 
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
{// alternative to ch->h->i[0].idx is is ch->mode_context_length
  return generate_tract_as_string (ch->context + (2 * ch->mode_context_id), ch->base, kmer_size, ch->h->i[0].idx); 
}

char*
context_histogram_generate_name (context_histogram_t ch, int kmer_size)
{
  return generate_name_from_flanking_contexts (ch->context + (2 * ch->mode_context_id), ch->base, kmer_size, ch->neg_strand); // assuming consecutive context[]
}

void
genomic_context_find_features (genomic_context_list_t genome)
{
  int i, j, n_features;
  gff3_fields *features;
  bwase_match_t match = new_bwase_match_t (genome->opt.reference_fasta_filename); // obtain ref seq names

  /* BWA was calculated in hopo_counter for each reference sequence (contig/genome) */
  for (i = 0; i < genome->n_hist; i++) {
    // find_gff3_fields_within_position (gff struct, name of reference genome exactly as in gff, position within ref genome, number of features found)
    j = genome->hist[i]->loc2d[1] + genome->opt.kmer_size;  // location of beginning of tract (after kmer-sized context)
    features = find_gff3_fields_within_position (genome->opt.gff, match->bns->anns[genome->hist[i]->loc2d[0]].name, j, &n_features);
    for (j = 0; j < n_features; j++) if (features[j].type.id != GFF3_TYPE_region) {
      genome->hist[i]->gffeature = features[j]; // store any feature (may be a gene, cds, ...); however ... 
      //printf ("DBG::%6d %6d | (%6d-%6d:%6d)  %s [%d]\n", i, j, features[j].start, features[j].end, match->m[i].position, features[j].attr_id.str, features[j].type.id);
      if (features[j].type.id == GFF3_TYPE_cds) {j = n_features; break; } // .. CDS have priority (overwrite previous and leave for loop)
    }
    if (features) free (features); 
  }
  del_bwase_match_t (match);
}

void
genomic_context_find_reference_location (genomic_context_list_t genome) // OBSOLETE (now it's done with hopo_counter)
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
genomic_context_merge_histograms_at_same_location (genomic_context_list_t genome) // a bit redundant with new_genomic_context_list() 
{
  context_histogram_t *new_h;
  int i, j = 0;

  new_h = (context_histogram_t*) biomcmc_malloc (genome->n_hist * sizeof (context_histogram_t));

  new_h[0] = genome->hist[0]; 
  for (i = 1; i < genome->n_hist; i++) {
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
    to->loc2d[0] = from->loc2d[0];
    to->loc2d[1] = from->loc2d[1];

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
      j = n1; // skip to next iteration (no need to compare with other to->context[] )
    } // if contexts are identical, do nothing

  del_context_histogram (from);
}
