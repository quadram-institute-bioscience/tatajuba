/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "genome_set.h"

const char *filename[] = {
  "selected_traits_unknown.csv",
  "selected_traits_annotated.csv",
  "tract_list.csv"
};

enum {FNAME_SELECTED_TRAITS_UNKNOWN, FNAME_SELECTED_TRAITS_ANNOTATED, FNAME_TRACT_LIST}; 
  
g_tract_vector_t new_g_tract_vector_from_genomic_context_list (genomic_context_list_t *genome, int n_genome);
void del_g_tract_vector (g_tract_vector_t tract);
void g_tract_vector_concatenate_tracts  (g_tract_vector_t tract, genomic_context_list_t *genome, int n_genome);
void update_g_tract_summary_from_context_histogram (g_tract_vector_t tract, int prev, int curr, int lev_distance, int n_genome);
void fill_g_tract_summary_tables (g_tract_s *this, context_histogram_t *concat, int prev, int curr);
void create_tract_in_reference_structure (genome_set_t g);
FILE * open_output_file (tatajuba_options_t opt, const char *file);
int lookup_bruteforce (char_vector haystack, const char *needle);

genome_set_t
new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt) 
{
  int64_t time0[2], time1[2];
  int i, j;
  hopo_counter hc;
  genome_set_t g = (genome_set_t) biomcmc_malloc (sizeof (struct genome_set_struct));
  g->n_genome = n_filenames;
  g->ref_counter = 1;
  g->ref_names = NULL;
  g->tract_ref = NULL; // set by create_tract_by_ref()
  double secs[2];
  g->secs[0] = g->secs[1] =  g->secs[2] = g->secs[3] = 0.;
  secs[0] = secs[1] = 0.;

  if (opt.paired_end) g->n_genome /= 2;
  
  g->genome = (genomic_context_list_t*) biomcmc_malloc (g->n_genome * sizeof (genomic_context_list_t));

  if (opt.paired_end) {
#ifdef _OPENMP
#pragma omp parallel for shared(g,opt,filenames) private(time0,time1,hc) schedule(dynamic) reduction(+:secs[:2])
#endif
    for (i = 0; i < g->n_genome; i++) {
      fprintf (stderr, "processing paired files %s and %s\n", filenames[2*i], filenames[2*i+1]);
      biomcmc_get_time (time0);
      hc = new_or_append_hopo_counter_from_file (NULL, filenames[2*i],   opt);
      hc = new_or_append_hopo_counter_from_file (hc,   filenames[2*i+1], opt);
      biomcmc_get_time (time1); secs[0] += biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];

      g->genome[i] = new_genomic_context_list (hc);
      del_hopo_counter (hc); hc = NULL;
      biomcmc_get_time (time1); secs[1] += biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];
    }
  }
  else {
#ifdef _OPENMP
#pragma omp parallel for shared(g,opt,filenames) private(time0,time1,hc) schedule(dynamic) reduction(+:secs[:2])
#endif
    for (i = 0; i < g->n_genome; i++) {
      fprintf (stderr, "processing file %s\n", filenames[i]);
      biomcmc_get_time (time0);
      hc = new_or_append_hopo_counter_from_file (NULL, filenames[i], opt);
      biomcmc_get_time (time1); secs[0] += biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];

      g->genome[i] = new_genomic_context_list (hc);
      del_hopo_counter (hc); hc = NULL;
      biomcmc_get_time (time1); secs[1] += biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];
    }
  }
  g->secs[0] = secs[0]; g->secs[1] = secs[1]; 
  biomcmc_get_time (time0);

  /* label each histogram with index of genome they belong to */
  for (i = 0; i < g->n_genome; i++) for (j = 0; j < g->genome[i]->n_hist; j++) g->genome[i]->hist[j]->index = i;
  /* generate comparisons between genomes */
  g->tract = new_g_tract_vector_from_genomic_context_list (g->genome, g->n_genome);
  biomcmc_get_time (time1); g->secs[2] = biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];

  create_tract_in_reference_structure (g);
  biomcmc_get_time (time1); g->secs[3] = biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];

  return g;
}

void
del_genome_set (genome_set_t g)
{
  int i;
  if (!g) return;
  if (--g->ref_counter) return;
  for (i = g->n_genome - 1; i >= 0; i--) del_genomic_context_list (g->genome[i]);
  if (g->genome) free (g->genome);
  if (g->tract_ref) { // do not free contig_name (which is just a pointer to char_vector ref_names)
    for (i = g->n_tract_ref-1; i >= 0; i--) if (g->tract_ref[i].seq)        free (g->tract_ref[i].seq); 
    for (i = g->n_tract_ref-1; i >= 0; i--) if (g->tract_ref[i].tract_name) free (g->tract_ref[i].tract_name); 
    free (g->tract_ref);
  }
  del_char_vector (g->ref_names);
  del_g_tract_vector (g->tract);
  free (g); 
}

// FIXME: aim is to deprecate tract_summary (i.e. recalculate after reference is added; not as we go, below) 
g_tract_vector_t
new_g_tract_vector_from_genomic_context_list (genomic_context_list_t *genome, int n_genome)
{
  int i, prev, lev_distance = 0, max_lev_distance = 0;
  bool is_same;
  g_tract_vector_t tract = (g_tract_vector_t) biomcmc_malloc (sizeof (struct g_tract_vector_struct));
  tract->summary = NULL;
  tract->concat = NULL; // big vector with all tracts from all genomes in linear order
  tract->n_summary = tract->n_concat = 0;
  g_tract_vector_concatenate_tracts  (tract, genome, n_genome); // updates tract->concat

  tract->concat[0]->tract_id = 0;
  for (prev=0, i = 1; i < tract->n_concat; i++) { // overlap() is true if tracts are the same, false if summary 
    is_same = context_histograms_overlap (tract->concat[i], tract->concat[prev], &lev_distance, genome[0]->opt);
    if (!is_same) {
      update_g_tract_summary_from_context_histogram (tract, prev, i, max_lev_distance, n_genome); // updates tract->summary
      tract->concat[i]->tract_id = tract->concat[prev]->tract_id + 1;
      prev = i;
      max_lev_distance = 0;
    }
    else {
      if (max_lev_distance < lev_distance) max_lev_distance = lev_distance; // tracts are similar enough 
      tract->concat[i]->tract_id = tract->concat[prev]->tract_id;
    }
  }
  // call context() just to calculate levenshtein distance 
  context_histograms_overlap (tract->concat[i-1], tract->concat[prev], &lev_distance, genome[0]->opt);
  if (max_lev_distance < lev_distance) max_lev_distance = lev_distance; // tracts are similar enough 
  update_g_tract_summary_from_context_histogram (tract, prev, i-1, max_lev_distance, n_genome); //  last block i == tract->n_concat

  return tract;
}

void
del_g_tract_vector (g_tract_vector_t tract)
{
  if (!tract) return;
  if (tract->concat) free (tract->concat);
  if (tract->summary) {
    for (int i = tract->n_summary - 1; i >= 0; i--) {
      if (tract->summary[i].d1) free (tract->summary[i].d1);
      if (tract->summary[i].tab0) free (tract->summary[i].tab0);
      if (tract->summary[i].genome_id) free (tract->summary[i].genome_id);
      del_context_histogram (tract->summary[i].example);
    }
    free (tract->summary);
  }
  free (tract);
}

void 
g_tract_vector_concatenate_tracts (g_tract_vector_t tract, genomic_context_list_t *genome, int n_genome)
{ 
  int i, j, k1, k2, dif, n_new_h;
  context_histogram_t *new_h = NULL;
  /** genomic_context_list is always sorted by position (with position=-1 at beginning, before ref_start) */
  /* brief explanation: each genome has a set of histograms in order; concat 'percolates' them with tracts from all genomes in one vector
   * genome A: a1 | a5 | a8  (i.e. locations 1,5,and 8) 
   * genome B: b2 | b5 | b6
   * concat: a1 | b2| a5 | b5| b6 | a8
   */
  /* 1. skip genomes without (genome-location) annotated histograms and start new histogram */ 
  for (i = 0; genome[i]->ref_start == genome[i]->n_hist; i++); // do nothing besides loop
  tract->n_concat = genome[i]->n_hist - genome[i]->ref_start;
  tract->concat   = (context_histogram_t*) biomcmc_malloc (sizeof (context_histogram_t) * tract->n_concat);
  for (j=0, k1 = genome[i]->ref_start; k1 < genome[i]->n_hist; j++, k1++) tract->concat[j] = genome[i]->hist[k1];
  new_h = tract->concat;
  n_new_h = tract->n_concat;

  /* 2. concat[] will be new_h and genome[i]->hist merged, in order (both are ordered) */
  for (i+=1; i < n_genome; i++) if (genome[i]->n_hist > genome[i]->ref_start) {
    tract->n_concat = n_new_h + genome[i]->n_hist - genome[i]->ref_start; // sum of lengths
    tract->concat = (context_histogram_t*) biomcmc_malloc (tract->n_concat * sizeof (context_histogram_t));
    for (j = 0, k1 = 0, k2 = genome[i]->ref_start; (k1 < n_new_h) && (k2 < genome[i]->n_hist); ) { 
      dif = new_h[k1]->location - genome[i]->hist[k2]->location;
      if (dif > 0) tract->concat[j++] = genome[i]->hist[k2++]; // k2 steps forward (if equal both step forward)
      else if (dif < 0) tract->concat[j++] = new_h[k1++];      // k1 steps forward
      else { // both are the same, both step forward and tract increases by two
        tract->concat[j++] = genome[i]->hist[k2++]; 
        tract->concat[j++] = new_h[k1++];
      }
    } // loop ends when one of them finishes, now we need to complete with the other
    for (; k2 < genome[i]->n_hist; k2++, j++) tract->concat[j] = genome[i]->hist[k2]; 
    for (; k1 < n_new_h; k1++, j++) tract->concat[j] = new_h[k1];
    free (new_h); // only place where we free is here (o.w. we only have pointers)
    new_h = tract->concat;
    n_new_h = tract->n_concat;
  }
}

void
update_g_tract_summary_from_context_histogram (g_tract_vector_t tract, int prev, int curr, int lev_distance, int n_genome)
{
  int i, j, k;
  double result[2];

  tract->n_summary++;
  tract->summary = (g_tract_s*) biomcmc_realloc ((g_tract_s*) tract->summary, tract->n_summary * sizeof (g_tract_s));
  g_tract_s *this = tract->summary + (tract->n_summary-1);

  this->location = tract->concat[prev]->location;
  this->d1 = NULL;
  this->lev_distance = lev_distance;
  this->id_in_concat = prev;
  this->n_genome_total = n_genome;
  this->n_genome_id = curr - prev;

  /* 1. indices of genomes where this tract is present */
  this->genome_id = (int*) biomcmc_malloc ((curr-prev) * sizeof (int));
  for (i = 0, j = prev; j < curr; j++, i++) this->genome_id[i] = tract->concat[j]->index;

  /* 2. example of a context_histogram_t (if we dont have access to parent genome_t */
  this->example = tract->concat[prev]; 
  this->example->ref_counter++;

  /* 3. table of values per genome (modal length, entropy, etc) allocated as 1D vector, with N_SUMMARY pointers */
  this->tab0 = (double*) biomcmc_malloc (this->n_genome_id * N_SUMMARY_TABLES * sizeof (double));
  for (i = 0; i < this->n_genome_id * N_SUMMARY_TABLES; i++) this->tab0[i] = 0.; 
  for (i = 0; i < N_SUMMARY_TABLES; i++) this->gentab[i] = this->tab0 + (this->n_genome_id * i); // pointers 
  fill_g_tract_summary_tables (this, tract->concat, prev, curr);

  /* 4. matrix of pairwise distances (1D vector sorted from highest to smallest distance) FIXME: may remove this*/
  this->n_dist = ((curr-prev) * (curr-prev-1))/2;
  if (!this->n_dist) return;
  this->d1 = (double*) biomcmc_malloc (2 * this->n_dist * sizeof (double)); // d1 and d2
  this->d2 = this->d1 + this->n_dist;

  for (k = 0, i = prev+1; i < curr; i++) for (j = prev; j < i; j++) {
    distance_between_context_histograms (tract->concat[i], tract->concat[j], result);
    this->d1[k]   = result[0]; 
    this->d2[k++] = result[1]; 
  }
  qsort (this->d1, this->n_dist, sizeof (double), compare_double_decreasing);
  qsort (this->d2, this->n_dist, sizeof (double), compare_double_decreasing); 

  return;
}

void  
fill_g_tract_summary_tables (g_tract_s *this, context_histogram_t *concat, int prev, int curr)
{
  int i1, i2, j;
  double x1, x2, **gentab = this->gentab;

  for (i2 = 0, i1 = prev; i1 < curr; i1++, i2++) { // for every genome i1 that has this tract  (i2 is new index)
    x1 = (double) (concat[i1]->h->i[0].freq)/(double)(concat[i1]->integral); // modal frequency
    gentab[0][i2] = x1;

    x2 = 0.; // weighted average tract length
    for (j = 0; j < concat[i1]->h->n; j++) if (concat[i1]->integral) 
      x2 += (double) (concat[i1]->h->i[j].freq * concat[i1]->h->i[j].idx)/(double)(concat[i1]->integral);
    gentab[1][i2] = x2;  

    // proportional coverage ("coverage" is the depth of most frequent kmer)
    gentab[2][i2] = (double)(concat[i1]->integral)/(double)(concat[i1]->coverage);
   
    // average coverage per context 
    gentab[3][i2] = (double)(concat[i1]->integral)/(double)(concat[i1]->n_context);

    x1 = 0.; // entropy
    for (j = 0; j < concat[i1]->h->n; j++) {
      if (concat[i1]->integral) x2 = (double) (concat[i1]->h->i[j].freq)/(double)(concat[i1]->integral);
      x1 += (x2 * log (x2)); 
    }
    gentab[4][i2] = - x1;
  }

  for (i1 = 0; i1 < N_SUMMARY_TABLES; i1++) {
    x1 = -FLT_MAX; x2 = FLT_MAX; // FLT and not DBL since we take negative 
    for (i2 = 0; i2 < this->n_genome_id; i2++) {
      if (x1 < gentab[i1][i2]) x1 = gentab[i1][i2]; // x1 will be max and x2 will be min
      if (x2 > gentab[i1][i2]) x2 = gentab[i1][i2];
    }
    if (x1 > DBL_MIN) this->reldiff[i1] = (x1 - x2)/*/x1*/;
    else this->reldiff[i1] = 0.;
  }
}

void
print_selected_g_tract_vector (genome_set_t g)
{
  int i, j, *cd_yes, *cd_no, n_yes=0, n_no=0;
  g_tract_s *t;
  gff3_fields gfi;
  bool to_print = false;
  FILE *fout = NULL;

  cd_yes = (int*) biomcmc_malloc (g->tract->n_summary * sizeof (int));
  cd_no  = (int*) biomcmc_malloc (g->tract->n_summary * sizeof (int));

  for (i = 0; i < g->tract->n_summary; i++) { 
    to_print = false;
    t = g->tract->summary + i;
    if (t->n_genome_id  < t->n_genome_total) to_print = true;
    if ((!to_print) && (t->lev_distance > 0)) to_print = true;
    if ((!to_print) && (t->reldiff[0] > 1e-6)) to_print = true;
    if ((!to_print) && (t->reldiff[1] > 1e-6)) to_print = true;
    if ((!to_print) && (t->reldiff[4] > 1e-6)) to_print = true;
    if (to_print) {
      if (gff3_fields_is_valid (t->example->gffeature)) cd_yes [n_yes++] = i; // GFF3 has info about this location 
      else cd_no[n_no++] = i;
    }
  } // for i in summary
  printf ("From %d tracts, %d interesting ones are annotated and %d interesting ones are not annotated\n", g->tract->n_summary, n_yes, n_no);

  fout = open_output_file (g->genome[0]->opt, filename[FNAME_SELECTED_TRAITS_UNKNOWN]);
  fprintf (fout, "tract    location    n_genomes    lev_distance | rd_frequency rd_avge_tract_length rd_coverage rd_context_covge rd_entropy\n"); 
  for (i = 0; i < n_no; i++) {
    t = g->tract->summary + cd_no[i];
    fprintf (fout, "%s\t%-8d %-5d %-5d | ", t->example->name, t->location, t->n_genome_id, t->lev_distance);
    for (j = 0; j < N_SUMMARY_TABLES; j++) fprintf (fout, "%8.6lf ", t->reldiff[j]);
    fprintf (fout, "\n");
    // if (t->n_dist > 0) printf ("%12e %12e\n", t->d1[0], t->d2[0]);
  }

  fclose (fout); fout = NULL;
  fout = open_output_file (g->genome[0]->opt, filename[FNAME_SELECTED_TRAITS_ANNOTATED]);
  fprintf (fout, "tract    location  [GFF3_info]  n_genomes    lev_distance | rd_frequency rd_avge_tract_length rd_coverage rd_context_covge rd_entropy\n"); 
  for (i = 0; i < n_yes; i++) {
    t = g->tract->summary + cd_yes[i]; 
    gfi = t->example->gffeature;
    fprintf (fout, "%s\t%-8d [%s %s] %-5d %-5d | ", t->example->name, t->location, gfi.type.str, gfi.attr_id.str, t->n_genome_id, t->lev_distance);
    for (j = 0; j < N_SUMMARY_TABLES; j++) fprintf (fout, "%8.6lf ", t->reldiff[j]);
    fprintf (fout, "\n");
    // if (t->n_dist > 0) printf ("%12e %12e\n", t->d1[0], t->d2[0]);
  }

  fclose (fout); 
  if (cd_yes) free (cd_yes);
  if (cd_no)  free (cd_no);
}

void
print_debug_g_tract_vector (genome_set_t g)
{
  int i,j;
  g_tract_s *t;
  printf ("tract location n_genomes lev_distance | relative_differences \n"); 
  for (i = 0; i < g->tract->n_summary; i++) { 
    t = g->tract->summary + i;
    printf ("%s\t%-8d %-5d %-5d | ", t->example->name, t->location, t->n_genome_id, t->lev_distance);
    for (j = 0; j < N_SUMMARY_TABLES; j++) printf ("%9.7e ", t->reldiff[j]);
    printf ("\n");
  }
}

FILE *
open_output_file (tatajuba_options_t opt, const char *file)
{
  char *filename = NULL;
  FILE *fout = NULL;
  size_t len = strlen (opt.outdir) + strlen (file) + 1;
  filename = biomcmc_malloc (sizeof (char) * len);
  strcpy (filename, opt.outdir); strcat (filename, file); filename[len-1] = '\0';
  fout = biomcmc_fopen (filename, "w");
  if (filename) free (filename);
  return fout;
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
    len = g->tract_ref[tid].max_length + 2 * g->genome[0]->opt.kmer_size + min_tract_size;
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

void
print_tract_list (genome_set_t g)
{
  FILE *fout = NULL;
  int tid;
  context_histogram_t hist;
  fout = open_output_file (g->genome[0]->opt, filename[FNAME_TRACT_LIST]);

  fprintf (fout, "tract_id\tcontig_name\tfeature_type\tfeature\tlocation_in_contig\tmax_tract_length\ttract\t\n");
  for (tid = 0; tid < g->n_tract_ref; tid++) {
    fprintf (fout, "tid_%06d\t", tid); 
    fprintf (fout, "%s\t", g->tract_ref[tid].contig_name); 
    hist = g->tract->concat[ g->tract_ref[tid].first_idx ]; 
    if (gff3_fields_is_valid (hist->gffeature))
      fprintf (fout, "%s\t%s\t", hist->gffeature.type.str, hist->gffeature.attr_id.str);
    else
      fprintf (fout, "nc\tunnanotated\t");
    fprintf (fout, "%d\t%d\t", g->tract_ref[tid].contig_location, g->tract_ref[tid].max_length);
    fprintf (fout, "%d\t%s\t", g->tract_ref[tid].tract_length, g->tract_ref[tid].tract_name);
    fprintf (fout, "%s\t", hist->name);

    fprintf (fout, "\n");
  }
  fclose (fout); fout = NULL;
}

/*! \brief hash lookup is based on strict identity of strings; if fails, we assume one of them is shorter (always the
 * case for BWA (before space) and our fasta (whole line) */
int
lookup_bruteforce (char_vector haystack, const char *needle)
{
  int i, j, minlen;
  size_t n_len = strlen (needle);
  for (i = 0; i < haystack->nstrings; i++) {
    minlen = BIOMCMC_MIN(n_len, haystack->nchars[i]);
    for (j = 0; (j < minlen) && (needle[j] == haystack->string[i][j]); j++);
    if (j == minlen) return i;
  }
  return -1;
}

