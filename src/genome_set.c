/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "genome_set.h"

const char *fixed_fname[] = {
  "per_sample_average_length.tsv",
  "per_sample_modal_frequency.tsv",
  "per_sample_proportional_coverage.tsv",
  "selected_tracts_unknown.tsv",
  "selected_tracts_annotated.tsv",
  "tract_list.tsv",
  "variable_tracts.bed"
};

int sample_print_precision[] = {2,2,5}; // how many decimals to print for each of the N_FNAME_SAMPLE files below

#define N_FNAME_SAMPLE 3
enum {FNAME_SAMPLE_AVGELENGTH, FNAME_SAMPLE_MODALFREQ, FNAME_SAMPLE_PROPCOV, FNAME_SELECTED_TRACTS_UNKNOWN, FNAME_SELECTED_TRACTS_ANNOTATED, FNAME_TRACT_LIST,
  FNAME_BEDFILE}; 

#define N_DESC_STATS 5
enum {DESC_STAT_avgelength, DESC_STAT_modalfreq, DESC_STAT_propcov, DESC_STAT_covpercontext, DESC_STAT_entropy}; 
  
void remove_empty_genomes (genome_set_t g);
void simplify_genome_names (genomic_context_list_t *genome, int n_genome);

g_tract_vector_t new_g_tract_vector_from_genomic_context_list (genomic_context_list_t *genome, int n_genome);
void del_g_tract_vector (g_tract_vector_t tract);
void g_tract_vector_concatenate_tracts  (g_tract_vector_t tract, genomic_context_list_t *genome, int n_genome);
void update_g_tract_summary_from_context_histogram (g_tract_vector_t tract, int prev, int curr, int lev_distance, int n_genome);
void fill_g_tract_summary_tables (g_tract_s *this, context_histogram_t *concat, int prev, int curr);
void create_tract_in_reference_structure (genome_set_t g);
FILE * open_output_file (tatajuba_options_t opt, const char *file);
void find_best_context_name_for_reference (tract_in_reference_s *ref_tid, char *dnacontig, size_t dnacontig_len, tatajuba_options_t  opt, context_histogram_t hist);

void describe_statistics_for_genome_set (genome_set_t g);
void initialise_files_descriptive_stats (genome_set_t g, FILE **fout);
bool update_descriptive_stats_for_this_trait (genome_set_t g, int prev, int curr, double *stats_per_hist, double *samples_per_trait);
void print_descriptive_stats_per_sample (genome_set_t g, FILE **fout, double *samples_per_trait, int tract_id);
void descriptive_stats_of_histogram (context_histogram_t concat, double *result);
double relative_difference_of_vector (double *vec, int n_vec);

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
#pragma omp parallel for shared(g,filenames,opt) private(time0,time1,hc) schedule(dynamic) reduction(+:secs[:2])
#endif
    for (i = 0; i < g->n_genome; i++) {
      fprintf (stderr, "processing paired files %s and %s\n", filenames[2*i], filenames[2*i+1]);
      biomcmc_get_time (time0);
      hc = new_or_append_hopo_counter_from_file (NULL, filenames[2*i],   opt);
      hc = new_or_append_hopo_counter_from_file (hc,   filenames[2*i+1], opt);
      biomcmc_get_time (time1); secs[0] += biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];
      g->genome[i] = new_genomic_context_list (hc); // returns null if no HT is mapped to reference
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

      g->genome[i] = new_genomic_context_list (hc); // returns null if no HT is mapped to reference
      del_hopo_counter (hc); hc = NULL;
      biomcmc_get_time (time1); secs[1] += biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];
    }
  }
  g->secs[0] = secs[0]; g->secs[1] = secs[1]; 
  biomcmc_get_time (time0);
  remove_empty_genomes (g); // exclude genomes without mapped HTs, realloc'ing and modifying n_genome 

  /* label each histogram with index of genome they belong to */
  for (i = 0; i < g->n_genome; i++) for (j = 0; j < g->genome[i]->n_hist; j++) g->genome[i]->hist[j]->index = i;
  /* simplify sample names, removing suffix and prefix */
  simplify_genome_names (g->genome, g->n_genome);
  /* generate comparisons between genomes */
  g->tract = new_g_tract_vector_from_genomic_context_list (g->genome, g->n_genome);
  biomcmc_get_time (time1); g->secs[2] = biomcmc_elapsed_time (time1, time0); time0[0] = time1[0]; time0[1] = time1[1];

  create_tract_in_reference_structure (g);
  describe_statistics_for_genome_set (g);
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
//    for (i = g->n_tract_ref-1; i >= 0; i--) if (g->tract_ref[i].seq)        free (g->tract_ref[i].seq); 
    for (i = g->n_tract_ref-1; i >= 0; i--) if (g->tract_ref[i].tract_name) free (g->tract_ref[i].tract_name); 
    free (g->tract_ref);
  }
  del_char_vector (g->ref_names);
  del_g_tract_vector (g->tract);
  free (g); 
}

void
remove_empty_genomes (genome_set_t g)
{
  int i, j;
  for (j = 0, i = 0; i < g->n_genome; i++) if (g->genome[i] != NULL) g->genome[j++] = g->genome[i];
  if (!j) { del_genome_set (g); biomcmc_error ("No sample has tracts mapped to reference genomel I can't proceed"); }
  if (j < i) biomcmc_warning ("Only %d out of %d genomes have mapped HTs and will thus be analysed", j, i);

  g->n_genome = j;
  g->genome = (genomic_context_list_t*) biomcmc_realloc ((genomic_context_list_t*) g->genome, g->n_genome * sizeof (genomic_context_list_t));
}

void 
simplify_genome_names (genomic_context_list_t *genome, int n_genome)
{
  int i, j;
  size_t n_i , n_j, n_small, pre, suf, min_pre = 0xffffff, min_suf = 0xffffff;
  if (n_genome < 2) biomcmc_warning ("Will not try to simplify sample names since only one sample available");
  for (i = 0; i < n_genome; i++) {
    n_i = strlen (genome[i]->name);
    for (j = 0; j < n_genome; j++) {
      n_small = n_j = strlen (genome[j]->name);
      if (n_small > n_i) n_small = n_i; // n_small = min (genome[i], genome[j])

      for (pre = 0; (pre < n_small) && (genome[i]->name[pre] == genome[j]->name[pre]); pre++);
      for (suf = 0; (suf < n_small) && (genome[i]->name[n_i - suf - 1] == genome[j]->name[n_j - suf - 1]); suf++);
      if (pre < min_pre) min_pre = pre;
      if (suf < min_suf) min_suf = suf;
    }
  }
  if ((min_pre == 0) && (min_suf == 0)) { 
    biomcmc_fprintf_colour (stderr, 0, 2, "Using original file names", " (could not find suffix or prefix in common)\n.");
    return; 
  }
  else {
    if (min_pre) biomcmc_fprintf_colour (stderr, 0, 2, "Modifying sample names: ", 
                                         "Removing prefix '%.*s' from sample names.\n", min_pre, genome[0]->name);
    if (min_suf) biomcmc_fprintf_colour (stderr, 0, 2, "Modifying sample names: ", 
                                         "Removing suffix '%s' from sample names.\n", genome[0]->name + strlen(genome[0]->name) - min_suf);

    for (i = 0; i < n_genome; i++) {
      n_j = strlen (genome[i]->name) - min_pre - min_suf; 
      for (pre = 0; pre < n_j; pre++) genome[i]->name[pre] = genome[i]->name[pre+min_pre];
      genome[i]->name[pre] = '\0';
      genome[i]->name = (char*) biomcmc_realloc ((char*) genome[i]->name, n_j + 1);
    }
  }
  return;
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
  tract->var_initial = NULL;
  tract->var_final = NULL;

  tract->n_summary = tract->n_concat = tract->n_var = 0;
  g_tract_vector_concatenate_tracts (tract, genome, n_genome); // updates tract->concat

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
  if (tract->var_initial) free (tract->var_initial);
  if (tract->var_final)   free (tract->var_final);
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
  /* 1. start new histogram */ 
  tract->n_concat = genome[0]->n_hist;
  tract->concat   = (context_histogram_t*) biomcmc_malloc (sizeof (context_histogram_t) * tract->n_concat);
  for (j=0; j < genome[0]->n_hist; j++) tract->concat[j] = genome[0]->hist[j];
  new_h = tract->concat;
  n_new_h = tract->n_concat;

  /* 2. concat[] will be new_h and genome[i]->hist merged, in order (both are ordered) */
  for (i = 1; i < n_genome; i++) {
    tract->n_concat = n_new_h + genome[i]->n_hist; // sum of lengths
    tract->concat = (context_histogram_t*) biomcmc_malloc (tract->n_concat * sizeof (context_histogram_t));
    for (j = 0, k1 = 0, k2 = 0; (k1 < n_new_h) && (k2 < genome[i]->n_hist); ) { 
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
print_selected_g_tract_vector (genome_set_t g) // called by main.c
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

  fout = open_output_file (g->genome[0]->opt, fixed_fname[FNAME_SELECTED_TRACTS_UNKNOWN]);
  fprintf (fout, "tract_id\tbegin_context\tn_genomes\tlev_distance\t|\trd_frequency\trd_avge_tract_length\trd_coverage\trd_context_covge\trd_entropy\n"); 
  for (i = 0; i < n_no; i++) {
    t = g->tract->summary + cd_no[i];
    //fprintf (fout, "%s\t%-8d %-5d %-5d | ", t->example->name, t->location, t->n_genome_id, t->lev_distance);
    fprintf (fout, "tid_%06d\t%8d\t%5d\t%5d\t|\t", cd_no[i], t->location, t->n_genome_id, t->lev_distance);
    for (j = 0; j < N_SUMMARY_TABLES; j++) fprintf (fout, "%8.6lf\t", t->reldiff[j]);
    fprintf (fout, "\n");
    // if (t->n_dist > 0) printf ("%12e %12e\n", t->d1[0], t->d2[0]);
  }

  fclose (fout); fout = NULL;
  fout = open_output_file (g->genome[0]->opt, fixed_fname[FNAME_SELECTED_TRACTS_ANNOTATED]);
  fprintf (fout, "tract_id\tGFF3_info\tbegin_context\tn_genomes\tlev_distance\t|\trd_frequency\trd_avge_tract_length\trd_coverage\trd_context_covge\trd_entropy\n"); 
  for (i = 0; i < n_yes; i++) {
    t = g->tract->summary + cd_yes[i]; 
    gfi = t->example->gffeature;
    //fprintf (fout, "%s\t%-8d %s %s %-5d %-5d | ", t->example->name, t->location, gfi.type.str, gfi.attr_id.str, t->n_genome_id, t->lev_distance);
    fprintf (fout, "tid_%06d\t%s\t%8d\t%5d\t%5d\t|\t", cd_yes[i], gfi.attr_id.str, t->location, t->n_genome_id, t->lev_distance);
    for (j = 0; j < N_SUMMARY_TABLES; j++) fprintf (fout, "%8.6lf\t", t->reldiff[j]);
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
    g->tract_ref[i].bwa_dist = 0xffe; // 12bits
    g->tract_ref[i].fasta_idx = g->tract_ref[i].ht_location = -1;
    g->tract_ref[i].contig_border[0] = 0xffff; g->tract_ref[i].contig_border[1] = -1; // limits given by BWA for HT+context, do not change afterwards 
    g->tract_ref[i].contig_name = NULL;
    g->tract_ref[i].tract_name = NULL;
  }

  /* 3. store information from context_histogram concat[] which will be used in copying reference sequence segments */
  for (i = 0; i < g->tract->n_concat; i++) {
    tid = g->tract->concat[i]->tract_id;
    if (!g->tract_ref[tid].contig_name) { // first time tid is met
      g->tract_ref[tid].contig_name = g->ref_names->string[ g->tract->concat[i]->loc2d[0] ]; // loc2d[0] is given by BWA index
      g->tract_ref[tid].contig_border[0] = g->tract->concat[i]->loc2d[1];
      g->tract_ref[tid].contig_border[1] = g->tract->concat[i]->loc2d[2];
      g->tract_ref[tid].bwa_dist   = g->tract->concat[i]->mismatches;
      g->tract_ref[tid].max_length = g->tract->concat[i]->mode_context_length;
      g->tract_ref[tid].concat_idx = i; 
    }
    else { // tid is present more than once: sequence in reference will be closest (fewer "mismatches", i.e. edit distance)   
      if (g->tract_ref[tid].contig_border[0] > g->tract->concat[i]->loc2d[1])
        g->tract_ref[tid].contig_border[0] = g->tract->concat[i]->loc2d[1]; // leftmost location

      if (g->tract_ref[tid].contig_border[1] < g->tract->concat[i]->loc2d[2])
        g->tract_ref[tid].contig_border[1] = g->tract->concat[i]->loc2d[2]; // rightmost  location

      if (g->tract_ref[tid].max_length < g->tract->concat[i]->mode_context_length) // longest HT
        g->tract_ref[tid].max_length = g->tract->concat[i]->mode_context_length;

      if (g->tract_ref[tid].bwa_dist > g->tract->concat[i]->mismatches) { // best HT across samples
        g->tract_ref[tid].bwa_dist   = g->tract->concat[i]->mismatches;
        g->tract_ref[tid].concat_idx = i; 
      }
    }
  }

  /* 4. copy tract info from fasta sequences */
  char_vector fasta = g->genome[0]->opt.gff->sequence;
  for (tid = 0; tid < g->n_tract_ref; tid++) {
    i = lookup_hashtable (g->genome[0]->opt.gff->seqname_hash, g->tract_ref[tid].contig_name); // gff struct can have missing seqs
    if (i < 0) biomcmc_error ("Contig/genome sequence %s not found in fasta file (check names in GFF and FASTA)", g->tract_ref[tid].contig_name);
    /* 4.1 create context+hopo name for reference as in samples */
    g->tract_ref[tid].fasta_idx = i; // index IN FASTA (embedded from gff), not from BWA ordering
    find_best_context_name_for_reference (&(g->tract_ref[tid]), fasta->string[i], fasta->nchars[i], g->genome[0]->opt, 
                                          g->tract->concat[g->tract_ref[tid].concat_idx]);
  }
  /* 5.  accumulate context_histogram from same ref.ht_location (i.e. merge tract_id) TODO */
}

void
find_best_context_name_for_reference (tract_in_reference_s *ref_tid, char *dnacontig, size_t dnacontig_len, tatajuba_options_t  opt, context_histogram_t hist)
{
  int extra_borders, min_tract_size, start_location, len, i, best_id, dist = 0xffff, best_dist = 0xffff; 
  bool neg_strand, need_monomer_search = true;
  hopo_counter hc = new_hopo_counter (opt.kmer_size);

  min_tract_size = opt.min_tract_size - 3; 
  if (min_tract_size < 2) min_tract_size = 2;
  extra_borders = opt.min_tract_size + 3; 
  start_location = ref_tid->contig_border[0] - extra_borders; // left shift (allow for mismatches) can even be longer than min tract length
  if (start_location < 0) start_location = 0;                 // since we handle spurious matches by chosing one with best distance
  len = ref_tid->contig_border[1] - ref_tid->contig_border[0] + 2 * extra_borders;
  if (len + start_location > (int) dnacontig_len) len = (int) dnacontig_len - start_location;

  // ht_location will point to beginning of homopolymer (instead of flanking region); for flanking region use contig_border[]
  update_hopo_counter_from_seq (hc, dnacontig + start_location, len, min_tract_size); 
  // if all homopolymers are from different base, or no homopolymer is found, then we assume reference might have a monomer
  for (i = 0; (i < hc->n_elem) && (need_monomer_search); i++) if (hist->base == hc->elem[i].base) need_monomer_search = false; 
  if (need_monomer_search) update_hopo_counter_from_seq_all_monomers (hc, dnacontig + start_location, len); // exhaustive monomers within kmers, if no HTs present

  if (!hc->n_elem) { // homopolymer not found; store the equivalent region from the reference (legacy code: should never happen since we get all monomers above)
    start_location = ref_tid->contig_border[0]; 
    len = ref_tid->contig_border[1] - ref_tid->contig_border[0];
    if (len + start_location > (int) dnacontig_len) len = (int) dnacontig_len - start_location; // never true but just as an assert
    ref_tid->tract_name = (char*) biomcmc_malloc (sizeof (char) * (len + 1));
    strncpy (ref_tid->tract_name, dnacontig + start_location, len);
    ref_tid->tract_name[len] = '\0';
    // points to beginning of where homopolymer should be (skips flanking region); since missing from reference, we assume it's midpoint
    ref_tid->ht_location = ref_tid->contig_border[0] + (int)(len/2);
    ref_tid->tract_length = 0;
    del_hopo_counter (hc);
    return;
  }
  // DEBUG // char *s1 = generate_name_from_flanking_contexts (hist->context + 2 * hist->mode_context_id, hist->base, opt.kmer_size, false);
  // DEBUG // printf("%s : %6d DEBUG \n", s1, start_location + opt.kmer_size); free (s1);
  best_id = 0; // in case the same HT (A/T or G/C) is not found
  for (i = 0; (i < hc->n_elem) && (best_dist > 0); i++) { // we could use bwa's edit distance (used to find hist), but it's safer to use same algo below as for sample HTs...
    dist = 0xfff;
    if (hist->base == hc->elem[i].base) dist = distance_between_context_kmer_pair (hc->elem[i].context, hist->context + 2 * hist->mode_context_id);
    if (dist < best_dist) { best_dist = dist; best_id = i; }
   /* 
    // BEGIN DEBUG
    s1 = generate_name_from_flanking_contexts (hc->elem[i].context, hc->elem[i].base, opt.kmer_size, false); 
    printf("%s :: score %6d\n", s1, dist); free(s1);
    */
  }
  ref_tid->tract_length = hc->elem[best_id].length;
  ref_tid->ht_location = start_location + hc->elem[best_id].read_offset + opt.kmer_size; // corrects read-based location by ref-based location; 
  if (hc->elem[best_id].canon_flag == 2) neg_strand = true; // contexts were flipped by tatajuba's canonical definition, we have them right to left
  else neg_strand = false; // stored order is same as reference genome (contexts are always stored in canonical)
  ref_tid->tract_name = generate_name_from_flanking_contexts (hc->elem[best_id].context, hc->elem[best_id].base, opt.kmer_size, neg_strand);
  del_hopo_counter (hc);
}

void
print_tract_list (genome_set_t g)
{
  FILE *fout = NULL;
  int tid;
  context_histogram_t hist;
  fout = open_output_file (g->genome[0]->opt, fixed_fname[FNAME_TRACT_LIST]);

  fprintf (fout, "tract_id\tcontig_name\tfeature_type\tfeature\tlocation_in_contig\tmax_tract_length\tref_tract_length\t\ttract\tref_tract\n");
  for (tid = 0; tid < g->n_tract_ref; tid++) {
    fprintf (fout, "tid_%06d\t", tid); 
    fprintf (fout, "%s\t", g->tract_ref[tid].contig_name); 
    hist = g->tract->concat[ g->tract_ref[tid].concat_idx ]; 
    if (gff3_fields_is_valid (hist->gffeature))
      fprintf (fout, "%s\t%s\t", hist->gffeature.type.str, hist->gffeature.attr_id.str);
    else
      fprintf (fout, "nc\tunannotated\t");
    fprintf (fout, "%d\t%d\t", g->tract_ref[tid].ht_location, g->tract_ref[tid].max_length);
    fprintf (fout, "%d\t", g->tract_ref[tid].tract_length);
    fprintf (fout, "%s\t", hist->name);
    fprintf (fout, "%s", g->tract_ref[tid].tract_name);
    fprintf (fout, "\n");
  }
  fclose (fout); fout = NULL;
}

/**   new functions   **/

void
describe_statistics_for_genome_set (genome_set_t g)
{
  int i, prev, tid;
  double *stats_per_hist, *samples_per_trait;
  FILE *fout[N_FNAME_SAMPLE], *fbed;
  bool to_print;
 
  for (i = 0; i < N_FNAME_SAMPLE; i++) fout[i] = NULL;
  stats_per_hist    = (double*) biomcmc_malloc (N_DESC_STATS * sizeof (double));
  samples_per_trait = (double*) biomcmc_malloc (N_DESC_STATS * g->n_genome * sizeof (double));

  initialise_files_descriptive_stats (g, fout);
  fbed = open_output_file (g->genome[0]->opt, fixed_fname[FNAME_BEDFILE]);
  //  fprintf (fbed, "chrom\tchromStart\tchromEnd\tname\n"); // BED files don't need the header?

  for (prev = 0, i = 1; i < g->tract->n_concat; i++) {
    if (g->tract->concat[prev]->tract_id != g->tract->concat[i]->tract_id) { // range  = prev, prev+1, ..., i-1  (doesnt include `i`)
      to_print = update_descriptive_stats_for_this_trait (g, prev, i, stats_per_hist, samples_per_trait); 
      if (to_print) {
        g->tract->var_initial = (int*) biomcmc_realloc ((int*) g->tract->var_initial, (g->tract->n_var + 1) * sizeof (int));
        g->tract->var_final   = (int*) biomcmc_realloc ((int*) g->tract->var_final,   (g->tract->n_var + 1) * sizeof (int));
        g->tract->var_initial[g->tract->n_var] = prev; 
        g->tract->var_final[g->tract->n_var++] = i;
        
        tid = g->tract->concat[prev]->tract_id;
        print_descriptive_stats_per_sample (g, fout, samples_per_trait, tid);
        // BED file
        fprintf(fbed, "%s\t%d\t%d\ttid_%06d\n", g->tract_ref[tid].contig_name,
                g->tract_ref[tid].contig_border[0] + g->genome[0]->opt.kmer_size, 
                g->tract_ref[tid].contig_border[1] - g->genome[0]->opt.kmer_size + 1, tid);
        //char *contig = g->genome[0]->opt.gff->sequence->string[g->tract_ref[tid].fasta_idx]; //DEBUG
        //printf("%5d : %5d %5d : ", tid, first, last); for (j=first; j < last; j++) {printf("%c", contig[j]);} printf(" %c\n", contig[j]);//DEBUG
        // END DEBUG
      }
      prev = i;
    }
  }
  to_print = update_descriptive_stats_for_this_trait (g, prev, i, stats_per_hist, samples_per_trait); // last trait
  if (to_print) {
    g->tract->var_initial = (int*) biomcmc_realloc ((int*) g->tract->var_initial, (g->tract->n_var + 1) * sizeof (int));
    g->tract->var_final   = (int*) biomcmc_realloc ((int*) g->tract->var_final,   (g->tract->n_var + 1) * sizeof (int));
    g->tract->var_initial[g->tract->n_var] = prev; 
    g->tract->var_final[g->tract->n_var++] = i;

    tid = g->tract->concat[prev]->tract_id;
    print_descriptive_stats_per_sample (g, fout, samples_per_trait, tid);
    // BED file
    fprintf(fbed, "%s\t%d\t%d\ttid_%06d\n", g->tract_ref[tid].contig_name,
            g->tract_ref[tid].contig_border[0] + g->genome[0]->opt.kmer_size, 
            g->tract_ref[tid].contig_border[1] - g->genome[0]->opt.kmer_size + 1, tid);
  }

  fclose(fbed);

  for (i = 0; i < N_FNAME_SAMPLE; i++) fclose (fout[i]);
  if (stats_per_hist) free (stats_per_hist);
  if (samples_per_trait) free (samples_per_trait);
  return;
}

void
initialise_files_descriptive_stats (genome_set_t g, FILE **fout)
{
  int i, j;
  for (j = 0; j < N_FNAME_SAMPLE; j++) {
    fout[j] = open_output_file (g->genome[0]->opt, fixed_fname[j]);
    fprintf (fout[j], "tract_id\tlocation\tfeature\t%s", "reference"); // cannot be ref genome name since there may be several contigs (user must check on tract_list.csv)
    for (i = 0; i < g->n_genome; i++) fprintf (fout[j], "\t%s", g->genome[i]->name);
    fprintf (fout[j], "\n");
  }
}

bool
update_descriptive_stats_for_this_trait (genome_set_t g, int prev, int curr, double *stats_per_hist, double *samples_per_trait)
{
  int i, j;
  double difference = 0.;
  for (i = 0; i < N_DESC_STATS * g->n_genome; i++) samples_per_trait[i] = 0.;
  for (i = prev; i < curr; i++) {
    descriptive_stats_of_histogram (g->tract->concat[i], stats_per_hist);
    for (j = 0; j < N_DESC_STATS; j++) samples_per_trait[g->tract->concat[i]->index + g->n_genome * j] = stats_per_hist[j];
  }
  if (curr - prev  < g->n_genome) return true; 
  // crude estimate of relative error
  difference  = relative_difference_of_vector (samples_per_trait + g->n_genome * 0, g->n_genome);
  difference += relative_difference_of_vector (samples_per_trait + g->n_genome * 1, g->n_genome);
  difference += relative_difference_of_vector (samples_per_trait + g->n_genome * 4, g->n_genome);
  if (difference > 1.e-5) return true;
  for (i = prev; i < curr; i++) if (g->tract_ref[ g->tract->concat[prev]->tract_id ].tract_length != g->tract->concat[i]->h->i[0].idx) return true;
  return false;
}

void
print_descriptive_stats_per_sample (genome_set_t g, FILE **fout, double *samples_per_trait, int tract_id)
{
  int i, j;
  context_histogram_t hist = g->tract->concat[ g->tract_ref[tract_id].concat_idx ]; // GFF is only in concat[], not tract_ref
  for (j = 0; j < N_FNAME_SAMPLE; j++) {
    fprintf (fout[j], "tid_%06d\t%d\t", tract_id, g->tract_ref[tract_id].ht_location);
    if (gff3_fields_is_valid (hist->gffeature)) fprintf (fout[j], "%s", hist->gffeature.attr_id.str);
    else                                        fprintf (fout[j], "%s", "unannotated");
  }

  // values for reference genome
  fprintf (fout[FNAME_SAMPLE_AVGELENGTH], "\t%.0lf", (double)(g->tract_ref[tract_id].tract_length)); // weighted length
  fprintf (fout[FNAME_SAMPLE_MODALFREQ], "\t%.0lf", 1.0); // modal freq = histogram bar length of modal length
  fprintf (fout[FNAME_SAMPLE_PROPCOV], "\t"); // tract coverage depth /  deepest coverage depth across all tracts

  for (j = 0; j < N_FNAME_SAMPLE; j++) {// [ [genome1, ..., ngenomes]d_1 , [genome1, ..., ngenomes]d_2, ...,[]d_N_STATS ]
    for (i = 0; i < g->n_genome; i++) {
      if ( samples_per_trait[j * g->n_genome + i] > 0.) fprintf (fout[j], "\t%.*lf", sample_print_precision[j], samples_per_trait[j * g->n_genome + i]); 
      else                                              fprintf (fout[j], "\t");
    }
    fprintf (fout[j], "\n");
  }
  return;
}

void
descriptive_stats_of_histogram (context_histogram_t concat, double *result)
{
  double x;
  int j;

  for (j = 0; j < N_DESC_STATS; j++) result[j] = 0.;

  for (j = 0; j < concat->h->n; j++) if (concat->integral) // weighted average tract length
    result[DESC_STAT_avgelength] += (double) (concat->h->i[j].freq * concat->h->i[j].idx)/(double)(concat->integral);

  result[DESC_STAT_modalfreq] = (double) (concat->h->i[0].freq)/(double)(concat->integral); // modal frequency

  // proportional coverage ("coverage" is the depth of most frequent kmer in fastq, an estimate of genomic coverage depth)
  result[DESC_STAT_propcov] = (double)(concat->integral)/(double)(concat->coverage);

  // average coverage per context 
  result[DESC_STAT_covpercontext] = (double)(concat->integral)/(double)(concat->n_context);

  for (j = 0; j < concat->h->n; j++) { // entropy
    if (concat->integral) x = (double) (concat->h->i[j].freq)/(double)(concat->integral);
    result[DESC_STAT_entropy] += (x * log (x)); 
  }
  result[DESC_STAT_entropy] *= -1.;
    
  //result[5] = (double) concat->h->i[0].idx;  // modal (=typical) tract length 

  return;
}

double
relative_difference_of_vector (double *vec, int n_vec)
{
  double x_max = -FLT_MAX, x_min = FLT_MAX; // FLT and not DBL since we take negative 
  int i;

  for (i = 0; i < n_vec; i++) {
    if (x_max < vec[i]) x_max = vec[i];
    if (x_min > vec[i]) x_min = vec[i];
  }
  if (x_max > DBL_MIN) return (x_max - x_min)/*/x1*/;
  else return 0;
}
