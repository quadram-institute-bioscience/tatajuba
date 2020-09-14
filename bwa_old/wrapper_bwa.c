/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file 
 *  \brief wrapper around BWA library developed by Heng Li
 *  https://github.com/lh3/bwa commit 13b5637 (dowloaded 2020.08.05)
 */

#include <unistd.h>  // access() returns zero in sucess, -1 if failure;  R_OK | W_OK | X_OK or exists: F_OK
#include "bwtindex.c"
#include "bwase.c"
#include "bwtaln.c"
#include "wrapper_bwa.h"

gap_opt_t * gap_opt_from_bwase_options_t (bwase_options_t bopt);
/*! \brief return the cigar string or assume it's a match of length p_len */ // Leo addition
char *alloc_cigar_string_bwaln (bwa_cigar_t *cigar, int n_cigar, int p_len);
char *alloc_cigar_string_bwmem (uint32_t *cigar, int n_cigar);

char*
save_bwa_index (const char *genome_filename, const char *suffix, char overwrite)
{
	int algo_type = BWTALGO_AUTO, block_size = 10000000;
  char *prefix;
  if (suffix) {
    prefix = (char*) biomcmc_malloc ((strlen(genome_filename) + strlen(suffix) + 2) * sizeof (char));
    strcpy (prefix, genome_filename); strcat (prefix, "."); strcat (prefix, suffix);
  }
  else {
    prefix = (char*) biomcmc_malloc ((strlen(genome_filename) + 1) * sizeof (char));
    strcpy (prefix, genome_filename);
  }
  if (overwrite == 0) { // verify if all index files exist
    char *idxname, *term[] = {".amb",".ann",".bwt",".pac",".sa"};
    char idx_file_exists = 1;
    int i;
    idxname = (char*) biomcmc_malloc ((strlen(prefix) + 5) * sizeof (char));
    for (i = 0; (i < 5) && idx_file_exists; i++) {
      strcpy(idxname, prefix); strcat(idxname, term[i]);
      if (access(idxname, F_OK) == -1) idx_file_exists = 0;
    }
    if (idxname) free (idxname);
    if (idx_file_exists) return prefix; // all files are present
  }
  bwa_idx_build (genome_filename, prefix, algo_type, block_size);
  return prefix;
}   

/** BWA-aln **/

int   // OLD
bwa_aln_bwase (const char *index_filename, char **seqname, char **dnaseq, char **qual, size_t *seq_len, int n_dnaseq, int n_occurrences, int **match_list, char sam_to_stdout)
{
  bwa_seq_t *seqs;
  int n_matches = 0;
  gap_opt_t *opt = gap_init_opt();
  char *prefix = save_bwa_index (index_filename, NULL, false);
  if (n_occurrences < 1) n_occurrences = 1;
  if (sam_to_stdout != 0) n_matches = -1; // bwa_sai2sam() will then print SAM to stdout instead of creating match_list[]

  // TODO: bwa works with chunks of 0x40000 (262k) query seqs 
  seqs = bwa_read_seq_from_vector (seqname, dnaseq, qual, seq_len, n_dnaseq, opt->trim_qual);
  seqs = bwa_aln_from_vector (prefix, seqs, n_dnaseq, opt);
  seqs = bwa_sai2sam_se_from_vector (prefix, seqs, n_dnaseq, n_occurrences, opt, match_list, &n_matches);
  bwa_free_read_seq (n_dnaseq, seqs);
  if (opt) free (opt);
  if (prefix) free (prefix);
  return n_matches; 
}

bwase_options_t
new_bwase_options_t (int level)
{ // this overwrites all gap_opt_t from bwa, so important to make them match here (have save defaults)
  bwase_options_t bopt;
  bopt.fnr = 0.04;    //-n NUM  missing prob under 0.02 err rate (bwa-aln allows integer for max_diff)
  bopt.max_gapo = 1;  //-o maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
  bopt.max_gape = -1; //-e maximum number of gap extensions, -1 for disabling long gaps [-1]\n"); -> becomes 6 if not changed
  bopt.s_mm = 3;      //-M mismatch penalty [%d]\n", opt->s_mm);
  bopt.s_gapo = 11;   //-O gap open penalty [%d]\n", opt->s_gapo);
  bopt.s_gape = 4;    //-E gap extension penalty [%d]\n", opt->s_gape);
  bopt.max_top2 = 30; //-R stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
  bopt.trim_qual = 0; //-q quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
  bopt.seed_len = 32; //-l seed length [%d]\n", opt->seed_len);
  bopt.max_seed_diff = 2;     //-k maximum differences in the seed [%d]\n", opt->max_seed_diff);
  bopt.indel_end_skip = 5;    //-i do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
  bopt.max_del_occ = 10;      //-d maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
  bopt.max_entries = 2000000; //-m maximum entries in the queue [%d]\n", opt->max_entries);
  bopt.logscaled = false;     //-L log-scaled gap penalty for long deletions\n");
  bopt.all_hits = false;      //-N non-iterative mode: search for all n-difference hits (slooow)\n");
  bopt.n_threads = 1;
  switch (level) {
    case 2: // short reads, loose but slow
      bopt.seed_len = 24;
      bopt.max_gapo = 3;
      bopt.max_gapo = 5;
      bopt.s_gape = 3;
      bopt.indel_end_skip = 2;
      bopt.max_seed_diff = 4;
      break;
    case 1: 
      bopt.seed_len = 35;
      bopt.max_gapo = 2;
      bopt.s_mm = 2;
      bopt.s_gapo = 10;
      bopt.s_gape = 3;
      bopt.fnr = 0.02;
      bopt.max_seed_diff = 6;
      bopt.logscaled = true;
      break;
    default: break; // do nothing 
  };
  return bopt;
}

gap_opt_t *
gap_opt_from_bwase_options_t (bwase_options_t bopt)
{         
  gap_opt_t *opt = gap_init_opt (); // with defaults from bwa-aln

  opt->max_gapo = bopt.max_gapo;
  opt->indel_end_skip = bopt.indel_end_skip;
  opt->max_del_occ = bopt.max_del_occ;
  opt->seed_len = bopt.seed_len; 
  opt->max_seed_diff = bopt.max_seed_diff;
  opt->max_entries = bopt.max_entries;
  opt->s_mm = bopt.s_mm;
  opt->s_gapo = bopt.s_gapo;
  opt->s_gape = bopt.s_gape;
  opt->max_top2 = bopt.max_top2;
  opt->trim_qual = bopt.trim_qual;
  opt->fnr = bopt.fnr;
  opt->n_threads = bopt.n_threads;

  if (bopt.max_gape > 0) {
    opt->max_gape = bopt.max_gape;
    opt->mode &= ~BWA_MODE_GAPE;  
  }
  else { // this is default behaviour on bwa-aln (with -1 from variable "opte")
    opt->max_gape = 6;
  }
  if (bopt.logscaled) opt->mode |= BWA_MODE_LOGGAP;
  if (bopt.all_hits)  { opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0xffff; } // bwa-aln even more agressive on max_top2
  return opt;
}

bwase_match_t
new_bwase_match_t (const char *index_filename)
{
  bwase_match_t match = (bwase_match_t) biomcmc_malloc (sizeof (struct bwase_match_struct));
  match->prefix = save_bwa_index (index_filename, NULL, false); // generates index files
  match->bns = bns_restore (match->prefix);
  match->m = NULL;
  match->n_m = 0;
  match->ref_counter = 1;
  return match;
}

void
del_bwase_match_t (bwase_match_t match)
{
  if (!match) return;
  if (--match->ref_counter) return;
  for (int i = 0; i < match->n_m; i++) if (match->m[i].cigar)  free (match->m[i].cigar); 
  if (match->m) free (match->m);
  if (match->bns) bns_destroy (match->bns);
  free (match);
  return;
}

bwase_match_t
new_bwase_match_from_bwa_and_char_vector (const char *index_filename, char_vector seqname, char_vector dnaseq, int n_occurrences, bwase_options_t bopt)
{
  int i, c, n_reads, chunks = 0x40000;
  bwa_seq_t *seqs;
  gap_opt_t *opt = gap_opt_from_bwase_options_t (bopt);

  if (seqname->nstrings != dnaseq->nstrings) biomcmc_error ("Vector length of read names do not match the one with reads in bwase_match()");
  bwase_match_t match = new_bwase_match_t (index_filename);

  for (c = 0; c < seqname->nstrings; c += chunks) { 
    n_reads = chunks;
    if ((c + n_reads) > seqname->nstrings) n_reads = seqname->nstrings - c; // last loop
    seqs = bwa_read_seq_from_vector (seqname->string + c, dnaseq->string + c, NULL, dnaseq->nchars + c, n_reads, opt->trim_qual);
    seqs = bwa_aln_from_vector (match->prefix, seqs, n_reads, opt);
    seqs = bwase_to_match_t (match, seqs, n_reads, n_occurrences, opt); // updates bwase_match_t
    bwa_free_read_seq (n_reads, seqs);
  }
  if (opt) free (opt);
  return match; 
}

bwa_seq_t *bwase_to_match_t (bwase_match_t match, bwa_seq_t *seqs, int n_dnaseq, int n_occ, gap_opt_t *opt)
{ // leo: prefix, sai (fn_sa), fastq (fn_fa); needs prefix since reads indices several times in different ways (cal_pac_pos())
  int i, j, m_aln = 0;
  int k, seqid, nn;
  long long tot_seqs = 0;
  bwt_aln1_t *aln = 0;
  bwa_seq_t *p;
  bwt_multi1_t *q;

  if (!match->bns) match->bns = bns_restore (match->prefix);
  bwase_initialize();
  srand48 (match->bns->seed);

  // copy alignment struct from seqs[i] to aln
  for (i = 0; i < n_dnaseq; i++) {
    p = seqs + i;
    if (p->n_aln > m_aln) {
      m_aln = p->n_aln;
      aln = (bwt_aln1_t*) biomcmc_realloc ((bwt_aln1_t*) aln, sizeof (bwt_aln1_t) * m_aln);
    }
    for (j = 0; j < p->n_aln; j++) *(aln + j) = *(p->aln + j); 
    bwa_aln2seq_core (p->n_aln, aln, p, 1, n_occ);
  }

  bwa_cal_pac_pos (match->bns, match->prefix, n_dnaseq, seqs, opt->max_diff, opt->fnr); // forward bwt will be destroyed here
  bwa_refine_gapped (match->bns, n_dnaseq, seqs, 0, 0); // last "0" -> reverse not needed , added by leo
  for (i = 0; i < n_dnaseq; ++i) if (seqs[i].type != BWA_TYPE_NO_MATCH) {
    p = seqs + i;
    j = pos_end (p) - p->pos; // j is the length of the reference in the alignment
    nn = bns_cnt_ambi (match->bns, p->pos, j, &seqid);  // get seqid of ref genome

    match->m = (bwase_elem_t*) biomcmc_realloc ((bwase_elem_t*) match->m, (match->n_m + 1) * sizeof (bwase_elem_t));

    match->m[match->n_m].query_id = i; // query number, name is  p->name
    match->m[match->n_m].ref_id   = seqid; // reference number , name is bns->anns[seqid].name)
    match->m[match->n_m].position = p->pos - match->bns->anns[seqid].offset; //  ZERO-based leftmost position
    match->m[match->n_m].gapo = p->n_gapo; // gap open  
    match->m[match->n_m].gape = p->n_gape; // gap extension
    match->m[match->n_m].mm   = p->n_mm;   // number of mismatches

    match->m[match->n_m].top1_hits = p->c1;  // shared between all hits for this read 
    match->m[match->n_m].top2_hits = p->c2;  // shared between all hits for this read 
    match->m[match->n_m].mapQ      = p->seQ; // shared between all hits for this read 

    match->m[match->n_m].neg_strand  = p->strand; // zero for positive strand, one for negative
    match->m[match->n_m].is_best_hit = 1;
    match->m[match->n_m].cigar = alloc_cigar_string_bwaln (p->cigar, p->n_cigar, p->len);
    match->n_m++;

    if (p->n_multi) for (k = 0; k < p->n_multi; ++k) { 
      q = p->multi + k;
      j = pos_end_multi (q, p->len) - q->pos;
      nn = bns_cnt_ambi (match->bns, q->pos, j, &seqid);
      
      match->m = (bwase_elem_t*) biomcmc_realloc ((bwase_elem_t*) match->m, (match->n_m + 1) * sizeof (bwase_elem_t));

      match->m[match->n_m].query_id = i; // query number, name is  p->name
      match->m[match->n_m].ref_id   = seqid; // reference number , name is bns->anns[seqid].name)
      match->m[match->n_m].position = q->pos - match->bns->anns[seqid].offset; //  ZERO-based leftmost position
      match->m[match->n_m].gapo = q->gapo; // gap open (added by leo) 
      match->m[match->n_m].gape = q->gape; // gap extension (added by leo)
      match->m[match->n_m].mm   = q->mm;   // number of mismatches (notice slighly different names, mm instead of n_mm etc.)

      match->m[match->n_m].top1_hits = p->c1;  // shared between all hits for this read (notice P and not Q)
      match->m[match->n_m].top2_hits = p->c2;  // shared between all hits for this read 
      match->m[match->n_m].mapQ      = p->seQ; // shared between all hits for this read 

      match->m[match->n_m].neg_strand  = q->strand; // zero for positive strand, one for negative
      match->m[match->n_m].is_best_hit = 0;
      match->m[match->n_m].cigar = alloc_cigar_string_bwaln (q->cigar, q->n_cigar, p->len);
      match->n_m++;
    } // multi
  }
  free(aln);
  return seqs;
}

char *
alloc_cigar_string_bwaln (bwa_cigar_t *cigar, int n_cigar, int p_len)
{
  int k;
  char *str;
  if (cigar) {
    str = (char*) biomcmc_malloc (sizeof (char) * n_cigar * 8); // overestimate
    str[0] = '\0';
    for (k = 0; k < n_cigar; ++k) sprintf (str, "%s%d%c", str, __cigar_len(cigar[k]), "MIDS"[__cigar_op(cigar[k])]);
    size_t len = strlen (str);
    str = (char*) biomcmc_realloc ((char*) str, len * sizeof (char));
    return str;
  } // if there is no cigar, then it's a perfect match 
  else {
    str = (char*) biomcmc_malloc (sizeof (char) * 8); // overestimate
    str[0] = '\0';
    sprintf (str, "%dM", p_len);
    size_t len = strlen (str);
    str = (char*) biomcmc_realloc ((char*) str, len * sizeof (char));
    return str;
  }
}

char *
bwase_match_ref_genome_name (bwase_match_t match, int i)
{
  return match->bns->anns[ match->m[i].ref_id ].name;
}

/** BWA-mem **/

bwmem_match_t
new_bwmem_match_t (const char *index_filename)
{
  bwmem_match_t match = (bwmem_match_t) biomcmc_malloc (sizeof (struct bwmem_match_struct));
  match->prefix = save_bwa_index (index_filename, NULL, false); // generates index files
  match->idx = bwa_idx_load(match->prefix, BWA_IDX_ALL);
  match->m = NULL;
  match->n_m = 0;
  match->ref_counter = 1;
  return match;
}

void
del_bwmem_match_t (bwmem_match_t match)
{
  if (!match) return;
  if (--match->ref_counter) return;
  for (int i = 0; i < match->n_m; i++) if (match->m[i].cigar)  free (match->m[i].cigar); 
  if (match->m) free (match->m);
  if (match->idx) bwa_idx_destroy (match->idx);
  free (match);
  return;
}

bwmem_match_t
new_bwmem_match_from_bwa_and_char_vector (const char *index_filename, char_vector dnaseq)
{
  mem_opt_t *opt = mem_opt_init ();
  bwmem_match_t match = new_bwmem_match_t (index_filename);
  for (int i = 0; i < dnaseq->nstrings; i++) update_bwmem_match (match, i, dnaseq->string[i], dnaseq->nchars[i], opt); 
  if (opt) free (opt);
  return match; 
}

void
update_bwmem_match (bwmem_match_t match, int query_id, char *seq, size_t len, mem_opt_t *opt)
{
  mem_alnreg_v ar;
  mem_aln_t a;
  bwmem_elem_t *p;

  ar = mem_align1 (opt, match->idx->bwt, match->idx->bns, match->idx->pac, len, seq); // all the hits
  if (!ar.n) return; // realloc (m, 0) is bad for business
  match->m = (bwmem_elem_t*) biomcmc_realloc ((bwmem_elem_t*) match->m, (match->n_m + ar.n) * sizeof (bwmem_elem_t));
  for (int i = 0; i < ar.n; ++i) { // traverse each hit
    a = mem_reg2aln (opt, match->idx->bns, match->idx->pac, len, seq, &ar.a[i]); // get forward-strand position and CIGAR
    p = match->m + match->n_m + i;
    p->query_id = query_id; 
    p->ref_id = a.rid;
    p->position = a.pos; 
    p->score = a.score;  // alternative would be ar.truesc
    p->rb = ar.a[i].rb; 
    p->re = ar.a[i].re;
    p->qb = ar.a[i].qb;
    p->qe = ar.a[i].qe;
    p->mapQ = a.mapq;
    p->edit_distance = a.NM;
    p->neg_strand = a.is_rev;
    p->is_primary = (ar.a[i].secondary < 0);
    p->cigar = alloc_cigar_string_bwmem (a.cigar, a.n_cigar);
  }
  match->n_m += ar.n;
  if (ar.a) free(ar.a); 
}

char *
alloc_cigar_string_bwmem (uint32_t *cigar, int n_cigar)
{
  int k;
  char *str;

  str = (char*) biomcmc_malloc (sizeof (char) * n_cigar * 12); // overestimate
  str[0] = '\0';
  for (k = 0; k < n_cigar; ++k) sprintf(str, "%s%d%c", str, cigar[k]>>4, "MIDSH"[cigar[k] & 0xf]);
  size_t len = strlen (str);
  str = (char*) biomcmc_realloc ((char*) str, len * sizeof (char));
  return str;
}
char *
bwmem_match_ref_genome_name (bwmem_match_t match, int i)
{
  return match->idx->bns->anns[ match->m[i].ref_id ].name;
}

