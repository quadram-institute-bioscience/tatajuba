/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "analyse_variable_tracts.h"

file_compress_t initialise_vcf_file (char *outdir, char *name);
void update_vcf_file_from_context_histogram_pilot (genome_set_t g, context_histogram_t concat, file_compress_t vcf);
void update_vcf_file_from_context_histogram (genome_set_t g, context_histogram_t concat, file_compress_t vcf);
int  get_next_ht_location_from_same_contig (genome_set_t g, context_histogram_t concat);
int find_ref_alt_ht_variants_from_strings (char *ref, char *alt, genome_set_t g, context_histogram_t concat);

void
generate_vcf_files (genome_set_t g)
{
  file_compress_t *vcf;
  int i,j;

  vcf = (file_compress_t*) biomcmc_malloc (g->n_genome * sizeof (file_compress_t));
#ifdef HAVE_ZLIB
  printf ("Will save VCF files with gzipped compression\n");
#else
  biomcmc_warning ("ZLIB library not installed, will sabe VCF files uncompressed\n");
#endif

  for (i = 0; i < g->n_genome; i++) vcf[i] = initialise_vcf_file (g->genome[0]->opt.outdir, g->genome[i]->name);

  for (i = 0; i < g->tract->n_var; i++) for (j = g->tract->var_initial[i]; j < g->tract->var_final[i]; j++) {
    update_vcf_file_from_context_histogram (g, g->tract->concat[j], vcf[ g->tract->concat[j]->index ]);
  }

  for (i = 0; i < g->n_genome; i++) biomcmc_close_compress (vcf[i]);
}

file_compress_t
initialise_vcf_file (char *outdir, char *name)
{
  file_compress_t vcf = NULL;
  size_t buffer_size = 8192;
  char *s = biomcmc_malloc (buffer_size * sizeof (char));

  memset (s, '\0', sizeof (char) * buffer_size); // strcat starts at first null char 
  strcpy (s, outdir);
  strcat (s, name);
  strcat (s, ".vcf");
#ifdef HAVE_ZLIB
  strcat (s, ".gz"); // by adding suffix create_compress() can guess the library to use
#endif
  vcf = biomcmc_create_compress_from_suffix (s); 

  memset (s, '\0', sizeof (char) * buffer_size);
  sprintf (s, "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##INFO=<ID=TID,Number=A,Type=String,Description=\"tract ID\">\n");
  biomcmc_write_compress (vcf, s);
  memset (s, '\0', sizeof (char) * buffer_size);
  sprintf (s, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",name);
  biomcmc_write_compress (vcf, s);

  if (s) free (s);
  return vcf;
}

void  // first version, superseded by one below
update_vcf_file_from_context_histogram_pilot (genome_set_t g, context_histogram_t concat, file_compress_t vcf)
{
  char *ref_contig, *query_sequence;
  size_t query_size, ref_size = concat->loc2d[2] - concat->loc2d[1] + 1; //last and first positions in ref_contig (which is whole chromosome/genome)
  int l[2], ht_limit = 0; // l[] = pref/suf lengths

  ref_contig = g->genome[0]->opt.gff->sequence->string[ g->tract_ref[ concat->tract_id ].fasta_idx ]; // tract_id=tid, fasta_idx=location in char_vector "sequence"
  query_sequence = generate_tract_as_string (concat->context + (2 * concat->mode_context_id), concat->base, 
                                             g->genome[0]->opt.kmer_size, concat->mode_context_length, concat->neg_strand);
  query_size = strlen (query_sequence);

  if (!common_prefix_suffix_lengths_from_strings (ref_contig + concat->loc2d[1], ref_size, query_sequence, query_size, l)) {
    if (query_sequence) free (query_sequence);
    return; // seqs are identical
  }
  ht_limit = g->genome[0]->opt.kmer_size + concat->mode_context_length + 1;
  if ((l[0] > ht_limit) || (l[1] > ht_limit)) {
    if (query_sequence) free (query_sequence);
    return; // modification is on context, not on HT (identical region spans HT in this sample)
  }

  printf ("DBG %.*s\n    %.*s %4d %4d tid_%06d\n", (int)ref_size, ref_contig+concat->loc2d[1], (int)query_size, query_sequence, l[0], l[1], concat->tract_id);

  size_t buffer_size = 8196;
  int i, j, last, pos; // pos=position in vcf file (one based)
  char *s = NULL, *altseq = NULL, *refseq = NULL;
  s = biomcmc_malloc (buffer_size * sizeof (char));
  altseq = biomcmc_malloc ((query_size + 2) * sizeof (char));
  refseq = biomcmc_malloc ((ref_size + 2) * sizeof (char));

  if (concat->loc2d[1] + l[0] < 1) { // if both are zero, start of chromosome is already different, then common base _after_ event must be included
    //ref: ABCDEF
    //alt:   CDEF  leads to ref:ABC alt:C starting at one
    pos = 1;
    last = concat->loc2d[2] - l[1] + 2; // adds one common base after event (loc2d[2] is last base _belonging_ to match)
    if (!l[1]) last--; // fringe case with no suffix in common (on top of no prefix in common!); maybe alignment is in middle and bwa predicted based on query length
//    printf ("DEBUG1::\t%d\t%d\t%d\t%lu\t%d\t%d\t", concat->loc2d[1], concat->loc2d[2], last, query_size, l[0], l[1]);
    for (i = 0; i < last; i++) refseq[i] = ref_contig[i]; 
    refseq[i] = '\0';
    last = (int)(query_size) - l[1] + 1; // adds one common base after event
//    printf ("%d\n", last);  // DEBUG
    for (i = 0; i < last; i++) altseq[i] = query_sequence[i]; 
    altseq[i] = '\0';
  }
  else {
    pos = concat->loc2d[1] + l[0]; // "+1" since vcf is one-based, "-1" since we start at common base
    last = concat->loc2d[2] - l[1] + 1; // even if no suffix in common, we just copy to the end (only case above needs rightmost common base)
//    printf ("DEBUG2::\t%d\t%d\t%d\t%lu\t%d\t%d\t(%d)\t", concat->loc2d[1], concat->loc2d[2], last, query_size, l[0], l[1], l[1]+l[0]);
    for (j=0, i = pos-1; i < last; j++, i++) refseq[j] = ref_contig[i]; 
    refseq[j] = '\0';
    last = (int)(query_size) - l[1];
//    printf ("%d\n", last); // DEBUG
    altseq[0] = ref_contig[pos-1]; // first must be the same
    for (j = 1, i = l[0]; i < last; j++, i++) altseq[j] = query_sequence[i]; 
    altseq[j] = '\0';
  }

  // FIXME: check if row is repeated (rare cases might have distinct tids for same location on same sample, but usually it's the same)
  // some that appear repeated are context changes spanning several HTs eg. two HTs  axcgggacTTTacac and axcGGGactttacac , where x differs
  memset (s, '\0', sizeof (char) * buffer_size);
  //example: "NC_017280 44  . G A . . TID=tid001   GT  1"
  sprintf (s, "%s\t%d\t.\t%s\t%s\t.\t.\tTID=tid_%06d\tGT\t1\n", g->tract_ref[ concat->tract_id ].contig_name, pos, refseq, altseq, concat->tract_id); 
  biomcmc_write_compress (vcf, s);

  if (query_sequence) free (query_sequence);
  if (s) free (s);
  if (refseq) free (refseq);
  if (altseq) free (altseq);
  return;
}

void
update_vcf_file_from_context_histogram (genome_set_t g, context_histogram_t concat, file_compress_t vcf)
{
  char *s = NULL, *ref_contig = NULL, *ref_sequence = NULL, *alt_sequence = NULL; // both _sequences will be modified to keep only region around HT which differ
  size_t buffer_size = 8196, ref_size = concat->loc2d[2] - concat->loc2d[1] + 1; //last and first positions, inclusive, in ref_contig (which is whole chromosome/genome)
  int position = -1, next_ht_length = 0; 

  // 1. create reference and query strings
  ref_contig = g->genome[0]->opt.gff->sequence->string[ g->tract_ref[ concat->tract_id ].fasta_idx ]; // tract_id=tid, fasta_idx=location in char_vector "sequence"
  ref_sequence = biomcmc_malloc ((ref_size + 1) * sizeof (char));
  memcpy (ref_sequence, ref_contig + concat->loc2d[1], ref_size);
  ref_sequence[ref_size] = '\0';
  // 1.1. create query string
  alt_sequence = generate_tract_as_string (concat->context + (2 * concat->mode_context_id), concat->base, 
                                             g->genome[0]->opt.kmer_size, concat->mode_context_length, concat->neg_strand);
  // 2. stop at next HT (as given by distance from end of string)
  next_ht_length = concat->loc2d[2] + 1 - get_next_ht_location_from_same_contig(g, concat); // how many bases, from end, overlap with next HT
  //printf("DEBUG:: %s\nDEBUG:: %s\ttid_%06d %4d BEFORE\n", ref_sequence, alt_sequence, concat->tract_id, next_ht_length);
  if (next_ht_length > 0) { // if negative then next HT starts after end of right context
    ref_size -= next_ht_length;
    ref_sequence[ref_size] = '\0'; // afterl null char is ignored
    alt_sequence[ strlen(alt_sequence) - next_ht_length ] = '\0';
  }

  position = find_ref_alt_ht_variants_from_strings (ref_sequence, alt_sequence, g, concat); // modifies strings
  if (position < 0) {
    if (ref_sequence) free (ref_sequence);
    if (alt_sequence) free (alt_sequence);
    return;
  }

  // FIXME: check if row is repeated (rare cases might have distinct tids for same location on same sample, but usually it's the same)
  // some that appear repeated are context changes spanning several HTs eg. two HTs  axcgggacTTTacac and axcGGGactttacac , where x differs
  s = biomcmc_malloc (buffer_size * sizeof (char));
  memset (s, '\0', sizeof (char) * buffer_size);
  //example: "NC_017280 44  . G A . . TID=tid001   GT  1"
  sprintf (s, "%s\t%d\t.\t%s\t%s\t.\t.\tTID=tid_%06d\tGT\t1\n", g->tract_ref[ concat->tract_id ].contig_name, position, ref_sequence, alt_sequence, concat->tract_id); 
  biomcmc_write_compress (vcf, s);

  if (ref_sequence) free (ref_sequence);
  if (alt_sequence) free (alt_sequence);
  if (s) free (s);
  return;
}

int 
get_next_ht_location_from_same_contig (genome_set_t g, context_histogram_t concat)
{
  int i;
  for (i = concat->tract_id + 1; 
       (i < g->tract->n_concat) && 
       (g->tract_ref[i].fasta_idx == g->tract_ref[concat->tract_id].fasta_idx) && 
       (g->tract_ref[i].ht_location == g->tract_ref[concat->tract_id].ht_location);
       i++);
  if((i < g->tract->n_concat) && (g->tract_ref[i].fasta_idx == g->tract_ref[concat->tract_id].fasta_idx)) return g->tract_ref[i].ht_location;
  return concat->loc2d[2] + 1; // default is to say it's after end of right context (remember that loc2d[2] is _inclusive_, i.e. last base _in_ context
}

int
find_ref_alt_ht_variants_from_strings (char *ref, char *alt, genome_set_t g, context_histogram_t concat)
{
  tract_in_reference_s tref = g->tract_ref[ concat->tract_id ];
  int ht_alt = g->genome[0]->opt.kmer_size, ht_ref = tref.ht_location - concat->loc2d[1];
  int len_ref = (int)(strlen (ref)) - ht_ref, len_alt = (int)(strlen (alt)) - ht_alt; // HT + right context
  int i, j, l[2];
  // ref: ABCDEF  | according to VCF specs this is a special case leading to ref:ABC alt:C starting at one
  // alt:   CDEF  | but we ignore modifications in left context (we only consider HTs which are kmer positions downstream of genome start)

  // if same HT with same length, then we skip it
  if ((ref[ht_ref] == alt[ht_alt]) && (tref.tract_length == concat->mode_context_length)) return -1; 
  // finds identical suffix and prefix (identical sequences were already skipped above, since they would have same HT size)
  if (!common_prefix_suffix_lengths_from_strings (ref + ht_ref, (size_t) len_ref, alt + ht_alt, (size_t) len_alt, l)) return -1;

  //printf("DEBUG:: %s  %s\nDEBUG:: %s  %s\ttid_%06d %4d %4d\n", ref, ref+ht_ref-1, alt, alt+ht_alt-1, concat->tract_id, l[0], l[1]);
  if (l[0] > 0) { // same base, thus an HT or a monomer in ref genome
    int start = ht_ref + l[0] - 1; // "-1" since we start at last base in common
    for (j = 0, i = start; j < len_ref - l[1]; j++, i++) ref[j] = ref[i];
    ref[j] = '\0';
    start = ht_alt + l[0] - 1; // "-1" since we start at last base in common
    for (j = 0, i = start; j < len_alt - l[1]; j++, i++) alt[j] = alt[i];
    alt[j] = '\0';
    return tref.ht_location + l[0]; // minus one since we start at common base, however vcf is one-based (thus plus one)
  }
  return -1;
}
