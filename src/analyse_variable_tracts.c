/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "analyse_variable_tracts.h"

file_compress_t initialise_vcf_file (char *outdir, char *name);
void update_vcf_file_from_context_histogram (genome_set_t g, context_histogram_t concat, file_compress_t vcf);

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

void
update_vcf_file_from_context_histogram (genome_set_t g, context_histogram_t concat, file_compress_t vcf)
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

#ifdef COMMENTED_OUT
void
update_vcf_file_from_context_histogram_new (genome_set_t g, context_histogram_t concat, file_compress_t vcf)
{
  char *s = NULL, *ref_contig = NULL, *ref_sequence = NULL, *alt_sequence = NULL; // both _sequences will be modified to keep only region around HT which differ
  size_t buffer_size = 8196, ref_size = concat->loc2d[2] - concat->loc2d[1] + 1; //last and first positions, inclusive, in ref_contig (which is whole chromosome/genome)
  int position = -1; 

  // create reference sequence
  ref_contig = g->genome[0]->opt.gff->sequence->string[ g->tract_ref[ concat->tract_id ].fasta_idx ]; // tract_id=tid, fasta_idx=location in char_vector "sequence"
  ref_sequence = biomcmc_malloc ((query_size + 1) * sizeof (char));
  memncpy (ref_sequence, ref_contig + concat->loc2d[1], ref_size);
  ref_sequence[ref_size] = '\0';
  // create query sequence
  alt_sequence = generate_tract_as_string (concat->context + (2 * concat->mode_context_id), concat->base, 
                                             g->genome[0]->opt.kmer_size, concat->mode_context_length, concat->neg_strand);

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
find_ref_alt_ht_variants_from_strings (char *ref, char *alt, genome_set_t g, context_histogram_t concat)
{

}
#endif
