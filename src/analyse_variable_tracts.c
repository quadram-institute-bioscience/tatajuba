/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "analyse_variable_tracts.h"

file_compress_t initialise_vcf_file (char *outdir, char *name);

void
generate_vcf_files (genome_set_t g)
{
  file_compress_t *vcf;
  int i;

  vcf = (file_compress_t*) biomcmc_malloc (g->n_genome * sizeof (file_compress_t));
#ifdef HAVE_ZLIB
  printf ("Will save VCF files with gzipped compression\n");
#else
  biomcmc_warning ("ZLIB library not installed, will sabe VCF files uncompressed\n");
#endif

  for (i = 0; i < g->n_genome; i++) vcf[i] = initialise_vcf_file (g->genome[0]->opt.outdir, g->genome[i]->name);


  for (i = 0; i < g->n_genome; i++) biomcmc_close_compress (vcf[i]);
}

file_compress_t
initialise_vcf_file (char *outdir, char *name)
{
  file_compress_t vcf = NULL;
  size_t l = strlen (name), buffer_size = 8192;
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
  sprintf (s, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",name);
  biomcmc_write_compress (vcf, s);

  if (s) free (s);
  return vcf;
}

/*
void
update_vcf_files_for_tid (genome_set_t g, int i1, int i2)
{
"NC_017280 44  . G A . PASS TID=tid001   GT  0"
}
*/
