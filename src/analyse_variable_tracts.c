/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "analyse_variable_tracts.h"

file_compress_t initialise_vcf_file (char *name);

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

  for (i = 0; i < g->n_genome; i++) vcf[i] = initialise_vcf_file (g->genome[i]->name);


  for (i = 0; i < g->n_genome; i++) biomcmc_close_compress (vcf[i]);
}

file_compress_t
initialise_vcf_file (char *name)
{
  file_compress_t vcf = NULL;
  size_t l = strlen (name), buffer_size = 8192;
  char *s = biomcmc_malloc (buffer_size * sizeof (char));
  strncpy (s, name, l);
  strncpy (s + l, ".vcf", 4);
  l += 4;
  s[l] = '\0';
#ifdef HAVE_ZLIB
  strncpy (s + l, ".gz", 3); // by adding suffix create_compress() can guess the library to use
  s[l+3] = '\0';
#endif
  vcf = biomcmc_create_compress_from_suffix (s); 

  //memset (s, '\0', sizeof (char) * buffer_size); // also strcat(s, tmp) to concat tmp to first null of s
  sprintf (s, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",name);
  biomcmc_write_compress (vcf, s);

  if (s) free (s);
  return vcf;
}
