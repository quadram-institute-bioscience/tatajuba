/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "hopo_counter.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);

static uint8_t dna_in_2_bits[256][2] = {{0xff}};
static char bit_2_dna[] = {'A', 'C', 'G', 'T'};
static void initialize_dna_to_bit_tables (void);
hopo_counter new_hopo_counter (void);
void update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int kmer_size, int min_hopo_size);

hopo_counter
new_hopo_counter_from_file (const char *filename)
{
  int i;
  gzFile fp = gzopen (filename, "r");
  kseq_t *seq = kseq_init (fp);
  hopo_counter hc = new_hopo_counter ();
  while ((i = kseq_read (seq)) >= 0) update_hopo_counter_from_seq (hc, seq->seq.s, seq->seq.l, 5, 3); // seq_l, kmer_size, min_homopol_size 
  kseq_destroy(seq); // other kseq_t parameters: seq->name.s, seq->seq.l, seq->qual.l
  gzclose(fp);
  return hc;
}

hopo_counter
new_hopo_counter (void)
{
  hopo_counter hc = (hopo_counter) biomcmc_malloc (sizeof (struct hopo_counter_struct));
  hc->ref_counter = 1;
  if (dna_in_2_bits[0][0] == 0xff) initialize_dna_to_bit_tables (); 

  return hc;
}

void
del_hopo_counter (hopo_counter hc)
{
  if (!hc) return;
  if (--hc->ref_counter) return;

  free (hc);
  return;
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

void
update_hopo_counter_from_seq (hopo_counter hc, char *seq, int seq_length, int kmer_size, int min_hopo_size)
{
  int i, j, k, *freq, count_same = 0, start_mono = -1;
  uint8_t context[2 * kmer_size]; 
  char prev_char = '$';

  printf ("\nDBG:: %s (%d)\n", seq, seq_length);
  count_same = 0;
  for (i = 0; i < (seq_length - kmer_size); i++) {
    if (seq[i] == prev_char) {
      count_same++;
      if ((count_same > min_hopo_size) && (start_mono >= kmer_size)) {
        while (i < (seq_length-kmer_size-1) && (seq[i+1] == prev_char)) { i++; count_same++; }
        for (j=start_mono; j <= i; j++) printf ("%c", seq[j]); //DEBUG
        if (dna_in_2_bits[prev_char][0] < dna_in_2_bits[prev_char][1]) { // A or C : forward strand
          k = 0;
          for (j = start_mono - kmer_size; j < start_mono; j++) context[k++] = dna_in_2_bits[ seq[j] ][0]; 
          for (j = i + 1; j <= i + kmer_size; j++) context[k++] = dna_in_2_bits[ seq[j] ][0]; 
          printf (" + ");
        }
        else if (dna_in_2_bits[prev_char][0] > dna_in_2_bits[prev_char][1]) { // T/U or G : reverse strand
          k = 0; // it would be easier to copy above, but with k--; however I wanna implement rolling hash in future
          for (j = i + kmer_size; j > i; j--) context[k++] = dna_in_2_bits[ seq[j] ][1]; 
          for (j = start_mono - 1; j >= start_mono - kmer_size; j--) context[k++] = dna_in_2_bits[ seq[j] ][1]; 
          printf (" - ");
        } // elsif (dna[0] > dna[1]); notice that if dna[0] == dna[1] then do nothing (not an unambiguous base)
        for (j=0; j < 2*kmer_size; j++) printf ("%c", bit_2_dna[ context[j] ]); //DEBUG
        printf (" ( %3d %3d)\n", i, count_same);
      } // if (count_same>2) [i.e. we found valid homopolymer]
    } else { 
      count_same = 1;
      prev_char = seq[i];
      start_mono = i;
    } // else (seq[i] == prev_char)
  } // for (i in seq[])
}
