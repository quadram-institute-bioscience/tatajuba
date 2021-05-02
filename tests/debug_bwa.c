#include <wrapper_bwa.h>
#include <biomcmc.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, *match_list = NULL;
  bwase_match_t match;
  bwase_options_t bopt = new_bwase_options_t (2);
  alignment aln;

  time0 = clock ();

  if (argc == 1) return TEST_SKIPPED;
  if (argc != 3) { fprintf (stderr, "usage: <ref.fa> <reads.fa>\n"); return 1; }
  aln = read_alignment_from_file (argv[2]);

  match = new_bwase_match_from_bwa_and_char_vector(argv[1], aln->taxlabel, aln->character, 5, bopt);
  for (i=0; i < match->n_m; i++) 
    printf ("0 %3d ]  read:%5d ref:%5d position:%5d mismatches:%4d gaps:%4d %4d  qual=%3d cigar=%s ref_name=%s\n", i, match->m[i].query_id, match->m[i].ref_id, match->m[i].position, 
            match->m[i].mm, match->m[i].gape, match->m[i].gapo, match->m[i].mapQ, match->m[i].cigar, bwase_match_ref_genome_name (match, i));
  
  del_bwase_match_t (match);
  bopt = new_bwase_options_t (1);
  match = new_bwase_match_from_bwa_and_char_vector(argv[1], aln->taxlabel, aln->character, 5, bopt);
  for (i=0; i < match->n_m; i++) 
    printf ("1 %3d ]  read:%5d ref:%5d position:%5d mismatches:%4d gaps:%4d %4d  qual=%3d cigar=%s ref_name=%s\n", i, match->m[i].query_id, match->m[i].ref_id, match->m[i].position, 
            match->m[i].mm, match->m[i].gape, match->m[i].gapo, match->m[i].mapQ, match->m[i].cigar, bwase_match_ref_genome_name (match, i));
  

  printf ("Now running again but printing SAM file to stdout (last parameter of bwa_aln_bwase())\n");

  bwa_aln_bwase (argv[1], aln->taxlabel->string, aln->character->string, NULL, aln->character->nchars, aln->ntax, 5, &match_list, 1);

  time1 = clock (); printf ("overall time: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  del_alignment (aln);
  del_bwase_match_t (match);
  return TEST_SKIPPED; 
}
