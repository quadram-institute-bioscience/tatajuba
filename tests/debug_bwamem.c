#include <wrapper_bwa.h>
#include <biomcmc.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int main (int argc, char **argv)
{
  clock_t time0, time1;
  int i;
  bwmem_match_t match;
  alignment aln;

  time0 = clock ();

  if (argc == 1) return TEST_SKIPPED;
  if (argc != 3) { fprintf (stderr, "usage: <ref.fa> <reads.fa>\n"); return 1; }
  aln = read_alignment_from_file (argv[2]);
  match = new_bwmem_match_from_bwa_and_char_vector(argv[1], aln->character);

  for (i=0; i < match->n_m; i++) 
    printf ("read:%5d ref:%5d position:%7ld score:%4d edist:%4d mapQ:%4d  primary:%1d cigar=%16s ref_name=%s\n", 
            match->m[i].query_id, match->m[i].ref_id, match->m[i].position, match->m[i].score,
            match->m[i].edit_distance, match->m[i].mapQ, match->m[i].is_primary, match->m[i].cigar, bwmem_match_ref_genome_name (match, i));

  time1 = clock (); printf ("overall time: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  del_alignment (aln);
  del_bwmem_match_t (match);
  return TEST_SKIPPED; 
}
