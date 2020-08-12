#include <wrapper_bwa.h>
#include <biomcmc.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int main (int argc, char **argv)
{
  clock_t time0, time1;
  alignment aln;

  time0 = clock ();

  if (argc == 1) return TEST_SKIPPED;
  if (argc != 3) { fprintf (stderr, "usage: <ref.fa> <reads.fa>\n"); return 1; }
  aln = read_alignment_from_file (argv[2]);
  bwa_aln_bwase (argv[1], aln->taxlabel->string, aln->character->string, NULL, aln->character->nchars, aln->ntax, 5);
  time1 = clock (); printf ("overall time: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  del_alignment (aln);
  return TEST_SKIPPED; 
}
