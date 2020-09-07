#include <wrapper_bwa.h>
#include <biomcmc.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

int main (int argc, char **argv)
{
  clock_t time0, time1;
  int i, n_matches, *match_list = NULL;
  alignment aln;

  time0 = clock ();

  if (argc == 1) return TEST_SKIPPED;
  if (argc != 3) { fprintf (stderr, "usage: <ref.fa> <reads.fa>\n"); return 1; }
  aln = read_alignment_from_file (argv[2]);
  n_matches = bwa_aln_bwase (argv[1], aln->taxlabel->string, aln->character->string, NULL, aln->character->nchars, aln->ntax, 5, &match_list, 0);

  for (i=0; i < n_matches; i++) printf ("read:%5d ref:%5d position:%5d mismathces:%5d gaps:%5d\n",
                                        match_list[5 * i],
                                        match_list[5 * i + 1],
                                        match_list[5 * i + 2],
                                        match_list[5 * i + 3],
                                        match_list[5 * i + 4]);
  printf ("Now running again but printing SAM file to stdout (last parameter of bwa_aln_bwase())\n");
  if (match_list) free (match_list);
  match_list = NULL;

  n_matches = bwa_aln_bwase (argv[1], aln->taxlabel->string, aln->character->string, NULL, aln->character->nchars, aln->ntax, 5, &match_list, 1);

  time1 = clock (); printf ("overall time: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  del_alignment (aln);
  return TEST_SKIPPED; 
}
