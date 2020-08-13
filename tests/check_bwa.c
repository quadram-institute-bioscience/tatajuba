#include <wrapper_bwa.h>
#include <biomcmc.h>
#include <check.h>

#define TEST_SUCCESS 0
#define TEST_FAILURE 1
#define TEST_SKIPPED 77
#define TEST_HARDERROR 99

#ifndef TEST_FILE_DIR
#define TEST_FILE_DIR "./files/"
#endif
// use it like memcpy(filename + prefix_size, "bacteria_riboprot.fasta", 23);
char filename[2048] = TEST_FILE_DIR; // now we can memcpy() file names _after_ prefix_size
size_t prefix_size = strlen(TEST_FILE_DIR); // all modifications to filename[] come after prefix_size

START_TEST(bwa_aln_bwase_function)
{
  clock_t time0, time1;
  alignment aln;

  time0 = clock ();
  memcpy (filename + prefix_size, "reads.fa", 9);
  aln = read_alignment_from_file (filename);
  memcpy(filename + prefix_size, "small.fasta", 11);
  bwa_aln_bwase (filename, aln->taxlabel->string, aln->character->string, NULL, aln->character->nchars, aln->ntax, 2);
  time1 = clock (); printf ("  time to read alignment: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  del_alignment (aln);
  if (false) ck_abort_msg ("dummy");
}
END_TEST

START_TEST(bwa_index_function)
{
  clock_t time0, time1;
  char *s;

  time0 = clock ();
  memcpy(filename + prefix_size, "small.fasta", 11);
  s = save_bwa_index (filename, "test", 1); free (s);
  time1 = clock (); printf ("  time to create index: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
  if (false) ck_abort_msg ("dummy");
}
END_TEST

Suite * bwa_suite(void)
{
  Suite *s;
  TCase *tc_case;

  s = suite_create("BWA bindings");
  tc_case = tcase_create("BWA simple");
  tcase_add_test(tc_case, bwa_aln_bwase_function);
  tcase_add_test(tc_case, bwa_index_function);
  suite_add_tcase(s, tc_case);
  return s;
}

int main(void)
{
  int number_failed;
  SRunner *sr;

  sr = srunner_create (bwa_suite());
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed > 0) ? TEST_FAILURE:TEST_SUCCESS;
}
