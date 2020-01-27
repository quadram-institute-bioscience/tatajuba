#include "tatajuba.h"

KSEQ_INIT(gzFile, gzread)

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_file *fasta;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .fasta   = arg_filen(NULL, NULL, NULL, 1, 1, "fasta file"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.fasta, params.end};
  params.argtable = argtable; 
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  arg_parse (argc, argv, params.argtable); // returns >0 if errors were found, but this info also on params.end->count
  print_usage (params, argv[0]);
  return params;
}

void
del_arg_parameters (arg_parameters params)
{
  if (params.help) free (params.help);
  if (params.version) free (params.version);
  if (params.fasta) free (params.fasta);
  if (params.end) free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return;

  if (params.end->count) {  // params.end holds error messages
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("Error when reading arguments from command line\n\n");
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("Just testing at the moment. Move along, nothing to see here!\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("I know you need help, but I'm useless\n");
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i;
  clock_t time0, time1;
  char_vector seqname = new_char_vector (1);
  char_vector dna = new_char_vector (1);
  char_vector align = NULL;

  time0 = clock ();
  arg_parameters params = get_parameters_from_argv (argc, argv);

  gzFile fp = gzopen((char*) params.fasta->filename[0], "r");
  kseq_t *seq = kseq_init(fp);
  while ((i = kseq_read(seq)) >= 0) {
    char_vector_add_string (seqname, seq->name.s);
    char_vector_add_string (dna, seq->seq.s);
  }
  kseq_destroy(seq);
  gzclose(fp);
  time1 = clock (); fprintf (stderr, "read in  %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); 

  align = kalign3_from_char_vector (dna);

  time1 = clock (); fprintf (stderr, "finished in  %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); 
  for (i= 0; i < align->nstrings; i++) printf (">%s\n%s\n", seqname->string[i], align->string[i]);

  del_char_vector (dna);
  del_char_vector (align);
  del_char_vector (seqname);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}

