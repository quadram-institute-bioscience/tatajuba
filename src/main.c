#include "hopo_counter.h"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_file *fastq;
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
    .fastq   = arg_filen(NULL, NULL, NULL, 1, 0xffff, "fasta or fastq file"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.fastq, params.end};
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
  if (params.fastq) free (params.fastq);
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
  hopo_counter hc;

  time0 = clock ();
  arg_parameters params = get_parameters_from_argv (argc, argv);

  for (i = 0; i < params.fastq->count; i++) {
    hc = new_hopo_counter_from_file (params.fastq->filename[i]);
//    printf ("\t%s\n", params.fastq->filename[i]);
    print_debug_hopo_counter (hc);
    del_hopo_counter (hc);
  }
  time1 = clock (); fprintf (stderr, "read in  %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); 

  del_arg_parameters (params);
  return EXIT_SUCCESS;
}

