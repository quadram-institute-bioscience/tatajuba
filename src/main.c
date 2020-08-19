#include "hopo_set.h"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_lit  *paired;
  struct arg_int  *kmer;
  struct arg_int  *minsize;
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
    .paired  = arg_litn("p","paired", 0, 1, "paired end (pairs of) files"),
    .kmer    = arg_int0("k","kmer","{2,...,15}", "kmer size flanking each side of homopolymer (default=8)"),
    .minsize = arg_int0("m","minsize","{1,...,32}", "minimum homopolymer tract length to be compared"),
    .fastq   = arg_filen(NULL, NULL, NULL, 1, 0xffff, "fasta or fastq file"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.paired, params.kmer, params.minsize, params.fastq, params.end};
  params.argtable = argtable; 
  params.kmer->ival[0]    = 8; // default values must be before parsing
  params.minsize->ival[0] = 3; // default values must be before parsing
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  arg_parse (argc, argv, params.argtable); // returns >0 if errors were found, but this info also on params.end->count
  print_usage (params, argv[0]);
  return params;
}

void
del_arg_parameters (arg_parameters params)
{
  if (params.help)    free (params.help);
  if (params.version) free (params.version);
  if (params.paired)  free (params.paired);
  if (params.kmer)    free (params.kmer);
  if (params.minsize) free (params.minsize);
  if (params.fastq)   free (params.fastq);
  if (params.end)     free (params.end);
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
  printf ("Compare histograms of homopolymeric tract lengths, within context.\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("'Context' is the pair of flanking k-mers. If you have paired end files, then their order should be strictly f1_R1, f1_R2, f2_R1 etc.\n");
    printf ("that is, R1 and R2 should be consecutive. \n");
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  hopo_set hs;
  distance_generator dg;

  time0 = clock ();
  arg_parameters params = get_parameters_from_argv (argc, argv);

  if (params.kmer->ival[0] < 2)  params.kmer->ival[0] = 2;
  if (params.kmer->ival[0] > 16) params.kmer->ival[0] = 16; 
  if (params.minsize->ival[0] < 1)  params.minsize->ival[0] = 1;
  if (params.minsize->ival[0] > 32) params.minsize->ival[0] = 32; 

  hs = new_hopo_set_from_files (params.fastq->filename, params.fastq->count, (bool) params.paired->count, params.kmer->ival[0], params.minsize->ival[0]);
  dg = new_distance_generator_from_hopo_set (hs);

  time1 = clock (); fprintf (stderr, "overall time: %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); time0 = time1; 
  fprintf (stderr, "internal timers::  %lf secs to read, %lf secs to normalise, and %lf secs to compare\n", hs->secs_read, hs->secs_finalise, hs->secs_comparison);

  del_distance_generator (dg);
  del_hopo_set (hs);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}

