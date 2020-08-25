#include "genome_set.h"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_lit  *paired;
  struct arg_int  *kmer;
  struct arg_int  *minsize;
  struct arg_int  *minread;
  struct arg_int  *maxdist;
  struct arg_int  *leven;
  struct arg_file *ref;
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
    .kmer    = arg_int0("k","kmer","{2,...,32}", "kmer size flanking each side of homopolymer (default=8)"),
    .minsize = arg_int0("m","minsize","{1,...,32}", "minimum homopolymer tract length to be compared"),
    .minread = arg_int0("i","minreads",NULL, "minimum number of reads for tract+context to be considered"),
    .maxdist = arg_int0("d","maxdist",NULL, "maximum distance between kmers of a flanking region to merge them into one context"),
    .leven   = arg_int0("l","leven",NULL, "levenshtein distance between flanking regions to merge them into one context (after ref genome mapping)"),
    .ref     = arg_file1("r", "ref", "<genome.fa|genome.fa.gz>", "reference genome file (bwa indices will be created if absent)"),
    .fastq   = arg_filen(NULL, NULL, "<fastq files>", 1, 0xffff, "fasta or fastq file with reads"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.paired, params.kmer, params.minsize, params.minread, params.maxdist, params.leven, params.ref, params.fastq, params.end};
  params.argtable = argtable; 
  params.kmer->ival[0]    = 16; // default values must be before parsing
  params.minsize->ival[0] = 4; 
  params.minread->ival[0] = 3;
  params.maxdist->ival[0] = 1;
  params.leven->ival[0] = -1;
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  if (arg_parse (argc, argv, params.argtable)) print_usage (params, argv[0]); // returns >0 if errors were found, info also on params.end->count
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
  if (params.minread) free (params.minread);
  if (params.maxdist) free (params.maxdist);
  if (params.leven)   free (params.leven);
  if (params.ref)     free (params.ref);
  if (params.fastq)   free (params.fastq);
  if (params.end)     free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return;

  if (params.end->count) {  // params.end holds error messages
    biomcmc_fprintf_colour (stdout, 0,1, "Error when reading arguments from command line:\n", NULL);
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("\n");
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("Compare histograms of homopolymeric tract lengths, within context.\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("'Context' is the pair of flanking k-mers. If you have paired end files, then their order should be strictly f1_R1, f1_R2, f2_R1 etc.\n");
    printf ("that is, R1 and R2 should be consecutive. A fasta reference genome must also be supplied, and bwa will create a series of index files. \n\n");
    printf ("The default values are as following:\nkmer=%3d\t minsize=%3d\t minread=%3d\t maxdist=%3d\n",
    params.kmer->ival[0], params.minsize->ival[0], params.minread->ival[0], params.maxdist->ival[0]);
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

tatajuba_options_t
get_options_from_argtable (arg_parameters params)
{
  tatajuba_options_t opt;
  int len = strlen (params.ref->filename[0]);
  opt.reference_genome_filename = (char*) biomcmc_malloc ((len+1) * sizeof (char));
  strcpy (opt.reference_genome_filename, params.ref->filename[0]);
  opt.paired_end = (params.paired->count? true: false);
  opt.kmer_size = params.kmer->ival[0];
  opt.min_tract_size = params.minsize->ival[0];
  opt.min_coverage = params.minread->ival[0]; 
  opt.max_distance_per_flank = params.maxdist->ival[0]; 
  opt.levenshtein_distance = params.leven->ival[0]; 

  if (opt.kmer_size < 2)  opt.kmer_size = 2; 
  if (opt.kmer_size > 32) opt.kmer_size = 32; 
  if (opt.min_tract_size < 1)  opt.min_tract_size = 1; 
  if (opt.min_tract_size > 32) opt.min_tract_size = 32; 
  if (opt.min_coverage < 0)      opt.min_coverage = 0; 
  if (opt.min_coverage > 0xffff) opt.min_coverage = 0xffff; 
  if (opt.max_distance_per_flank < 0) opt.max_distance_per_flank = 0; 
  if (opt.max_distance_per_flank > opt.kmer_size/2) opt.max_distance_per_flank = opt.kmer_size/2; 
  if (opt.levenshtein_distance < 0) opt.levenshtein_distance = opt.max_distance_per_flank + 1; 
  return opt;
}

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  genome_set_t g;

  time0 = clock ();
  arg_parameters params = get_parameters_from_argv (argc, argv);
  tatajuba_options_t opt = get_options_from_argtable (params);
  print_tatajuba_options (opt);

  g = new_genome_set_from_files (params.fastq->filename, params.fastq->count, opt); 
  print_selected_g_tract_vector (g);

  time1 = clock (); 
  biomcmc_fprintf_colour (stderr, 0,2, "Internal timer::", "%9lf secs to read and generate initial histograms\n", g->secs[0]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal timer::", "%9lf secs to merge and map histograms\n", g->secs[1]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal timer::", "%9lf secs to compare across genomes\n", g->secs[2]);
  biomcmc_fprintf_colour (stderr, 0,2, "Overall timer ::", "%9lf secs\n\n", (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); 
  biomcmc_fprintf_fortune (stderr);

  del_genome_set (g);
  del_arg_parameters (params);
  if (opt.reference_genome_filename) free(opt.reference_genome_filename);
  return EXIT_SUCCESS;
}

