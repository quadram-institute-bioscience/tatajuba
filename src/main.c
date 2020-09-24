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
  struct arg_int  *threads;
  struct arg_file *gff;
  struct arg_file *fna;
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
    .threads = arg_int0("t","nthreads",NULL, "suggested number of threads (default is to let system decide; I may not honour your suggestion btw)"),
    .gff     = arg_file1("g", "gff", "<genome.gff3|genome.gff3.gz>", "reference genome file in GFF3, preferencially with sequence"),
    .fna     = arg_file0(NULL, "fasta", "<genome.fna|genome.fna.gz>", "reference genome file in fasta format, if absent from GFF3"),
    .fastq   = arg_filen(NULL, NULL, "<fastq files>", 1, 0xffff, "fastq file with reads (weirdly, fasta also possible as long as contains all reads and not only contigs)"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.paired, params.kmer, params.minsize, params.minread, params.maxdist, params.leven, 
    params.threads, params.gff, params.fna, params.fastq, params.end};
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
  if (params.threads) free (params.threads);
  if (params.gff)     free (params.gff);
  if (params.fna)     free (params.fna);
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
    printf ("  that is, R1 and R2 should be consecutive. A GFF3 reference genome must also be supplied, and bwa will create a series of index files\n");
    printf ("  if tatajuba can find the DNA sequences at end of GFF3 file, after pragma '##FASTA'\n\n");
    printf ("You can also supply a fasta file, e.g. when you download GFF3 from NCBI it may lack the DNA genome at the end --- then downloading the\n");
    printf ("  associated fna may work (untested). The internal library updates the GFF3 structure with provided fasta sequences (so a fasta file may\n");
    printf ("  overwrite DNA sequences present in the GFF3 file. If you modify the fasta file please delete the index files generated by the BWA library\n");
    printf ("  so that they can be regenerated with the updated information.\n\n");
    printf ("Notice that tatajuba creates files with same prefix and in same location as the GFF3 file, which may overwrite existing ones.\n");
    printf ("  On the other hand, as suggested above, you can recreate all the generated files by deleting them and running tatajuba again.\n\n");
    printf ("The default values are as following:\nkmer=%3d\t minsize=%3d\t minread=%3d\t maxdist=%3d\n",
    params.kmer->ival[0], params.minsize->ival[0], params.minread->ival[0], params.maxdist->ival[0]);
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

tatajuba_options_t
get_options_from_argtable (arg_parameters params)
{
  tatajuba_options_t opt;

  opt.gff = read_gff3_from_file (params.gff->filename[0]);
  if (params.fna->count) {
   size_t len = strlen (params.fna->filename[0]);
   opt.reference_fasta_filename = (char*) biomcmc_malloc ((len+1) * sizeof (char));
   strcpy (opt.reference_fasta_filename, params.fna->filename[0]);
  } 
  else opt.reference_fasta_filename = save_fasta_from_gff3 (opt.gff, NULL, false); // NULL is custom name; false=dont overwrite if file exists 
  if (!opt.reference_fasta_filename) {
    del_gff3_t (opt.gff);
    biomcmc_error ("No fasta provided and GFF3 file doesn't contain sequences\n You must provide a fasta file with "
                   "reference genome sequence(s) that match the GFF3 features, or you should find a GFF3 file with a '##FASTA' section at the end.\n"); 
 
  }

  opt.paired_end = (params.paired->count? true: false);
  opt.kmer_size = params.kmer->ival[0];
  opt.min_tract_size = params.minsize->ival[0];
  opt.min_coverage = params.minread->ival[0]; 
  opt.max_distance_per_flank = params.maxdist->ival[0]; 
  opt.levenshtein_distance = params.leven->ival[0]; 
  if (params.threads->count) {
    opt.n_threads = params.threads->ival[0];
    if (opt.n_threads < 1) opt.n_threads = 1;
#ifdef _OPENMP
    omp_set_num_threads (opt.n_threads); // try to fix n_threads to input number 
#endif
  }
#ifdef _OPENMP
  opt.n_threads = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
#else
  opt.n_threads = 0; // compiled without openMP support (e.g. --disable-openmp)
#endif
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

void
del_tatajuba_options (tatajuba_options_t opt)
{
  if (opt.reference_fasta_filename) free(opt.reference_fasta_filename);
  del_gff3_t (opt.gff);
}

int
main (int argc, char **argv)
{
  int64_t time0[2];
  genome_set_t g;

  biomcmc_get_time (time0);
  arg_parameters params = get_parameters_from_argv (argc, argv);
  tatajuba_options_t opt = get_options_from_argtable (params);
  print_tatajuba_options (opt);
  fprintf (stderr, "Read GFF3 reference genome in %15lf secs\n\n", biomcmc_update_elapsed_time (time0)); 


  g = new_genome_set_from_files (params.fastq->filename, params.fastq->count, opt); 
  print_selected_g_tract_vector (g);

  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to read and generate initial histograms\n", g->secs[0]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to merge and map histograms\n", g->secs[1]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to compare across genomes\n", g->secs[2]);
  biomcmc_fprintf_colour (stderr, 0,2, "Non-threaded timing      ::", " %15lf secs\n\n", biomcmc_update_elapsed_time (time0)); 
  biomcmc_fprintf_fortune (stderr);

  del_genome_set (g);
  del_arg_parameters (params);
  del_tatajuba_options (opt);
  return EXIT_SUCCESS;
}

