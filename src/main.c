#ifndef _tatajuba_main_c_
#define _tatajuba_main_c_
#include "analyse_variable_tracts.h"
#endif

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_lit  *paired;
  struct arg_lit  *keepbias;
  struct arg_lit  *vcf;
  struct arg_int  *kmer;
  struct arg_int  *minsize;
  struct arg_int  *minread;
  struct arg_int  *maxdist;
  struct arg_int  *leven;
  struct arg_int  *threads;
  struct arg_file *gff;
  struct arg_file *fna;
  struct arg_file *fastq;
  struct arg_file *outdir;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);
char* string_for_outdir (const char *dname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help     = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version  = arg_litn("v","version",0, 1, "print version and exit"),
    .paired   = arg_litn("p","paired", 0, 1, "paired end (pairs of) files"),
    .keepbias = arg_litn("b","keep_bias", 0, 1, "keep biased tracts, i.e. present only in reverse or only in forward strains (default=remove)"),
    .vcf      = arg_litn("V","vcf", 0, 1, "generate VCF files for each sample, around the HT regions (EXPERIMENTAL) (default=not to save)"),
    .kmer     = arg_int0("k","kmer","{2,...,32}", "kmer size flanking each side of homopolymer (default=25)"),
    .minsize  = arg_int0("m","minsize","{1,...,32}", "minimum homopolymer tract length to be compared (default=4)"),
    .minread  = arg_int0("i","minreads",NULL, "minimum number of reads for tract+context to be considered (default=5)"),
    .maxdist  = arg_int0("d","maxdist",NULL, "maximum distance between kmers of a flanking region to merge them into one context (default=1)"),
    .leven    = arg_int0("l","leven",NULL, "levenshtein distance between flanking regions to merge them into one context (after ref genome mapping)"),
    .threads  = arg_int0("t","nthreads",NULL, "suggested number of threads (default is to let system decide; I may not honour your suggestion btw)"),
    .gff      = arg_file1("g","gff", "<genome.gff3|genome.gff3.gz>", "reference genome file in GFF3, preferencially with sequence"),
    .fna      = arg_file0("f","fasta", "<genome.fna|genome.fna.gz>", "reference genome file in fasta format, if absent from GFF3"),
    .fastq    = arg_filen(NULL, NULL, "<fastq files>", 1, 0xffff, "fastq file with reads (weirdly, fasta also possible as long as contains all reads and not only contigs)"),
    .outdir   = arg_file0("o", "outdir", NULL, "output directory, or 'random' for generating random dir name (default=current dir '.')"),
    .end      = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.paired, params.keepbias, params.vcf, params.kmer, params.minsize, params.minread, params.maxdist, params.leven, 
    params.threads, params.gff, params.fna, params.fastq, params.outdir, params.end};
  params.argtable = argtable; 
  params.kmer->ival[0]    = 25; // default values must be before parsing
  params.minsize->ival[0] = 4;
  params.minread->ival[0] = 5;
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
  if (params.keepbias)free (params.keepbias);
  if (params.vcf)     free (params.vcf);
  if (params.kmer)    free (params.kmer);
  if (params.minsize) free (params.minsize);
  if (params.minread) free (params.minread);
  if (params.maxdist) free (params.maxdist);
  if (params.leven)   free (params.leven);
  if (params.threads) free (params.threads);
  if (params.gff)     free (params.gff);
  if (params.fna)     free (params.fna);
  if (params.fastq)   free (params.fastq);
  if (params.outdir)  free (params.outdir);
  if (params.end)     free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return;

  if (params.end->count && (!params.help->count)) {  // params.end holds error messages
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
    printf ("  associated fna may work. The internal library updates the GFF3 structure with provided fasta sequences (so a fasta file may\n");
    printf ("  overwrite DNA sequences present in the GFF3 file. If you modify the fasta file please delete the index files generated by the BWA library\n");
    printf ("  so that they can be regenerated with the updated information.\n\n");
    printf ("Notice that tatajuba creates files with same prefix and in same location as the GFF3 file, which may overwrite existing ones.\n");
    printf ("  On the other hand, as suggested above, you can recreate all the generated files by deleting them and running tatajuba again.\n\n");
    printf ("The generation of the VCF files is still experimental, but seems to be accepted by snpEff. \n\n");
    printf ("The default values are as following:\nkmer=%3d\t minsize=%3d\t minread=%3d\t maxdist=%3d\n",
    params.kmer->ival[0], params.minsize->ival[0], params.minread->ival[0], params.maxdist->ival[0]);
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

tatajuba_options_t
get_options_from_argtable (arg_parameters params)
{
  tatajuba_options_t opt;

  opt.outdir = NULL;
  if (params.outdir->count) opt.outdir = string_for_outdir (params.outdir->filename[0]);
  else opt.outdir = string_for_outdir ("./");

  opt.gff = read_gff3_from_file (params.gff->filename[0]);
  if (params.fna->count) {
   size_t len = strlen (params.fna->filename[0]);
   opt.reference_fasta_filename = (char*) biomcmc_malloc ((len+1) * sizeof (char));
   strcpy (opt.reference_fasta_filename, params.fna->filename[0]);
   alignment aln = read_fasta_alignment_from_file (opt.reference_fasta_filename, false); 
   add_fasta_to_gff3 (opt.gff, aln->taxlabel, aln->character);
   del_alignment (aln);
  } 
  else opt.reference_fasta_filename = save_fasta_from_gff3 (opt.gff, NULL, false); // NULL is custom name; false=dont overwrite if file exists 
  if (!opt.reference_fasta_filename) {
    del_gff3_t (opt.gff);
    biomcmc_error ("No fasta provided and GFF3 file doesn't contain sequences\n You must provide a fasta file with "
                   "reference genome sequence(s) that match the GFF3 features, or you should find a GFF3 file with a '##FASTA' section at the end.\n"); 
  }

  // create index here to make sure it is done only once (not inside threads)
  char *s = save_bwa_index (opt.reference_fasta_filename, NULL, false); // returns file name with suffix (NULL in our case)
  if (s) free (s);

  opt.paired_end = (params.paired->count? true: false);
  opt.save_vcf   = (params.vcf->count? true: false);
  opt.remove_biased = (params.keepbias->count? false: true); // if set, then do not remove bias; default is to remove
  opt.n_samples = (params.paired->count? params.fastq->count/2: params.fastq->count);
  if (opt.n_samples < 2) {
    biomcmc_warning ("More than one sample is needed, since this program is based on differences between samples; Proceed at your own peril.");
  }
  opt.kmer_size = params.kmer->ival[0];
  opt.min_tract_size = params.minsize->ival[0];
  opt.min_coverage = params.minread->ival[0]; 
  opt.max_distance_per_flank = params.maxdist->ival[0]; 
  opt.levenshtein_distance = params.leven->ival[0]; 

#ifdef _OPENMP
  opt.n_threads = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
#else
  opt.n_threads = 0; // compiled without openMP support (e.g. --disable-openmp)
  biomcmc_warning ("Program compiled without multithread support");
#endif

#ifdef _OPENMP
  if (params.threads->count) {
    opt.n_threads = params.threads->ival[0];
    if (opt.n_threads < 1) opt.n_threads = 1;
  }
  if (opt.n_samples < opt.n_threads) { 
    opt.n_threads = opt.n_samples; 
    biomcmc_warning ("Decreasing number of threads to match number of samples");
  }
  omp_set_num_threads (opt.n_threads); // try to fix n_threads to input number 
  opt.n_threads = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
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
  if (opt.outdir) free (opt.outdir); 
  return;
}

char*
string_for_outdir (const char *dname)
{
  char *s = NULL;
  if (!dname) {
    s = biomcmc_malloc (sizeof (char) * 3);
    s[0] = '\0';
    strcpy (s, "./");
    s[2] = '\0';
  }
  else if ((strstr (dname, "random") != NULL) && (strlen (dname) < 8)) { // "randomiser" is fine 
    uint32_t r3 = biomcmc_rng_get () & 0xffff;
    uint32_t r2 = getpid  () & 0xff; // PID is usually small number
    uint32_t r1 = getppid () & 0xff; // parent PID
    s = biomcmc_malloc (sizeof (char) * 13);
    s[0] = '\0';
    sprintf (s, "out%x%x%x/", r1, r2, r3);
    s[12] = '\0';
  }
  else {
    size_t len = strlen (dname);
    if (dname[len-1] != '/') len++; // added in any case
    s = biomcmc_malloc (sizeof (char) * (len+1));
    s[0] = '\0';
    strcpy (s, dname);
    s[len-1] = '/';
    s[len]   = '\0';
  }

  if ((strstr (s, "./") != NULL) && (strlen (s) < 4)) return s; // no need to mkdir()

  if (mkdir((const char *) s, 0777)) { // zero on success, negative if error
    if (errno == EEXIST) biomcmc_warning ("File '%s' already exists; I'll assume it's a directory. Contents will be overwritten", s);
    else biomcmc_error ("Problem trying to create directory '%s'; error description:\n%s", s, strerror (errno));
  }
  return s;
}

int
main (int argc, char **argv)
{
  int64_t time0[2];
  genome_set_t g;

  biomcmc_get_time (time0);
  biomcmc_random_number_init (0);
  arg_parameters params = get_parameters_from_argv (argc, argv);
  tatajuba_options_t opt = get_options_from_argtable (params);
  print_tatajuba_options (opt);
  fprintf (stderr, "Read GFF3 reference genome in %15lf secs\n\n", biomcmc_update_elapsed_time (time0)); 

  g = new_genome_set_from_files (params.fastq->filename, params.fastq->count, opt); 
  print_selected_g_tract_vector (g);
  if (opt.save_vcf) generate_vcf_files (g);
  print_tract_list (g);

  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to read and generate initial histograms\n", g->secs[0]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to merge and map histograms\n", g->secs[1]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to compare across sample genomes\n", g->secs[2]);
  biomcmc_fprintf_colour (stderr, 0,2, "Internal (threaded) timer::", " %15lf secs to compare with reference\n", g->secs[3]);
  biomcmc_fprintf_colour (stderr, 0,2, "Non-threaded timing      ::", " %15lf secs\n\n", biomcmc_update_elapsed_time (time0)); 
  biomcmc_fprintf_fortune (stderr);

  del_genome_set (g);
  del_arg_parameters (params);
  del_tatajuba_options (opt);
  biomcmc_random_number_finalize ();
  return EXIT_SUCCESS;
}

