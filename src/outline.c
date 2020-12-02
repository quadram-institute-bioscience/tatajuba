
genome_set_t new_genome_set_from_files (const char **filenames, int n_filenames, tatajuba_options_t opt) {
  hc = new_or_append_hopo_counter_from_file (hc,   filenames[2*i+1], opt);
  g->genome[i] = new_genomic_context_list (hc);
  g->tract = new_g_tract_vector_from_genomic_context_list (g->genome, g->n_genome);
  create_tract_in_reference_structure (g);
}


hopo_counter new_or_append_hopo_counter_from_file (hopo_counter hc, const char *filename, tatajuba_options_t opt) {
  update_hopo_counter_from_seq (hc_local, seq->seq.s, seq->seq.l, opt.min_tract_size); 
  --> add_kmer_to_hopo_counter (hc, context, hopo_base_int, count_same);
}

genomic_context_list_t new_genomic_context_list (hopo_counter hc)  { 
  finalise_hopo_counter (hc);
