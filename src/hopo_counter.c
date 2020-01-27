/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "hopo_counter.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

hopo_counter
new_hopo_counter_from_file (char *filename)
{
  int i;
  gzFile fp = gzopen (filename, "r");
  kseq_t *seq = kseq_init (fp);
  hopo_counter hc = new_hopo_counter ();
  while ((i = kseq_read (seq)) >= 0) update_hopo_counter_from_seq (hc, seq->seq.s); // seq->name.s, seq->seq.l, seq->qual.l
  kseq_destroy(seq);
  gzclose(fp);
  return hc;
}

hopo_coutner
new_hopo_counter (void)
{
  hopo_counter hc = (hopo_counter) biomcmc_malloc (sizeof (struct hopo_counter_struct));
  hc->ref_counter = 1;
  return hc;
}

void
del_hopo_counter (hopo_counter hc)
{
  if (!hc) return;
  if (--hc->ref_counter) return;


  free (hc);
  return;
}
