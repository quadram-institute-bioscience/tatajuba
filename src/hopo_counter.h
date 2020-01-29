/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

/*! \file
 *  \brief homopolymer counter
 */

#ifndef _hopo_counter_h_
#define _hopo_counter_h_

#include <kalign.h>

typedef struct hopo_counter_struct* hopo_counter;

struct hopo_counter_struct
{
  int ref_counter;
};

hopo_counter new_hopo_counter_from_file (const char *filename);
void del_hopo_counter (hopo_counter hc);


#endif
