/*
 *  YASS 1.15
 *  Copyright (C) 2004-...
 *  the YASS team
 *  Laurent Noe, Gregory Kucherov, Mikhail Roytberg,
 *  Steven Corroy, Antoine De Monte, Christophe Valmir.
 *
 *  laurent.noe|<A>|univ-lille.fr
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the CeCILL License as published by
 *  the CEA-CNRS-INRIA; either version 2 of the License, or (at your
 *  option) any later version, and the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This software contains code derived from the GNU libavl library.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */

#ifndef __ALIGN_H_
#define __ALIGN_H_
#include "tuple.h"
#include "threads.h"

/*
 * align tuplelists that are no longer extended (this create MA),
 * otherwise remove tuples that do not lead to interesting alignments
 */

long int AlignTuples(tuplelist * tuple_list,
                 char *data_query, long int datasize_query,
                 char *data_text, long int datasize_text,
                 Feature *feature
                 );

/*
 * select tuplelists that are no longer extended and call previous function
 */

long int AlignAndFree(
                  char *data_query, long int datasize_query,
                  char *data_text, long int datasize_text,
                  long int force,
                  Feature *feature
                  );

/*
 * previous function, but adapted for one thread
 */

#ifdef THREAD_ASSEMBLE_ALIGN
#if defined(WIN32) || defined(WIN64)
DWORD WINAPI thread_work_align(PVOID fvoid);
#else
void * thread_work_align(void * fvoid);
#endif
#else
void * thread_work_align(void * fvoid);
#endif

#endif
