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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* 1) include utils macro */
#include "util.h"
/* 2) include global variables */
#include "global_var.h"
/* 3) the current file defs */
#include "assemble.h"

/* 4) other files */
#include "kword.h"
#include "tuple.h"
#include "prdyn.h"
#include "proba.h"
#include "align.h"
#include "display.h"

/*
 * diagonal localy good hash table
 */

#ifdef  CACHE
#define HASH_SHIFT          (7)
#define HASH_DIAG(diagonal) ((diagonal)>>(HASH_SHIFT))
#endif

/*
 * Modulo access to the last tuple
 */

#ifdef MEMOPT
  /* query table index (little more time consuming but memory economic ) */
  #define MEMOPT_SETMOD   { long int i=0; memopt_modmask = 0; \
        for (i=0; (                                           \
                   (ALIGN_EVERY_NB_ITER + 1 +                 \
                    keysize_query       + 1 +                 \
                    gp_rho_stat         + 1 +                 \
                    gp_delta_stat       + 1                   \
                    )                                         \
                    >= memopt_modmask)                        \
                  ;i++)                                       \
                 memopt_modmask |=  (1 << i);                 \
  }

  #define MEMOPT_MOD      &(memopt_modmask)
  #define MEMOPT_TABSIZE   (memopt_modmask+1)
#else
  /* full table index (memory consuming) */
  #define MEMOPT_SETMOD
  #define MEMOPT_MOD
  #define MEMOPT_TABSIZE   (keysize_query + keysize_text)
#endif

/*
 * Shift tables used in assemble algorithms
 */


long int initialise_deltashift()
{
  long int i, s=1;
  if (gv_delta_shift == NULL) {
    gv_delta_shift = (long int * ) MALLOC ((2+2*gp_delta_stat)*sizeof(long int));
    ASSERT(gv_delta_shift, "initialise_deltashift");

    for (i = 0; i <= 1 + 2 * gp_delta_stat; i++) {
      gv_delta_shift[i] = s * (i+1) / 2;
      s *= -1;
    }
  }
  return 0;
}

/*
 *  definition des fonctions de la table hash ou tableaux
 */
#define INIT_TAB(table, type, size, set)    {                                 \
     table =  (type) MALLOC(size);                                            \
     ASSERT(table, errmessage);                                               \
     memset(table, set, size);                                                \
}

#define GET_TAB_MIN(table, index, minValue) ((table)[(index) MEMOPT_MOD])
#define GET_TAB(table, index)               ((table)[(index) MEMOPT_MOD])
#define EXIST_TAB(table, index)             ((table)[(index) MEMOPT_MOD])
#define PUT_TAB(table, index, value)        ((table)[(index) MEMOPT_MOD] = value)
#define DEL_TAB(table, index, value)        ((table)[(index) MEMOPT_MOD] = value)
#define FREE_TAB(table, sizeelement)        {FREE(table, (MEMOPT_TABSIZE)*(sizeelement));}
#define RESET_TAB(table, size, set)         {memset(table, set, (size));}




long int Assemble_Single(/*in*/ char * data, /*in*/ long int datasize, Feature * f) {

  /* main variables
   * for the linked list computation
   */
  long int * keylist[MAX_SEED];
  long int   keylist_size[MAX_SEED];
  long int   keysize_query  = datasize;
#ifdef MEMOPT
  long int   memopt_modmask = 0;
#else
  long int   keysize_text   = 0;
#endif

  long int * last_tuple_pos_with_diag = NULL;
  tuple   ** last_tuple_ptr_with_diag = NULL;

#ifdef CACHE
  long int * last_tuple_pos_with_diag_hash = NULL;
#endif

#ifdef STATS
  f->last_clock = clock();
#endif

  /* to keep lists of tuple found */
  if (!f->first_tl)
    f->last_tl = f->first_tl = CreateTupleList();

  /* [0] left correction factor set to zero */
  f->left_correction =  0;

  MEMOPT_SETMOD;

  /* [1] key linked list creation */
  {
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++)
      CreateKeyList(data, datasize, seed, &(keylist[seed]), &(keylist_size[seed]), f);
  }

  /* [2] allocations */
  INIT_TAB(last_tuple_pos_with_diag, long int*, (MEMOPT_TABSIZE)*sizeof(long int), 0xff);
  INIT_TAB(last_tuple_ptr_with_diag, tuple**, (MEMOPT_TABSIZE)*sizeof(tuple*), 0x00);

#ifdef CACHE
  last_tuple_pos_with_diag_hash = (long int *) MALLOC((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
  ASSERT(last_tuple_pos_with_diag_hash, "Assemble");
  memset(last_tuple_pos_with_diag_hash, 0x7f, (1 + HASH_DIAG(MEMOPT_TABSIZE) * sizeof(long int));
#endif

  STATS_ADD_CLOCK(f, clock_pre);

  /* [3] algorithm to create linked tuples */

  /* for each pos */
  for (f->i_current = 0; f->i_current < datasize; f->i_current++) {
    /* for each seed */
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++) {

      if (f->i_current <= datasize - gp_seeds_span[seed]) {
        long int i_current_end;
        long int code = 0;
        KEY(code, i_current_end, data, f->i_current, seed); /* we take the code of the kword 'w' located at f->i_current position : i_current_end is the last pos of the seed */
        if (code >= 0) {

#ifdef KEYLISTCOMPRESS
          long int index      = f->keycode_first[seed][code];  /* search the first occurrence of 'w' on the sequence */
          long int index_end  = f->keycode_first[seed][code+1]; /* warning : undefined on "fast" first implemetation */
          long int i_previous = -1;
          if (index >= 0) {
            i_previous = keylist[seed][index];
          }

          /* for all the occurrences of 'w' before i_current_end */
          while (i_previous < i_current_end - gp_distdiag && i_previous >= 0 /* for "fast" only */ && index != index_end /* for "full" only */) {
#else
          long int index      = f->keycode_first[seed][code];  /* search the first occurrence of 'w' on the sequence */
          long int i_previous = index;

          /* consider all the occurrences of 'w' before i_current */
          while (i_previous < i_current_end - gp_distdiag && index >= 0) {
#endif

            long int diagonal = i_current_end - i_previous;

#ifdef PREFETCH
            __builtin_prefetch(&GET_TAB(last_tuple_pos_with_diag, diagonal), 1, 3);
            __builtin_prefetch(&GET_TAB(last_tuple_ptr_with_diag, diagonal), 1, 3);
#endif

            /* consider all the interesting diagonals by adding indel bias  0 +1 -1 +2...
             * until we reach a new tuple or the upper indel bound */

            long int observed_diagonal = diagonal;
            long int delta_diagonal = 0;

#ifdef CACHE
            {
              long int i_hash     = HASH_DIAG(MAX((diagonal - gp_delta_stat), 0)MEMOPT_MOD);
              long int i_hash_end = HASH_DIAG(MIN((diagonal + gp_delta_stat), MEMOPT_TABSIZE)MEMOPT_MOD);
              do {
                  long int i_previous_hash = last_tuple_pos_with_diag_hash[i_hash];
                  if ((i_previous_hash >= i_current_end - gp_rho_stat) &&
                      (i_previous_hash <= i_current_end) &&
                      (i_previous_hash - observed_diagonal <= i_current_end - diagonal) &&
                      (i_previous_hash >= 0)
                      )
                    goto diag_start;
                  if (i_hash  == i_hash_end)
                    goto maj_last_tuple_ptr_with_diag;
                  i_hash++;
                  i_hash %= HASH_DIAG(MEMOPT_TABSIZE);
                }  while (1);

            }
          diag_start:
#endif

            while (delta_diagonal < 2 * gp_delta_stat + 1) {

              /* search for interesting tuples */
              long int j = GET_TAB_MIN(last_tuple_pos_with_diag, observed_diagonal, i_current_end - gp_rho_stat);

              /* if tuples found are respecting criteria (neighboors on diagonal and on sequences) */
              if (j >= 0                                                           &&
#ifdef OLDSELECTION
                  j                     <= i_current_end                           &&
#endif
                  j                     >= i_current_end            - gp_rho_stat  &&
#ifdef OLDSELECTION
                  j - observed_diagonal <= i_current_end - diagonal                &&
#endif
                  j - observed_diagonal >= i_current_end - diagonal - gp_rho_stat   ) {

                tuple * t, * t_;

                /*---------------------------------------------------------------------*/

                /* (1) possibly create one instance of tuple [j]
                 * if this one was not in any list.
                 */

                if (EXIST_TAB(last_tuple_ptr_with_diag, observed_diagonal) == 0) {

                  STATS_NB_CHAINS_BUILT_INC(f);

                  /* create a tuple [j] with a small seed and update lists */
                  CREATETUPLE(t, j, observed_diagonal, gp_seeds_span_min);
                  /* create a tuple list element */
                  CREATETUPLELIST(f->last_tl->next, t, NULL);
                  f->last_tl = f->last_tl->next;

                  PUT_TAB(last_tuple_ptr_with_diag, observed_diagonal, t);
                } else {

                  /* (2) we search for the last tuple [j] instanciation
                   */
                  t = (tuple *)GET_TAB(last_tuple_ptr_with_diag, observed_diagonal);
                }

                /*
                 * (3) iterate tuples and add this new one to the most adapted one
                 */

                do {
                  long int advance;
                  if (diagonal == t->diagonal  && ((advance = i_current_end - t->occurrence) <= gp_seeds_span[seed] + 1)) {

                    /* (3a) possible fusion with this tuple (no instanciation)
                     *
                     * [inter-tuple overlap]
                     */

                    if (advance >= 0) {
                      t->leftsize          += advance;
                      t->occurrence         = i_current_end;
                      goto no_maj_last_tuple_ptr_with_diag;
                    } else {
                      goto next_key;
                    }
                  }

                  /* list search for the most fitted */
                  if (!(t->next))
                    break;
                  t = t->next;

                } while (1);

                /*
                 * (3d) any fusion has not been possible
                 *
                 * we instanciate the tuple and add this one at the end of the chain.
                 *
                 */

                CREATETUPLE(t_, i_current_end, diagonal, gp_seeds_span[seed]);
                t->next = t_;
                PUT_TAB(last_tuple_ptr_with_diag, diagonal, t_);

                /*
                 * (4) dont reset last_tuple_ptr_with_diag[diagonal] to NULL
                 * because it has been modified in (3)
                 */
                goto no_maj_last_tuple_ptr_with_diag;

              }

              /*---------------------------------------------------------------------*/

              /* diagonal walk d, d+1, d-1, ... */
              delta_diagonal++;
              observed_diagonal = diagonal + gv_delta_shift[delta_diagonal];
              observed_diagonal = MAX(MIN(observed_diagonal, f->i_current), 0);

              /* border effects avoided*/

            } /* while : diagonal walk */


            if (gp_hitcriterion == 1) {
              SINGLEHITDIAGONAL(data, datasize, data, datasize);
            }

            /* reset last tuple instancied to NULL (because not instance is created)
             */

#ifdef CACHE
          maj_last_tuple_ptr_with_diag:
#endif
            DEL_TAB(last_tuple_ptr_with_diag, diagonal, NULL);

          no_maj_last_tuple_ptr_with_diag:

            /*
             * update last tuple position  with diagonal = 'diagonal'
             */
            PUT_TAB(last_tuple_pos_with_diag, diagonal, i_current_end);

#ifdef CACHE
            last_tuple_pos_with_diag_hash[HASH_DIAG(diagonal MEMOPT_MOD)] = i_current_end;
#endif
            /* iterate 'w' code list before i_current_end position */
          next_key:
            STATS_NB_SEEDS_INC(f);

#ifdef KEYLISTCOMPRESS
            index ++;
            i_previous = keylist[seed][index];
#else
            i_previous = index = keylist[seed][index];
#endif
          } /* while i_previous */
        } /* if (code >= 0) */
      } /* if (f->i_current <= datasize - gp_seeds_span[seed]) */
    } /* for each seed */


    if (!(f->i_current % ALIGN_EVERY_NB_ITER)) {
#ifdef TRACE
      Display_Progress(f->i_current, (keysize_query), f);
#endif

      STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
      DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
      if (f->thread_align) {
        WAIT_THREAD(f->thread_align);
        f->thread_align = 0;
      }
      if (f->i_current) {
        /* dont start at very first time :
         * - avoid too small sequences being threaded
         * - give some elements to be processed to the work_align
         */
        CREATE_THREAD(f->thread_align, thread_work_align, f);
      }
#else
      AlignAndFree(
                   data, datasize,
                   data, datasize,
                   0,
                   f);
#endif

      STATS_ADD_CLOCK(f, clock_align);

    }
  } /* for each pos */

#ifdef TRACE
  Display_Progress(f->i_current, (keysize_query), f);
#endif

  STATS_ADD_CLOCK(f, clock_chain);

  /* aligning remaining tuples */
#ifdef DEBUG_ASSEMBLE
  DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
  if (f->thread_align) {
    WAIT_THREAD(f->thread_align);
    f->thread_align = 0;
  }
#endif

  /*
   * CREATE_THREAD(f->thread_align, thread_work_align, f);
   */

  AlignAndFree(
               data, datasize,
               data, datasize,
               1,
               f);

  STATS_ADD_CLOCK(f, clock_align);

#ifdef CACHE
  FREE(last_tuple_pos_with_diag_hash, (1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif

  FREE_TAB(last_tuple_pos_with_diag, sizeof(long int));
  FREE_TAB(last_tuple_ptr_with_diag, sizeof(tuple*));
  {
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++) {
      FREE(keylist[seed], keylist_size[seed]);
      /* FREE(first[seed],   first_size[seed]);*/
    }
  }

  STATS_ADD_CLOCK(f, clock_pre);

  return 0;
}


long int Assemble_SingleRev(/*in*/   char * datarev, /*in*/  char * data, /*in*/ long int datasize, Feature  *f) {

  /* main variables
   * for the linked list computation
   */
  long int * keylist[MAX_SEED];
  long int   keylist_size[MAX_SEED];
  long int   keysize_query  = datasize;
#ifdef MEMOPT
  long int   memopt_modmask = 0;
#else
  long int   keysize_text   = keysize_query;
#endif

  long int * last_tuple_pos_with_diag = NULL;
  tuple   ** last_tuple_ptr_with_diag = NULL;

#ifdef CACHE
  long int * last_tuple_pos_with_diag_hash = NULL;
#endif

#ifdef STATS
  f->last_clock = clock();
#endif

  /* to keep lists of tuple found */
  if (!f->first_tl)
    f->last_tl = f->first_tl = CreateTupleList();

  /* [0] left correction factor set to keysize_query */
  f->left_correction = keysize_query;

  MEMOPT_SETMOD;

  /* [1] key linked list creation */
  {
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++)
      CreateKeyList(datarev, datasize, seed, &(keylist[seed]), &(keylist_size[seed]), f);
  }

  /* [2] allocations */
  INIT_TAB(last_tuple_pos_with_diag, long int *, (MEMOPT_TABSIZE) * sizeof(long int), 0xff);
  INIT_TAB(last_tuple_ptr_with_diag, tuple **, (MEMOPT_TABSIZE) * sizeof(tuple*), 0x00);

#ifdef CACHE
  last_tuple_pos_with_diag_hash = (long int *) MALLOC((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int)));
  ASSERT(last_tuple_pos_with_diag_hash, "Assemble");
  memset(last_tuple_pos_with_diag_hash, 0x7f, (1 + HASH_DIAG(MEMOPT_TABSIZE) * sizeof(long int)));
#endif

  STATS_ADD_CLOCK(f, clock_pre);

  /* [3] algorithm to create linked tuples */

  /* for each pos */
  for (f->i_current = 0; f->i_current < datasize; f->i_current++) {
    /* for each seed */
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++) {

      if (f->i_current <= datasize - gp_seeds_span[seed]) {
        long int i_current_end;
        long int code = 0;
        KEY(code, i_current_end, data, f->i_current, seed); /* we take the code of the kword 'w' located at f->i_current position : i_current_end is the last pos of the seed */
        if (code >= 0) {

#ifdef KEYLISTCOMPRESS
          long int index      = f->keycode_first[seed][code];  /* search the first occurrence of 'w' on the sequence */
          long int index_end  = f->keycode_first[seed][code+1]; /* warning : undefined on "fast" first implemetation */
          long int i_previous = -1;
          if (index >= 0) {
            i_previous = keylist[seed][index];
          }

          /* for all the occurrences of 'w' before i_current_end */
          while (i_previous < keysize_query - i_current_end && i_previous >= 0 /* for "fast" only */ && index != index_end /* for "full" only */) {
#else
          long int index      = f->keycode_first[seed][code];  /* search the first occurrence of 'w' on the sequence */
          long int i_previous = index;

          /* consider all the occurrences of 'w' before f->i_current */
          while (i_previous < keysize_query - i_current_end  && index >= 0) {
#endif

            long int diagonal = i_current_end - i_previous + keysize_query;

#ifdef PREFETCH
            __builtin_prefetch(&GET_TAB(last_tuple_pos_with_diag, diagonal), 1, 3);
            __builtin_prefetch(&GET_TAB(last_tuple_ptr_with_diag, diagonal), 1, 3);
#endif

            /* consider all the interesting diagonals by adding indel bias  0 +1 -1 +2...
             * until we reach a new tuple or the upper indel bound */

            long int observed_diagonal = diagonal;
            long int delta_diagonal = 0;

#ifdef CACHE
            {
              long int i_hash     = HASH_DIAG(MAX((diagonal - gp_delta_stat), 0)MEMOPT_MOD);
              long int i_hash_end = HASH_DIAG(MIN((diagonal + gp_delta_stat), MEMOPT_TABSIZE)MEMOPT_MOD);
              do {
                  long int i_previous_hash = last_tuple_pos_with_diag_hash[i_hash];
                  if ((i_previous_hash >= i_current_end - gp_rho_stat) &&
                      (i_previous_hash <= i_current_end) &&
                      (i_previous_hash - observed_diagonal <= i_current_end - diagonal) &&
                      (i_previous_hash >= 0)
                      )
                    goto diag_start;
                  if (i_hash  == i_hash_end)
                    goto maj_last_tuple_ptr_with_diag;
                  i_hash++;
                  i_hash %= HASH_DIAG(MEMOPT_TABSIZE);
                }  while (1);

            }
          diag_start:
#endif

            while (delta_diagonal < 2 * gp_delta_stat + 1) {

              /* search for interesting tuples */
              long int j = GET_TAB_MIN(last_tuple_pos_with_diag, observed_diagonal, i_current_end - gp_rho_stat);

              /* if tuples found are respecting criteria (neighboors on diagonal and on sequences) */
              if (j >= 0                                                           &&
#ifdef OLDSELECTION
                  j                     <= i_current_end                           &&
#endif
                  j                     >= i_current_end            - gp_rho_stat  &&
#ifdef OLDSELECTION
                  j - observed_diagonal <= i_current_end - diagonal                &&
#endif
                  j - observed_diagonal >= i_current_end - diagonal - gp_rho_stat   ) {

                tuple * t, * t_;

                /*---------------------------------------------------------------------*/

                /* (1) possibly create one instance of tuple [j]
                 * if this one was not in any list.
                 */

                if (EXIST_TAB(last_tuple_ptr_with_diag, observed_diagonal) == 0) {

                  STATS_NB_CHAINS_BUILT_INC(f);

                  /* create a tuple [j] with a small seed and update lists */
                  CREATETUPLE(t, j, observed_diagonal, gp_seeds_span_min);
                  /* create a tuple list element */
                  CREATETUPLELIST(f->last_tl->next, t, NULL);
                  f->last_tl = f->last_tl->next;

                  PUT_TAB(last_tuple_ptr_with_diag, observed_diagonal, t);
                } else {

                  /* (2) we search for the last tuple [j] instanciation
                   */
                  t = (tuple *)GET_TAB(last_tuple_ptr_with_diag, observed_diagonal);
                }

                /*
                 * (3) iterate tuples and add this new one to the most adapted one
                 */

                do {
                  long int advance;
                  if (diagonal == t->diagonal  && ((advance = i_current_end - t->occurrence) <= gp_seeds_span[seed] + 1)) {

                    /* (3a) possible fusion with this tuple (no instanciation)
                     *
                     * [inter-tuple overlap]
                     */

                    if (advance >= 0) {
                      t->leftsize          += advance;
                      t->occurrence         = i_current_end;
                      goto no_maj_last_tuple_ptr_with_diag;
                    } else {
                      goto next_key;
                    }
                  }

                  /* list search for the most fitted */
                  if (!(t->next))
                    break;
                  t = t->next;

                } while (1);

                /*
                 * (3d) any fusion has not been possible
                 *
                 * we instanciate the tuple and add this one at the end of the chain.
                 *
                 */

                CREATETUPLE(t_, i_current_end, diagonal, gp_seeds_span[seed]);
                t->next = t_;
                PUT_TAB(last_tuple_ptr_with_diag, diagonal, t_);

                /*
                 * (4) dont reset last_tuple_ptr_with_diag[diagonal] to NULL
                 * because it has been modified in (3)
                 */
                goto no_maj_last_tuple_ptr_with_diag;

              }

              /*---------------------------------------------------------------------*/

              /* diagonal walk d, d+1, d-1, ... */
              delta_diagonal++;
              observed_diagonal = diagonal + gv_delta_shift[delta_diagonal];
              observed_diagonal = MAX(MIN(observed_diagonal, f->i_current+keysize_query), i_current_end);
              /* border effects avoided*/

            } /* while : diagonal walk */


            if (gp_hitcriterion == 1) {
              SINGLEHITDIAGONAL(datarev, datasize, data, datasize);
            }

            /* reset last tuple instancied to NULL (because not instance is created)
             */

#ifdef CACHE
          maj_last_tuple_ptr_with_diag:
#endif
            DEL_TAB(last_tuple_ptr_with_diag, diagonal, NULL);

          no_maj_last_tuple_ptr_with_diag:

            /*
             * update last tuple position  with diagonal = 'diagonal'
             */
            PUT_TAB(last_tuple_pos_with_diag, diagonal, i_current_end);

#ifdef CACHE
            last_tuple_pos_with_diag_hash[HASH_DIAG(diagonal MEMOPT_MOD)] = i_current_end;
#endif
            /* iterate 'w' code list before i_current_end position */
          next_key:
            STATS_NB_SEEDS_INC(f);

#ifdef KEYLISTCOMPRESS
            index ++;
            i_previous = keylist[seed][index];
#else
            i_previous = index = keylist[seed][index];
#endif
          } /* while i_previous */
        } /* if (code >= 0) */
      } /* if (f->i_current <= datasize - gp_seeds_span[seed]) */
    } /* for each seed */


    if (!(f->i_current % ALIGN_EVERY_NB_ITER)) {
#ifdef TRACE
      Display_Progress(f->i_current, (keysize_query), f);
#endif

      STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
      DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
      if (f->thread_align) {
        WAIT_THREAD(f->thread_align);
        f->thread_align = 0;
      }
      if (f->i_current) {
        /* dont start at very first time :
         * - avoid too small sequences being threaded
         * - give some elements to be processed to the work_align
         */
        CREATE_THREAD(f->thread_align, thread_work_align, f);
      }
#else
      AlignAndFree(
                   datarev, datasize,
                   data,    datasize,
                   0,
                   f);
#endif

      STATS_ADD_CLOCK(f, clock_align);

    }
  } /* for each pos */

#ifdef TRACE
  Display_Progress(f->i_current, (keysize_query), f);
#endif

  STATS_ADD_CLOCK(f, clock_chain);

  /* aligning remaining tuples */
#ifdef DEBUG_ASSEMBLE
  DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
  if (f->thread_align) {
    WAIT_THREAD(f->thread_align);
    f->thread_align = 0;
  }
#endif

  /*
   * CREATE_THREAD(f->thread_align, thread_work_align, f);
   */

  AlignAndFree(
               datarev, datasize,
               data,    datasize,
               1,
               f);

  STATS_ADD_CLOCK(f, clock_align);

#ifdef CACHE
  FREE(last_tuple_pos_with_diag_hash, (1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif

  FREE_TAB(last_tuple_pos_with_diag, sizeof(long int));
  FREE_TAB(last_tuple_ptr_with_diag, sizeof(tuple*));
  {
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++) {
      FREE(keylist[seed], keylist_size[seed]);
      /* FREE(first[seed],   first_size[seed]); */
    }
  }

  STATS_ADD_CLOCK(f, clock_pre);

  return 0;
}


long int MultiAssemble_Double(char * data_query, long int datasize_query,
                              char * data_text,  long int datasize_text,
                              long int nbchunks_text, char ** chunkname_text,
                              long int * chunksize_text,
                              long int * chunkstrt_text,
                              Feature  *f) {
  /* used for chunk selection */
  long int i_base                     = 0;
  long int max_chunk_datasize_text    = 0;
  char *  chunk_data_text        = data_text;

  /*  main variables
   *  for the linked list computation
   */

  long int * keylist_query[MAX_SEED];
  long int   keylist_query_size[MAX_SEED];
  long int   keysize_query = datasize_query;
#ifdef MEMOPT
  long int   memopt_modmask = 0;
#else
  long int   keysize_text = 0;
#endif

  long int * last_tuple_pos_with_diag = NULL;
  tuple **   last_tuple_ptr_with_diag = NULL;
  /* >> */
  long int * last_chunk_diag_used     = NULL;
  long int   nb_last_chunk_diag_used  = 0;
  /* << */

#ifdef CACHE
  long int *last_tuple_pos_with_diag_hash;
#endif

#ifdef STATS
  f->last_clock = clock();
#endif

  /* to keep lists of tuple found */
  if (!f->first_tl)
    f->last_tl = f->first_tl = CreateTupleList();

  {
    long int i_chunk;
    for (i_chunk = 0; i_chunk < nbchunks_text; i_chunk++)
      max_chunk_datasize_text = MAX(max_chunk_datasize_text, chunksize_text[i_chunk]);
  }

#ifndef MEMOPT
  keysize_text = max_chunk_datasize_text;
#endif

  /* [0] left correction factor set to keysize_query */
  f->left_correction =  keysize_query;

  MEMOPT_SETMOD;

  /* [1] key linked list creation */
  {
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++)
      CreateKeyList(data_query, datasize_query, seed, &(keylist_query[seed]), &(keylist_query_size[seed]), f);
  }

  /* [2] allocations */
  INIT_TAB(last_tuple_pos_with_diag, long int *, (MEMOPT_TABSIZE) * sizeof(long int), 0xff);
  INIT_TAB(last_tuple_ptr_with_diag, tuple **, (MEMOPT_TABSIZE) * sizeof(tuple*), 0x00);

#ifdef CACHE
  last_tuple_pos_with_diag_hash = (long int *) MALLOC((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
  ASSERT(last_tuple_pos_with_diag_hash, "Assemble");
  memset(last_tuple_pos_with_diag_hash, 0xff, ((1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif

  /* >> */
  last_chunk_diag_used  = (long int *) MALLOC(MAX_DIAG_POS_KEPT_DATACHUNK * sizeof(long int));
  ASSERT(last_chunk_diag_used, "Assemble");
  memset(last_chunk_diag_used, 0x00, (MAX_DIAG_POS_KEPT_DATACHUNK * sizeof(long int)));
  /* << */

  /* for each chunk element */
  for (f->i_chunk = 0; f->i_chunk < nbchunks_text; f->i_chunk++) {

    /*
     * [3] clean used tables  and
     *     and set parameters
     */

    long int chunk_datasize_text  =  chunksize_text[f->i_chunk];
    long int chunk_keysize_text   =  chunksize_text[f->i_chunk];
    chunk_data_text          =  data_text + chunkstrt_text[f->i_chunk];
    i_base                   =  chunkstrt_text[f->i_chunk];

    /* >> */
    if (nb_last_chunk_diag_used < MAX_DIAG_POS_KEPT_DATACHUNK) {
      int i;

      /* sort before for cache efficiency */
      if (nb_last_chunk_diag_used >= 256)
        qsort(last_chunk_diag_used, nb_last_chunk_diag_used,  sizeof(long int), long_int_cmp);

      for (i = 0; i < nb_last_chunk_diag_used; i++) {
        DEL_TAB(last_tuple_pos_with_diag, last_chunk_diag_used[i], -1);
        DEL_TAB(last_tuple_ptr_with_diag, last_chunk_diag_used[i], NULL);
      }
      nb_last_chunk_diag_used = 0;
    } else {
      nb_last_chunk_diag_used = 0;
      RESET_TAB(last_tuple_pos_with_diag, (MEMOPT_TABSIZE)*sizeof(long int), 0xff);
      RESET_TAB(last_tuple_ptr_with_diag, (MEMOPT_TABSIZE)*sizeof(tuple*), 0x00);
    }
    /* << */

    STATS_ADD_CLOCK(f, clock_pre);

    /* [4] algorithm to create linked tuples */
#ifdef DEBUG_ASSEMBLE
    fprintf(stdout, "> f->i_chunk : %ld, chunk_keysize_text = %ld \n", f->i_chunk, chunk_keysize_text);fflush(NULL);
#endif
    /* for each pos  */
    for (f->i_current = 0; f->i_current < chunk_keysize_text; f->i_current++) {
      /* for each seed */
      long int seed;
      for (seed = 0; seed < gp_nb_seeds; seed++) {

        if (f->i_current <= chunk_keysize_text - gp_seeds_span[seed]) {
          long int i_current_end;
          long int code = 0;
          KEY(code, i_current_end, chunk_data_text, f->i_current, seed); /* we take the code of the kword 'w' located at f->i_current position : i_current_end is the last pos of the seed */
          if (code >= 0) {
#ifdef KEYLISTCOMPRESS
            long int index      = f->keycode_first[seed][code];  /* search the first occurrence of 'w' on the sequence */
            long int index_end  = f->keycode_first[seed][code+1]; /* warning : undefined on "fast" first implemetation */
            long int i_previous = -1;
            if (index >= 0) {
              i_previous = keylist_query[seed][index];
            }

            /* for all occurrences of 'w' on the "Query" sequence */
            while (i_previous >= 0 /* for "fast" only */ && index != index_end /* for "full" only */) {
#else
            long int index      = f->keycode_first[seed][code];  /* search the first occurrence of 'w' on the sequence */
            long int i_previous = index;

            /* for each occurrence of 'w' on the "Query" sequence */
            while (index >= 0) {
#endif

              long int diagonal = i_current_end - i_previous + keysize_query;

#ifdef PREFETCH
            __builtin_prefetch(&GET_TAB(last_tuple_pos_with_diag, diagonal), 1, 3);
            __builtin_prefetch(&GET_TAB(last_tuple_ptr_with_diag, diagonal), 1, 3);
#endif

            /* consider all the interesting diagonals by adding indel bias  0 +1 -1 +2...
             * until we reach a new tuple or the upper indel bound */

              long int observed_diagonal = diagonal;
              long int delta_diagonal = 0;

#ifdef CACHE
              {
                long int i_hash     = HASH_DIAG(MAX((diagonal - gp_delta_stat), 0)MEMOPT_MOD);
                long int i_hash_end = HASH_DIAG(MIN((diagonal + gp_delta_stat), MEMOPT_TABSIZE)MEMOPT_MOD);
                do {
                    long int i_previous_hash = last_tuple_pos_with_diag_hash[i_hash];
                    if ((i_previous_hash >= i_current_end - gp_rho_stat) &&
                        (i_previous_hash <= i_current_end) &&
                        (i_previous_hash - observed_diagonal <= i_current_end - diagonal) &&
                        (i_previous_hash >= 0)
                        )
                      goto diag_start;
                    if (i_hash  == i_hash_end)
                      goto maj_last_tuple_ptr_with_diag;
                    i_hash++;
                    i_hash %= HASH_DIAG(MEMOPT_TABSIZE);
                  }  while (1);

              }
            diag_start:
#endif

              while (delta_diagonal < 2 * gp_delta_stat + 1) {

                /* search for interesting tuples */
                long int j = GET_TAB_MIN(last_tuple_pos_with_diag, observed_diagonal, i_current_end - gp_rho_stat);

                /* if tuples found are respecting criteria (neighboors on diagonal and on sequences) */
                if (j >= 0                                                           &&
#ifdef OLDSELECTION
                    j                     <= i_current_end                           &&
#endif
                    j                     >= i_current_end            - gp_rho_stat  &&
#ifdef OLDSELECTION
                    j - observed_diagonal <= i_current_end - diagonal                &&
#endif
                    j - observed_diagonal >= i_current_end - diagonal - gp_rho_stat   ) {

                  tuple * t, * t_;

                  /*---------------------------------------------------------------------*/

                  /* (1) possibly create one instance of tuple [j]
                   * if this one was not in any list.
                   */

                  if (EXIST_TAB(last_tuple_ptr_with_diag, observed_diagonal) == 0) {

                    STATS_NB_CHAINS_BUILT_INC(f);

                    /* create a tuple [j] with a small seed and update lists */
                    CREATETUPLE(t, j, observed_diagonal, gp_seeds_span_min);
                    /* create a tuple list element */
                    CREATETUPLELIST(f->last_tl->next, t, NULL);
                    f->last_tl = f->last_tl->next;

                    PUT_TAB(last_tuple_ptr_with_diag, observed_diagonal, t);

                  } else {

                    /* (2) we search for the last tuple [j] instanciation
                     */
                    t = (tuple *)GET_TAB(last_tuple_ptr_with_diag, observed_diagonal);
                  }

                  /*
                   * (3) iterate tuples and add this new one to the most adapted one
                   */

                  do {
                    long int advance;
                    if (diagonal == t->diagonal  && ((advance = i_current_end - t->occurrence) <= gp_seeds_span[seed] + 1)) {

                      /* (3a) possible fusion with this tuple (no instanciation)
                       *
                       * [inter-tuple overlap]
                       */

                      if (advance >= 0) {
                        t->leftsize          += advance;
                        t->occurrence         = i_current_end;
                        goto no_maj_last_tuple_ptr_with_diag;
                      } else {
                        goto next_key;
                      }
                    }

                    /* list search for the most fitted */
                    if (!(t->next))
                      break;
                    t = t->next;

                  } while (1);

                  /*
                   * (3d) any fusion has not been possible
                   *
                   * we instanciate the tuple and add this one at the end of the chain.
                   *
                   */

                  CREATETUPLE(t_, i_current_end, diagonal, gp_seeds_span[seed]);
                  t->next = t_;
                  PUT_TAB(last_tuple_ptr_with_diag, diagonal, t_);

                  /*
                   * (4) dont reset last_tuple_ptr_with_diag[diagonal] to NULL
                   * because it has been modified in (3)
                   */
                  goto no_maj_last_tuple_ptr_with_diag;

                }

                /*---------------------------------------------------------------------*/

                /* diagonal walk d, d+1, d-1, ... */
                delta_diagonal++;
                observed_diagonal = diagonal + gv_delta_shift[delta_diagonal];
                observed_diagonal = MAX(MIN(observed_diagonal, f->i_current+keysize_query), i_current_end);
                /* border effects avoided*/

              } /* while : diagonal walk */


              if (gp_hitcriterion == 1) {
                SINGLEHITDIAGONAL_MULTI(data_query, datasize_query, chunk_data_text, chunk_datasize_text);
              }

              /* reset last tuple instancied to NULL (because not instance is created)
               */

#ifdef CACHE
            maj_last_tuple_ptr_with_diag:
#endif
              DEL_TAB(last_tuple_ptr_with_diag, diagonal, NULL);

            no_maj_last_tuple_ptr_with_diag:

              /*
               * update last tuple position  with diagonal = 'diagonal'
               */
              PUT_TAB(last_tuple_pos_with_diag, diagonal, i_current_end);
              /* >> */
              if (nb_last_chunk_diag_used < MAX_DIAG_POS_KEPT_DATACHUNK)
                last_chunk_diag_used[nb_last_chunk_diag_used++] = diagonal MEMOPT_MOD;
              /* << */

#ifdef CACHE
              last_tuple_pos_with_diag_hash[HASH_DIAG(diagonal MEMOPT_MOD)] = i_current_end;
#endif
              /* iterate 'w' code list before f->i_current position */
            next_key:
              STATS_NB_SEEDS_INC(f);

#ifdef KEYLISTCOMPRESS
              index ++;
              i_previous = keylist_query[seed][index];
#else
              i_previous = index = keylist_query[seed][index];
#endif
            } /* while (index != indexfin) */
          } /* if (code >= 0) */
        } /* if (f->i_current <= chunk_keysize_text - gp_seeds_span[seed]) */
      } /* for each seed */


      if (!(f->i_current % ALIGN_EVERY_NB_ITER)) {
#ifdef TRACE
        Display_Progress(f->i_current+i_base, (datasize_text), f);
#endif

        STATS_ADD_CLOCK(f, clock_chain);

#ifdef DEBUG_ASSEMBLE
        DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
        if (f->thread_align) {
          WAIT_THREAD(f->thread_align);
          f->thread_align = 0;
        }
        if (f->i_current) {
          /* dont start at very first time :
           * - avoid too small sequences being threaded
           * - give some elements to be processed to the work_align
           */
          CREATE_THREAD(f->thread_align, thread_work_align, f);
        }
#else
        AlignAndFree(
                     data_query,        datasize_query,
                     chunk_data_text,   chunk_datasize_text,
                     0,
                     f);
#endif

        STATS_ADD_CLOCK(f, clock_align);

      }
    } /* for each pos */

#ifdef TRACE
    Display_Progress(f->i_current+i_base, (datasize_text), f);
#endif

    STATS_ADD_CLOCK(f, clock_chain);

    /* aligning remaining tuples */
#ifdef DEBUG_ASSEMBLE
    DisplayListTupleList(f->first_tl);
#endif


#ifdef THREAD_ASSEMBLE_ALIGN
    if (f->thread_align) {
      WAIT_THREAD(f->thread_align);
      f->thread_align = 0;
    }
#endif

    /*
     * CREATE_THREAD(f->thread_align, thread_work_align, f);
     */

    AlignAndFree(
                 data_query     ,  datasize_query,
                 chunk_data_text,  chunk_datasize_text,
                 1,
                 f);

    STATS_ADD_CLOCK(f, clock_align);

  } /* end of [4] (chunk loop) */

  /* >> */
  FREE(last_chunk_diag_used, MAX_DIAG_POS_KEPT_DATACHUNK * sizeof(long int));
  /* << */

#ifdef CACHE
  FREE(last_tuple_pos_with_diag_hash, (1 + HASH_DIAG(MEMOPT_TABSIZE)) * sizeof(long int));
#endif

  FREE_TAB(last_tuple_pos_with_diag, sizeof(long int));
  FREE_TAB(last_tuple_ptr_with_diag, sizeof(tuple*));
  {
    long int seed;
    for (seed = 0; seed < gp_nb_seeds; seed++) {
      FREE(keylist_query[seed], keylist_query_size[seed]);
      /* FREE(first_query[seed],   first_query_size[seed]); */
    }
  }

  STATS_ADD_CLOCK(f, clock_pre);
  return 0;
}



