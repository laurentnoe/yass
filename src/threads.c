/*
 *  YASS 1.15
 *  Copyright (C) 2004-2016
 *  the YASS team
 *  Laurent Noe, Gregory Kucherov, Mikhail Roytberg, 
 *  Steven Corroy, Antoine De Monte, Christophe Valmir.
 *
 *  laurent.noe|<A>|univ-lille1.fr
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

/*
 *  YASS 1.15
 *  Copyright (C) 2004-2016
 *  the YASS team
 *  Laurent Noe, Gregory Kucherov, Mikhail Roytberg, 
 *  Steven Corroy, Antoine De Monte, Christophe Valmir.
 *
 *  laurent.noe|<A>|univ-lille1.fr
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
#include "global_var.h"
#include "threads.h"
#include "prdyn.h"

Feature * gv_feature[MAX_QUERY_CHUNK_THREADS][2];

/*
 *  Affichage du grouping
 */

int             regrouping_begin = 0;
int             regrouping_end   = 0;
int             filtering_begin  = 0;
int             filtering_end    = 0;
int             sorting_begin    = 0;
int             sorting_end      = 0;



#ifdef THREAD_QUERY_CHUNK
#if defined(WIN32) || defined(WIN64)
HANDLE merge_ma_mutex;
HANDLE query_chunk_mutex;
HANDLE gv_threads[MAX_QUERY_CHUNK_THREADS];
#else
pthread_mutex_t merge_ma_mutex;
pthread_mutex_t query_chunk_mutex;
pthread_t gv_threads[MAX_QUERY_CHUNK_THREADS];
#endif
#else
int merge_ma_mutex = 0;
int query_chunk_mutex = 0;
#endif



#ifdef THREAD_FORWARD_REVERSE
#if defined(WIN32) || defined(WIN64)

/* WINDOWS multi proc */
HANDLE regrouping_mutex_begin;
HANDLE regrouping_mutex_end  ;
HANDLE filtering_mutex_begin ;
HANDLE filtering_mutex_end   ;
HANDLE sorting_mutex_begin   ;
HANDLE sorting_mutex_end     ;
#else

/* UNIX multi proc*/
pthread_mutex_t regrouping_mutex_begin;
pthread_mutex_t regrouping_mutex_end  ;
pthread_mutex_t filtering_mutex_begin ;
pthread_mutex_t filtering_mutex_end   ;
pthread_mutex_t sorting_mutex_begin   ;
pthread_mutex_t sorting_mutex_end     ;

#endif

#else

/* Mono Proc*/
int regrouping_mutex_begin      = 0;
int regrouping_mutex_end        = 0;
int filtering_mutex_begin       = 0;
int filtering_mutex_end         = 0;
int sorting_mutex_begin         = 0;
int sorting_mutex_end           = 0;
#endif

void AllocInitFeature(Feature ** p_feature) {

  Feature * feature = (Feature *) MALLOC(sizeof(Feature));
  ASSERT(feature, InitFeature);

  feature->first_tl    = NULL;
  feature->last_tl     = NULL;

  feature->first_MA    = NULL;
  feature->last_MA     = NULL;

  feature->i_current   = 0;
  feature->i_chunk     = 0;

  feature->left_correction = 0;

  feature->MAminscore      = 0;
  {
    int i;
    for (i = 0; i < MAX_SEED; i++) {
      feature->keycode_first[i]  = NULL;
      feature->keycode_bucket[i] = NULL;
      feature->keycode_count[i]  = NULL;

      feature->keycode_first_list_pos[i] = MALLOC(MAX_FIRST_POS_KEPT_QUERYCHUNK * sizeof(long int));
      ASSERT(feature->keycode_first_list_pos[i],AllocInitFeature);
      feature->keycode_first_list_pos_size[i] = 0;
    }
  }

#ifdef STATS
  feature->clock_pre        = 0;
  feature->clock_align      = 0;
  feature->clock_chain      = 0;
  feature->clock_post       = 0;
  feature->nb_keys_removed  = 0;
  feature->nb_seeds         = 0;
  feature->nb_single_tests  = 0;
  feature->nb_single_hits   = 0;
  feature->nb_chains_built  = 0;
  feature->nb_ma            = 0;
  feature->nb_postprocessed_grouping_tests = 0;
  feature->nb_postprocessed_ma             = 0;
  feature->nb_postprocessed_grouping_links = 0;
#endif

#ifdef THREAD_FORWARD_REVERSE
  feature->thread_assemble = 0;
#endif

#ifdef THREAD_ASSEMBLE_ALIGN
  feature->thread_align    = 0;
#endif

  feature->last_point      = 1;

  /* mem allocations (separate paging of features ) */
  initialise_alignment(4 * MEGA, feature);
  CreateCountMA(feature);
  feature->nb_pair_of_triplets = lint_directtable(64,64,0);
  (*p_feature) = feature;
}



long int CreateCountMA(Feature * f)
{
  f->MAcount = (long int *) MALLOC(sizeof(long int)*8*sizeof(long int));
  ASSERT(f->MAcount, CreateCountMA);
  memset(f->MAcount, '\0', sizeof(long int)*8*sizeof(long int));
  return 0;
}


long int MinScoreOnCountMA(Feature * f, long int score)
{
  long int i  = 0;
  long int pi = 1;
  while (score > 0) {
    score >>= 1;
    i++;
    pi    <<= 1;
  }

  f->MAcount[i]++;

  if (f->MAcount[i] > gp_nbmaxlines && f->MAminscore < pi - 1) {
    fprintf(stderr,"\n");
    _WARNING("-O limit reached => increasing score threshold");
    f->MAminscore = pi - 1;
    fprintf(stderr,"  the score threshold is now = %ld\n", pi - 1);
  }

  return 0;
}







