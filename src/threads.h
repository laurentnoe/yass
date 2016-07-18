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

#ifndef __THREADS_H_
#define __THREADS_H_
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "tuple.h"
#include "util.h"
#include "list.h"

/*
 *    Include windows or linux for a thread
 *
 */

#ifdef THREAD_ASSEMBLE_ALIGN
#define THREAD
#endif

#ifdef THREAD_FORWARD_REVERSE
#define THREAD
#endif

#ifdef THREAD_QUERY_CHUNK
#define THREAD
#endif



#ifdef THREAD
#if defined(WIN32) || defined(WIN64)
#include <windows.h>
#else
#include <pthread.h>
#endif
#endif



#ifdef CHOOSERBTREE
#include "red_black.h"
#else
#ifdef CHOOSEAVLTREE
#include "avl.h"
#else
#error No tree definition : uncomment either "CHOOSERBTREE" or "CHOOSEAVLTREE" in the "util.h" file
#endif
#endif



/*
 * forward-backward thread features
 */
typedef struct _FEATURE_ {
  char * chunk_query;
  long int chunk_query_size;
  long int    reverse;     /* on the complementary strain */
  long int    left_correction;
  long int  * buffer01;    /* buffer for prdyn */
  long int    last_point;  /* count the last dot displayed */
  long int  * MAcount;     /* count number of MA with score between 2^i and < 2^(i+1) */
  long int    MAminscore;  /* minimal score to be reached */

  /* index part */
  long int * keycode_first[MAX_SEED]; /* keep a pointed to the "first" for each seed (otherwise set to NULL) */
  long int * keycode_bucket[MAX_SEED]; /* keep a pointed to the "bucket" for each seed (otherwise set to NULL) */
  long int * keycode_count[MAX_SEED]; /* keep a pointed to the "bucket" for each seed (otherwise set to NULL) */

  /* index part : used for rewriting the "first/bucket/count" tables by indexing when sparse ... */
  long int * keycode_first_list_pos[MAX_SEED];
  long int   keycode_first_list_pos_size[MAX_SEED];

  tuplelist *first_tl;
  tuplelist *last_tl;

  MA *first_MA;
  MA *last_MA;

  volatile long int i_current; /* position on the current chunk */
  long int i_chunk;   /* current chunk number */
  long int j_chunk;   /* current query chunk number */

#ifdef THREAD_FORWARD_REVERSE
#if defined(WIN32) || defined(WIN64)
  HANDLE    thread_assemble;
#else
  pthread_t thread_assemble;
#endif
#endif

#ifdef THREAD_ASSEMBLE_ALIGN
#if defined(WIN32) || defined(WIN64)
 HANDLE    thread_align;
#else
 pthread_t thread_align;
#endif
#endif

  /* B) Statistics */
#ifdef STATS
  clock_t last_clock;  /* lask clock memorized by the thread */

  /* a) time */
  clock_t clock_pre;
  clock_t clock_chain;
  clock_t clock_align;
  clock_t clock_post;

  /* b) numbers */

  /* b.1) pre-process */
  long int    nb_keys_removed;
  /* b.2) assemble */
  long int    nb_seeds;
  long int    nb_single_tests;
  long int    nb_single_hits;
  long int    nb_chains_built;
  /* b.3) align */
  long int    nb_ma;
  /* b.4) post-process */
  long int    nb_postprocessed_grouping_tests;
  long int    nb_postprocessed_ma;
  long int    nb_postprocessed_grouping_links;
#endif

  /*
   * grouping variables
   */

  /* The windows size used to regroup aligments */
  int windows;
  /* The current position (start of the window)*/
  int win_position;

  /* The global queue of MA, ordered by right pos end*/
  struct _queue_MA Q;

  /* We use a tree which can be a red black tree or an avl tree*/
#ifdef CHOOSERBTREE
  struct _rb_tree T;
#else
#ifdef CHOOSEAVLTREE
  struct _avl_tree T;
#else
#error No tree definition : uncomment either "CHOOSERBTREE" or "CHOOSEAVLTREE" inside the "util.h" file
#endif
#endif

  /* used by complexity filters : triplet counters */
  long int nb_triplet_count1[64];
  long int nb_triplet_count2[64];
  long int nb_non_mutated_triplet_count[64];
  long int ** nb_pair_of_triplets;

} Feature;


extern Feature * gv_feature[MAX_QUERY_CHUNK_THREADS][2];

/*********************
 *
 * Thread definition
 *
 ***********************/

/*
 *  Create Thread
 */

#ifdef THREAD
#if defined(WIN32) || defined(WIN64)
#define CREATE_THREAD(thread,function,parameters)                                                     \
if ((thread=CreateThread(NULL,0,function, (PVOID) parameters,0, (DWORD *)&gv_thread_result))==NULL) { \
      fprintf (stderr,"* Error : \"cant create thread\"\n");                                          \
      fprintf (stderr,"please use mono-threaded yass\n");                                             \
      exit(1);                                                                                        \
   }

#else
#define CREATE_THREAD(thread,function,parameters)                   \
  if (pthread_create(&thread,NULL,function,(void *)parameters)) {   \
      fprintf (stderr,"* Error : \"cant create thread\"\n");        \
      fprintf (stderr,"please use mono-threaded yass\n");           \
      exit(1);                                                      \
   }
#endif
#else
#define CREATE_THREAD(thread,function,parameters)    { function(parameters); }
#endif

/*
 *  End Thread
 */

#ifdef THREAD
#if defined(WIN32) || defined(WIN64)
#define END_THREAD() {return 0;}
#else
#define END_THREAD() {pthread_exit(NULL);}
#endif
#else
#define END_THREAD() {return NULL;}
#endif

/*
 *  Wait Thread
 */

#ifdef THREAD
#if defined(WIN32) || defined(WIN64)
#define WAIT_THREAD(THREAD)  { WaitForSingleObject(THREAD, INFINITE); }
#else
#define WAIT_THREAD(THREAD)  { pthread_join (THREAD, (void *) &gv_thread_result); }
#endif
#else
#define WAIT_THREAD(THREAD)  ;
#endif



/*
 * Threads for query chunk
 */

#ifdef THREAD_QUERY_CHUNK
#if defined(WIN32) || defined(WIN64)
/* Windows multi cpu */
extern HANDLE merge_ma_mutex        ;
extern HANDLE query_chunk_mutex     ;
#define INIT_QUERY_MUTEX() {                                       \
 query_chunk_mutex = CreateMutex(NULL,FALSE,"query_chunk_mutex");  \
 merge_ma_mutex = CreateMutex(NULL,FALSE,"merge_ma_mutex");        \
}
#define LOCK(NAMEMUTEX)   {WaitForSingleObject(NAMEMUTEX,INFINITE);}
#define UNLOCK(NAMEMUTEX) {ReleaseMutex(NAMEMUTEX);}

#else
/* Unix multi cpu */
extern pthread_mutex_t merge_ma_mutex        ;
extern pthread_mutex_t query_chunk_mutex     ;
#define INIT_QUERY_MUTEX() {                    \
 pthread_mutex_init(&query_chunk_mutex,NULL);   \
 pthread_mutex_init(&merge_ma_mutex,NULL);      \
}

#define LOCK(NAMEMUTEX)   {pthread_mutex_lock(&NAMEMUTEX);}
#define UNLOCK(NAMEMUTEX) {pthread_mutex_unlock(&NAMEMUTEX);}

#endif

#else
/* Mono cpu */
extern int merge_ma_mutex     ;
extern int query_chunk_mutex  ;
#define INIT_QUERY_MUTEX()    ;
#define LOCK(NAMEMUTEX)       ;
#define UNLOCK(NAMEMUTEX)     ;

#endif


/*
 *  Print
 */

extern int             regrouping_begin ;
extern int             regrouping_end   ;
extern int             filtering_begin  ;
extern int             filtering_end    ;
extern int             sorting_begin    ;
extern int             sorting_end      ;

#ifdef THREAD_FORWARD_REVERSE
#if defined(WIN32) || defined(WIN64)

/* WINDOWS multi cpu */
extern HANDLE regrouping_mutex_begin;
extern HANDLE regrouping_mutex_end  ;
extern HANDLE filtering_mutex_begin ;
extern HANDLE filtering_mutex_end   ;
extern HANDLE sorting_mutex_begin   ;
extern HANDLE sorting_mutex_end     ;

#define INIT_MUTEX() {                                                     \
 regrouping_mutex_begin = CreateMutex(NULL,FALSE,"regroup_mutex_begin");   \
 regrouping_mutex_end   = CreateMutex(NULL,FALSE,"regroup_mutex_end");     \
 filtering_mutex_begin  = CreateMutex(NULL,FALSE,"filtering_mutex_begin"); \
 filtering_mutex_end    = CreateMutex(NULL,FALSE,"filtering_mutex_end");   \
 sorting_mutex_begin    = CreateMutex(NULL,FALSE,"sorting_mutex_begin");   \
 sorting_mutex_end      = CreateMutex(NULL,FALSE,"sorting_mutex_end");     \
}


#define DISPLAY_BEGIN(NAMEMUTEX,NOMVAR,STR) {                             \
  if (gp_selection_fasta) {                                               \
     WaitForSingleObject(NAMEMUTEX,INFINITE);                             \
     NOMVAR++;                                                            \
     if (NOMVAR >= gp_reverse) {fprintf(stderr,"%s",STR);fflush(NULL);}   \
     ReleaseMutex(NAMEMUTEX);                                             \
  }                                                                       \
}

#define DISPLAY_END(NAMEMUTEX,NOMVAR)  {                                  \
  if (gp_selection_fasta) {                                               \
     WaitForSingleObject(NAMEMUTEX,INFINITE);                             \
     NOMVAR++;                                                            \
     if (NOMVAR>=gp_reverse) {fprintf(stderr,"finished\n");fflush(NULL);}  \
     ReleaseMutex(NAMEMUTEX);                                             \
  }                                                                       \
}

#else

/* UNIX multi cpu */
extern pthread_mutex_t regrouping_mutex_begin;
extern pthread_mutex_t regrouping_mutex_end  ;
extern pthread_mutex_t filtering_mutex_begin ;
extern pthread_mutex_t filtering_mutex_end   ;
extern pthread_mutex_t sorting_mutex_begin   ;
extern pthread_mutex_t sorting_mutex_end     ;


#define INIT_MUTEX()  \
 pthread_mutex_init(&regrouping_mutex_begin,NULL); \
 pthread_mutex_init(&regrouping_mutex_end,NULL);   \
 pthread_mutex_init(&filtering_mutex_begin,NULL);  \
 pthread_mutex_init(&filtering_mutex_end,NULL);    \
 pthread_mutex_init(&sorting_mutex_begin,NULL);    \
 pthread_mutex_init(&sorting_mutex_end,NULL);      \


#define DISPLAY_BEGIN(NAMEMUTEX,NOMVAR,STR) {                             \
  if (gp_selection_fasta) {                                               \
    pthread_mutex_lock(&NAMEMUTEX);                                       \
    NOMVAR++;                                                             \
    if (NOMVAR >= gp_reverse) { fprintf(stderr,"%s",STR);fflush(NULL); }  \
    pthread_mutex_unlock(&NAMEMUTEX);                                     \
  }                                                                       \
}

#define DISPLAY_END(NAMEMUTEX,NOMVAR) {                                   \
  if (gp_selection_fasta) {                                               \
    pthread_mutex_lock(&NAMEMUTEX);                                       \
    NOMVAR++;                                                             \
    if (NOMVAR>=gp_reverse) {fprintf(stderr,"finished\n");fflush(NULL);}  \
    pthread_mutex_unlock(&NAMEMUTEX);                                     \
  }                                                                       \
}

#endif
#else

/* Mono cpu */
extern int regrouping_mutex_begin;
extern int regrouping_mutex_end  ;
extern int filtering_mutex_begin;
extern int filtering_mutex_end  ;
extern int sorting_mutex_begin;
extern int sorting_mutex_end  ;


#define INIT_MUTEX()  ;

#define DISPLAY_BEGIN(NAMEMUT,NOMVAR,STRING) {                            \
  if (gp_selection_fasta) {                                               \
    NOMVAR = NOMVAR + 1;                                                  \
    if ((NOMVAR)==2) {fprintf(stderr,"%s",(STRING));fflush(NULL);}        \
  }                                                                       \
}

#define DISPLAY_END(NAMEMUT,NOMVAR) {                                     \
    if (gp_selection_fasta) {                                             \
      NOMVAR++;                                                           \
      if (NOMVAR==2) {fprintf(stderr,"finished\n");fflush(NULL);}         \
    }                                                                     \
  }

#endif

void AllocInitFeature(Feature ** p_feature);

long int CreateCountMA(Feature * f);

long int MinScoreOnCountMA(Feature * f, long int score);




#define RESETFEATURE(f,query_chunk_minscore,query_chunk_data,query_chunk_size,query_chunk_reverse,query_chunk_nb) { \
      (f)->first_MA         = NULL;                                                                                 \
      (f)->last_MA          = NULL;                                                                                 \
      (f)->i_current        = 0;                                                                                    \
      (f)->i_chunk          = 0;                                                                                    \
      (f)->left_correction  = 0;                                                                                    \
      (f)->last_point       = 0;                                                                                    \
      (f)->MAminscore       = query_chunk_minscore;                                                                 \
      (f)->chunk_query      = query_chunk_data;                                                                     \
      (f)->chunk_query_size = query_chunk_size;                                                                     \
      (f)->reverse          = query_chunk_reverse;                                                                  \
      (f)->j_chunk          = query_chunk_nb;                                                                       \
}

#ifdef STATS

#define STATS_ADD_CLOCK(f,variable)             { \
   clock_t _current_clock = clock();              \
   f->variable  += _current_clock - f->last_clock;\
   f->last_clock = _current_clock;                \
}

#define STATS_NB_SINGLE_TESTS(f)                     {f->nb_single_tests++;}
#define STATS_NB_SINGLE_HITS(f)                      {f->nb_single_hits++;}
#define STATS_NB_KEYS_REMOVED_INC(f)                 {f->nb_keys_removed++;}
#define STATS_NB_SEEDS_INC(f)                        {f->nb_seeds++;}
#define STATS_NB_CHAINS_BUILT_INC(f)                 {f->nb_chains_built++;}
#define STATS_NB_MA_INC(f)                           {f->nb_ma++;}
#define STATS_NB_POSTPROCESSED_GROUPING_TESTS_INC(f) {f->nb_postprocessed_grouping_tests++;}
#define STATS_NB_POSTPROCESSED_MA_INC(f)             {f->nb_postprocessed_ma++;}
#define STATS_NB_POSTPROCESSED_GROUPING_LINKS_INC(f) {f->nb_postprocessed_grouping_links++;}

#else

#define STATS_ADD_CLOCK(variable)
#define STATS_NB_SINGLE_TESTS(f)
#define STATS_NB_SINGLE_HITS(f)
#define STATS_NB_KEYS_REMOVED_INC(f)
#define STATS_NB_SEEDS_INC(f)
#define STATS_NB_CHAINS_BUILT_INC(f)
#define STATS_NB_MA_INC(f)
#define STATS_NB_POSTPROCESSED_GROUPING_TESTS_INC(f)
#define STATS_NB_POSTPROCESSED_MAX_INC(f)
#define STATS_NB_POSTPROCESSED_GROUPING_LINKS_INC(f)

#endif


#ifdef THREAD_QUERY_CHUNK
#if defined(WIN32) || defined(WIN64)
extern HANDLE gv_threads[MAX_QUERY_CHUNK_THREADS];
#else
extern pthread_t gv_threads[MAX_QUERY_CHUNK_THREADS];
#endif
#endif

#endif
