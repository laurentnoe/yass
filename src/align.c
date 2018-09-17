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

/* 1) include utils macro */
#include "util.h"
/* 2) include global variables */
#include "global_var.h"
/* 3) the current file defs */
#include "align.h"
/* 4) other files */
#include "kword.h"
#include "tuple.h"
#include "prdyn.h"
#include "display.h"
#include "util.h"


#define RESOLVE_OVERLAPS(pta,ta,tb)                                               \
if (TEL_POS(tb) >= TBL_POS(ta) && TER_POS(tb) >= TBR_POS(ta)) {                   \
  if (TBL_POS(tb) < TEL_POS(ta) || TBR_POS(tb) < TER_POS(ta)) {                   \
    long int tdelta = MAX(TEL_POS(ta) - TBL_POS(tb), TER_POS(ta) - TBR_POS(tb));  \
    (tb)->leftsize -= MAX(tdelta,gp_seeds_span_min);                              \
    if ((tb)->leftsize < 0)                                                       \
      (tb)->leftsize = 0;                                                         \
    (ta)->leftsize   -= MAX(tdelta,gp_seeds_span_min);                            \
    (ta)->occurrence -= MAX(tdelta,gp_seeds_span_min);                            \
    if ((ta)->leftsize < 0) {                                                     \
      (ta)->occurrence -= (ta)->leftsize;                                         \
      (ta)->leftsize = 0;                                                         \
    }                                                                             \
  }                                                                               \
} else {                                                                          \
  _WARNING("negative gap size in RESOLVE_OVERLAPS");                              \
  DisplayTuple(ta);DisplayTuple(tb);fflush(NULL);                                 \
  if (TSIZE(ta) <= TSIZE(tb) && pta) {                                            \
    pta->next = ta->next;                                                         \
    FREE(ta,sizeof(tuple));                                                       \
    ta = pta->next;                                                               \
  } else {                                                                        \
    tuple * __t__ = ta->next;                                                     \
    ta->next  = ta->next->next;                                                   \
    FREE(__t__,sizeof(tuple));                                                    \
  }                                                                               \
  goto full_restart;                                                              \
}

#define MEMORISE_IN_DISPLAY_OR_ERASE(prev_first_ta,first_ta,last_ta,t_score)                      { \
      long int left_score      = 0,                                                                 \
               left_pos_begin  = TBL_POS(first_ta),                                                 \
               right_pos_begin = TBR_POS(first_ta),                                                 \
               right_score     = 0,                                                                 \
               left_pos_end    = TEL_POS(last_ta),                                                  \
               right_pos_end   = TER_POS(last_ta);                                                  \
      long int tbl_space   = MIN(left_pos_begin,gp_border);                                         \
      long int tbr_space   = MIN(right_pos_begin,gp_border);                                        \
      long int tel_space   = MIN(MAX(datasize_query -  left_pos_end,0),gp_border);                  \
      long int ter_space   = MIN(MAX(datasize_text  - right_pos_end,0),gp_border);                  \
      long int tp_score    = t_score;                                                               \
      long int left_flag   = (tbl_space == gp_border && tbr_space == gp_border);                    \
      long int right_flag  = (tel_space == gp_border && ter_space == gp_border);                    \
                                                                                                    \
      if (left_flag && right_flag) {                                                                \
        left_alignment_SG_Lz(                                                                       \
                             data_query, TBL_POS(first_ta),                                         \
                             tbl_space,                                                             \
                             data_text,  TBR_POS(first_ta),                                         \
                             tbr_space,                                                             \
                             gp_xdrop, &left_score, &left_pos_begin, &right_pos_begin, feature);    \
        right_alignment_SG_Lz(                                                                      \
                              data_query, TEL_POS(last_ta),                                         \
                              tel_space,                                                            \
                              data_text,  TER_POS(last_ta),                                         \
                              ter_space,                                                            \
                              gp_xdrop, &right_score, &left_pos_end, &right_pos_end, feature);      \
      } else {                                                                                      \
        left_alignment_SG_Border(                                                                   \
                          data_query, TBL_POS(first_ta),                                            \
                          tbl_space,                                                                \
                          data_text,  TBR_POS(first_ta),                                            \
                          tbr_space,                                                                \
                          gp_xdrop, &left_score, &left_pos_begin, &right_pos_begin, feature);       \
        right_alignment_SG_Border(                                                                  \
                           data_query, TEL_POS(last_ta),                                            \
                           tel_space,                                                               \
                           data_text,  TER_POS(last_ta),                                            \
                           ter_space,                                                               \
                           gp_xdrop, &right_score, &left_pos_end, &right_pos_end, feature);         \
      }                                                                                             \
      tp_score += left_score + right_score;                                                         \
      MinScoreOnCountMA(feature,tp_score);                                                          \
      if (tp_score >  feature->TUminscore) {                                                        \
         double entropy_query = Entropy3mer(data_query, left_pos_begin,                             \
                                            left_pos_end-left_pos_begin);                           \
         if (entropy_query > gp_entropy_min) {                                                      \
            double entropy_text = Entropy3mer(data_text, right_pos_begin,                           \
                                              right_pos_end-right_pos_begin);                       \
            if (entropy_text > gp_entropy_min) {                                                    \
                                                                                                    \
                  if (left_flag && right_flag) {                                                    \
                    left_score  = 0,                                                                \
                    left_pos_begin  = TBL_POS(first_ta),                                            \
                    right_pos_begin = TBR_POS(first_ta),                                            \
                    left_alignment_SG_Border(                                                       \
                          data_query, TBL_POS(first_ta),                                            \
                          tbl_space,                                                                \
                          data_text,  TBR_POS(first_ta),                                            \
                          tbr_space,                                                                \
                          gp_xdrop, &left_score, &left_pos_begin, &right_pos_begin, feature);       \
                    t_score += left_score;                                                          \
                    right_score     = 0,                                                            \
                    left_pos_end    = TEL_POS(last_ta),                                             \
                    right_pos_end   = TER_POS(last_ta);                                             \
                    right_alignment_SG_Border(                                                      \
                           data_query, TEL_POS(last_ta),                                            \
                           tel_space,                                                               \
                           data_text,  TER_POS(last_ta),                                            \
                           ter_space,                                                               \
                           gp_xdrop, &right_score, &left_pos_end, &right_pos_end, feature);         \
                    t_score += right_score;                                                         \
                  }                                                                                 \
                  {                                                                                 \
              /*MA * ma =*/   CreateMA(left_pos_begin,                                              \
                                       right_pos_begin,                                             \
                                       left_pos_end,                                                \
                                       right_pos_end,                                               \
                                       left_score,                                                  \
                                       right_score,                                                 \
                                       tuple_list,                                                  \
                                       prev_first_ta,                                               \
                                       first_ta,                                                    \
                                       last_ta,                                                     \
                                       t_score,                                                     \
                                       MIN(entropy_query,entropy_text),                             \
                                       feature);                                                    \
                    /*DisplayMA(ma);*/                                                              \
                    /*display_alignment_SG_on_MA(data_query,data_text,ma);*/                        \
                  }                                                                                 \
               } else {                                                                             \
                  CleanTupleList(tuple_list, prev_first_ta, first_ta, last_ta);                     \
               }                                                                                    \
           } else {                                                                                 \
              CleanTupleList(tuple_list, prev_first_ta, first_ta, last_ta);                         \
           }                                                                                        \
         } else {                                                                                   \
            CleanTupleList(tuple_list, prev_first_ta, first_ta, last_ta);                           \
         }                                                                                          \
  }                                                                                                 \




/***********************************************************/


#ifdef INLINE
inline
#endif
long int AlignTuples( tuplelist * tuple_list,
                      char * data_query, long int datasize_query,
                      char * data_text,  long int datasize_text,
                      Feature *feature)
{
  long int left_correction = feature->left_correction;

  /* "first" and "last" tuples in the current alignment, "prev_last_tuple_aligned" helps to edit, and "t" helps to navigate in tuples */
  tuple * first_tuple_aligned, * last_tuple_aligned, * prev_last_tuple_aligned, * t;
  long int prev_total_score, last_seedscore, total_score;

  /* set a restarting point when many "tuple overlaps" remain unsolved */
 full_restart:

  /* (1/2) First tuple alignment */
  first_tuple_aligned      =  tuple_list->first_tuple;
  last_tuple_aligned       =  tuple_list->first_tuple;
  prev_last_tuple_aligned  =  NULL;
  t                        =  tuple_list->first_tuple->next;
  prev_total_score         =  0;

#ifdef PREFETCH
  __builtin_prefetch(&data_query[TBL_POS(tuple_list->first_tuple)],0);
  __builtin_prefetch(&data_text[TBR_POS(tuple_list->first_tuple)],0);
  __builtin_prefetch(t,0);
#endif

  last_seedscore =  alignment_SG_Strait ( data_query, TBL_POS(tuple_list->first_tuple),
                                          data_text , TBR_POS(tuple_list->first_tuple),
                                          TSIZE(tuple_list->first_tuple),
                                          feature);
  total_score    =  last_seedscore;

#ifdef DEBUG_ALIGNTUPLES
  printf("# AlignTuples :\n");
  DisplayTupleList(tuple_list);
  fflush(NULL);
  printf("# -- first tuple : seedscore = %ld (total_score = %ld) :\n",last_seedscore,total_score);
  DisplayTuple(last_tuple_aligned);
  fflush(NULL);
#endif

  /* find a good start : remove bad scores from the begining first, computing only the seed score : non resolved (zero lenght) case are removed here if its a full_restart */
  while (last_seedscore <= 0) {

#ifdef DEBUG_ALIGNTUPLES
    printf("# -- [[<<delete first tuple>>]] :\n");
    DisplayTuple(tuple_list->first_tuple);
    fflush(NULL);
#endif
    t                       = tuple_list->first_tuple;
    tuple_list->first_tuple = tuple_list->first_tuple->next;
    FREE(t,sizeof(tuple));

    if (tuple_list->first_tuple) {
      /* retry on the next one */
      first_tuple_aligned = last_tuple_aligned = tuple_list->first_tuple;
      t                   = tuple_list->first_tuple->next;
      last_seedscore      = alignment_SG_Strait ( data_query, TBL_POS(tuple_list->first_tuple),
                                                  data_text , TBR_POS(tuple_list->first_tuple),
                                                  TSIZE(tuple_list->first_tuple),
                                                  feature);
      total_score         = last_seedscore;
    } else {
      /* all seedscore are with bad score ... and have been removed */
      return 0;
    }
#ifdef DEBUG_ALIGNTUPLES
    printf("# -- NEW first tuple : seedscore = %ld (total_score = %ld) :\n",last_seedscore,total_score);
    DisplayTuple(last_tuple_aligned);
    fflush(NULL);
#endif
  }

  /* (2/2) Chaining tuple alignment */
  while (t) {

#ifdef PREFETCH
    __builtin_prefetch(&data_query[TBL_POS(t)],0);
    __builtin_prefetch(&data_text[TBR_POS(t)],0);
    __builtin_prefetch(t->next,0);
#endif


    long int left_gap  = TGAP_L( last_tuple_aligned, t);
    long int right_gap = TGAP_R( last_tuple_aligned, t);
    long int seedscore = alignment_SG_Strait ( data_query, TBL_POS(t),
                                               data_text , TBR_POS(t),
                                               TSIZE(t),
                                               feature);
    long int score = 0;



    /* I) compute scores */
    if (left_gap >= 0 && right_gap >= 0) {
      /* non-overlapped case */
      score = alignment_SG_Border ( data_query, TEL_POS(last_tuple_aligned),  left_gap,
                                    data_text , TER_POS(last_tuple_aligned), right_gap,
                                    gp_xdrop+seedscore,
                                    feature);
    } else {
      /* overlapped case */
      long int max_neg_gap = MIN(left_gap,right_gap);
      score =
        gp_cost_max_substitution_matrix * max_neg_gap
        +
        gp_cost_gap_opened
        +
        gp_cost_gap_continued * ABS(left_gap - right_gap);
    }

#ifdef DEBUG_ALIGNTUPLES
    printf("# -- current tuple [t] : seedscore = %ld, score = %ld (total_score = %ld) :\n",seedscore,score,total_score);
    DisplayTuple(t);
    fflush(NULL);
#endif


    /* II) try "increasing scores" by adding "t" to the end of "last_tuple_aligned" */
    {
      long int dscore = score + seedscore;
      if (total_score + dscore > MAX(total_score,seedscore)) {

        /* a) check with next tuples if some "tnext" can added to the end of "last_tuple_aligned"  */
        tuple *tnext = t->next;
        tuple *tnext_prev = t;
        long int l_next_left_gap  = 0;
        long int l_next_right_gap = 0;

#ifdef DEBUG_ALIGNTUPLES
        printf("   - (II) increasing scores (total_score + dscore > MAX(total_score,seedscore))\n");
        printf("     -> tnext loop\n");
        fflush(NULL);
#endif
        while (tnext &&
               (l_next_left_gap  = TGAP_L( last_tuple_aligned, tnext)) <= gp_rho_stat &&
               (l_next_right_gap = TGAP_R( last_tuple_aligned, tnext)) <= gp_rho_stat) {


          /* 1) compute between "last_tuple_aligned" and "tnext" */
          long int next_seedscore   = alignment_SG_Strait (data_query, TBL_POS(tnext),
                                                           data_text , TBR_POS(tnext),
                                                           TSIZE(tnext),
                                                           feature);
          if (/*X*/next_seedscore > seedscore) {

            long int l_next_score;
            if (l_next_left_gap >= 0 && l_next_right_gap >= 0) {
              /* non-overlapped case */
              l_next_score = alignment_SG_Border ( data_query, TEL_POS(last_tuple_aligned), l_next_left_gap,
                                                   data_text , TER_POS(last_tuple_aligned), l_next_right_gap,
                                                   gp_xdrop + next_seedscore,
                                                   feature);
            } else {
              /* overlapped case */
              long int l_next_max_neg_gap = MIN(l_next_left_gap,l_next_right_gap);
              l_next_score =
                gp_cost_max_substitution_matrix * l_next_max_neg_gap
                +
                gp_cost_gap_opened
                +
                gp_cost_gap_continued * ABS(l_next_left_gap - l_next_right_gap);
            }

            {
              long int l_next_dscore =  next_seedscore + l_next_score;

              if (l_next_dscore > 0) {

                /* 2) compute between "t" and "tnext" */
                {
                  long int t_next_left_gap  = TGAP_L(t, tnext);
                  long int t_next_right_gap = TGAP_R(t, tnext);
                  long int t_next_score     = 0;

                  if (t_next_left_gap >= 0 && t_next_right_gap >= 0) {
                    /* non-overlapped case */
                    t_next_score = alignment_SG_Border ( data_query, TEL_POS(t), t_next_left_gap,
                                                         data_text , TER_POS(t), t_next_right_gap,
                                                         gp_xdrop + next_seedscore,
                                                         feature);
                  } else {
                    /* overlapped case */
                    long int t_next_max_neg_gap = MIN(t_next_left_gap,t_next_right_gap);
                    t_next_score =
                      gp_cost_max_substitution_matrix * t_next_max_neg_gap
                      +
                      gp_cost_gap_opened
                      +
                      gp_cost_gap_continued * ABS(t_next_left_gap - t_next_right_gap);
                  }


                  /*
                   *   last           t        ...       tnext
                   *  [====] ----- [=====] ------------ [=====]
                   *         <--->         <---------->
                   *         score <-----> t_next_score <----->
                   *               seedscore            next_seedscore
                   *         <-----------> <------------------>
                   *            dscore        t_next_dscore
                   *         <------------------------>
                   *                 l_next_score
                   *         <-------------------------------->
                   *                 l_next_dscore
                   */

#ifdef DEBUG_ALIGNTUPLES
                  printf("   - next tuple :  next_seedscore = %ld, t_next_score = %ld (l_next_score = %ld) :\n",next_seedscore,t_next_score,l_next_score);
                  DisplayTuple(tnext);
                  fflush(NULL);
#endif

                  /* 3) check if this is interesting to "jump" to t_next */
                  {
                    long int t_next_dscore =  t_next_score + next_seedscore;
                    if (l_next_dscore > dscore + t_next_dscore) {
                      /* remove from "last_tuple_aligned" (non included) so from "t" to "prev_next" ("tnext" non included) */
                      RESOLVE_OVERLAPS(prev_last_tuple_aligned,last_tuple_aligned,tnext);
                      CleanTupleList(tuple_list, last_tuple_aligned, t, tnext_prev);
                      t                        = tnext;

#ifdef DEBUG_ALIGNTUPLES
                      printf("     *** selected as t : l_next_dscore (%ld) > dscore (%ld) + t_next_dscore (%ld) ***\n",
                             l_next_dscore,dscore,t_next_dscore);
                      printf("     <- end of tnext loop [nt]\n");
                      fflush(NULL);
#endif
                      goto nt;
                    } else {
#ifdef DEBUG_ALIGNTUPLES
                      printf("     ### not selected as t ###\n");
                      fflush(NULL);
#endif
                    }
                  }
                }
              } /* if (l_next_score > 0) */
            } /* local block */
          } /* if (next_seedscore > seedscore) */
          tnext_prev = tnext;
          tnext = tnext->next;
        } /* while */

          /* No tnext improves this so keep link "last -> t" */
        RESOLVE_OVERLAPS(prev_last_tuple_aligned,last_tuple_aligned,t);

#ifdef DEBUG_ALIGNTUPLES
        printf("#<-> :  linking last_tuple_aligned to t ::\n");
        DisplayTuple(last_tuple_aligned);
        DisplayTuple(t);
        fflush(NULL);
#endif

        prev_total_score        = total_score;
        total_score            += dscore;
        prev_last_tuple_aligned = last_tuple_aligned;
        last_tuple_aligned      = t;
        t = t ->next;
        goto nt;
      } else {
#ifdef DEBUG_ALIGNTUPLES
        printf("# -- (III) non increasing scores\n");
        fflush(NULL);
#endif

        /* III) non increasing score : greedy strategy with only one step further or backward is enought (fast and "almost" accurate) */

        /* 1) PREVIOUS ONE : does not "greedy" increase between "last_tuple_aligned" and "t" : try with "prev_last_tuple_aligned" and "t" */

        if (prev_last_tuple_aligned) {
#ifdef DEBUG_ALIGNTUPLES
          printf("# -- (prev) :\n");
          DisplayTuple(prev_last_tuple_aligned);
          fflush(NULL);
#endif

          long int prev_left_gap    =  TGAP_L(prev_last_tuple_aligned, t);
          long int prev_right_gap   =  TGAP_R(prev_last_tuple_aligned, t);
          if (prev_left_gap <= gp_rho_stat && prev_right_gap <= gp_rho_stat ) {

            long int prev_score = 0;
            if (prev_left_gap >= 0 && prev_right_gap >= 0) {
              /* non-overlapped case */
              prev_score = alignment_SG_Border ( data_query, TEL_POS(prev_last_tuple_aligned), prev_left_gap,
                                                 data_text , TER_POS(prev_last_tuple_aligned), prev_right_gap,
                                                 gp_xdrop + seedscore,
                                                 feature);
            } else {
              /* overlapped case */
              long int prev_max_neg_gap = MIN(prev_left_gap,prev_right_gap);
              prev_score =
                gp_cost_max_substitution_matrix * prev_max_neg_gap
                +
                gp_cost_gap_opened
                +
                gp_cost_gap_continued * ABS(prev_left_gap - prev_right_gap);
            }
            {
              long int prev_dscore  = prev_score + seedscore;
              if (prev_total_score + prev_dscore >= MAX(total_score,seedscore)) {
                /* jumping last_tuple_aligned "gready" increases the score (remove it) */
                tuple * t_dummy = NULL;
                RESOLVE_OVERLAPS(t_dummy,prev_last_tuple_aligned,t);

#ifdef DEBUG_ALIGNTUPLES
                printf("#<*> :  linking {prev}_last_tuple_aligned to [t] :\n");
                DisplayTuple(t);
                fflush(NULL);
#endif

                FREE(last_tuple_aligned,sizeof(tuple));
                prev_last_tuple_aligned->next = t;

                last_tuple_aligned            = t;
                total_score                   = prev_total_score + prev_dscore;
                t = t->next;
                goto nt;
              }
            }
          }
        }
        {
          tuple *tnext = t->next;

          /* 2) NEXT ONE : does not "gready" increase between "last_tuple_aligned" and "t" neither
           * between "prev_last_tuple_aligned" and "t" : try with "last_tuple_aligned" and "t"
           */
          if (tnext) {
            long int next_left_gap  =  TGAP_L( last_tuple_aligned, tnext);
            long int next_right_gap =  TGAP_R( last_tuple_aligned, tnext);

#ifdef DEBUG_ALIGNTUPLES
            printf("   - (next) :\n");
            DisplayTuple(tnext);
            fflush(NULL);
#endif

            if (next_left_gap <= gp_rho_stat && next_right_gap <= gp_rho_stat ) {
              long int next_seedscore  = alignment_SG_Strait (data_query, TBL_POS(tnext),
                                                              data_text , TBR_POS(tnext),
                                                              TSIZE(tnext),
                                                              feature);

              long int next_score      = 0;
              if (next_left_gap >= 0 && next_right_gap >= 0) {
                /* non-overlapped case */
                next_score = alignment_SG_Border ( data_query, TEL_POS(last_tuple_aligned), next_left_gap,
                                                   data_text , TER_POS(last_tuple_aligned), next_right_gap,
                                                   gp_xdrop + seedscore,
                                                   feature);
              } else {
                /* overlapped case */
                long int next_max_neg_gap = MIN(next_left_gap,next_right_gap);
                next_score =
                  gp_cost_max_substitution_matrix * next_max_neg_gap
                  +
                  gp_cost_gap_opened
                  +
                  gp_cost_gap_continued * ABS(next_left_gap - next_right_gap);
              }
              {
                long int next_dscore  = next_score + next_seedscore;
                if (total_score + next_dscore >=  MAX(total_score,next_seedscore)) {
                  /* jumping t "gready "increases" the score (remove it)*/
                  RESOLVE_OVERLAPS(prev_last_tuple_aligned,last_tuple_aligned,tnext);

#ifdef DEBUG_ALIGNTUPLES
                  printf("#<*> :  linking last_tuple_aligned to {tnext} :\n");
                  DisplayTuple(tnext);
                  fflush(NULL);
#endif


                  FREE(t,sizeof(tuple));
                  last_tuple_aligned->next = tnext;

                  prev_last_tuple_aligned  = last_tuple_aligned;
                  last_tuple_aligned       = tnext;
                  prev_total_score         = total_score;
                  total_score             += next_dscore;
                  last_seedscore           = next_seedscore;
                  t = tnext->next;
                  goto nt;
                }
              }
            }
          }

#ifdef DEBUG_ALIGNTUPLES
          printf("# -- (delete)\n");
          fflush(NULL);
#endif

          /* 3) keep the best score of the two if this cannot be solved
           *    the less interesting will be saved without extension
           *    if its score is > TUminscore.
           */

          if (seedscore < total_score) {
#ifdef DEBUG_ALIGNTUPLES
            printf("# -- (seedscore < total_score)\n");
            fflush(NULL);
#endif
            if (seedscore > feature->TUminscore) {
              /* FIXME : No left/right extends here */
#ifdef DEBUG_ALIGNTUPLES
              printf("# -- (mem)\n");
              fflush(NULL);
#endif
              MEMORISE_IN_DISPLAY_OR_ERASE(last_tuple_aligned,t,t,seedscore);
            } else {
#ifdef DEBUG_ALIGNTUPLES
              printf("# -- [[<<delete t>>]] :\n");
              DisplayTuple(t);
              fflush(NULL);
#endif
              last_tuple_aligned->next = tnext;
              FREE(t,sizeof(tuple));
            }
            t = tnext;
            goto nt;
          } else {
#ifdef DEBUG_ALIGNTUPLES
            printf("# -- (seedscore > total_score)\n");
            fflush(NULL);
#endif

            if (total_score > feature->TUminscore) {
#ifdef DEBUG_ALIGNTUPLES
              printf("# -- (mem)\n");
              fflush(NULL);
#endif
              /* FIXME : No left/right extends here */
              MEMORISE_IN_DISPLAY_OR_ERASE(NULL,first_tuple_aligned,last_tuple_aligned,total_score);
            } else {
#ifdef DEBUG_ALIGNTUPLES
              printf("# -- [[<<delete fta->lta !!>>]]:\n");
              DisplayTuple(first_tuple_aligned);
              DisplayTuple(last_tuple_aligned);
              fflush(NULL);
#endif
              CleanTupleList(tuple_list,NULL,first_tuple_aligned,last_tuple_aligned);
            }
            first_tuple_aligned      = t;
            last_tuple_aligned       = t;
            total_score              = seedscore;
            last_seedscore           = seedscore;
            prev_last_tuple_aligned  = NULL;
            t = tnext;
            goto nt;
          }
        }
      }  /* (II) or (III) */
    }
  nt:;
  } /* while (t) */
  MEMORISE_IN_DISPLAY_OR_ERASE(NULL,first_tuple_aligned,last_tuple_aligned,total_score);
#ifdef DEBUG_ALIGNTUPLES
  printf("# ------------------------------------------- \n\n");
  fflush(NULL);
#endif

  return 0;
}







long int AlignAndFree(
                 char * data_query, long int datasize_query,
                 char * data_text,  long int datasize_text,
                 long int force,
                 Feature *feature) {

  /* list of lists manipulating vars */
  tuplelist * tl_ = NULL , *tl_prev_ = NULL, *tl_last_;
  tuple     * t_  = NULL;
  long int       countreturn = 5;
  /* statistical vars */
  long int number_of_tuples = 0;


  /* head lists selection and processing (alignement) :
   *
   *  Those lists are the older ones, so if i_list < gv_i_current - threshold,
   *  those lists will not be used by assemble
   *  we can process those and free some memory space
   */


  /*[a] tuplelists start  */
  long int i_current  = feature->i_current;
  tl_prev_ = tl_ = feature->first_tl;
  tl_last_ =       feature->last_tl;

  while ((force && (tl_ != NULL)) || (!force && (tl_ != tl_last_))) {

#ifdef PREFETCH
    __builtin_prefetch(tl_->next,0);
#endif

    t_ = tl_->first_tuple;
    number_of_tuples = 0;

    /* (1) list used by ASSEMBLE and ALIGN (beware of reading/writing) */
    while (t_ != NULL) {

#ifdef PREFETCH
      __builtin_prefetch(t_->next,0);
#endif

      /* if list goes further, we do not consider it
       * and we stop this algorithm :
       * - all the following list have also chances to
       *   have their last element inside consideration
       *   zone
       */

      if (!force)
        if (t_->occurrence >= i_current - gp_border) {
          countreturn--;
          if (!countreturn)
            return 0;
          else
            goto nexttuplelist;
        }

      number_of_tuples++;
      t_ = t_->next;
    }/* (1) */

    /* (2) specific consideration of each tuple list */

    /* lists iteration has not given any element in the active update
     * zone. Those alignements have been "completed" for the assemble
     * algorithm and thus can be aligned
     */

    if (number_of_tuples <= 0) {

      /* (2.0) The first list is always empty (do nothing
       * and jump to the next list)
       */
    nexttuplelist:
      tl_prev_ = tl_;
      tl_      = tl_->next;

    } else {

      /* (2.1) alignements of one of several tuples
       */

      AlignTuples(tl_,
                  data_query, datasize_query,
                  data_text , datasize_text,
                  feature);

      if (!force)
        FreeTupleList(&tl_,&tl_prev_,&(tl_last_)); /* however tl_last_ must never be reached */
      else
        FreeTupleList(&tl_,&tl_prev_,&(feature->last_tl));
    }
  }/* [a] tuple list iteration */

  return 0;
}




/*
 *   Function to be lauched by threads
 */

#ifdef THREAD_ASSEMBLE_ALIGN
#if defined(WIN32) || defined(WIN64)
DWORD WINAPI thread_work_align(PVOID fvoid)
#else
     void * thread_work_align(void * fvoid)
#endif

#else
     void * thread_work_align(void * fvoid)
#endif
{
  Feature * f = (Feature *) fvoid;

  /*
   * [] Align work done here
   */

  AlignAndFree(
               f->chunk_query,
               f->chunk_query_size,
               gp_text + gp_chunkstrt_text[f->i_chunk],
               gp_chunksize_text[f->i_chunk],
               0,
               f);

#ifdef THREAD_ASSEMBLE_ALIGN
  END_THREAD();
#else
  return 0;
#endif
}
