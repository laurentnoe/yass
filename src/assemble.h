/*
 *  YASS 1.16
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

#ifndef __ASSEMBLE_H_
#define __ASSEMBLE_H_

#include "tuple.h"
#include "threads.h"

long int initialise_deltashift();
long int Assemble_Single( /*in */ char *data, /*in */ long int datasize,
                    Feature  *f
                    );
long int Assemble_SingleRev( /*in */ char *datarev, /*in */ char *data, /*in */long int datasize,
                       Feature  *f
                       );


long int Assemble_Double( /*in */ char *data_query, /*in */ long int datasize_query,
                          /*in */ char *data_text,  /*in */ long int datasize_text,
                            long int reverse_repeat,
                          /*out*/ MA **first_MA , /*out*/ MA **last_MA);

long int MultiAssemble_Double( /*in */ char *data_query,        /*in */ long int datasize_query,
                         /*in */ char *data_text,         /*in */ long int datasize_text,
                         /*in */ long int nbchunks_text,  /*in */ char **chunkname_text,
                         /*in */ long int *chunksize_text,/*in */ long int *chunkstrt_text,
                         Feature  *f);




#define SINGLEHITDIAGONAL(data_query,data_query_size,data_text,data_text_size) {                      \
  long int j_last = GET_TAB_MIN(last_tuple_pos_with_diag,diagonal,f->i_current - gp_rho_stat);        \
  if (j_last < 0 || j_last < f->i_current - gp_rho_stat) {                                            \
    long int max_right_score = 0;                                                                     \
    long int max_left_score  = 0;                                                                     \
    long int max_i = i_current_end;                                                                   \
    long int min_i = i_current_end;                                                                   \
    {                                                                                                 \
      long int i = i_current_end;                                                                     \
      long int j = i_previous;                                                                        \
      long int right_score = 0;                                                                       \
      while (i < data_text_size && j < data_query_size) {                                             \
        right_score += gp_substitution_matrix[(long int)((data_query)[j])][(long int)(data_text)[i]]; \
        if (right_score > max_right_score) {                                                          \
          max_i           = i+1;                                                                      \
          max_right_score = right_score;                                                              \
        } else                                                                                        \
          if (right_score < -gp_xdrop-max_right_score)                                                \
            break;                                                                                    \
        i++;j++;                                                                                      \
      }                                                                                               \
    }                                                                                                 \
    {                                                                                                 \
      long int i = i_current_end;                                                                     \
      long int j = i_previous;                                                                        \
      long int left_score = 0;                                                                        \
      while (j > 0 && i > 0) {                                                                        \
        i--;j--;                                                                                      \
        left_score += gp_substitution_matrix[(long int)((data_query)[j])][(long int)(data_text)[i]];  \
        if (left_score > max_left_score) {                                                            \
          min_i          = i;                                                                         \
          max_left_score = left_score;                                                                \
        } else                                                                                        \
          if (left_score < -gp_xdrop-max_left_score)                                                  \
            break;                                                                                    \
      }                                                                                               \
    }                                                                                                 \
    STATS_NB_SINGLE_TESTS(f);                                                                         \
    if (max_left_score + max_right_score >  f->TUminscore) {                                          \
      tuple * t;                                                                                      \
      STATS_NB_SINGLE_HITS(f);                                                                        \
      STATS_NB_CHAINS_BUILT_INC(f);                                                                   \
      /* create a tuple */                                                                            \
      CREATETUPLE(t,max_i,diagonal,max_i-min_i);                                                      \
      PUT_TAB(last_tuple_ptr_with_diag,diagonal,t);                                                   \
      PUT_TAB(last_tuple_pos_with_diag,diagonal,max_i);                                               \
      /* create a tuple list element */                                                               \
      CREATETUPLELIST(f->last_tl->next,t,NULL);                                                       \
      f->last_tl = f->last_tl->next;                                                                  \
      goto next_key;                                                                                  \
    }                                                                                                 \
  } /* j_last < i_current_end - gp_rho_stat */                                                        \
}


#define SINGLEHITDIAGONAL_MULTI(data_query,data_query_size,data_text,data_text_size) {                \
  long int j_last = GET_TAB_MIN(last_tuple_pos_with_diag,diagonal,i_current_end - gp_rho_stat);       \
  if (j_last < 0 || j_last < i_current_end - gp_rho_stat ) {                                          \
    long int max_right_score = 0;                                                                     \
    long int max_left_score  = 0;                                                                     \
    long int max_i = i_current_end;                                                                   \
    long int min_i = i_current_end;                                                                   \
    {                                                                                                 \
      long int i = i_current_end;                                                                     \
      long int j = i_previous;                                                                        \
      long int right_score = 0;                                                                       \
      while (i < data_text_size && j < data_query_size) {                                             \
        right_score += gp_substitution_matrix[(long int)((data_query)[j])][(long int)(data_text)[i]]; \
        if (right_score > max_right_score) {                                                          \
          max_i           = i+1;                                                                      \
          max_right_score = right_score;                                                              \
        } else                                                                                        \
          if (right_score < -gp_xdrop-max_right_score)                                                \
            break;                                                                                    \
        i++;j++;                                                                                      \
      }                                                                                               \
    }                                                                                                 \
    {                                                                                                 \
      long int i = i_current_end;                                                                     \
      long int j = i_previous;                                                                        \
      long int left_score = 0;                                                                        \
      while (j > 0 && i > 0) {                                                                        \
        i--;j--;                                                                                      \
        left_score += gp_substitution_matrix[(long int)((data_query)[j])][(long int)(data_text)[i]];  \
        if (left_score > max_left_score) {                                                            \
          min_i          = i;                                                                         \
          max_left_score = left_score;                                                                \
        } else                                                                                        \
          if (left_score < -gp_xdrop-max_left_score)                                                  \
            break;                                                                                    \
      }                                                                                               \
    }                                                                                                 \
    STATS_NB_SINGLE_TESTS(f);                                                                         \
    if (max_left_score + max_right_score >  f->TUminscore) {                                          \
      tuple * t;                                                                                      \
      STATS_NB_SINGLE_HITS(f);                                                                        \
      STATS_NB_CHAINS_BUILT_INC(f);                                                                   \
      /* create a tuple */                                                                            \
      CREATETUPLE(t,max_i,diagonal,max_i-min_i);                                                      \
      PUT_TAB(last_tuple_ptr_with_diag,diagonal,t);                                                   \
      PUT_TAB(last_tuple_pos_with_diag,diagonal,max_i);                                               \
      /* create a tuple list element */                                                               \
      CREATETUPLELIST(f->last_tl->next,t,NULL);                                                       \
      f->last_tl = f->last_tl->next;                                                                  \
      /* >> */                                                                                        \
      if (nb_last_chunk_diag_used < MAX_DIAG_POS_KEPT_DATACHUNK)                                      \
        last_chunk_diag_used[nb_last_chunk_diag_used++] = diagonal;                                   \
      /* << */                                                                                        \
      goto next_key;                                                                                  \
    }                                                                                                 \
  } /* j_last < i_current_end - gp_rho_stat */                                                        \
}


#endif
