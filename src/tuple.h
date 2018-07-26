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

#ifndef __TUPLE_H_
#define __TUPLE_H_

/*
 * It represent an ungapped alignment it term of
 * "diagonal" and "occurrence" on the right sequence
 *
 * - "left" is the "lenght" of already matching hits on the
 * same diagonal (a tuple is always on one single diagonal).
 */

typedef struct _tuple {
  struct _tuple *next;
  long int occurrence;
  long int diagonal;
  long int leftsize;
} tuple;

typedef struct _tuplelist {
  struct _tuplelist *next;
  struct _tuple *first_tuple;
} tuplelist;


/*
 * MA : "Memorized Alignment"
 *       represents an alignment in term
 *       of left/rigth (query/text) sequences and
 *       begining-ending position on each sequence
 *
 *       In order to quicky align sequences,
 *       chain of tuples is kept as a "backbone"
 *       for the alignment
 */

typedef struct _MA {

  struct _MA * next;
  struct _tuple *first_tuple;

  /* position of memorized alignment */
  long int left_pos_begin;
  long int right_pos_begin;
  long int left_pos_end;
  long int right_pos_end;

  /* blastscore */
  long int left_blastscore;
  long int right_blastscore;
  long int blastscore;
  float entropy;
  float mutual;

  /* stats */
  long int trinomial_count1[3];
  long int transindels[7]; /*
                      * #transition, #transversions , #indels,
                      * #indels query, #indel text,
                      * #block indels query, #block indels text
                      *
                      */

  long int trinomial_count2[3];

  /* chunk used on the first sequence (query) */
  long int j_chunk;

  /* chunk used on the second sequence (text) */
  long int i_chunk;
  /* reverse (boolean) */
  unsigned char reverse;

} MA;



/*
 *      LEFT SEQUENCE : QUERY                               RIGHT SEQUENCE : TEXT
 *
 *
 *
 *
 *  TBL_POS(t1)   TEL_POS(t1)                            TBR_POS(t1)  TER_POS(t1)
 *     |            |                                       |             |
 *     |            |                                       |             |
 *-----[=========#=]*----[====#=]---------------------------[==========#=]*-[====#=]-------------
 *               |            |                                        |         |
 *               |            |                                        |         |
 *          TL_POS(t1)  TL_POS(t2)                              TR_POS(t1)    TR_POS(t2)
 *
 *
 *
 *
 *
 *
 *        TSIZE(t1)
 *     <----------->
 *     .           .
 *     .           .
 *-----[=========#=]-----[====#=]---------------------------[==========#=]--[====#=]-------------
 *     .           .     .      .                           .            .  .      .
 *     .           .     .      .                           .            .  .      .
 *     .           .     .      .                           .            .  .      .
 *     .           <----->      .                           .            .  .      .
 *     .         TGAP_L(t1,t2)  .                           .           ->--<-     .
 *     .                        .                           .     TGAP_R(t1,t2)    .
 *     .                        .                           .                      .
 *     <------------------------>                           .                      .
 *            TSEG_L(t1,t2)                                 <---------------------->
 *                                                               TSEG_R(t1,t2)
 */


#include "global_var.h"
#include "threads.h"

/* position on the text (R as right) */
#define TR_POS(t) ((t)->occurrence)

/* position on the query (L as left)
 *   TL_POS should be modified according to "left_correction"
 *   whith can be either 0 or keysize
 *
 */
#define TL_POS1(t) ((t)->occurrence - (t)->diagonal)
#define TL_POS(t)  ((t)->occurrence - (t)->diagonal + left_correction)


/* position begin "B" and end "E" on each Left and Rigth seed */

#define TBR_POS(t)  (TR_POS(t) - (t)->leftsize)
#define TBL_POS(t)  (TL_POS(t) - (t)->leftsize)

/* End position + 1 !!! */
#define TER_POS(t)  (TR_POS(t))
#define TEL_POS(t)  (TL_POS(t))


/* length */

/* length of one segment (t1<t2) between 2 seeds (including seeds size) */
#define TSEG_L(t1,t2) (TL_POS1(t2) - TL_POS1(t1) + (t1)->leftsize)
#define TSEG_R(t1,t2) (TR_POS(t2)  - TR_POS(t1)  + (t1)->leftsize)

/* length of one gap     (t1<t2) between 2 seeds (excluding seeds size) */
#define TGAP_L(t1,t2) (TL_POS1(t2) - TL_POS1(t1) - (t2)->leftsize)
#define TGAP_R(t1,t2) (TR_POS(t2)  - TR_POS(t1)  - (t2)->leftsize)

/* length of one "unsized" gap  (t1<t2) between 2 seeds (seeds size) */
#define TUSGAP_L(t1,t2) (TL_POS1(t2) - TL_POS1(t1))
#define TUSGAP_R(t1,t2) (TR_POS(t2)  - TR_POS(t1))

/* full seed size */
#define TSIZE(t)  ((t)->leftsize)


/*
 * [1] tuple/tuplelists/ma (memorized alignments) management
 */


/* Build a tuple */
tuple     * CreateTuple(long int occurrence, long int diagonal, long int size);

#define CREATETUPLE(__tuple__,__occurrence__,__diagonal__,__leftsize__) { \
  tuple * __tuple_tmp__ = (tuple *) MALLOC(sizeof(tuple));                \
  ASSERT(__tuple_tmp__, CreateTuple);                                     \
  __tuple_tmp__->occurrence = __occurrence__;                             \
  __tuple_tmp__->diagonal   = __diagonal__;                               \
  __tuple_tmp__->leftsize   = __leftsize__;                               \
  __tuple_tmp__->next       = NULL;                                       \
  __tuple__                 = __tuple_tmp__;                              \
}

/* Build an empty tuple list */
tuplelist * CreateTupleList();

#define CREATETUPLELIST(_tl,_t,_ttl) {                                      \
  tuplelist * _tl_tmp = (tuplelist *) MALLOC((unsigned) sizeof(tuplelist)); \
  ASSERT(_tl_tmp, CreateTupleList);                                         \
  _tl_tmp->first_tuple = _t;                                                \
  _tl_tmp->next        = _ttl;                                              \
  _tl                  = _tl_tmp;                                           \
}


/* Build a Memorized alignment */
MA        * CreateMA(long int left_pos_begin,  long int right_pos_begin,
                     long int left_pos_end,    long int right_pos_end,
                     long int left_blastscore, long int right_blastscore,
                     tuplelist * tuple_list,
                     tuple * prev_first_tuple_aligned,
                     tuple * first_tuple_aligned,
                     tuple * last_tuple_aligned,
                     long int blastscore, double entropy,
                     Feature * feature
                 );
/* Free a single tuple */
void FreeTuple(tuple * t);
/* Free a list of tuples (based on next chain) */
void FreeTuples(tuple * t);
/* Free a complete MA and its tuple list */
void FreeMA(MA * ma);

/* Free a tuplelist pointed out by p_tl (and its list of tuple),
 * you have to pass "previous" tuplelist element and also a "last_tl"
 * ---
 * Specific function only used in "align"
 */
void FreeTupleList(tuplelist ** p_tl, tuplelist ** p_tl_prev,tuplelist **last_tl);
/* Cut and paste the interesting chain of "tl" ("first,last") into "ma" */
void KeepTupleList(tuplelist * tl, tuple * prev_tfirst, tuple * tfirst, tuple * tlast, MA * ma);
/* Delete the non-interesting chain of "tl"  ("first,last") */
void CleanTupleList(tuplelist * tl, tuple * prev_tfirst, tuple * tfirst, tuple * tlast);
/* Left correction is applied */
void LeftCorrection_MA(MA * ma, long int left_correction);



/*
 * [2] Sorting/filtering functions
 */

/* sorting functions : for one MA */
typedef long int (SortCrit) (MA *);
long int SortCriterionScore(MA * ma);
long int SortCriterionEntropy(MA * ma);
long int SortCriterionMutual(MA * ma);
long int SortCriterionScoreWithEntropy(MA * ma);
long int SortCriterionQueryBegin(MA * ma);
long int SortCriterionTextBegin(MA * ma);
long int SortCriterionPercentIdentityAlign(MA * ma);
long int SortCriterionPercentIdentityQuery(MA * ma);
long int SortCriterionPercentIdentityText(MA * ma);



/* sorting functions : gives blocks consistency */
typedef long int (SortBlocksCrit) (MA *);
long int SortBlocksCriterionQueryNumber(MA * ma);
long int SortBlocksCriterionTextNumber(MA * ma);
long int SortBlocksCriterionQueryTextNumber(MA * ma);
long int SortBlocksCriterionTextQueryNumber(MA * ma);

/* recursive function */

void ListSort_MAList(MA ** p_firstMA, MA ** p_lastMA, long int mask /* 0x40000000 */, long int flag_sortblocks);
/* main sorting function called by one thread */
void Sort_MAList(Feature * f);

/* merge results provided by several threads */

/* merge the current list with "gv_first/last_MA" (global) */
void MergeSort_MAList(MA * first_MA, MA * last_MA);
/* merge two reverse and forward features : does not affect stats ... */
void MergeSort_forward_reverse(Feature * f1, Feature * f2, MA ** p_first_MA, MA ** p_last_MA);

/* filtering functions : for one MA */
long int EntropyAndScoreFilter_MA(Feature * feature, MA * ma, char *dataquery,  char *datachunktext);
/* main filtering function called by one thread */
long int EntropyAndScoreFilter_MAList(Feature * feature, char *dataquery, char *datatext, long int *datatext_start);



/*
 * [3] Display (provided for debugging)
 */
long int DisplayListTupleList(tuplelist * tuplelist);
long int DisplayTupleList(tuplelist * tuplelist);
long int DisplayListTuple(tuplelist * tuple);
long int DisplayTuple(tuple * tuple);
long int DisplayMA(MA * ma);
long int DisplayListMA (MA * ma);

#endif
