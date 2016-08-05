
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_var.h"
#include "tuple.h"
#include "prdyn.h"
#include "proba.h"
#include "kword.h"


#ifdef INLINE
inline
#endif
tuple * CreateTuple(long int occurrence, long int diagonal, long int size)
{

    tuple *t = (tuple *) MALLOC(sizeof(tuple));
    ASSERT(t, CreateTuple);
    t->occurrence = occurrence;
    t->diagonal   = diagonal;
    t->leftsize   = size;
    t->next       = NULL;
    return t;
}


#ifdef INLINE
inline
#endif
tuplelist * CreateTupleList()
{
  tuplelist *tl = (tuplelist *) MALLOC(sizeof(tuplelist));
  ASSERT(tl, CreateTupleList);
  tl->first_tuple = NULL;
  tl->next = NULL;
  return tl;
}



#ifdef INLINE
inline
#endif
void KeepTupleList(tuplelist * tl, tuple * t_prevfirst, tuple * tfirst, tuple * tlast, MA * ma)
{
  tuple * tlast_next = tlast->next;
  if (t_prevfirst)
    t_prevfirst->next  = tlast_next;
  else
    tl->first_tuple    = tlast_next;

  tlast->next = NULL;

  ma->first_tuple = tfirst;
}

#ifdef INLINE
inline
#endif
void CleanTupleList(tuplelist * tl, tuple * t_prevfirst, tuple * tfirst, tuple * tlast)
{
  tuple * tlast_next = tlast->next;
  tuple * t = tfirst;

  if (t_prevfirst)
    t_prevfirst->next  = tlast_next;
  else
    tl->first_tuple    = tlast_next;

  while (t != tlast_next) {
    tuple * t_prev = t;
    t = t->next;
    FREE(t_prev,sizeof(tuple));
  }
}



#ifdef INLINE
inline
#endif
void FreeTupleList(tuplelist ** p_tl, tuplelist ** p_tl_prev, tuplelist **last_tl) {
    tuple *t, *t_prev;

    if (*last_tl == *(p_tl))
        *last_tl = *(p_tl_prev);

    t = (*(p_tl))->first_tuple;
    while (t != NULL) {
        t_prev = t;
        t = t->next;
        FREE(t_prev,sizeof(tuple));
    }

    (*(p_tl_prev))->next = (*(p_tl))->next;
    FREE(*(p_tl),sizeof(tuplelist));

    *(p_tl) = (*(p_tl_prev))->next;
}



#ifdef INLINE
inline
#endif
void FreeTuples(tuple * t)
{
  tuple * t_ =  t, * t_prev_;
  while (t_ != NULL) {
    t_prev_ = t_;
    t_ = t_->next;
    FREE(t_prev_,sizeof(tuple));
  }
}



#ifdef INLINE
inline
#endif
MA * CreateMA(long int left_pos_begin,  long int right_pos_begin,
              long int left_pos_end,    long int right_pos_end,
              long int left_blastscore, long int right_blastscore,
              tuplelist * tuple_list,
              tuple * prev_first_tuple_aligned,
              tuple * first_tuple_aligned,
              tuple * last_tuple_aligned,
              long int blastscore, double entropy,
              Feature * feature
                  )
{

    MA * ma = NULL;

    ma = (MA *) MALLOC(sizeof(MA));
    ASSERT(ma, CreateMA);

    /* 1) initialise ma */
    ma->next = NULL;
    ma->left_pos_begin   = left_pos_begin;
    ma->right_pos_begin  = right_pos_begin;
    ma->left_pos_end     = left_pos_end;
    ma->right_pos_end    = right_pos_end;
    ma->left_blastscore  = left_blastscore;
    ma->right_blastscore = right_blastscore;
    ma->blastscore       = blastscore;
    ma->entropy          = entropy;

    /* 2) keep the seeds chain inside MA */
    KeepTupleList(tuple_list, prev_first_tuple_aligned, first_tuple_aligned, last_tuple_aligned, ma);

    ma->reverse = (char) feature->reverse;
    ma->i_chunk = feature->i_chunk;
    ma->j_chunk = feature->j_chunk;

    /* correction due to the algorithm "Assemble_Double" over "Assemble_Single" */

    /* left_correction take effect on TL_POS(t): we can do alignments before calling the "left_correction"
     * function on each MA, while left_correction global variable is != 0, or either after the "left_correction()"
     * call provided that left_correction global variable IS = 0.
     */


    LeftCorrection_MA(ma,feature->left_correction);
    /* 3) add  "ma" in the MA_linked_list */
    if (feature->first_MA == NULL) {
        feature->first_MA = ma;
    } else {
        feature->last_MA->next = ma;
    }
    feature->last_MA = ma;

    STATS_NB_MA_INC(feature);

    return ma;
}

#ifdef INLINE
inline
#endif
void FreeMA(MA * ma)
{
    tuple *t_, *t_prev_;
    t_ = ma->first_tuple;
    while (t_) {
        t_prev_ = t_;
        t_ = t_->next;
        FREE(t_prev_,sizeof(tuple));
    }
    FREE(ma,sizeof(MA));
}

#ifdef INLINE
inline
#endif
void LeftCorrection_MA(MA * ma,long int left_correction)
{
    tuple *t = ma->first_tuple;
    while (t) {
        t->diagonal -= left_correction;
        t = t->next;
    }
}

/*

long int QuickSort_MAList(MA ** p_firstMA, MA ** p_lastMA, SortFunct * f )
{


  MA * _ma_top    = *(p_firstMA);
  MA * _ma_1      = NULL;
  MA * _ma_1_top  = NULL;
  MA * _ma_0      = NULL;
  MA * _ma_0_top  = NULL;
  MA * _ma        = NULL;

  if (_ma_top) {

    _ma = _ma_top->next;

    while (_ma) {
      if (f(_ma,_ma_top)) {
        if (!_ma_1)
          _ma_1 = _ma;
        else
          _ma_1_top->next = _ma;
        _ma_1_top = _ma;
      } else {
        if (!_ma_0)
          _ma_0 = _ma;
        else
          _ma_0_top->next = _ma;
        _ma_0_top = _ma;
      }
      _ma = _ma->next;
    }


    if (_ma_0_top) {
      _ma_0_top->next = NULL;
      if (_ma_0 != _ma_0_top)
        QuickSort_MAList(&_ma_0,&_ma_0_top,f);
    }

    if (_ma_1_top) {
      _ma_1_top->next = NULL;
      if (_ma_1 != _ma_1_top)
        QuickSort_MAList(&_ma_1,&_ma_1_top,f);
    }




    _ma_top->next =_ma_0;
    if (_ma_1_top) {
      _ma_1_top->next = _ma_top;
      *(p_firstMA) = _ma_1;
    } else {
      *(p_firstMA) = _ma_top;
    }

    if (p_lastMA) {
      if (_ma_0_top)
        *(p_lastMA) = _ma_0_top;
      else
        *(p_lastMA) = _ma_top;
    }

    if (*(p_lastMA))
      (*(p_lastMA))->next = NULL;
  }
  return 0;
}

*/


/* ============================================================================ */

long int SortCriterionScore(MA * ma)            { return (ma->blastscore);  }
long int SortCriterionEntropy(MA * ma)          { return (long int)((ma->entropy)*1000);}
long int SortCriterionMutual(MA * ma)           { return (long int)((ma->mutual)*1000);}
long int SortCriterionScoreWithEntropy(MA * ma) { return (long int)((ma->entropy) * (ma->blastscore));}
long int SortCriterionQueryBegin(MA * ma)       { return (long int)((ma->reverse) ?
                                                                    (ma->left_pos_begin):
                                                                    (gp_chunksize_query[ma->j_chunk] - ( ma->left_pos_end)));}
long int SortCriterionTextBegin(MA * ma)        { return (long int) (0x7fffffff - ma->right_pos_begin);}
long int SortCriterionPercentIdentityAlign(MA * ma)  {
  return (long int)
    1000
    -
    (ma->transindels[0] +  ma->transindels[1] +  ma->transindels[2])  /*transitions+transversion+indels */
    *
    1000
    /
    (
     ma->transindels[3]
     +
     ma->right_pos_end
     -
     ma->right_pos_begin
     ); /* alignment full length */
}

long int SortCriterionPercentIdentityQuery(MA * ma) {
  return (long int)
    1000
    -
    (ma->transindels[0] +  ma->transindels[1] +  ma->transindels[3])  /* transitions+transversion+indels on query */
    *
    1000
    /
    (
     ABS(ma->left_pos_end - ma->left_pos_begin) /* query length */
    );
}

long int SortCriterionPercentIdentityText(MA * ma) {
  return (long int)
    1000
    -
    (ma->transindels[0] +  ma->transindels[1] +  ma->transindels[4])  /* transitions+transversion+indels on text */
    *
    1000
    /
    (
     ma->right_pos_end - ma->right_pos_begin /* text length */
    );
}


long int SortBlocksCriterionQueryNumber(MA * ma)        { return (long int) (ma->j_chunk); }
long int SortBlocksCriterionTextNumber(MA * ma)         { return (long int) (ma->i_chunk); }
long int SortBlocksCriterionQueryTextNumber(MA * ma)    { return (long int) (ma->j_chunk)+(ma->i_chunk)*gp_nbchunks_text; }
long int SortBlocksCriterionTextQueryNumber(MA * ma)    { return (long int) (ma->i_chunk)+(ma->j_chunk)*gp_nbchunks_query; }

/*
 * ListSort
 */

void ListSort_MAList(MA ** p_firstMA, MA ** p_lastMA, long int mask /* 0x40000000 */, long int flag_sortblocks)
{

    MA *ma_0 = NULL, *ma_1 = NULL;
    MA *ma_0_top = NULL, *ma_1_top = NULL;

    MA *_ma = *(p_firstMA);

    while (_ma) {
        if ((gp_sortcriterion_func(_ma)) & mask) {
            if (!ma_1)
                ma_1 = _ma;
            else
                ma_1_top->next = _ma;

            /* post-block sorting  (fast forward in sorted list) */
            if (flag_sortblocks) {
              while (_ma->next && (gp_sortblockscriterion_func(_ma->next) == gp_sortblockscriterion_func(_ma))) {
                _ma = _ma->next;
              }
            }
            ma_1_top = _ma;
        } else {
            if (!ma_0)
                ma_0 = _ma;
            else
                ma_0_top->next = _ma;

            /* post-block sorting (fast forward in sorted list) */
            if (flag_sortblocks) {
              while (_ma->next && (gp_sortblockscriterion_func(_ma->next) == gp_sortblockscriterion_func(_ma))) {
                _ma = _ma->next;
              }
            }
            ma_0_top = _ma;
        }
        _ma = _ma->next;
    }

    mask >>= 1;

    if (ma_0_top) {
        ma_0_top->next = NULL;
        if (mask > 0)
          ListSort_MAList(&ma_0, &ma_0_top, mask, flag_sortblocks);
    }
    if (ma_1_top) {
        ma_1_top->next = NULL;
        if (mask > 0)
          ListSort_MAList(&ma_1, &ma_1_top, mask, flag_sortblocks);
    }


    /* chain the two sublists into one */
    if (ma_1_top) {
        *(p_firstMA) = ma_1;
        if (ma_0_top) {
            ma_1_top->next = ma_0;
            *(p_lastMA) = ma_0_top;
        } else {
            *(p_lastMA) = ma_1_top;
        }
    } else {
        *(p_firstMA) = ma_0;
        *(p_lastMA)  = ma_0_top;

    }
    /* lastMA->next is set to NULL */
    if (*(p_lastMA))
      (*(p_lastMA))->next = NULL;
}





/*
 * main function used to sort MA
 */

void Sort_MAList(Feature * f) {

  if (f->first_MA) {

    if (gp_sortcriterion >= 2 * NBPOTENTIALCRITERIA) {

      /* [A] Sort with chunk selection */
      MA * _ma_after_end    = f->first_MA;
      MA * _ma_before_start = NULL;

      while (_ma_after_end) {
        MA *       _ma_end   = _ma_after_end;
        MA *       _ma_start = _ma_after_end;
        long int   i_chunk   = _ma_after_end->i_chunk;

        /* search sublist 's' with same chunk */
        _ma_after_end = _ma_after_end->next;
        while (_ma_after_end && (_ma_after_end->i_chunk == i_chunk)) {
          _ma_end        = _ma_after_end;
          _ma_after_end  = _ma_after_end->next;
        }

        /* sort this sublist if more than one element */
        if (_ma_start != _ma_end) {
          _ma_end->next = NULL;
          ListSort_MAList(&_ma_start,&_ma_end,0x40000000,FALSE);
          _ma_end->next = _ma_after_end;
        }

        /* insert the sublist once sorted */
        if (!_ma_before_start)
          f->first_MA            = _ma_start;
        else
          _ma_before_start->next = _ma_start;
        _ma_before_start = _ma_end;

        f->last_MA = _ma_end;
      }

    } else {

      /* [B] Sort without chunk selection */
      ListSort_MAList(&(f->first_MA),&(f->last_MA),0x40000000,FALSE);
    }
  }
}


/*
 * Merge with "gv_first/last_MA" (global)
 */

void MergeSort_MAList(MA * first_MA, MA * last_MA) {

  LOCK(merge_ma_mutex);

  if (first_MA && gv_first_MA) {

    if ( /* [FIXME] */
        (gp_sortblockscriterion == 1) || (gp_sortblockscriterion == 3) ||
        (gp_sortblockscriterion == 4) || (gp_sortblockscriterion >= 6)
       ) {

      /*
       * mergesort by query chunks : must be fast if the order of threads is preserved ...
       */

      if (gv_last_MA->j_chunk <  first_MA->j_chunk) { /* chain at the end */
        gv_last_MA->next =  first_MA;
        gv_last_MA       =  last_MA;
      } else if (gv_first_MA->j_chunk >  last_MA->j_chunk) { /* chain at the beginning */
        last_MA->next = gv_first_MA;
        gv_first_MA   = first_MA;
      } else {  /* chain at the middle (more costly) */


        MA *current_MA   = gv_first_MA;
        MA *previous_MA  = gv_first_MA;

        while (current_MA && current_MA->j_chunk < first_MA->j_chunk) {
          previous_MA  = current_MA;
          current_MA   = current_MA->next;
        }

        if (last_MA->next != NULL) {
          _WARNING("MergeSort_MALists() : \"query chunk sorted\" list is corrupted");
        }

        last_MA->next = current_MA;
        if (current_MA == gv_first_MA)
          gv_first_MA       = first_MA;
        else
          previous_MA->next = first_MA;
      }
    } else if ((gp_sortblockscriterion == 2) || (gp_sortblockscriterion == 5)) {

      /*
       * mergesort by text chunks : if alone, must be slower than the previous "query chunk" one
       */

      MA * _ma_1        = first_MA;
      MA * _ma_2        = gv_first_MA;
      MA * _ma_12_end   = NULL;

      if ( (_ma_1->i_chunk < _ma_2->i_chunk) ||
           (_ma_1->i_chunk == _ma_2->i_chunk && gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2)) ) {
        gv_first_MA = _ma_12_end = _ma_1;
        _ma_1       = _ma_1->next;
      } else {
        gv_first_MA = _ma_12_end = _ma_2;
        _ma_2       = _ma_2->next;
      }

      while (_ma_1 && _ma_2) {
        if ( (_ma_1->i_chunk < _ma_2->i_chunk) ||
             (_ma_1->i_chunk == _ma_2->i_chunk && gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2)) ) {
          _ma_12_end->next = _ma_1;
          _ma_12_end       = _ma_1;
          _ma_1 = _ma_1->next;
        } else {
          _ma_12_end->next = _ma_2;
          _ma_12_end       = _ma_2;
          _ma_2 = _ma_2->next;
        }
      }

      if (!_ma_1) {
        _ma_12_end->next = _ma_2;
        _ma_12_end       = gv_last_MA;
      } else {
        if (!_ma_2) {
          _ma_12_end->next = _ma_1;
          _ma_12_end       = last_MA;
        } else {
          _WARNING("MergeSort_MALists() : \"text chunk sorted\" list is corrupted");
        }
      }
      gv_last_MA = _ma_12_end;

    } else {

      /*
       * mergesort without any chunk consideration (only score, position ...) : slower
       */

      MA * _ma_1 = first_MA;
      MA * _ma_2 = gv_first_MA;
      MA * _ma_12_end   = NULL;


      if (gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2)) {
        gv_first_MA = _ma_12_end = _ma_1;
        _ma_1       = _ma_1->next;
      } else {
        gv_first_MA = _ma_12_end = _ma_2;
        _ma_2       = _ma_2->next;
      }

      while (_ma_1 && _ma_2) {
        if (gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2)) {
          _ma_12_end->next = _ma_1;
          _ma_12_end       = _ma_1;
          _ma_1    = _ma_1->next;
        } else {
          _ma_12_end->next = _ma_2;
          _ma_12_end       = _ma_2;
          _ma_2    = _ma_2->next;
        }
      }

      if (!_ma_1) {
        _ma_12_end->next = _ma_2;
        _ma_12_end       = gv_last_MA;
      } else {
        if (!_ma_2) {
          _ma_12_end->next = _ma_1;
          _ma_12_end       = last_MA;
        } else {
          _WARNING("MergeSort_MALists() : \"non chunk sorted\" list is corrupted");
        }
      }
      gv_last_MA = _ma_12_end;
    }
  } else {
    if (first_MA) {
      gv_first_MA = first_MA;
      gv_last_MA  =  last_MA;
    }
  }

  UNLOCK(merge_ma_mutex);
}



void MergeSort_forward_reverse(Feature * f1, Feature * f2, MA ** p_first_MA, MA ** p_last_MA) {

  MA * first_MA = NULL;
  MA * last_MA  = NULL;

  if (f1->first_MA && f2->first_MA) {

    if (gp_sortcriterion >= 2 * NBPOTENTIALCRITERIA) {

      /* [A] MergeSort with text chunk selection */

      MA * _ma_1 = f1->first_MA;
      MA * _ma_2 = f2->first_MA;

      if ( (_ma_1->i_chunk < _ma_2->i_chunk) ||
           ((_ma_1->i_chunk == _ma_2->i_chunk) && (gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2))) ) {
        first_MA = last_MA = _ma_1;
        _ma_1    = _ma_1->next;
      } else {
        first_MA = last_MA = _ma_2;
        _ma_2    = _ma_2->next;
      }

      while (_ma_1 && _ma_2) {
        if ( (_ma_1->i_chunk < _ma_2->i_chunk) ||
             ((_ma_1->i_chunk == _ma_2->i_chunk) && (gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2))) ) {
          last_MA->next = _ma_1;
          last_MA       = _ma_1;
          _ma_1    = _ma_1->next;
        } else {
          last_MA->next = _ma_2;
          last_MA       = _ma_2;
          _ma_2    = _ma_2->next;
        }
      }

      if (!_ma_1) {
        last_MA->next = _ma_2;
        last_MA = f2->last_MA;
      } else {
        if (!_ma_2) {
          last_MA->next = _ma_1;
          last_MA = f1->last_MA;
        } else {
          _WARNING("MergeSort_forward_reverse() : \"text chunk sorted\" list is corrupted");
          last_MA->next = NULL;
        }
      }


    } else {

      /* [B] MergeSort without text chunk selection */

      MA * _ma_1 = f1->first_MA;
      MA * _ma_2 = f2->first_MA;

      if (gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2)) {
        first_MA = last_MA = _ma_1;
        _ma_1    = _ma_1->next;
      } else {
        first_MA = last_MA = _ma_2;
        _ma_2    = _ma_2->next;
      }

      while (_ma_1 && _ma_2) {
        if (gp_sortcriterion_func(_ma_1) > gp_sortcriterion_func(_ma_2)) {
          last_MA->next = _ma_1;
          last_MA       = _ma_1;
          _ma_1    = _ma_1->next;
        } else {
          last_MA->next = _ma_2;
          last_MA       = _ma_2;
          _ma_2    = _ma_2->next;
        }
      }

      if (!_ma_1) {
        last_MA->next = _ma_2;
        last_MA = f2->last_MA;
      } else {
        if (!_ma_2) {
          last_MA->next = _ma_1;
          last_MA = f1->last_MA;
        } else {
          _WARNING("MergeSort_forward_reverse() : \"non text chunk sorted\" list is corrupted");
          last_MA->next = NULL;
        }
      }
    }
  } else {
    if (f1->first_MA) {
      first_MA = f1->first_MA;
      last_MA  = f1->last_MA;
    } else {
      first_MA = f2->first_MA;
      last_MA  = f2->last_MA;
    }
  }

  /* set first_MA and last_MA */
  *p_first_MA = first_MA;
  *p_last_MA  = last_MA;
  return;
}


/*
 * return true if the filter does apply on the current MA
 * and thus the MA has to be removed
 */

long int EntropyAndScoreFilter_MA(Feature * f, MA * ma, char *dataquery,  char *datachunktext)
{

  /* Score filter */
  if (ma->blastscore <= f->MAminscore)
    return 1;

  /* Entropy filter */
  memset(ma->trinomial_count1, '\0', 3 * sizeof(long int));
  memset(ma->trinomial_count2, '\0', 3 * sizeof(long int));
  memset(ma->transindels,      '\0', 7 * sizeof(long int));

  memset(f->nb_triplet_count1,            '\0', 64 * sizeof(long int));
  memset(f->nb_triplet_count2,            '\0', 64 * sizeof(long int));
  memset(f->nb_non_mutated_triplet_count, '\0', 64 * sizeof(long int));
  memset(f->nb_pair_of_triplets[0],       '\0', 64 * 64 * sizeof(long int));

  alignment_SG_stats_on_MA(dataquery, datachunktext, 0/*left_correction*/,
                           ma, f , 1 /* full or only straits alignments */);

  ma->entropy =  (float)
                 MIN(
                  ma->entropy,
                  MIN3(entropyTriplet(f->nb_non_mutated_triplet_count),
                       entropyTriplet(f->nb_triplet_count1),
                       entropyTriplet(f->nb_triplet_count2))
            );

  ma->mutual  =  (float)
                 mutualInformationTriplet(f->nb_pair_of_triplets);

  return ma->entropy < gp_entropy_min ;
}



/*
 * Entropy filter : remove low complexity regions
 */

long int EntropyAndScoreFilter_MAList(Feature * f, char *dataquery,
                                      char *datatext, long int *datatext_start)
{
    MA * ma = NULL, * ma_prev = NULL;

    ma = f->first_MA;

    /* (0) clean-up if nothing */
    if (!ma) {
      f->last_MA = NULL;
      return 0;
    }

    /* (1) find the first element OK and set it to *(p_firstMA) */
    while (EntropyAndScoreFilter_MA (f, ma, dataquery, datatext + datatext_start[ma->i_chunk])) {
      MA * to_free_ma = ma;
      ma = ma->next;
      FreeMA(to_free_ma);
      if (!ma) {
        f->first_MA = f->last_MA = NULL;
        return 0;
      }
    }
    f->first_MA = f->last_MA = ma;

    ma_prev = ma;
    ma = ma_prev->next;

    /* continue with the others MAs */
    while (ma) {
      if (EntropyAndScoreFilter_MA(f, ma, dataquery, datatext + datatext_start[ma->i_chunk])) {
          MA * to_free_ma = ma;
          ma = ma->next;
          ma_prev->next = ma;
          FreeMA(to_free_ma);
        } else {
          ma_prev = ma;
          ma = ma->next;
        }
    }
    ma_prev->next = NULL;
    f->last_MA = ma_prev;

    return 0;
}






/*
 * Display for debuging effort
 */

long int DisplayTuple(tuple * tuple)
{
  long int left_correction = 0;
  printf("\t tuple : [(TBL:%9ld,TEL:%9ld),(TBR:%9ld,TBL:%9ld)] \t diag : %9ld \t size : %9ld \n", TBL_POS(tuple), TEL_POS(tuple), TBR_POS(tuple), TER_POS(tuple), tuple->diagonal , TSIZE(tuple));
  return 0;
}


long int DisplayListTupleList(tuplelist * tuplelist)
{

  printf("**display tuple list**\n");
  while (tuplelist) {
    printf("\t--tl--\n");
    DisplayTupleList(tuplelist);
    tuplelist = tuplelist->next;
  }
  return 0;
}


long int DisplayTupleList(tuplelist * tuplelist)
{
  tuple *tuple;
  if (tuplelist->first_tuple) {
    tuple = tuplelist->first_tuple;
    while (tuple) {
      DisplayTuple(tuple);
      tuple = tuple->next;
    }
    printf("\n");
  }
  return 0;
}

long int DisplayGap(tuple * t1, tuple * t2)
{
  printf("\t intergaps {%ld,%ld}\n",TGAP_L(t1,t2),TGAP_R(t1,t2));
  return 0;
}

long int DisplayMA(MA * ma)
{

  tuple * tuple = ma->first_tuple;

  while (tuple) {
    printf("\t");DisplayTuple(tuple);
    if (tuple->next) {
      printf("\t");DisplayGap(tuple,tuple->next);
    }
    tuple = tuple->next;
  }
  printf("\t = {MA : (%ld,%ld), (%ld,%ld) chunk:%ld, score:%ld, reverse:%d}\n", ma->left_pos_begin, ma->left_pos_end, ma->right_pos_begin, ma->right_pos_end, ma->i_chunk, ma->blastscore, ma->reverse);
  printf("\n");
  return 0;
}

long int DisplayListMA (MA * ma) {
  while (ma) {
    DisplayMA(ma);
    ma = ma -> next;
  }
  return 0;
}

