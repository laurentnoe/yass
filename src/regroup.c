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
#include <math.h>

#ifdef DEBUG_REGROUP
#include <assert.h>
#endif

#include "global_var.h"
#include "tuple.h"
#include "regroup.h"


/*********************************
 *  RedBlack-AVL tree functions  *
 *********************************/


/*
 * TREE_TYPE(_node * smallest (TREE_TYPE(_tree) *T)
 *
 * Return the smallest element af a tree
 */

TREE_TYPE(_node) * smallest(TREE_TYPE(_tree) * T)
{
  TREE_TYPE(_node) * y = T->root;
  TREE_TYPE(_node) * s = NULL;

  while (y != NULL) {
    s = y;
    y = y->TREE_TYPE(_link[0]);

  }
  return s;
}



/*
 * void init_tree(TREE_TYPE(_tree) * tree)
 *
 * Reset tree data
 */

void init_tree(TREE_TYPE(_tree) * tree)
{
  tree->root = NULL;
  tree->k    = 0;
  memset(tree->pa,'\0',TREE_MAX_HEIGHT * sizeof(TREE_TYPE(_node *)));
  memset(tree->da,'\0',TREE_MAX_HEIGHT * sizeof(unsigned char));
  tree->to_delete = NULL;
}





/*
 * list_MA * tree_q_delete (tree_type *tree, MA *ma)
 *
 * Delete a node from the tree
 */

list_MA *tree_q_delete(TREE_TYPE(_tree) * tree, MA * ma)
{

  tree_data *data = NULL;
  list_MA *list_MA_curr = NULL;
#ifdef DEBUG_REGROUP
  assert(tree != NULL && ma != NULL);
#endif
  /* we search the element in the tree */
  data = TREE_TYPE(_search) (tree, DIAGB(ma));

  if (data) {
    /* when found we delete it from its queue in the tree */
    list_MA_curr = queue_extract_MA(&data->queue, ma);

    /* after that, if the queue is empty we delete the node */
    if (data->queue.first == NULL) {
      free_data(data);
      TREE_TYPE(_delete_rebalancing) (tree);
    }

    /* and free the  */

  } else {
    _WARNING("tree_q_delete() searched element not found : data empty");
  }
  return list_MA_curr;
}





/*
 * void display_tree_list (TREE_TYPE(_tree *T))
 *
 * Display the ordered list of nodes of a tree
 */

void display_tree_list(TREE_TYPE(_tree) * T)
{

  TREE_TYPE(_node) * y = smallest(T);
  tree_data *data = NULL;

  list_MA *l;
  long int i = 0;

  if (y) {
    data = &(y->data);
  } else {
    printf("+++ Empty List ! +++\n");
  }

  while (data) {
    l = data->queue.first;
    i = 0;
    while (l) {
      i++;
      l = l->next;
    }
    printf("(dist=%ld, elems=%ld)", data->distance, i);
    data = data->list[1];
  }
  printf("\n");

}


/*
 * void clear_tree_and_queue(Feature *feature)
 *
 * Delete all the "tree" and "queue" structures allocated and pointed by "feature"
 */

void clear_tree_and_queue(Feature *feature)
{

  list_MA * list_MA_l;
  list_MA * list_MA_pac;

  while (feature->T.root) {
    list_MA_l = (feature->T.root->data).queue.first;

    list_MA_pac = queue_extract_MA(&feature->Q, list_MA_l->ma);
    if (list_MA_pac) {
      FREE(list_MA_pac,sizeof(list_MA));
    } else {
      _WARNING("clear_tree() : missing prod_and_cust in Q ");
    }

    list_MA_pac = tree_q_delete(&feature->T, list_MA_l->ma);
    if (list_MA_pac) {
      if (list_MA_pac->customer) {
      FREE(list_MA_pac->customer,sizeof(table_MA));
      }
      if (list_MA_pac->producer) {
      FREE(list_MA_pac->producer,sizeof(table_MA));
      }
      FREE(list_MA_pac,sizeof(list_MA));
    } else {
      _WARNING("clear_tree() : missing prod_and_cust in T ");
    }
  }
}


/*
 * Compute the minimal "prmin" dyn table.
 */


#define DYNSIZE()                                                                                    { \
 if ( TUSGAP_L(t_current_MA_p,t_target_MA_c) >= 0 && TUSGAP_R(t_current_MA_p,t_target_MA_c) >= 0 ) {   \
     /* a) between "t_current_MA_p" and "t_target_MA_c" */                                             \
     if (TGAP_L(t_current_MA_p,t_target_MA_c) <= 0 || TGAP_R(t_current_MA_p,t_target_MA_c) <= 0 ) {    \
       prmin = 0;                                                                                      \
       t_current_best = t_current_MA_p;                                                                \
       t_target_best  = t_target_MA_c;                                                                 \
     } else {                                                                                          \
          /* b)  both gaps are > 0 */                                                                  \
           if (TGAP_L(t_current_MA_p,t_target_MA_c) *                                                  \
               TGAP_R(t_current_MA_p,t_target_MA_c) < prmin) {                                         \
               prmin = TGAP_L(t_current_MA_p,t_target_MA_c) *                                          \
                       TGAP_R(t_current_MA_p,t_target_MA_c);                                           \
             t_current_best = t_current_MA_p;                                                        \
             t_target_best  = t_target_MA_c;                                                         \
           }                                                                                           \
     }                                                                                                 \
  }                                                                                           \
}



/***************************************************************************/



/*
 * void window_deleter (long int pos);
 *
 * Delete all the MA which are no longer in the windows at position "pos"
 */

#ifdef INLINE
inline
#endif
void window_deleter(long int pos, Feature *feature)
{

  list_MA *itself;
  MA *current_MA;
  MA *target_MA;
  long int left_correction=0;

  tuple *t_current_MA,   *t_target_MA;
  tuple *t_current_MA_c, *t_current_MA_p, *t_target_MA_c;

  list_MA *list_MA_left;
  list_MA *list_MA_right;
  list_MA *list_MA_pac;

  list_MA *first_customer = NULL;
  list_MA *first_producer = NULL;

  /* data and data_sn are temporary variable */
  tree_data *data;
  tree_data *data_sn;

  list_MA *liste;

  long int i;
  long int dscore = 0;


  /* while elements are in the window */
  while ((feature->Q.first != NULL) && (feature->Q.first->ma->right_pos_end < pos)) {

    current_MA = feature->Q.first->ma;

    /* we delete the element from the tree and the global queue */
    list_MA_left = tree_q_delete(&feature->T, current_MA);
    FREE(queue_pop_MA(&feature->Q),sizeof(list_MA));

    /* we search for a customer of the deleted MA which best producer is the deleted MA */
    if (list_MA_left) {
      long int canbegrouped = 1;
    nextone:
      /* TEST */


      /* (0)  we check if they can be grouped (0) */
      if (list_MA_left->customer->index == 0) {
      goto cannotbegrouped;
      }

      /* (1) we check if they can be grouped (classify) */
      first_customer =   list_MA_left->customer->table[0];
      first_producer = first_customer->producer->table[0];
      if (first_producer != list_MA_left) {
      goto cannotbegrouped;
      }


      /* (2) we check if they can be grouped (score) */
      target_MA      = first_customer->ma;
      dscore = DIAG_COST(DIAG_LENGTH(current_MA, target_MA)) + INDEL_COST(INDEL_LENGTH(current_MA,target_MA));
#ifdef DEBUG_REGROUP
      if (target_MA->blastscore + dscore < 0 || current_MA->blastscore + dscore < 0 ) {
      _WARNING("window_deleter : scoring mismatch ");
      goto cannotbegrouped;
      }
#endif


      /* (3) we check if they can be grouped (tuple) and also link their own tuplelists between if possible */
      t_current_MA        = current_MA->first_tuple;
      t_target_MA         = target_MA->first_tuple;
#ifdef DEBUG_REGROUP
      assert(t_current_MA);
      assert(t_target_MA);
#endif

      /* (3.1) avoid overlaps : cross the first tuple list until first second tuple list is reached */
      t_current_MA_c      = t_current_MA;
      t_current_MA_p      = NULL;
      while ( (t_current_MA_c) &&
            TEL_POS(t_current_MA_c) <= TL_POS(t_target_MA) &&
            TER_POS(t_current_MA_c) <= TR_POS(t_target_MA)) {
      t_current_MA_p = t_current_MA_c;
      t_current_MA_c = t_current_MA_c->next;
      }

      /* (3.2) avoid overlaps : cross both tuple lists */
      if (t_current_MA_p) {
      long int prmin = MAX_WINDOW_REGROUP*4;
      tuple *t_current_best = NULL, *t_target_best = NULL;
      t_target_MA_c      = t_target_MA;
      /* a) Memorisez la dynsize */
      DYNSIZE();
      /* b) parcours conjoint */
      while (t_current_MA_p) {
        while ((t_target_MA_c) &&  !(TL_POS(t_target_MA_c) >= TEL_POS(t_current_MA_p) &&  TR_POS( t_target_MA_c) >= TER_POS(t_current_MA_p))) {
          t_target_MA_c    =  t_target_MA_c->next;
        }
        if (t_target_MA_c) {
          /* a) Memorisez la dynsize */
          DYNSIZE();
        }/* if (t_target_MA_c) */
        t_current_MA_p = t_current_MA_p->next;
      } /* while (t_current_MA_p) */

      /* Adjust and Chain */
      if (t_current_best && t_target_best) {
        while (TEL_POS(t_current_best) > TBL_POS(t_target_best) || TER_POS(t_current_best) > TBR_POS(t_target_best)) {
          t_target_best->leftsize--;
#ifdef DEBUG_REGROUP
          assert( t_target_best->leftsize >= 0 );
#endif
        }
        t_current_best->next =  t_target_best;
        /* FIXME : tuple cleaning */
        /*
        printf("s1=%ld,s2=%ld %ld->%ld(%ld), %ld->%ld(%ld)\n", current_MA->blastscore, target_MA->blastscore,
             TEL_POS(t_current_best),TBL_POS(t_target_best),TBL_POS(t_target_best)-TEL_POS(t_current_best),
             TER_POS(t_current_best),TBR_POS(t_target_best),TBR_POS(t_target_best)-TER_POS(t_current_best));
        fflush(NULL);
        */
      } else {
#ifdef DEBUG_REGROUP
        _WARNING("window_deleter : unable to join tuples du to tuple constraints : jump a");
        fflush(NULL);
        DisplayMA(current_MA);
        DisplayMA(target_MA);
#endif
        table_MA_delete(first_producer->customer, first_customer);
        table_MA_delete(first_customer->producer, list_MA_left);
        goto nextone;
      }

      } else {
      /* (3.2) avoid overlaps : cross the second tuple list because there is only one tuple in the first list */
      t_target_MA_c      = t_target_MA;
      while ((t_target_MA_c) &&  (TL_POS(t_target_MA_c) < TEL_POS(t_current_MA) ||  TR_POS( t_target_MA_c) < TER_POS(t_current_MA))) {
        t_target_MA_c    =  t_target_MA_c->next;
      }

      if (t_target_MA_c) {
#ifdef DEBUG_REGROUP
        assert(TEL_POS(t_current_MA) <  TL_POS(t_target_MA_c) && TER_POS(t_current_MA) <  TR_POS(t_target_MA_c));
#endif
        while (TEL_POS(t_current_MA) > TBL_POS(t_target_MA_c) || TER_POS(t_current_MA) > TBR_POS(t_target_MA_c)) {
          t_target_MA_c->leftsize--;
#ifdef DEBUG_REGROUP
          assert( t_target_MA_c->leftsize >= 0 );
#endif
        }
        /* FIXME : tuple cleaning */
        t_current_MA->next = t_target_MA_c;
      } else {
        /* (3.3) nothing can be grouped : but try the next one*/
#ifdef DEBUG_REGROUP
        _WARNING("window_deleter : unable to join tuples du to tuple constraints : jump b");
        DisplayMA(current_MA);
        DisplayMA(target_MA);
#endif
        table_MA_delete(first_producer->customer, first_customer);
        table_MA_delete(first_customer->producer, list_MA_left);
        goto nextone;
      }
      }



      /* if MAs can be grouped together (note this is always true here )*/
      if (canbegrouped) {

      STATS_NB_POSTPROCESSED_GROUPING_LINKS_INC(feature);

      /*we delete the grouped MA from the tree */
      list_MA_right = tree_q_delete(&feature->T, target_MA);

#ifdef DEBUG_REGROUP
      assert(list_MA_right->ma == target_MA);

      if (list_MA_right) {
#endif

        /* we delete the grouped MA from the global queue */
        list_MA_pac = queue_extract_MA(&feature->Q, target_MA);
#ifdef DEBUG_REGROUP
        if (list_MA_pac) {
#endif
          if (list_MA_pac->customer) {
             _WARNING("window_deleter : list_MA_pac->customer is not empty");
          }
          if (list_MA_pac->producer) {
             _WARNING("window_deleter : list_MA_pac->customer is not empty");
          }

          FREE(list_MA_pac,sizeof(list_MA));


#ifdef DEBUG_REGROUP
        } else {
          _WARNING("window_deleter :  missing prod_and_cust in target_MA");
        }
#endif



        /* update alignment positions and scores */

        target_MA->first_tuple    = NULL;
        current_MA->blastscore   += target_MA->blastscore + dscore;
        current_MA->left_pos_end  = target_MA->left_pos_end;
        current_MA->right_pos_end = target_MA->right_pos_end;


        /* we signal to other customer from the left MA that it won't be no more a producer */
        for (i = 1; i < list_MA_left->customer->index; i++) {
          list_MA *i_customer =
            list_MA_left->customer->table[i];
          table_MA_delete(i_customer->producer,
                      list_MA_left);
        }

        /* we signal to other producer from the right MA that it won't be no more a customer */
        for (i = 1; i < list_MA_right->producer->index; i++) {
          list_MA *i_producer =
            list_MA_right->producer->table[i];
          table_MA_delete(i_producer->customer,
                      list_MA_right);
        }

        /* and free the unused structures */
#ifdef DEBUG_REGROUP
        memset(list_MA_left->customer,  0, sizeof(table_MA));
        memset(list_MA_right->producer, 0, sizeof(table_MA));
#endif
        FREE(list_MA_left->customer,  sizeof(table_MA));
        FREE(list_MA_right->producer, sizeof(table_MA));


        /* then we push it in order for it to be the last of the list */
        itself = queue_push_and_sort_MA(&feature->Q, current_MA, NULL);

        /* and then we insert it into the tree */
        data_sn = TREE_TYPE(_find) (&feature->T, DIAGB(current_MA));
        data    = TREE_TYPE(_insert) (&feature->T, DIAGB(current_MA));

        if (data->queue.first != NULL) {
          liste = queue_push_and_sort_MA(&data->queue, current_MA, itself);
        } else {
          liste = new_data(data, data_sn, itself);
        }

        /* we give back to liste the producer of the left MA and the customer of the right MA */
        liste->customer = list_MA_right->customer;
        liste->producer = list_MA_left->producer;

        /* we signal to the customer of list that its adress has changed */
        for (i = 0; i < liste->customer->index; i++) {
          list_MA *i_customer = liste->customer->table[i];
          table_MA_delete(i_customer->producer, list_MA_right);
          table_MA_insert(i_customer->producer, liste        ,current_MA->blastscore, i_customer);
        }

        /* we try to find new producers for list */
        regroup(liste,feature);

#ifdef DEBUG_REGROUP
        memset(list_MA_left, 0, sizeof(list_MA));
        memset(list_MA_right, 0, sizeof(list_MA));
#endif
        FREE(list_MA_left,  sizeof(list_MA));
        FREE(list_MA_right, sizeof(list_MA));

#ifdef DEBUG_REGROUP

      } else {
        _WARNING("window_deleter : missing list_MA_right in T ");
        goto cannotbegrouped;
      }


#endif

      } else {
      cannotbegrouped:
      /* If the MA can't be regrouped we signal to its customer that it can't be no more a producer */
      for (i = 0; i < list_MA_left->customer->index; i++) {
        list_MA * i_customer = list_MA_left->customer->table[i];
        table_MA_delete(i_customer->producer, list_MA_left);
      }

#ifdef DEBUG_REGROUP
      memset(list_MA_left->customer, 0, sizeof(table_MA));
      memset(list_MA_left->producer, 0, sizeof(table_MA));
      memset(list_MA_left, 0, sizeof(list_MA));
#endif
      /* and free associated structures */
      FREE(list_MA_left->producer, sizeof(table_MA));
      FREE(list_MA_left->customer, sizeof(table_MA));
      FREE(list_MA_left, sizeof(list_MA));

      }
    } else {
      _WARNING("window_deleter : missing list_MA_left in T ");
    }
  }
}


/*
 * void regroup(list_MA * liste,  Feature * feature);
 *
 * Test with blascore and distance criterions if a MA can be linked with an
 * other MA and proceed if possible
 */

#ifdef INLINE
inline
#endif
void regroup(list_MA * liste,  Feature * feature)
  {

  /* y is the current node */
  tree_data *y = NULL;
  /* y_plus is y->right */
  tree_data *y_plus = NULL;
  /* best is the nearest node from liste->ma between y and y_plus */
  tree_data *tree_data_best = NULL;

  /* d is the value of the distance to diagonale of the begin of liste->ma */
  long int d;

  /* min and max are the bound between which we can make a regroupement */
  long int min;
  long int max;

  /* if hit == 1 then we make a regroupement */
  long int hit = 0;

  /* the list_MA which is regrouped with liste */
  list_MA *l = NULL;

  /* First we choose which data will be the most efficient to start our comparisons */
  y = TREE_TYPE(_find) (&feature->T, DIAGB(liste->ma));
  d = DIAGE(liste->ma);

  min = d - INDEL_BOUND(liste->ma);
  max = d + INDEL_BOUND(liste->ma);

  y_plus          = y->list[1];
  tree_data_best  = y;

  /* At this point we start the search of an MA which can be regroup with liste->ma */

  /* While potential candidates do exist */

  while (((y != NULL) || (y_plus != NULL)) && (!hit)) {

    l = tree_data_best->queue.last;

    /* For each element in the node queue which are valuable for linking */

    while ( (l)  && (!hit)) {
      /*
       * get first tuple from liste->ma
       * and last tuple from l->ma
       */
      tuple * liste_ma_tuple_begin = liste->ma->first_tuple;
      tuple * l_ma_tuple_end       = l->ma->first_tuple;
      while (l_ma_tuple_end->next)
        l_ma_tuple_end = l_ma_tuple_end->next;

#ifdef DEBUG_REGROUP
            printf("liste->ma : (%ld,%ld),(%ld,%ld)  l->ma :(%ld,%ld),(%ld,%ld)\n",
               liste->ma->left_pos_begin, liste->ma->right_pos_begin,
               liste->ma->left_pos_end,   liste->ma->right_pos_end,
               l->ma->left_pos_begin,     l->ma->right_pos_begin,
               l->ma->left_pos_end,       l->ma->right_pos_end );
#endif


      /*********************************************************************************************
       *  ATTENTION CECI EST PROVISOIRE ET N'EST UTILE QUE POUR EVITER LES BUGS DE L'OPTION -d 1   *
       *  On ne regroupe pas les MA qui aurait une trop grosse ouverture de GAP                    *
       *   la memoire risquerait d'exploser quand on recherche le meilleur chemin en les 2 MAs     *
       *   au moment de l'alignement                                                               *
       *********************************************************************************************/
       if (l->ma == liste->ma) {
       hit = 0;
       } else if (
              ((ABS(liste->ma->right_pos_begin - l->ma->right_pos_end) + 1) > MAX_WINDOW_REGROUP / (1+2*gp_delta_stat))
              ||
              ((ABS(liste->ma->left_pos_begin  - l->ma->left_pos_end ) + 1) > MAX_WINDOW_REGROUP / (1+2*gp_delta_stat))
              ||
              (
               (
                MIN(
                  TGAP_L(l_ma_tuple_end,liste_ma_tuple_begin),
                  TGAP_R(l_ma_tuple_end,liste_ma_tuple_begin)
                  )
               )
               *
               (
                (
                 ABS(TGAP_R(l_ma_tuple_end,liste_ma_tuple_begin) - TGAP_L(l_ma_tuple_end,liste_ma_tuple_begin))
                 )
                + 2*gp_delta_stat + 1
               )
               )
              >  MAX_WINDOW_REGROUP
              )
       {
         hit = 0;
#ifdef DEBUG_REGROUP
         printf(" * HIT AVOIDED : left gap = %ld , right gap = %ld , product = %ld\n",
              /*
                (ABS(liste->ma->right_pos_begin - l->ma->right_pos_end) + 1),
                (ABS(liste->ma->left_pos_begin  - l->ma->left_pos_end ) + 1),
                (ABS(liste->ma->right_pos_begin - l->ma->right_pos_end) + 1) *
                (ABS(liste->ma->left_pos_begin  - l->ma->left_pos_end ) + 1)
              */
              (ABS(TGAP_L(l_ma_tuple_end,liste_ma_tuple_begin)) + 1),
              (ABS(TGAP_R(l_ma_tuple_end,liste_ma_tuple_begin)) + 1),
              (ABS(TGAP_L(l_ma_tuple_end,liste_ma_tuple_begin)) + 1) *
              (ABS(TGAP_R(l_ma_tuple_end,liste_ma_tuple_begin)) + 1)
              );
#endif

       } else {
         hit = GOOD_LINK(l,liste);

#ifdef DEBUG_REGROUP
         if (hit) {
           printf(" * HIT\n");
           printf(" -> gap length %ld,%ld\n",liste->ma->right_pos_begin - l->ma->right_pos_end+1,liste->ma->left_pos_begin - l->ma->left_pos_end+1);
           printf(" -> scores : ma_before : %ld, ma_after : %ld, between : %ld\n",
                liste->ma->blastscore,
                l->ma->blastscore,
                DIAG_COST(DIAG_LENGTH(l->ma, liste->ma)) + INDEL_COST(INDEL_LENGTH(l->ma,liste->ma)));
           DisplayMA(l->ma);DisplayMA(liste->ma);
         }
#endif
       }


      /*********************************************************************************************/



      /* if the two MA can be regrouped */
      if (hit) {
      long int score = SCORE_ALIGN(l,liste);

      /* we insert the new MA in the customer table of the old MA */
      if (table_MA_insert(l->customer, liste, score, l)) {
        /* if this insertion has succeeded we insert the old MA in the producer table of the new MA */
        if (!table_MA_insert(liste->producer, l, score, liste)) {
          /* if this insertin has failed we delete the new MA from the customer table of the old MA */
          table_MA_delete(l->customer, liste);
        }
      }
      hit = 0;
      }

      STATS_NB_POSTPROCESSED_GROUPING_TESTS_INC(feature);

      l = l->prev;

    } /* while (l) ... */



    if (!hit) {

      /* we change y or y_plus depending on who was the best */
      if (tree_data_best == y) {
        y = y->list[0];
      } else {
        y_plus = y_plus->list[1];
      }

      /* we test if we are under the bounds */
      if ((y != NULL) && (y->distance < min)) {
      y = NULL;
      }

      if ((y_plus != NULL) && (y_plus->distance > max)) {
      y_plus = NULL;
      }


      /* then we decide who is the best between y and y_plus (if we go on the left or on the right of the list) */
      if (y != NULL) {
        tree_data_best = y;
        if ((y_plus != NULL)
            && (y_plus->distance - d) < (d - y->distance)) {
          tree_data_best = y_plus;
        }
      } else {
        if (y_plus != NULL) {
          tree_data_best = y_plus;
        }
      }
    }/* if (hit) */


  }
}




/*
 * MA * main_regroup (MA * first_MA, Feature *feature);
 *
 */

#ifdef INLINE
inline
#endif
MA * main_regroup(MA * first_MA, Feature *feature) {

  /* the current MA */
  MA * current_MA = first_MA, * prev_MA = NULL; /* FIXME : */

  long int i_chunk = first_MA->i_chunk;


  tree_data *data;
  tree_data *data_sn;            /*sn = smallest nearest */

  list_MA *itself;
  list_MA *l;

  /* while we can make new regroupement */

  current_MA = first_MA;

  /* init structures */
  init_tree(&feature->T);

  feature->Q.first = NULL;
  feature->Q.last  = NULL;

  feature->win_position = 0;

  /* for each MA which has the same i_chunk */

  while (current_MA && current_MA->i_chunk == i_chunk) {

    /* we insert it inside the global queue */

    itself = queue_push_and_sort_MA(&feature->Q, current_MA, NULL);

    /* we insert it inside the tree */

    data_sn = TREE_TYPE(_find)    (&feature->T, DIAGB(current_MA));
    data    = TREE_TYPE(_insert)  (&feature->T, DIAGB(current_MA));

    /* we fill the tree_data field */

    if (data->queue.first != NULL) {

      l = queue_push_and_sort_MA(&data->queue, current_MA, itself);
      l->customer = (table_MA *)MALLOC(sizeof(table_MA));
      ASSERT(l->customer,"main_regroup");

      l->producer = (table_MA *)MALLOC(sizeof(table_MA));
      ASSERT(l->producer,"main_regroup");

      table_MA_init(l->customer);
      table_MA_init(l->producer);

    } else {

      l = new_data(data, data_sn, itself);
    }

    /* we try to regroup the MA */
    regroup(l,feature);


    /* if the new MA is out of the window we move it */
    if (current_MA->right_pos_begin - feature->windows > feature->win_position) {
      feature->win_position = current_MA->right_pos_begin - feature->windows;
      /* then we delete the elements which are now out of the window */
      window_deleter(feature->win_position,feature);
    }

    current_MA = current_MA->next;
    STATS_NB_POSTPROCESSED_MA_INC(feature);
  } /* while (ma && i_chunk) */

    /* end of the process : we push the window on the right border */
  feature->win_position = gp_textsize;
  window_deleter(feature->win_position,feature);

  clear_tree_and_queue(feature);

  /* At this point we delete all the MA which have been grouped : their tuplelist is empty */

  current_MA = first_MA;

  while (current_MA && current_MA->i_chunk == i_chunk) {
    if (current_MA->first_tuple == NULL) {
      MA * to_free_MA = current_MA;
      current_MA      = current_MA->next;
      if (prev_MA) {
        prev_MA->next = current_MA;
      } else {
        _WARNING("main_regroup() : prev_MA undefined"); /* FIXME : this test was not here before */
        feature->last_MA  =  current_MA; /* FIXME : dangerous on the first tuple */
      }
      FREE(to_free_MA,sizeof(MA));
    } else {
      feature->last_MA    =  prev_MA  = current_MA;
      current_MA    = current_MA->next;
    }
  }

  /* return first_MA of the next i_chunk */
  return current_MA;
}




void Regroup_MAList(Feature * feature) {

  MA * ma =  feature->first_MA, *chunkma = feature->first_MA;

  if (gp_win_mul > 1.0) {
    /* foreach chunk */
    while (ma) {
      feature->windows = gp_win_min;
      /* foreach window size */
      while (feature->windows <= gp_win_max) {
        chunkma = main_regroup(ma, feature);
        {
          long int new_windows = (long int)(feature->windows * gp_win_mul);
          if (feature->windows == new_windows)
            feature->windows++;
          else
            feature->windows = new_windows;
        }
      }
      ma = chunkma;
    }
  }
  feature->last_MA->next = NULL;
}

