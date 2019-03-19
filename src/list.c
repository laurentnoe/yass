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
#include <math.h>
#include <errno.h>

#include "global_var.h"
#include "tuple.h"
#include "red_black.h"


/*
 * void display_queue (queue_MA *q)
 *
 * Display a queue
 */

void display_queue(queue_MA * q) {
  list_MA *l = q->first;
  while (l) {
    printf("|(%ld,%ld,adr=%p,blast=%ld)",
           DIAGB(l->ma),
           l->ma->right_pos_end,
           (void *) l,
           l->ma->blastscore);
    l = l->next;
  }
  printf("|X\n");
}

void display_table(table_MA * t) {
  long int i;
  printf("------------------\n");
  for (i = 0; i < CUSTOMER_SIZE; i++) {
    if (t->table[i]) {
      printf("| %p|%ld\n", (void *) t->table[i], t->table_blastscore[i]);
    } else {
      printf("|  NULL  |    0\n");
    }
  }
  printf("------------------\n");
}

/*
 * void new_data (tree_data *data, tree_data *data_sn, list_MA *itself)
 *
 * Fill a tree_data field (and link the list field)
 */

#ifdef INLINE
inline
#endif
list_MA * new_data(tree_data * data, tree_data * data_sn,
                   list_MA * itself) {

  tree_data * next;
  list_MA   * first;

  if (!data) {
    return NULL;
  }

  /* create a new list_MA */
  first = (list_MA *)MALLOC(sizeof(list_MA));
  ASSERT(first, "new_data()");
  memset(first, 0, sizeof(list_MA));
  first->ma     = itself->ma;
  first->itself = itself;
  first->next = NULL;
  first->prev = NULL;

  /* and a new tree_data  */
  data->distance = DIAGB(itself->ma);
  data->list[0] = NULL;
  data->list[1] = NULL;
  data->queue.first = first;
  data->queue.last  = first;

  /* After that we link the lists field of the tree data structure */
  if (data_sn) {  /* "sn" means "smallest nearest" */

    /* if a smaller element than data exists */
    next = data_sn->list[1];

    if (data_sn->distance > data->distance) {
      data->list[0] = NULL;
      data->list[1] = data_sn;
      data_sn->list[0] = data;
    } else {

      /* if data is the smallest element */
      data->list[0] = data_sn;
      data_sn->list[1] = data;
      data->list[1] = next;
      if (next != NULL) {
        next->list[0] = data;
      }
    }
  }
  return first;
}


/*
 * void free_data(tree_data *data)
 *
 * Unlink the list fields of a tree_data and free the queue
 */

#ifdef INLINE
inline
#endif
void free_data(tree_data * data) {
  /* jump the linked list over this "tree_data" */
  tree_data *prev = data->list[0];
  tree_data *next = data->list[1];
  if (prev != NULL) {
    prev->list[1] = next;
  }
  if (next != NULL) {
    next->list[0] = prev;
  }
  /* erase any link of this "tree_data", but dont "free" it */
  data->list[0] = NULL;
  data->list[1] = NULL;
}


/*
 * list_MA * queue_push_and_sort_MA (queue_MA *q, MA *item, list_MA *itself)
 *
 * Put a MA on the top of a queue and keep the queue sorted
 */

#ifdef INLINE
inline
#endif
list_MA * queue_push_and_sort_MA(queue_MA * q,
                                 MA * item,
                                 list_MA * itself) {
  list_MA *next_list_MA = NULL;
  list_MA *curr_list_MA = (q->last);

  list_MA *new_list_MA  = (list_MA *) MALLOC(sizeof(list_MA));
  ASSERT(new_list_MA, "queue_push_and_sort_MA");


  new_list_MA->ma     = item;
  new_list_MA->next   = NULL;
  new_list_MA->prev   = NULL;
  new_list_MA->customer = NULL; /* unused in the queue */
  new_list_MA->producer = NULL; /* unused in the queue */
  new_list_MA->itself = itself;


  /* Insert "new_list_MA" at the right position in the queue */
  while (curr_list_MA) {
    if (curr_list_MA->ma->right_pos_end <= item->right_pos_end)
      break;
    curr_list_MA = curr_list_MA->prev;
  }

  if (curr_list_MA) {
    next_list_MA = curr_list_MA->next;
    curr_list_MA->next = new_list_MA;
  } else {
    next_list_MA = q->first;
    q->first = new_list_MA;
  }

  new_list_MA->prev = curr_list_MA;
  new_list_MA->next = next_list_MA;

  if (next_list_MA) {
    next_list_MA->prev = new_list_MA;
  } else {
    q->last = new_list_MA;
  }

  /* return the newly created element */
  return new_list_MA;
}




#ifdef INLINE
inline
#endif
list_MA * queue_insert_and_sort_MA(queue_MA * q,
                                   MA * item,
                                   list_MA * itself) {
  list_MA *new_list_MA  = (list_MA *) MALLOC(sizeof(list_MA));
  ASSERT(new_list_MA, "queue_insert_and_sort_MA");

  new_list_MA->ma     = item;
  new_list_MA->next   = NULL;
  new_list_MA->prev   = NULL;
  new_list_MA->customer = NULL; /* unused in the queue */
  new_list_MA->producer = NULL; /* unused in the queue */
  new_list_MA->itself = itself;

  /* Insert "new_list_MA" easily when empty */
  if (!(q->first) && !(q->last)) {
    q->first = q->last = new_list_MA;
    return new_list_MA;
  } else {
    /* Insert "new_list_MA" at the correct position in the queue, by checking Both directions here */
    list_MA *curr_list_MA_from_begin = (q->first);
    list_MA *curr_list_MA_from_end   = (q->last);
    while (curr_list_MA_from_end && curr_list_MA_from_begin) {
      /* (a) from the end */
      if ((curr_list_MA_from_end->ma->right_pos_end < item->right_pos_end) || ((curr_list_MA_from_end->ma->right_pos_end == item->right_pos_end) && (curr_list_MA_from_end->ma->left_pos_end < item->left_pos_end))) {
        list_MA * tmp_list_MA = curr_list_MA_from_end->next;
        curr_list_MA_from_end->next = new_list_MA;
        new_list_MA->prev = curr_list_MA_from_end;
        if (tmp_list_MA) {
          new_list_MA->next = tmp_list_MA;
          tmp_list_MA->prev = new_list_MA;
        } else
          q->last = new_list_MA;
        return new_list_MA;
      }
      /* (b) from the beginning */
      if ((curr_list_MA_from_begin->ma->right_pos_end > item->right_pos_end) || ((curr_list_MA_from_begin->ma->right_pos_end == item->right_pos_end) && (curr_list_MA_from_begin->ma->left_pos_end >= item->left_pos_end))) {
        list_MA * tmp_list_MA = curr_list_MA_from_begin->prev;
        curr_list_MA_from_begin->prev = new_list_MA;
        new_list_MA->next = curr_list_MA_from_begin;
        if (tmp_list_MA) {
          new_list_MA->prev = tmp_list_MA;
          tmp_list_MA->next = new_list_MA;
        } else
          q->first = new_list_MA;
        return new_list_MA;
      }
      curr_list_MA_from_begin = curr_list_MA_from_begin->next;
      curr_list_MA_from_end   = curr_list_MA_from_end->prev;
    }
  }

  /* if the list is correctly set, would never reach this step */
  _WARNING("queue_insert_and_sort_MA() : badly formed queue");
  return NULL;
}


/*
 * list_MA * queue_pop_MA (queue_MA *q)
 *
 * Delete the first MA of a queue
 */

#ifdef INLINE
inline
#endif
list_MA * queue_pop_MA(queue_MA * q) {
  list_MA *curr_list_MA = q->first;
  if (curr_list_MA) {
    if (q->first == q->last) {
      q->first = NULL;
      q->last  = NULL;
    } else {
      q->first = curr_list_MA->next;
      if (q->first)
        q->first->prev = NULL;
    }
    /* separate the list_MA from any link ... */
    curr_list_MA->next = NULL;
    curr_list_MA->prev = NULL;
  }
  return curr_list_MA;
}


/*
 * list_MA * queue_extract_MA (queue_MA *q, list_MA *l, list_MA ** p_list_MA_hint)
 *
 * extract a MA of a queue, by searching both ends or by following a hint
 */

#ifdef INLINE
inline
#endif
list_MA * queue_extract_MA(queue_MA * q, MA * item, list_MA ** p_list_MA_hint) {

  list_MA *curr_list_MA = NULL;

  /* (a) search for the MA "item" in the queue "q" */
  if (p_list_MA_hint && (*p_list_MA_hint)) {
    /* (a.1) use a possible "hint" to follow close elements in the queue "q" */
    list_MA *curr_list_MA_from_hint_prev = *p_list_MA_hint;
    list_MA *curr_list_MA_from_hint_next = *p_list_MA_hint;
    while (curr_list_MA_from_hint_next || curr_list_MA_from_hint_prev) {
      if (curr_list_MA_from_hint_next) {
        if (curr_list_MA_from_hint_next->ma != item)
          curr_list_MA_from_hint_next = curr_list_MA_from_hint_next->next;
        else {
          curr_list_MA = curr_list_MA_from_hint_next;
          break;
        }
      }
      if (curr_list_MA_from_hint_prev) {
        if (curr_list_MA_from_hint_prev->ma != item)
          curr_list_MA_from_hint_prev = curr_list_MA_from_hint_prev->prev;
        else {
          curr_list_MA = curr_list_MA_from_hint_prev;
          break;
        }
      }
    }
    if (!curr_list_MA) {
      _WARNING("queue_extract_MA() : missing MA from hint in the queue, or invalid queue");
      return NULL;
    }
  } else {
    /* (a.2) no "hint", so follow queue "q" both ends */
    list_MA *curr_list_MA_from_begin = q->first;
    list_MA *curr_list_MA_from_end   = q->last;

    if (!curr_list_MA_from_begin || !curr_list_MA_from_end) {
      _WARNING("queue_extract_MA() : empty queue, missing MA in the queue");
      return NULL;
    }

    while (curr_list_MA_from_begin->ma != item && curr_list_MA_from_end->ma != item) {
      curr_list_MA_from_end   = curr_list_MA_from_end->prev;
      curr_list_MA_from_begin = curr_list_MA_from_begin->next;
      if (!curr_list_MA_from_end || !curr_list_MA_from_begin) {
        _WARNING("queue_extract_MA() : missing MA in the queue or invalid queue");
        return NULL;
      }
    }

    if (curr_list_MA_from_begin && curr_list_MA_from_begin->ma == item)
      curr_list_MA = curr_list_MA_from_begin;
    else
      if (curr_list_MA_from_end && curr_list_MA_from_end->ma == item)
        curr_list_MA = curr_list_MA_from_end;
      else {
        _WARNING("queue_extract_MA() : missing MA in the queue with invalid queue");
        return NULL;
      }
  }


  /* (b) remove the "item" from the queue "q" */
  list_MA *prev_list_MA  = curr_list_MA->prev;
  list_MA *next_list_MA  = curr_list_MA->next;
  curr_list_MA->prev = NULL;
  curr_list_MA->next = NULL;

  if (prev_list_MA)
    prev_list_MA->next = next_list_MA;
  else {
    q->first = next_list_MA;
    if (next_list_MA)
      next_list_MA->prev = NULL;
  }
  if (next_list_MA)
    next_list_MA->prev = prev_list_MA;
  else {
    q->last = prev_list_MA;
    if (prev_list_MA)
      prev_list_MA->next = NULL;
  }

  /* (c) put the "prev" or "next" to the "hint" if any variable was set to */
  if (p_list_MA_hint) {
    if (prev_list_MA)
      *p_list_MA_hint = prev_list_MA;
    else
      if (next_list_MA)
        *p_list_MA_hint = next_list_MA;
      else
	*p_list_MA_hint = NULL;
  }
  return curr_list_MA;
}


#ifdef INLINE
inline
#endif
void table_MA_init(table_MA * T)
{
  long int i;
  T->size = 0;
  for (i = 0; i < CUSTOMER_SIZE; i++) {
    T->table[i]            = NULL;
    T->table_blastscore[i] = 0;
  }
}


#ifdef INLINE
inline
#endif
long int table_MA_insert(table_MA * T, list_MA * insert, long int dscore, list_MA * l) {

  long int curr = 0;            /* used to sort the table */

  /* if the table is already full */
  if (T->size >= CUSTOMER_SIZE) {
    /* if our MA has to be inserted */
    if (dscore > T->table_blastscore[CUSTOMER_SIZE - 1]) {
      /* we delete the MA with the lowest score */
      list_MA *old_list = T->table[CUSTOMER_SIZE - 1];
      T->table           [CUSTOMER_SIZE - 1] = insert;
      T->table_blastscore[CUSTOMER_SIZE - 1] = dscore;
      curr = CUSTOMER_SIZE - 1;
      /* we signal to other MA the deletion */
      if (l->producer == T)
        table_MA_delete(old_list->customer, l);
      if (l->customer == T)
        table_MA_delete(old_list->producer, l);
    } else {
      return 0;
    }
  } else {
    /* if the table is not full */
    T->table           [T->size] = insert;
    T->table_blastscore[T->size] = dscore;
    curr = T->size;
    T->size++;
  }


  /* table sort */
  while (curr > 0 && T->table_blastscore[curr] > T->table_blastscore[curr - 1]) {
    list_MA *old_list   = T->table[curr];
    long int old_dscore = T->table_blastscore[curr];

    T->table[curr]            = T->table[curr - 1];
    T->table_blastscore[curr] = T->table_blastscore[curr - 1];

    T->table[curr - 1]            = old_list;
    T->table_blastscore[curr - 1] = old_dscore;

    curr--;
  }

  return 1;
}



#ifdef INLINE
inline
#endif
long int table_MA_delete(table_MA * T, list_MA * suppr) {
  long int i = 0;
  while (i < T->size) {
    if (T->table[i] == suppr)
      break;
    i++;
  }

  if (i == T->size) {
    _WARNING("table_MA_delete() element not found");
    return 0;
  }

  T->size--;

  for (; i < T->size; i++) {
    T->table           [i] = T->table[i + 1];
    T->table_blastscore[i] = T->table_blastscore[i + 1] ;
  }

  T->table           [T->size] = NULL;
  T->table_blastscore[T->size] = 0;
  return 1;
}
