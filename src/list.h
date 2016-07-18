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

#ifndef  __LIST_H__
#define  __LIST_H__

#define LPB(x)            ((x)->left_pos_begin)
#define RPB(x)            ((x)->right_pos_begin)
#define LPE(x)            ((x)->left_pos_end)
#define RPE(x)            ((x)->right_pos_end)
#define DIAGE(x)          (LPB(x)-RPB(x))
#define DIAGB(x)          (LPE(x)-RPE(x))

#define COMPARE(a,b)  ((a)<(b)?(0):(1))

#define CUSTOMER_SIZE 4


typedef struct _table_MA {
    long int index;
    long int table_blastscore[CUSTOMER_SIZE];
    struct _list_MA *table[CUSTOMER_SIZE];
} table_MA;


/*
* List of MA structure, itself represents a pointer on the same MA
* in another list
*/

typedef struct _list_MA {

    MA *ma;
    struct _list_MA *next;
    struct _list_MA *prev;
    struct _list_MA *itself; /* the same but inside the queue from the tree */
    struct _table_MA *customer;
    struct _table_MA *producer;

} list_MA;


/*Queue of MA*/

typedef struct _queue_MA {

    struct _list_MA *first;
    struct _list_MA *last;

} queue_MA;


/*
* The data structure inside the tree, it contains the distance field which
* represents the distance to the MA diagonale, a queue and two links which
* point on the smaller and the bigger MA
*/

typedef struct _tree_data {
    long int distance;
    struct _queue_MA   queue;
    struct _tree_data *list[2];      /* list[0] : pointer to smaller neigboor node,
                                      * list[1] : pointer to bigger neigboor node
                                      */
} tree_data;


/*
 * void display_queue (queue_MA *q)
 *
 * Display a queue
 */

void display_queue(queue_MA * q);


/*
* list_MA * void new_data (tree_data *data, tree_data *data_sn, list_MA *itself)
*
* Fill a tree_data field (malloc a queue and link the list fiels
*/

list_MA *new_data(tree_data * data, tree_data * data_sn, list_MA * itself);


/*
* void free_data(tree_data *data)
*
* Unlink the list fields of a tree_data and free the queue
*/

void free_data(tree_data * data);


/*
* list_MA * queue_push_and_sort_MA (queue_MA *q, MA *item, list_MA *itself)
*
* Put a MA on the top of a queue and keep the queue sorted
*/

list_MA *queue_push_and_sort_MA(queue_MA * q, MA * item, list_MA * itself);


/*
* list_MA * queue_pop_MA (queue_MA *q)
*
* Delete the first MA of a queue
*/

list_MA *queue_pop_MA(queue_MA * q);

/*
* list_MA * queue_extract_MA (queue_MA *q, list_MA *l)
*
* Delete a list_MA from a queue
*/

list_MA * queue_extract_MA(queue_MA * q, MA * item);


/*
* void table_MA_init (table_MA *T)
*
* Initializes a MA table
*/

void table_MA_init(table_MA * T);


/*
* long int  table_MA_insert (table_MA *T, list_MA *ma, long int blastscore, list_MA *l)
*
* Insert a MA in a table of MA
*/

long int table_MA_insert(table_MA * T, list_MA * ma, long int blastscore,
                         list_MA * l);


/*
* int  table_MA_delete (table_MA *T, list_MA *ma);
*
* Delete a MA from a table of MA
*/

long int table_MA_delete(table_MA * T, list_MA * ma);

void display_table(table_MA * t);

#endif
