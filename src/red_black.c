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

#include "global_var.h"
#include "tuple.h"
#include "red_black.h"

/* Print k spaces */

void rb_print_espace(long int k)
{

    while (k != 0) {

      printf("|__");
      k--;

    }
}

/* Display a tree */

void rb_display_tree(rb_node * n, long int h)
{

    long int i = 0;
    list_MA *l;

    if (n != NULL) {

      l = (n->data).queue.first;
      rb_display_tree(n->rb_link[0], h + 1);
      rb_print_espace(h);
      while (l) {
      i++;
      l = l->next;
      }
      printf("(dist=%ld, elems=%ld)", (n->data).distance, i);
      l = (n->data).queue.first;
      while (l) {
      printf("|%ld", l->ma->right_pos_end);
      l = l->next;
      }
      printf("\n");
      rb_display_tree(n->rb_link[1], h + 1);

    } else {

      rb_print_espace(h);
      printf("X\n");

    }
}


/*
 * tree_data * rb_find (rb_tree *T, long int distance)
 *
 * Search for a MA distance in the tree, if it is not present it returns the just smaller value
 */

tree_data *rb_find(rb_tree * T, long int distance)
{

    rb_node *y = T->root;
    rb_node *sn = NULL;
    rb_node *prec = NULL;

    while (y != NULL) {

      /*if the value is present we return the data */
      if ((y->data).distance == distance) {
          return &(y->data);
      }

      /*if not we search for the just smaller value */
      if ((y->data).distance < distance) {

          sn = y;
          y = y->rb_link[1];

      } else {

          prec = y;
          y = y->rb_link[0];

      }
    }

    /*if the value is the smallest */
    if (sn == NULL) {
      if (prec != NULL) {
          /*we return prec which is the just bigger value */
          return &(prec->data);
      } else {
          /*if the tree is empty we return NULL */
          return NULL;
      }
    } else {
      /*we return the just smaller value */
      return &(sn->data);
    }

}


/*
 * rb_node * new_rb_node(rb_tree *tree)
 *
 * Allocate a new node and insert it in the tree
 */

rb_node *new_rb_node(rb_tree * tree)
{

    rb_node *n;

    /*The allocation */

    n = (rb_node *) MALLOC(sizeof(rb_node));
    ASSERT(n, "new_rb_node");

    n->rb_link[0] = NULL;
    n->rb_link[1] = NULL;
    n->rb_color   = RB_RED;

    /*The insertion */

    if (tree->root == NULL) {
      tree->root = n;
    } else {
      tree->pa[tree->k - 1]->rb_link[tree->da[tree->k - 1]] = n;
    }

    return n;

}


/*
 * tree_data * rb_search (rb_tree *tree, long int distance)
 *
 * Search a distance in the tree and fill two table : one with the path and the other with the nodes constituting the path
 */

tree_data *rb_search(rb_tree * tree, long int distance)
{

    rb_node *p;
    long int cmp;

    tree->pa[0] = NULL;
    tree->da[0] = 0;
    tree->k = 1;

    for (p = tree->root; p != NULL; p = p->rb_link[tree->da[tree->k - 1]]) {

      if (distance == (p->data).distance) {
          /* if the distance already exists, then we return the data and put it in to_delete */
          tree->to_delete = p;
          return &(p->data);
      }

      cmp = COMPARE(distance, (p->data).distance);
      tree->pa[tree->k] = p;
      tree->da[tree->k++] = cmp > 0;

      if (tree->k >= TREE_MAX_HEIGHT) {
          _ERROR("maximum tree heigh exceeded");
      }
    }

    return NULL;

}


/*
 * void rb_insert_rebalancing (rb_tree *tree)
 *
 * Rebalance the tree after insertion
 *
 * Taken from the GNU libavl 2.0.1
 */

void rb_insert_rebalancing(rb_tree * tree)
{

    rb_node *y;
    rb_node *x;

    while (tree->k >= 3 && tree->pa[tree->k - 1]->rb_color == RB_RED) {

      if (tree->da[tree->k - 2] == 0) {

          y = tree->pa[tree->k - 2]->rb_link[1];

          if (y != NULL && y->rb_color == RB_RED) {

            tree->pa[tree->k - 1]->rb_color = y->rb_color = RB_BLACK;
            tree->pa[tree->k - 2]->rb_color = RB_RED;
            tree->k -= 2;
          } else {


            if (tree->da[tree->k - 1] == 0)
                y = tree->pa[tree->k - 1];

            else {
                x = tree->pa[tree->k - 1];
                y = x->rb_link[1];
                x->rb_link[1] = y->rb_link[0];
                y->rb_link[0] = x;
                tree->pa[tree->k - 2]->rb_link[0] = y;

            }

            x = tree->pa[tree->k - 2];
            x->rb_color = RB_RED;
            y->rb_color = RB_BLACK;

            x->rb_link[0] = y->rb_link[1];
            y->rb_link[1] = x;
            if (tree->pa[tree->k - 3] != NULL) {
                tree->pa[tree->k - 3]->rb_link[tree->da[tree->k - 3]] =
                  y;
            } else {
                tree->root = y;
            }

            break;
          }

      } else {

          y = tree->pa[tree->k - 2]->rb_link[0];
          if (y != NULL && y->rb_color == RB_RED) {
            tree->pa[tree->k - 1]->rb_color = y->rb_color = RB_BLACK;
            tree->pa[tree->k - 2]->rb_color = RB_RED;
            tree->k -= 2;
          } else {


            if (tree->da[tree->k - 1] == 1) {

                y = tree->pa[tree->k - 1];

            } else {

                x = tree->pa[tree->k - 1];
                y = x->rb_link[0];
                x->rb_link[0] = y->rb_link[1];
                y->rb_link[1] = x;
                tree->pa[tree->k - 2]->rb_link[1] = y;

            }


            x = tree->pa[tree->k - 2];
            x->rb_color = RB_RED;
            y->rb_color = RB_BLACK;
            x->rb_link[1] = y->rb_link[0];
            y->rb_link[0] = x;

            if (tree->pa[tree->k - 3] != NULL) {
                tree->pa[tree->k - 3]->rb_link[tree->da[tree->k - 3]] =
                  y;
            } else {
                tree->root = y;
            }

            break;

          }

      }

    }

    tree->root->rb_color = RB_BLACK;

}


/*
* tree_data * rb_insert (rb_tree *tree, long int distance)
*
* Insert a node in the tree
*/

tree_data *rb_insert(rb_tree * tree, long int distance)
{

    rb_node *n;
    tree_data *data;

    /* We search for the right place to insert */
    data = rb_search(tree, distance);

    /* If the node doesn't exists */
    if (!data) {

      /* we create the node */
      n = new_rb_node(tree);
      (n->data).queue.first = NULL;
      (n->data).queue.last  = NULL;

      /* then we rebalance the tree */
      rb_insert_rebalancing(tree);
      return &(n->data);

    }

    return data;

}


/*
 * void rb_delete_rebalancing (rb_tree *tree)
 *
 * Rebalance the tree after deletion
 *
 * Taken from the GNU libavl 2.0.1
 */

void rb_delete_rebalancing(rb_tree * tree)
{
    long int j = 0;
    color t = RB_BLACK;
    rb_node *p = NULL;            /* The node to delete, or a node tree->part way to it. */
    rb_node *r = NULL;
    rb_node *s = NULL;
    rb_node *x = NULL;
    rb_node *w = NULL;
    rb_node *y = NULL;

    p = tree->to_delete;

    if (p->rb_link[1] == NULL) {

      if (tree->pa[tree->k - 1] != NULL) {
          tree->pa[tree->k - 1]->rb_link[tree->da[tree->k - 1]] =
            p->rb_link[0];
      } else {
          tree->root = p->rb_link[0];
      }

    } else {

      r = p->rb_link[1];
      if (r->rb_link[0] == NULL) {
          r->rb_link[0] = p->rb_link[0];
          t = r->rb_color;
          r->rb_color = p->rb_color;
          p->rb_color = t;

          if (tree->pa[tree->k - 1] != NULL) {
            tree->pa[tree->k - 1]->rb_link[tree->da[tree->k - 1]] = r;
          } else {
            tree->root = r;
          }


          tree->da[tree->k] = 1;
          tree->pa[tree->k++] = r;

      } else {
          j = tree->k++;

          while (1) {

            tree->da[tree->k] = 0;
            tree->pa[tree->k++] = r;
            s = r->rb_link[0];

            if (s->rb_link[0] == NULL)
                break;
            r = s;
          }

          tree->da[j] = 1;
          tree->pa[j] = s;

          if (tree->pa[j - 1] != NULL) {
            tree->pa[j - 1]->rb_link[tree->da[j - 1]] = s;
          } else {
            tree->root = s;
          }

          s->rb_link[0] = p->rb_link[0];
          r->rb_link[0] = s->rb_link[1];
          s->rb_link[1] = p->rb_link[1];

          t = s->rb_color;
          s->rb_color = p->rb_color;
          p->rb_color = t;

      }
    }


    if (p->rb_color == RB_BLACK) {

      while (1) {

          if (tree->pa[tree->k - 1] != NULL) {
            x = tree->pa[tree->k - 1]->rb_link[tree->da[tree->k - 1]];
          } else {
            x = tree->root;
          }


          if (x != NULL && x->rb_color == RB_RED) {
            x->rb_color = RB_BLACK;
            break;
          }

          if (tree->k < 2)
            break;

          if (tree->da[tree->k - 1] == 0) {

            w = tree->pa[tree->k - 1]->rb_link[1];

            if ((w != NULL) && (w->rb_color == RB_RED)) {


                w->rb_color = RB_BLACK;
                tree->pa[tree->k - 1]->rb_color = RB_RED;

                tree->pa[tree->k - 1]->rb_link[1] = w->rb_link[0];
                w->rb_link[0] = tree->pa[tree->k - 1];


                if (tree->pa[tree->k - 2] != NULL) {
                  tree->pa[tree->k -
                         2]->rb_link[tree->da[tree->k - 2]] = w;
                } else {
                  tree->root = w;
                }

                tree->pa[tree->k] = tree->pa[tree->k - 1];
                tree->da[tree->k] = 0;
                tree->pa[tree->k - 1] = w;
                tree->k++;

                w = tree->pa[tree->k - 1]->rb_link[1];

            }


            if (w == NULL) {
                break;
            }


            if ((w->rb_link[0] == NULL
                 || w->rb_link[0]->rb_color == RB_BLACK)
                && (w->rb_link[1] == NULL
                  || w->rb_link[1]->rb_color == RB_BLACK)) {
                w->rb_color = RB_RED;
            } else {
                if (w->rb_link[1] == NULL
                  || w->rb_link[1]->rb_color == RB_BLACK) {
                  y = w->rb_link[0];
                  y->rb_color = RB_BLACK;
                  w->rb_color = RB_RED;
                  w->rb_link[0] = y->rb_link[1];
                  y->rb_link[1] = w;
                  w = tree->pa[tree->k - 1]->rb_link[1] = y;
                }

                w->rb_color = tree->pa[tree->k - 1]->rb_color;
                tree->pa[tree->k - 1]->rb_color = RB_BLACK;
                w->rb_link[1]->rb_color = RB_BLACK;

                tree->pa[tree->k - 1]->rb_link[1] = w->rb_link[0];
                w->rb_link[0] = tree->pa[tree->k - 1];

                if (tree->pa[tree->k - 2] != NULL) {
                  tree->pa[tree->k -
                         2]->rb_link[tree->da[tree->k - 2]] = w;
                } else {
                  tree->root = w;
                }
                break;
            }

          } else {
            w = tree->pa[tree->k - 1]->rb_link[0];

            if ((w != NULL) && (w->rb_color == RB_RED)) {
                w->rb_color = RB_BLACK;
                tree->pa[tree->k - 1]->rb_color = RB_RED;

                tree->pa[tree->k - 1]->rb_link[0] = w->rb_link[1];
                w->rb_link[1] = tree->pa[tree->k - 1];

                if (tree->pa[tree->k - 2] != NULL) {
                  tree->pa[tree->k -
                         2]->rb_link[tree->da[tree->k - 2]] = w;
                } else {
                  tree->root = w;
                }

                tree->pa[tree->k] = tree->pa[tree->k - 1];
                tree->da[tree->k] = 1;
                tree->pa[tree->k - 1] = w;
                tree->k++;

                w = tree->pa[tree->k - 1]->rb_link[0];

            }

            if (w == NULL) {
                break;
            }

            if ((w->rb_link[0] == NULL
                 || w->rb_link[0]->rb_color == RB_BLACK)
                && (w->rb_link[1] == NULL
                  || w->rb_link[1]->rb_color == RB_BLACK)) {
                w->rb_color = RB_RED;
            } else {
                if (w->rb_link[0] == NULL
                  || w->rb_link[0]->rb_color == RB_BLACK) {
                  y = w->rb_link[1];
                  y->rb_color = RB_BLACK;
                  w->rb_color = RB_RED;
                  w->rb_link[1] = y->rb_link[0];
                  y->rb_link[0] = w;
                  w = tree->pa[tree->k - 1]->rb_link[0] = y;

                }

                w->rb_color = tree->pa[tree->k - 1]->rb_color;
                tree->pa[tree->k - 1]->rb_color = RB_BLACK;
                w->rb_link[0]->rb_color = RB_BLACK;

                tree->pa[tree->k - 1]->rb_link[0] = w->rb_link[1];

                w->rb_link[1] = tree->pa[tree->k - 1];

                if (tree->pa[tree->k - 2] != NULL) {
                  tree->pa[tree->k -
                         2]->rb_link[tree->da[tree->k - 2]] = w;
                } else {
                  tree->root = w;
                }
            }
          }
          tree->k--;
      }
    }

    memset(p, 0x00, sizeof(rb_node));
    FREE(p,sizeof(rb_node));
    tree->to_delete = NULL;

}
