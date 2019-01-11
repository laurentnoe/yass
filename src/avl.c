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
#ifdef DEBUG_AVL
#include <assert.h>
#endif

#include "global_var.h"
#include "avl.h"

/* Print k spaces */
void avl_print_espace(long int k)
{
    while (k != 0) {
      printf("|__");
      k--;
    }
}

/* Display a tree */

void avl_display_tree(avl_node * n, long int h)
{

    long int i = 0;
    list_MA *l;

    if (n != NULL) {

      l = (n->data).queue.first;
      avl_display_tree(n->avl_link[0], h + 1);
      avl_print_espace(h);
      while (l) {
      printf("%p ", (void *)l);
      i++;
      l = l->next;
      }
      printf("\n");
      avl_print_espace(h);
      printf("(dist=%ld, elems=%ld, bal = %d)", (n->data).distance, i,
           n->avl_balance);
      l = (n->data).queue.first;
      while (l) {
      printf("|%ld", l->ma->right_pos_end);
      l = l->next;
      }
      printf("\n");
      avl_display_tree(n->avl_link[1], h + 1);

    } else {

      avl_print_espace(h);
      printf("X\n");

    }
}

long int avl_check_tree(avl_node * n)
{

    long int h1;
    long int h2;

    if (n == NULL) {
      return 0;
    }

    else {

      h1 = avl_check_tree(n->avl_link[0]);
      h2 = avl_check_tree(n->avl_link[1]);

      if (((h2 - h1) < -1) || ((h2 - h1) > 1)
          || ((h2 - h1) != n->avl_balance)) {

        printf(" ******* ERREUR D'ARBRE ********* %lx %lx\n", h1, h2);
        avl_display_tree(n, 0);
        exit(0);
          return 0;

      } else {

          return 1 + MAX(h1, h2);

      }

    }

}


/*
* tree_data * avl_find (avl_tree *T, int distance)
*
* Search for a MA in the tree, if it is not present it returns the just smaller value
*/

tree_data *avl_find(avl_tree * T, long int distance)
{

    avl_node *y = T->root;
    avl_node *sn = NULL;
    avl_node *prec = NULL;

    while (y != NULL) {

      /*if the value is present we return the data */
      if ((y->data).distance == distance) {
          return &(y->data);
      }

      /*if not we search for the just smaller value */
      if ((y->data).distance < distance) {

          sn = y;
          y = y->avl_link[1];

      } else {

          prec = y;
          y = y->avl_link[0];

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
 * avl_node * new_avl_node(avl_tree *tree)
 *
 * Allocate a new node
 */

avl_node *new_avl_node(avl_tree * tree)
{

    avl_node *n;

    /*The allocation */

    n = (avl_node *)MALLOC(sizeof(avl_node));
    ASSERT(n, "new_avl_node");

    n->avl_link[0] = NULL;
    n->avl_link[1] = NULL;
    n->avl_balance = 0;

    return n;
}



/*
* tree_data * avl_insert (avl_tree *tree, long int distance)
*
* Insert a MA in the tree
*
* Taken from the GNU libavl 2.0.1
*/

tree_data *avl_insert(avl_tree * tree, long int distance)
{

    avl_node *y;
    avl_node *z;            /* Top node to uptree->date balance factor, and parent. */
    avl_node *p;
    avl_node *q;            /* Iterator, and parent. */
    avl_node *n;            /* Newly inserted node. */
    avl_node *w;            /* New root of rebalanced subtree. */
    avl_node *x = NULL;

    unsigned char dir;                  /* Direction to descend. */
    long int cmp;

    tree->k = 0;
    tree->pa[tree->k] = 0;

    z = NULL;
    y = tree->root;
    dir = 0;

    /* search */
    for (q = z, p = y; p != NULL; q = p, p = p->avl_link[dir]) {

      cmp = COMPARE(distance, (p->data).distance);

      if (distance == (p->data).distance)
          return &(p->data);



      if (p->avl_balance != 0) {
          z = q;
          y = p;
          tree->k = 0;
      }
      tree->da[tree->k++] = dir = cmp > 0;

      if (tree->k >= TREE_MAX_HEIGHT) {
        _ERROR("maximal avl tree heigh exceeded, (try with -w 0 parameter)\n");
      }
    }



    /* insert */
    if (q != NULL) {
      n = q->avl_link[dir] = new_avl_node(tree);
    } else {
      n = tree->root = new_avl_node(tree);
    }
    (n->data).queue.first = NULL;
    (n->data).queue.last  = NULL;


    /* rebalance */
    if (y == NULL)
      return &(n->data);

    for (p = y, tree->k = 0; p != n;
       p = p->avl_link[tree->da[tree->k]], tree->k++) {
      if (tree->da[tree->k] == 0) {
          p->avl_balance--;
      } else {
          p->avl_balance++;
      }
    }

    if (y->avl_balance == -2) {

      x = y->avl_link[0];

      if (x->avl_balance == -1) {

          w = x;
          y->avl_link[0] = x->avl_link[1];
          x->avl_link[1] = y;
          x->avl_balance = y->avl_balance = 0;

      } else {
#ifdef DEBUG_AVL
          assert(x->avl_balance == +1);
#endif
          w = x->avl_link[1];
          x->avl_link[1] = w->avl_link[0];
          w->avl_link[0] = x;
          y->avl_link[0] = w->avl_link[1];
          w->avl_link[1] = y;

          if (w->avl_balance == -1) {

            x->avl_balance = 0, y->avl_balance = +1;

          } else {

            if (w->avl_balance == 0) {

                x->avl_balance = y->avl_balance = 0;

            } else {      /* w->avl_balance == +1 */

                x->avl_balance = -1, y->avl_balance = 0;

            }

          }

          w->avl_balance = 0;

      }

    } else {

      if (y->avl_balance == +2) {

          x = y->avl_link[1];

          if (x->avl_balance == +1) {

            w = x;
            y->avl_link[1] = x->avl_link[0];
            x->avl_link[0] = y;
            x->avl_balance = y->avl_balance = 0;

          } else {

#ifdef DEBUG_AVL
            assert(x->avl_balance == -1);
#endif
            w = x->avl_link[0];
            x->avl_link[0] = w->avl_link[1];
            w->avl_link[1] = x;
            y->avl_link[1] = w->avl_link[0];
            w->avl_link[0] = y;

            if (w->avl_balance == +1) {

                x->avl_balance = 0, y->avl_balance = -1;

            } else {

                if (w->avl_balance == 0) {

                  x->avl_balance = y->avl_balance = 0;

                } else {      /* w->avl_balance == -1 */

                  x->avl_balance = +1, y->avl_balance = 0;

                }

            }

            w->avl_balance = 0;

          }

      } else {

          return &(n->data);

      }

    }

    if (z != NULL) {
      z->avl_link[y != z->avl_link[0]] = w;
    } else {
      tree->root = w;
    }

    return &(n->data);

}



/*
 *   tree_data * rb_search (rb_tree *tree, long int distance)
 *
 * Search a distance in the tree and fill two table : one with the path and the other with the nodes constituting the path
 */


tree_data *avl_search(avl_tree * tree, long int distance)
{

    avl_node *p = NULL;
    long int cmp = 0;
    tree->k = 1;

    for (p = tree->root; p != NULL; p = p->avl_link[tree->da[tree->k - 1]]) {

      if (distance == (p->data).distance) {
          /* if the distance already exists, then we return the data and put it in to_delete */
          tree->to_delete = p;
          return &(p->data);
      }

      cmp = COMPARE(distance, (p->data).distance);
      tree->pa[tree->k] = p;
      tree->da[tree->k++] = cmp > 0;

      if (tree->k >= TREE_MAX_HEIGHT) {
        _ERROR("maximal avltree heigh exceeded, (try with -w 0 parameter)");
      }
    }

    return NULL;

}



/*
 * void avl_delete_rebalancing (avl_tree *tree)
 *
 * Rebalance the tree after deletion
 *
 * Taken from the GNU libavl 2.0.1
 */


void avl_delete_rebalancing(avl_tree * tree)
{
    long int j = 0;
    avl_node *p = NULL;            /* Traverses tree to find node to delete. */
    avl_node *r = NULL;
    avl_node *s = NULL;
    avl_node *y = NULL;
    avl_node *x = NULL;
    avl_node *w = NULL;

    p = tree->to_delete;

    if (p->avl_link[1] == NULL) {

      if (tree->pa[tree->k - 1] != NULL) {
          tree->pa[tree->k - 1]->avl_link[tree->da[tree->k - 1]] =
            p->avl_link[0];
      } else {
          tree->root = p->avl_link[0];
      }
    } else {

      r = p->avl_link[1];
      if (r->avl_link[0] == NULL) {

          r->avl_link[0] = p->avl_link[0];
          r->avl_balance = p->avl_balance;

          if (tree->pa[tree->k - 1] != NULL) {
            tree->pa[tree->k - 1]->avl_link[tree->da[tree->k - 1]] = r;
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
            s = r->avl_link[0];

            if (s->avl_link[0] == NULL)
                break;
            r = s;
          }

          s->avl_link[0] = p->avl_link[0];
          r->avl_link[0] = s->avl_link[1];
          s->avl_link[1] = p->avl_link[1];
          s->avl_balance = p->avl_balance;

          if (tree->pa[j - 1] != NULL) {
            tree->pa[j - 1]->avl_link[tree->da[j - 1]] = s;
          } else {
            tree->root = s;
          }
          tree->da[j] = 1;
          tree->pa[j] = s;

      }
    }

    memset(p, 0x00, sizeof(avl_node));
    FREE(p,sizeof(avl_node));
    tree->to_delete = NULL;
#ifdef DEBUG_AVL
    assert(tree->k > 0);
#endif

    while (--tree->k > 0) {

      y = tree->pa[tree->k];

      if (tree->da[tree->k] == 0) {

          y->avl_balance++;
          if (y->avl_balance == +1) {
            break;
          } else {
            if (y->avl_balance == +2) {
                x = y->avl_link[1];
#ifdef DEBUG_AVL
                assert(x != NULL);
#endif
                if (x->avl_balance == -1) {
#ifdef DEBUG_AVL
                  assert(x->avl_balance == -1);
#endif
                  w = x->avl_link[0];
                  x->avl_link[0] = w->avl_link[1];
                  w->avl_link[1] = x;
                  y->avl_link[1] = w->avl_link[0];
                  w->avl_link[0] = y;

                  if (w->avl_balance == +1) {
                      x->avl_balance = 0, y->avl_balance = -1;
                  } else {

                      if (w->avl_balance == 0) {
                        x->avl_balance = y->avl_balance = 0;
                      } else {      /* w->avl_balance == -1 */
                        x->avl_balance = +1, y->avl_balance = 0;
                      }
                  }
                  w->avl_balance = 0;

                  if (tree->pa[tree->k - 1] != NULL) {
                      tree->pa[tree->k -
                             1]->avl_link[tree->da[tree->k - 1]] =
                        w;
                  } else {
                      tree->root = w;
                  }
                } else {
                  y->avl_link[1] = x->avl_link[0];
                  x->avl_link[0] = y;

                  if (tree->pa[tree->k - 1] != NULL) {
                      tree->pa[tree->k -
                             1]->avl_link[tree->da[tree->k - 1]] =
                        x;
                  } else {
                      tree->root = x;
                  }

                  if (x->avl_balance == 0) {
                      x->avl_balance = -1;
                      y->avl_balance = +1;
                      break;
                  } else {
                      x->avl_balance = y->avl_balance = 0;
                  }
                }
            }
          }
      } else {            /*da[k] != 0 */
          y->avl_balance--;
          if (y->avl_balance == -1) {
            break;
          } else {
            if (y->avl_balance == -2) {
                x = y->avl_link[0];
#ifdef DEBUG_AVL
                assert(x != NULL);
#endif
                if (x->avl_balance == +1) {

#ifdef DEBUG_AVL
                  assert(x->avl_balance == +1);
#endif
                  w = x->avl_link[1];
                  x->avl_link[1] = w->avl_link[0];
                  w->avl_link[0] = x;
                  y->avl_link[0] = w->avl_link[1];
                  w->avl_link[1] = y;

                  if (w->avl_balance == -1) {
                      x->avl_balance = 0, y->avl_balance = +1;
                  } else {
                      if (w->avl_balance == 0) {
                        x->avl_balance = y->avl_balance = 0;
                      } else {      /* w->avl_balance == +1 */
                        x->avl_balance = -1, y->avl_balance = 0;
                      }
                  }

                  w->avl_balance = 0;
                  if (tree->pa[tree->k - 1] != NULL) {
                      tree->pa[tree->k -
                             1]->avl_link[tree->da[tree->k - 1]] =
                        w;
                  } else {
                      tree->root = w;
                  }
                } else {
                  y->avl_link[0] = x->avl_link[1];
                  x->avl_link[1] = y;

                  if (tree->pa[tree->k - 1] != NULL) {
                      tree->pa[tree->k -
                             1]->avl_link[tree->da[tree->k - 1]] =
                        x;
                  } else {
                      tree->root = x;
                  }
                  if (x->avl_balance == 0) {
                      x->avl_balance = +1;
                      y->avl_balance = -1;
                      break;
                  } else {
                      x->avl_balance = y->avl_balance = 0;
                  }
                }
            }
          }
      }
    }
}
