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

#ifndef __RB_H__
#define __RB_H__

#include "list.h"
#include "util.h"

/* color of a node, can be red or black */
typedef enum {
    RB_BLACK,                  /* Black. */
    RB_RED                  /* Red. */
} color;

/* A red-black tree node */
typedef struct _rb_node {
    color rb_color;
    struct _rb_node *rb_link[2];
    tree_data data;
} rb_node;


/* The tree structure, it contains a root, the current position in the tree, the path to the
 * current node and all the precedent nodes and a node to delete
 */
typedef struct _rb_tree {
    struct _rb_node *root;
    long int k;
    struct _rb_node *pa[TREE_MAX_HEIGHT];      /* Nodes on stack. */
    unsigned char da[TREE_MAX_HEIGHT];      /* Directions moved from stack nodes. */
    struct _rb_node *to_delete;
} rb_tree;


/*
 * void display_tree (rb_node *n, long int h)
 *
 * Display a tree
 */

void rb_display_tree(rb_node * n, long int h);


/*
 * tree_data * rb_find (rb_tree *T, long int distance)
 *
 * Search for a MA distance in the tree, if it is not present it returns the just smaller value
 */

tree_data *rb_find(rb_tree * T, long int distance);

/*
 * tree_data * rb_search (rb_tree *tree, long int distance)
 *
 * Search a distance in the tree and fill two table : one with the path and the other with the nodes constituting the path
 */

tree_data *rb_search(rb_tree * tree, long int distance);

/*
 * tree_data *rb_insert(rb_tree * tree, long int distance);
 *
 * Insert a MA in a tree
 */

tree_data *rb_insert(rb_tree * tree, long int distance);


/*
 * void rb_delete_rebalancing (rb_tree *tree)
 *
 * Rebalance the tree after deletion
 *
 * Taken from the GNU libavl 2.0.1
 */

void rb_delete_rebalancing(rb_tree * tree);



#endif
