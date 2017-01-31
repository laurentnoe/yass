/*
 *  YASS 1.15
 *  Copyright (C) 2004-2017
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

#ifndef __AVL_H__
#define __AVL_H__

#include "list.h"
#include "util.h"

/* An AVL tree node. */
typedef struct _avl_node {
    signed char avl_balance;      /* Balance factor. */
    struct _avl_node *avl_link[2];      /* Subtrees. */
    tree_data data;
} avl_node;

/* The tree structure, it contains a root, the current position in the tree, the path to the
 * current node and all the precedent nodes and a node to delete
 */
typedef struct _avl_tree {
    struct _avl_node *root;
    long int k;
    struct _avl_node *pa[TREE_MAX_HEIGHT];      /* Nodes on stack. */
    unsigned char da[TREE_MAX_HEIGHT];      /* Directions moved from stack nodes. */
    struct _avl_node *to_delete;
} avl_tree;


/*
 * void display_tree (avl_node *n, long int h)
 *
 * Display a tree
 */

void avl_display_tree(avl_node * n, long int h);


/*
 * tree_data *avl_find(avl_tree * T, long int distance);
 *
 * Search for a MA in the tree, if it is not present it returns the just smaller value
 */

tree_data *avl_find(avl_tree * T, long int distance);


/*
 * tree_data * avl_search (avl_tree *tree, long int distance)
 *
 * Search a distance in the tree and fill two table : one with the path and the other with the nodes constituting the path
 */

tree_data *avl_search(avl_tree * tree, long int distance);

/*
 * tree_data *avl_insert(avl_tree * tree, long int distance);
 *
 * Insert a MA in a tree
 */

tree_data *avl_insert(avl_tree * tree, long int distance);

/*
 * void avl_delete_rebalancing (avl_tree *tree)
 *
 * Rebalance the tree after deletion
 *
 * Taken from the GNU libavl 2.0.1
 */

void avl_delete_rebalancing(avl_tree * tree);


#endif
