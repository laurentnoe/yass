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

#ifndef __REGROUP_H__
#define __REGROUP_H__

#define DIAG_LENGTH(map,mai)    (MIN((mai)->right_pos_begin - (map)->right_pos_end, (mai)->left_pos_begin-(map)->left_pos_end))
#define INDEL_LENGTH(dp,di)     (ABS((DIAGE(di))-(DIAGB(dp)))) /*FIXME : need to be simplified */
#define DIAG_COST(l)            (-ABS(gp_cost_max_substitution_matrix * l))
#define INDEL_COST(i)           (gp_cost_gap_opened+(i)*gp_cost_gap_continued)
#define GOOD_LINK(l1,l2)        ((l1)->ma->blastscore + (l2)->ma->blastscore + DIAG_COST(DIAG_LENGTH((l1)->ma, (l2)->ma)) + INDEL_COST(INDEL_LENGTH((l1)->ma,(l2)->ma)) > MAX(((l1)->ma->blastscore),((l2)->ma->blastscore)))
#define SCORE_LINK(l1,l2)      (DIAG_COST(DIAG_LENGTH((l1)->ma, (l2)->ma)) + INDEL_COST(INDEL_LENGTH((l1)->ma,(l2)->ma)))
#define SCORE_ALIGN(l1,l2)      ((l1)->ma->blastscore + (l2)->ma->blastscore + DIAG_COST(DIAG_LENGTH((l1)->ma, (l2)->ma)) + INDEL_COST(INDEL_LENGTH((l1)->ma,(l2)->ma)))
#define INDEL_BOUND(x)          ((((x)->blastscore) / -gp_cost_gap_continued) + (gp_cost_gap_opened/gp_cost_gap_continued) + 1)


#ifdef CHOOSERBTREE
#define TREE_TYPE(x) rb##x
#include "red_black.h"
#else
#ifdef CHOOSEAVLTREE
#define TREE_TYPE(x) avl##x
#include "avl.h"
#else
#error No tree definition : uncomment either "CHOOSERBTREE" or "CHOOSEAVLTREE" in the "util.h" file
#endif
#endif


#include "global_var.h"
#include "threads.h"
#include "list.h"

/* Element use in statistical functions */
typedef struct _elem_stat {
    float size;
    long int weight;
} elem;


/*
 *
 *
 *
 *
 *
 *
 *             MA1                          MA2
 * |----------------------|.....|--------------------------|
 *  * list_MA_left :
 *     - customer :
 *        [ , , , ]
 *         |__________________________^
 *     - producer :
 *        [ , , , ]
 *
 *
 */


/*
 * void window_deleter(long int pos,  Feature *feature);
 *
 * Delete all the MA which are no longer in the windows at position "pos"
 */

void window_deleter(long int pos,  Feature *feature);


/*
 * void regroup(list_MA * listes, Feature *feature);
 *
 * Test with blascore and distance criterions if a MA can be link with an
 * other MA and proceed if possible
 */

void regroup(list_MA * listes, Feature *feature);


/*
 * Proceed the global regroupement
 */

MA * main_regroup(MA * first_MA, Feature *feature);


/*
 * void Regroup_MAList(Feature *feature);
 *
 * Proceed multi regroupement
 */

void Regroup_MAList(Feature *feature);

#endif
