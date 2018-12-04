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

#ifndef  __GLOBAL_VAR_H_
#define  __GLOBAL_VAR_H_
#include <time.h>
#include "util.h"
#include "tuple.h"



/********************
 * GLOBALVARIABLES  *
 ********************/


/* [1]  MODIFIED WHEN STARTING PROGRAM
 *========================================*/

/* parameters and flags */
extern char *gp_motifs[MAX_SEED];
extern long int gp_seeds_span[MAX_SEED];
extern long int gp_seeds_bitweight[MAX_SEED];
extern long int gp_seeds_span_min;
extern long int gp_seeds_span_max;
extern long int gp_nb_seeds;

extern long int gp_mutations_percent;
extern double gp_mutations;
extern long int gp_indels_percent;
extern double gp_indels;
extern long int gp_alpha_percent;
extern double gp_alpha;

extern double gp_entropy_min;
extern double gp_expectation_value;

extern long int gp_hitcriterion;
extern long int gp_reverse;
extern long int gp_display;
extern long int gp_nbmaxlines;
extern long int gp_distdiag;
extern long int gp_lowercase;
extern double gp_t;

/* statisticals parameters */
extern long int gp_rho_stat;
extern long int gp_delta_stat;
extern long int gp_border;

/* post-processing */
extern long int gp_win_min;
extern long int gp_win_max;
extern double gp_win_mul;

/* sort function index */
extern long int          gp_sortcriterion;
extern SortCrit *        gp_sortcriterion_func;
extern long int          gp_sortblockscriterion;
extern SortBlocksCrit *  gp_sortblockscriterion_func;

/* files number et filenames */
extern long int gp_nbfiles;
extern char gp_files[2][16384];
extern long int gp_selection_fasta;


/* cost */
extern long int SUBMATRIXTABLE[NBMATRICES][4];
extern long int INDELSTABLE[NBMATRICES][2];
extern long int gp_cost_gap_opened;
extern long int gp_cost_gap_continued;
extern long int gp_costs[4];
extern long int ** gp_substitution_matrix;
extern long int gp_cost_max_substitution_matrix;
extern long int gp_matrix;
extern long int gp_adhoc_matrix;
extern double gp_k_blast;
extern double gp_lambda_blast;
extern long int gp_xdrop;

/* data tables and their sizes  */
extern char *gp_query;
extern char *gp_query_rev;
extern long int gp_querysize;

extern char *gp_text;
extern long int gp_textsize;

/* first file chunks */
extern long int     gp_nbchunks_query;
extern char ** gp_chunkname_query;
extern long int   * gp_chunksize_query;
extern long int   * gp_chunkstrt_query;

/* second file chunks */
extern long int     gp_nbchunks_text;
extern char ** gp_chunkname_text;
extern long int   * gp_chunksize_text;
extern long int   * gp_chunkstrt_text;

/* statistics on each sequence */
extern long int       gp_nb_letters[2][4];
extern double         gp_freq_letters[2][4];
extern long int       gp_nb_triplets[2][64];
extern double         gp_freq_background[4][4];
extern double         gp_freq_tripletbackground[64][64];

/* lockup/code reverals and other stuff */
extern char lookup[32];
extern long int  backcode[32];
extern char complement[32];
extern long int  unindexable[32];
extern SortCrit* sortcriteria[NBSORTCRITERIA];
extern SortBlocksCrit*  sortblockscriteria[NBSORTBLOCKSCRITERIA];
#ifdef TRACE
extern long int * gp_dots;
#endif

#ifdef MEM_ALLOCATED
extern unsigned long gv_mem_allocated;
extern unsigned long gp_max_mem_allocated;
#endif


/* [2]   MODIFIED WHEN RUNNING ...
 *================================*/
extern FILE *     gv_outstream;
extern long int   gv_thread_result;
extern long int * gv_delta_shift;

/* statistics on time spent */
#ifdef STATS
extern time_t gv_time_spent;
#endif

extern long int gv_chunk_nb;
extern long int gv_chunk_nb_end;
extern long int gv_thread_num[MAX_QUERY_CHUNK_THREADS];
extern MA *     gv_first_MA;
extern MA *     gv_last_MA;

extern long int gv_last_print_is_a_dot;
#endif
