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

#ifndef __PRDYN_H_
#define __PRDYN_H_

#include "global_var.h"
#include "tuple.h"
#include "threads.h"


long int initialise_alignment(long int size,
                              Feature *feature);

/* I) alignment methods:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 * WARNING : use the preallocated table provided by "initialise_alignment"
 */


/* I.a) Extensions ("left" or "right") methods :
 */

long int left_alignment_SG(char *data1, long int pos1, long int len1,
                           char *data2, long int pos2, long int len2, long int XDROP,
                           /*out */ long int *p_gainscore,
                           /*out */ long int *p_last_pos1,
                           /*out */ long int *p_last_pos2,
                           Feature *feature);

long int right_alignment_SG(char *data1, long int pos1, long int len1,
                            char *data2, long int pos2, long int len2, long int XDROP,
                            /*out */ long int *p_gainscore,
                            /*out */ long int *p_last_pos1,
                            /*out */ long int *p_last_pos2,
                            Feature *feature);

long int left_alignment_SG_Diag (char * data1, long int pos1, long int len1,
                                 char * data2, long int pos2, long int len2,
                                 long int XDROP,
                                 /*out*/ long int * p_gainscore,
                                 /*out*/ long int * p_last_pos1,
                                 /*out*/ long int * p_last_pos2,
                                 Feature *feature);

long int right_alignment_SG_Diag (char * data1, long int pos1, long int len1,
                                  char * data2, long int pos2, long int len2,
                                  long int XDROP,
                                  /*out*/ long int * p_gainscore,
                                  /*out*/ long int * p_last_pos1,
                                  /*out*/ long int * p_last_pos2,
                                  Feature *feature);

long int left_alignment_SG_Lz(char *data1, long int pos1, long int len1,
                              char *data2, long int pos2, long int len2, long int XDROP,
                              /*out */ long int *p_gainscore,
                              /*out */ long int *p_last_pos1,
                              /*out */ long int *p_last_pos2,
                              Feature *feature);

long int right_alignment_SG_Lz(char *data1, long int pos1, long int len1,
                               char *data2, long int pos2, long int len2, long int XDROP,
                               /*out */ long int *p_gainscore,
                               /*out */ long int *p_last_pos1,
                               /*out */ long int *p_last_pos2,
                               Feature *feature);

long int left_alignment_SG_Border(char *data1, long int pos1, long int len1,
                                  char *data2, long int pos2, long int len2, long int XDROP,
                                  /*out */ long int *p_gainscore,
                                  /*out */ long int *p_last_pos1,
                                  /*out */ long int *p_last_pos2,
                                  Feature *feature);

long int right_alignment_SG_Border(char *data1, long int pos1, long int len1,
                                   char *data2, long int pos2, long int len2, long int XDROP,
                                   /*out */ long int *p_gainscore,
                                   /*out */ long int *p_last_pos1,
                                   /*out */ long int *p_last_pos2,
                                   Feature *feature);

/* I.b) Alignment methods :
 */

long int alignment_SG(char *data1, long int pos1, long int len1,
                      char *data2, long int pos2, long int len2,
                      Feature *feature);

long int alignment_SG_Strait(char *data1, long int pos1,
                             char *data2, long int pos2, long int len,
                             Feature *feature);

long int alignment_SG_DROP(char *data1, long int pos1, long int len1,
                           char *data2, long int pos2, long int len2, long int XDROP,
                           Feature *feature);

long int alignment_SG_DROP_opt(char *data1, long int pos1, long int len1,
                               char *data2, long int pos2, long int len2, long int XDROP,
                               Feature *feature);

long int alignment_SG_Lz(char *data1, long int pos1, long int len1,
                         char *data2, long int pos2, long int len2, long int XDROP,
                         Feature *feature);

long int alignment_SG_Border(char *data1, long int pos1, long int len1,
                             char *data2, long int pos2, long int len2, long int XDROP,
                             Feature *feature);

long int alignment_SG_on_MA(char *data1, char *data2, MA * ma,long int left_correction,
                            Feature *feature);


/*-----------------------------------------------------------------*/
/*
 * II) Statistical functions:
 *  - (compute transition/transversion ratio)
 *  - (compute triplet mutation frequency)
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */


/* II.a) Extensions ("left" or "right") methods :
 */


long int left_alignment_SG_stats_Border(char * data1, long int pos1, long int len1,
                                        char * data2, long int pos2, long int len2,
                                        unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                        MA   * ma,
                                        Feature * feature);

long int right_alignment_SG_stats_Border(char * data1, long int pos1, long int len1,
                                         char * data2, long int pos2, long int len2,
                                         unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                         MA   * ma,
                                         Feature * feature);

/* II.b) Alignment methods :
 */

long int alignment_SG_stats(char * data1, long int pos1, long int len1,
                            char * data2, long int pos2, long int len2,
                            unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                            MA   * ma,
                            Feature * feature);

long int alignment_SG_stats_Border(char * data1, long int pos1, long int len1,
                                   char * data2, long int pos2, long int len2,
                                   unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                   MA   * ma,
                                   Feature * feature);

long int alignment_SG_stats_Strait(char * data1,  long int pos1, char * data2,  long int pos2, long int len,
                                   unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                   MA * ma,
                                   Feature * feature);

long int alignment_SG_stats_on_MA(char * data1, char * data2,
                                  long int left_correction,
                                  MA * ma,
                                  Feature * feature,
                                  long int full);

/*-----------------------------------------------------------------*/
/*
 * III) PSL output functions:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */


/* III.a) Extensions ("left" or "right") methods :
 */


long int left_alignment_SG_PSL_Border(char * data1, long int pos1, long int len1,
                                      char * data2, long int pos2, long int len2,
                                      /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                      /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                      long int choice);

long int right_alignment_SG_PSL_Border(char * data1, long int pos1, long int len1,
                                       char * data2, long int pos2, long int len2,
                                       /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                       /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                       long int choice);

/* III.b) Alignment methods :
 */

long int alignment_SG_PSL(char * data1, long int pos1, long int len1,
                          char * data2, long int pos2, long int len2,
                          /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                          /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                          long int choice);


long int alignment_SG_PSL_Border(char * data1, long int pos1, long int len1,
                                 char * data2, long int pos2, long int len2,
                                 /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                 /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                 long int choice);


long int alignment_SG_PSL_Strait(char * data1, long int pos1,
                                 char * data2, long int pos2, long int len,
                                 /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                 /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                 long int choice);

long int alignment_SG_PSL_on_MA(char * data1, char * data2,
                                long int left_correction,
                                MA * ma,
                                long int choice);

/*-----------------------------------------------------------------*/

/* IV) YASS output functions:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */



/* IV.a) Extensions ("left" or "right") methods :
 */

long int display_left_alignment_SG_Border_noflush(long int querychunk, long int reverse,
                                                  char * data1, long int pos1, long int len1,
                                                  char * data2, long int pos2, long int len2);

long int display_right_alignment_SG_Border_noflush(long int querychunk, long int reverse,
                                                   char * data1, long int pos1, long int len1,
                                                   char * data2, long int pos2, long int len2);


/* IV.b) Alignment methods :
 */

long int display_alignment_SG_Border_noflush(long int querychunk, long int reverse,
                                             char * data1, long int pos1, long int len1,
                                             char * data2, long int pos2, long int len2);

long int display_alignment_SG_noflush(long int querychunk, long int reverse,
                                      char *data1, long int pos1, long int len1,
                                      char *data2, long int pos2, long int len2);

long int display_alignment_SG_Strait_noflush(long int querychunk, long int reverse,
                                             char *data1, long int pos1,
                                             char *data2, long int pos2, long int len);

long int display_alignment_SG(long int querychunk, long int reverse,
                              char *data1, long int pos1, long int len1,
                              char *data2, long int pos2, long int len2);


long int display_alignment_SG_on_MA(char *data1, char *data2, MA * ma);

/*-----------------------------------------------------------------*/

/* V) Score debugging functions:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */

/* V.a) Extensions ("left" or "right") methods :
 */


/* V.b) Alignment methods :
 */

long int alignment_SG_score_on_MAs(char * data1, char * data2,MA * firstMa,
                                   long int left_correction);

long int alignment_SG_score_on_MA(char * data1, char * data2,
                                  MA * ma,long int left_correction);

long int alignment_SG_score_Strait(char * data1,  long int pos1,
                                   char * data2,  long int pos2, long int len,
                                   /*out*/long int    *p_score);

long int alignment_SG_score(char * data1, long int pos1, long int len1,
                            char * data2, long int pos2, long int len2,
                            /*out*/long int    *p_score);


#endif
