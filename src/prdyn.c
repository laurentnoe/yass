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
#include "global_var.h"
#include "util.h"
#include "tuple.h"
#include "prdyn.h"

#define  CHAR_EQU '|'
#define  CHAR_SS  ':'
#define  CHAR_SV  '.'
#define  CHAR_DEL ' '
#define  CHAR_INS ' '
#define  CHAR_DUM '-'

#define TNORM(a) (((a)%16) == 7?(3):((a)%16))

/*
 * Run this method before any use of "quick" alignment algorithm
 * (it allocates dynamic programming table)
 *
 */

#ifdef INLINE
inline
#endif
long int initialise_alignment(long int size,
                              Feature *feature) {
  feature->buffer01 = (long int *)MALLOC((size+1)* sizeof(long int));
  ASSERT(feature->buffer01,initialise_aligment);
  return 0;
}



/* I) alignment methods:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 * WARNING : use the preallocated table provided by "initialise_alignment"
 */

#define I2(i,j) (feature->buffer01[(i)+ (((j)&1)?(len1+1):(0))])
#define M2(i,j) (feature->buffer01[(i)+ (((j)&1)?(len1+1):(0)) + ((len1+1)*2)])
#define S2(i,j) (gp_substitution_matrix[(long int)(data1[pos1+(i)-1])][ ((long int)(data2[pos2+(j)-1]))])


/* I.a) Extensions ("left" or "right") methods :
 */


#ifdef INLINE
inline
#endif
long int left_alignment_SG(char * data1, long int pos1, long int len1,
                           char * data2, long int pos2, long int len2,
                           long int XDROP,
                           /*out*/ long int * p_gainscore,
                           /*out*/ long int * p_last_pos1,
                           /*out*/ long int * p_last_pos2,
                           Feature *feature
                           )
{
  long int i,j;
  long int max = 0;

  *p_gainscore = 0;
  *p_last_pos1 = pos1;
  *p_last_pos2 = pos2;

  M2(0,0) = 0;
  I2(0,0) = gp_cost_gap_opened -  gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M2(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I2(i,0) = I2(i-1,0) +  gp_cost_gap_continued;

  for (j = 1; j <= len2; j++) {

    M2(0,j) =  -INFINITY_INT;
    I2(0,j) =  I2(0,j-1)  + gp_cost_gap_continued;
    max     =  M2(0,j);

    for (i = 1; i <= len1; i++) {
      I2(i,j) = MAX (
                     MAX( M2(i,j-1) + gp_cost_gap_opened    , M2(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I2(i,j-1) + gp_cost_gap_continued ,  I2(i-1,j) + gp_cost_gap_continued )
                     );
      M2(i,j) = MAX ( M2(i-1,j-1) + S2(-i+1,-j+1), I2(i-1,j-1) + S2(-i+1,-j+1) );
      max     = MAX(max,M2(i,j));
      if (max > *p_gainscore) {
        /* [SFS] : single file stop */
        if (data1 + pos1 - i == data2 + pos2 - j)
          return max;

        *p_gainscore = max;
        *p_last_pos1 = pos1 - i;
        *p_last_pos2 = pos2 - j;
      }
    }
    if (max < -XDROP)
      return max;
  }
  return max;
}



#ifdef INLINE
inline
#endif
long int right_alignment_SG(char * data1, long int pos1, long int len1,
                            char * data2, long int pos2, long int len2,
                            long int XDROP,
                            /*out*/ long int * p_gainscore,
                            /*out*/ long int * p_last_pos1,
                            /*out*/ long int * p_last_pos2,
                            Feature *feature
                            )
{
  long int i,j;
  long int max = 0;

  *p_gainscore = 0;
  *p_last_pos1 = pos1;
  *p_last_pos2 = pos2;

  M2(0,0) = 0;
  I2(0,0) = gp_cost_gap_opened -  gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M2(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I2(i,0) = I2(i-1,0) +  gp_cost_gap_continued;

  for (j = 1; j <= len2; j++) {

    M2(0,j) =  -INFINITY_INT;
    I2(0,j) =  I2(0,j-1)  + gp_cost_gap_continued;
    max     =  M2(0,j);

    for (i = 1; i <= len1; i++) {
      I2(i,j) = MAX (
                     MAX( M2(i,j-1) + gp_cost_gap_opened    , M2(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I2(i,j-1) + gp_cost_gap_continued ,  I2(i-1,j) + gp_cost_gap_continued )
                     );
      M2(i,j) = MAX ( M2(i-1,j-1) + S2(i,j), I2(i-1,j-1) + S2(i,j) );
      max     = MAX(max,M2(i,j));
      if (max > *p_gainscore) {
        /* [SFS] : single file stop */
        if (data1 + pos1 + i == data2 + pos2 + j)
          return max;

        *p_gainscore = max;
        *p_last_pos1 = pos1 + i;
        *p_last_pos2 = pos2 + j;
      }
    }
    if (max < -XDROP)
      return max;
  }
  return max;
}


#define D_I2(k,i) (feature->buffer01[(i) + ((k)&1?(m_len+1):0) ])
#define D_M2(k,i) (feature->buffer01[(i) + ((k)&1?(m_len+1):0) + 2*(m_len+1) ])
#define D_S2(i,j) (gp_substitution_matrix[(long int)(data1[pos1+(i)])][(long int)(data2[pos2+(j)])])


#ifdef INLINE
inline
#endif
long int left_alignment_SG_Diag (char * data1, long int pos1, long int len1,
                                 char * data2, long int pos2, long int len2,
                                 long int XDROP,
                                 /*out*/ long int * p_gainscore,
                                 /*out*/ long int * p_last_pos1,
                                 /*out*/ long int * p_last_pos2,
                                 Feature *feature)
{
  long int i,j;
  long int xdrop_current,xdrop_prev,xdrop_current2,xdrop_best=-INFINITY_INT;
  long int d_i_prev,d_i;
  long int m_i_prev,m_i;
  long int m_len = MIN(len1,len2);

  /* (1) initialisation */
  for (i = 0; i < m_len; i++) {
    D_I2(0,i) = -INFINITY_INT;
    D_I2(1,i) = -INFINITY_INT;
    D_M2(0,i) = -INFINITY_INT;
    D_M2(1,i) = -INFINITY_INT;
  }

  D_M2(0,0) = 0;
  *(p_gainscore) = 0;
  *(p_last_pos1) = pos1;
  *(p_last_pos2) = pos2;

  xdrop_current2 = 0;
  xdrop_current  = -INFINITY_INT;
  xdrop_best     = 0;

  /* (2) compute (stop if Xdrop) */
  for (j = 0; j < m_len; j++) {

    d_i = -INFINITY_INT;
    m_i = -INFINITY_INT;

    xdrop_prev    = xdrop_current;
    xdrop_current = -INFINITY_INT;

    for (i = 0; i <= j+1; i++) {

      d_i_prev = d_i;
      m_i_prev = m_i;

      if (i<j) {
        d_i      = D_I2(j+1,i);
        m_i      = D_M2(j+1,i);
      } else {
        d_i      = -INFINITY_INT;
        m_i      = -INFINITY_INT;
      }

      /* compute  d_i(j+1) */
      D_I2(j+1,i) = MAX(
                        MAX( ((i) ? D_M2(j,i-1) : -INFINITY_INT ), D_M2(j,i)) + gp_cost_gap_opened
                        ,
                        MAX( ((i) ? D_I2(j,i-1) : -INFINITY_INT ), D_I2(j,i)) + gp_cost_gap_continued
                        );

      /*  compute  m_i(j+1) */
      D_M2(j+1,i) = MAX (
                         m_i_prev  + D_S2(-i,(-j+i-1)),
                         d_i_prev  + D_S2(-i,(-j+i-1))
                         );


      xdrop_current = MAX(xdrop_current,D_M2(j+1,i));

    }/* i */
    if (j&1)
      xdrop_current2 = MAX(xdrop_current,xdrop_prev);
    if (xdrop_current2 <= -XDROP) {
      /* set the parameters beeing returned */
      return xdrop_current2;
    } else {
      if (xdrop_current > xdrop_best) {
        xdrop_best     = xdrop_current;
        *(p_gainscore) = xdrop_best;
        for (i = 0; i <= j+1; i++)
          if (D_M2(j+1,i) == xdrop_best) {
            *(p_last_pos1) = pos1 - i;
            *(p_last_pos2) = pos2 - j + i  - 1;
            break;
          }
      }
    }
  }
  return 0;
}


#ifdef INLINE
inline
#endif
long int right_alignment_SG_Diag(char * data1, long int pos1, long int len1,
                                 char * data2, long int pos2, long int len2,
                                 long int XDROP,
                                 /*out*/ long int * p_gainscore,
                                 /*out*/ long int * p_last_pos1,
                                 /*out*/ long int * p_last_pos2,
                                 Feature *feature)
{
  long int i,j;
  long int xdrop_current,xdrop_prev,xdrop_current2,xdrop_best=-INFINITY_INT;
  long int d_i_prev,d_i;
  long int m_i_prev,m_i;
  long int m_len = MIN(len1,len2);

  /* (1) initialisation */
  for (i = 0; i < m_len; i++) {
    D_I2(0,i) = -INFINITY_INT;
    D_I2(1,i) = -INFINITY_INT;
    D_M2(0,i) = -INFINITY_INT;
    D_M2(1,i) = -INFINITY_INT;
  }

  D_M2(0,0) = 0;
  *(p_gainscore) = 0;
  *(p_last_pos1) = pos1;
  *(p_last_pos2) = pos2;

  xdrop_current2 = 0;
  xdrop_current  = -INFINITY_INT;
  xdrop_best     = 0;

  /* (2) compute (stop if Xdrop) */
  for (j = 0; j < m_len; j++) {/* for each diagonal */

    d_i = -INFINITY_INT;
    m_i = -INFINITY_INT;

    xdrop_prev    = xdrop_current;
    xdrop_current = -INFINITY_INT;

    for (i = 0;i <= j+1; i++) {

      d_i_prev = d_i;
      m_i_prev = m_i;

      if (i<j) {
        d_i      = D_I2(j+1/*j-1*/,i);
        m_i      = D_M2(j+1/*j-1*/,i);
      } else {
        d_i      = -INFINITY_INT;
        m_i      = -INFINITY_INT;
      }

      /* compute d_i(j+1) */
      D_I2(j+1,i) = MAX(
                        MAX( ((i) ? D_M2(j,i-1) : -INFINITY_INT ), D_M2(j,i)) + gp_cost_gap_opened
                        ,
                        MAX( ((i) ? D_I2(j,i-1) : -INFINITY_INT ), D_I2(j,i)) + gp_cost_gap_continued
                        );

      /* compute m_i(j+1) */
      D_M2(j+1,i) = MAX (
                         m_i_prev  + D_S2(i,(j-i+1)),
                         d_i_prev  + D_S2(i,(j-i+1))
                         );


      xdrop_current = MAX(xdrop_current,D_M2(j+1,i));

    }/* i */
    if (j&1)
      xdrop_current2 = MAX(xdrop_current,xdrop_prev);
    if (xdrop_current2 <= -XDROP) {
      /* set the parameters being returned */
      return xdrop_current2;
    } else {
      if (xdrop_current > xdrop_best) {
        xdrop_best = xdrop_current;
        *(p_gainscore) = xdrop_best;
        for (i = 0; i <= j+1; i++)
          if (D_M2(j+1,i) == xdrop_best) {
            *(p_last_pos1) = pos1 + i  +  1;
            *(p_last_pos2) = pos2 + j - i + 2;
            break;
          }
      }
    }
  }
  return 0;
}

#ifdef INLINE
inline
#endif
long int left_alignment_SG_Lz(char * data1, long int pos1, long int len1,
                              char * data2, long int pos2, long int len2,
                              long int XDROP,
                              /*out*/ long int * p_gainscore,
                              /*out*/ long int * p_last_pos1,
                              /*out*/ long int * p_last_pos2,
                              Feature *feature) {
  long int len = MIN(len1,len2) - gp_delta_stat;
  long int last_indels = 0, indels = 0;
  long int last_i      = 0, i      = 0;
  long int gainscore = 0;

  while (i < len) {
    if (indels != 0)
      gainscore = *p_gainscore +  gp_cost_gap_opened + ABS(indels - last_indels) * gp_cost_gap_continued;

    for (i = last_i+1; i < len; i++) {
      gainscore += D_S2(-i,-i-indels);

      if (gainscore < -XDROP)
        break;

      /* [SFS] : single file stop */
      if (data1 + pos1 - i == data2 + pos2 - i - indels)
        break;

      if (gainscore > *p_gainscore) {

        *p_gainscore  = gainscore;
        *p_last_pos1  = pos1-i;
        *p_last_pos2  = pos2-i-indels;
        last_i        = i;
        last_indels   = indels;
      }
    }/* for i */

    if (indels > 0) {
      indels = -indels;
    } else {
      indels = -indels+1;
    }
    if (indels > gp_delta_stat)
      break;
  }/* while */

  return 1;
}



#ifdef INLINE
inline
#endif
long int right_alignment_SG_Lz(char * data1, long int pos1, long int len1,
                               char * data2, long int pos2, long int len2,
                               long int XDROP,
                               /*out*/ long int * p_gainscore,
                               /*out*/ long int * p_last_pos1,
                               /*out*/ long int * p_last_pos2,
                               Feature *feature) {
  long int len = MIN(len1,len2) - gp_delta_stat;
  long int last_indels = 0, indels = 0;
  long int last_i      = 0, i      = 0;
  long int gainscore = 0;

  while (i < len) {
    if (indels != 0)
      gainscore = *p_gainscore +  gp_cost_gap_opened + ABS(indels - last_indels) * gp_cost_gap_continued;

    for (i = last_i; i < len; i++) {
      gainscore += D_S2(i,i+indels);

      if (gainscore < -XDROP)
        break;

      /* [SFS] : single file stop */
      if (data1 + pos1 + i == data2 + pos2 + i + indels)
        break;

      if (gainscore > *p_gainscore) {

        *p_gainscore  = gainscore;
        *p_last_pos1  = pos1+i+1;
        *p_last_pos2  = pos2+i+indels+1;
        last_i        = i;
        last_indels   = indels;
      }
    }/* for i */

    if (indels > 0) {
      indels = -indels;
    } else {
      indels = -indels+1;
    }
    if (indels > gp_delta_stat)
      break;
  }/* while */

  return 1;
}



#ifdef INLINE
inline
#endif
long int left_alignment_SG_Border(char * data1, long int pos1, long int len1,
                                  char * data2, long int pos2, long int len2,
                                  long int XDROP,
                                  /*out*/ long int * p_gainscore,
                                  /*out*/ long int * p_last_pos1,
                                  /*out*/ long int * p_last_pos2,
                                  Feature *feature) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;
  long int row2       = 2*row1;

  /* dp tables */
  long int * bm0 = (long int *)feature->buffer01;
  long int * bm1 = (long int *)feature->buffer01+row1;
  long int * bi0 = (long int *)feature->buffer01+row2;
  long int * bi1 = (long int *)feature->buffer01+row2+row1;

  /* maximal and minimal positions taken into account */
  long int jmin = 0, jmax = row1-1;
  {
    long int j;
    for (j = 0; j < row1; j++) {
      bm1[j] = -INFINITY_INT;
      bi1[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
      bm0[j] = -INFINITY_INT;
      bi0[j] = -INFINITY_INT;
    }
    bm1[bandwidth] = 0;
    bi1[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 1; i <= len1; i++) {

      /* switch tables */
      {
        long int * bmt = bm0; long int * bit = bi0;
        bm0 = bm1; bi0 = bi1;
        bm1 = bmt; bi1 = bit;
      }

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(jmin,i+bandwidth-len1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i-bandwidth+j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i-bandwidth+j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 - i - bandwidth + j;
            *p_last_pos2 = pos2 - i;
          }
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(jmax,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 - i;
            *p_last_pos2 = pos2 - i + bandwidth - j;
          }
        }
      }

      /* MinMax evaluations */
      while (jmin < bandwidth
             && bm1[jmin] < -XDROP
             && bm0[jmin] < -XDROP)
        {
          bm0[jmin] = -INFINITY_INT;
          bm1[jmin] = -INFINITY_INT;
          bi0[jmin] = -INFINITY_INT;
          bi1[jmin] = -INFINITY_INT;
          jmin++;
        }
      while (jmax > bandwidth
             && bm1[jmax] < -XDROP
             && bm0[jmax] < -XDROP)
        {
          bm0[jmax] = -INFINITY_INT;
          bm1[jmax] = -INFINITY_INT;
          bi0[jmax] = -INFINITY_INT;
          bi1[jmax] = -INFINITY_INT;
          jmax--;
        }
      if (jmax <= jmin)
        return -XDROP;

      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>jmin; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmin < bandwidth)
        bi1[jmin] = MAX(bi1[jmin+1] + gp_cost_gap_continued, bm1[jmin+1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<jmax; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmax > bandwidth)
        bi1[jmax] = MAX(bi1[jmax-1] + gp_cost_gap_continued, bm1[jmax-1] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 1; i <= len2; i++) {

      /* switch tables */
      {
        long int * bmt = bm0; long int * bit = bi0;
        bm0 = bm1; bi0 = bi1;
        bm1 = bmt; bi1 = bit;
      }

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(jmin,i+bandwidth-len2); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 - i;
            *p_last_pos2 = pos2 - i - bandwidth + j;
          }
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(jmax,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i+bandwidth-j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i+bandwidth-j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 - i + bandwidth - j;
            *p_last_pos2 = pos2 - i;
          }
        }
      }

      /* MinMax evaluations */
      while (jmin < bandwidth
             && bm1[jmin] < -XDROP
             && bm0[jmin] < -XDROP)
        {
          bm0[jmin] = -INFINITY_INT;
          bm1[jmin] = -INFINITY_INT;
          bi0[jmin] = -INFINITY_INT;
          bi1[jmin] = -INFINITY_INT;
          jmin++;
        }
      while (jmax > bandwidth
             && bm1[jmax] < -XDROP
             && bm0[jmax] < -XDROP)
        {
          bm0[jmax] = -INFINITY_INT;
          bm1[jmax] = -INFINITY_INT;
          bi0[jmax] = -INFINITY_INT;
          bi1[jmax] = -INFINITY_INT;
          jmax--;
        }
      if (jmax <= jmin)
        return -XDROP;

      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>jmin; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmin < bandwidth)
        bi1[jmin] = MAX(bi1[jmin+1] + gp_cost_gap_continued, bm1[jmin+1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<jmax; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmax > bandwidth)
        bi1[jmax] = MAX(bi1[jmax-1] + gp_cost_gap_continued, bm1[jmax-1] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */
  return bm1[bandwidth];
}



#ifdef INLINE
inline
#endif
long int right_alignment_SG_Border(char * data1, long int pos1, long int len1,
                                   char * data2, long int pos2, long int len2,
                                   long int XDROP,
                                   /*out*/ long int * p_gainscore,
                                   /*out*/ long int * p_last_pos1,
                                   /*out*/ long int * p_last_pos2,
                                   Feature *feature) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;
  long int row2       = 2*row1;

  /* dp tables */
  long int * bm0 = (long int *)feature->buffer01;
  long int * bm1 = (long int *)feature->buffer01+row1;
  long int * bi0 = (long int *)feature->buffer01+row2;
  long int * bi1 = (long int *)feature->buffer01+row2+row1;

  /* maximal and minimal positions taken into account */
  long int jmin = 0, jmax = row1-1;
  {
    long int j;
    for (j = 0; j < row1; j++) {
      bm1[j] = -INFINITY_INT;
      bi1[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
      bm0[j] = -INFINITY_INT;
      bi0[j] = -INFINITY_INT;
    }
    bm1[bandwidth] = 0;
    bi1[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      /* switch tables */
      {
        long int * bmt = bm0; long int * bit = bi0;
        bm0 = bm1; bi0 = bi1;
        bm1 = bmt; bi1 = bit;
      }

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(jmin,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 + i + bandwidth - j + 1;
            *p_last_pos2 = pos2 + i + 1;
          }
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(jmax,len2+bandwidth-i-1); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 + i + 1;
            *p_last_pos2 = pos2 + i - bandwidth + j + 1;
          }
        }
      }

      /* MinMax evaluations */
      while (jmin < bandwidth
             && bm1[jmin] < -XDROP
             && bm0[jmin] < -XDROP)
        {
          bm0[jmin] = -INFINITY_INT;
          bm1[jmin] = -INFINITY_INT;
          bi0[jmin] = -INFINITY_INT;
          bi1[jmin] = -INFINITY_INT;
          jmin++;
        }
      while (jmax > bandwidth
             && bm1[jmax] < -XDROP
             && bm0[jmax] < -XDROP)
        {
          bm0[jmax] = -INFINITY_INT;
          bm1[jmax] = -INFINITY_INT;
          bi0[jmax] = -INFINITY_INT;
          bi1[jmax] = -INFINITY_INT;
          jmax--;
        }
      if (jmax <= jmin)
        return -XDROP;

      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>jmin; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmin < bandwidth)
        bi1[jmin] = MAX(bi1[jmin+1] + gp_cost_gap_continued, bm1[jmin+1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<jmax; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmax > bandwidth)
        bi1[jmax] = MAX(bi1[jmax-1] + gp_cost_gap_continued, bm1[jmax-1] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      /* switch tables */
      {
        long int * bmt = bm0; long int * bit = bi0;
        bm0 = bm1; bi0 = bi1;
        bm1 = bmt; bi1 = bit;
      }

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(jmin,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 + i + 1;
            *p_last_pos2 = pos2 + i + bandwidth - j + 1;
          }
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(jmax,len1+bandwidth-i-1); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
          if (bm1[j] > *p_gainscore) {
            *p_gainscore = bm1[j];
            *p_last_pos1 = pos1 + i - bandwidth + j + 1;
            *p_last_pos2 = pos2 + i + 1;
          }
        }
      }

      /* MinMax evaluations */
      while (jmin < bandwidth
             && bm1[jmin] < -XDROP
             && bm0[jmin] < -XDROP)
        {
          bm0[jmin] = -INFINITY_INT;
          bm1[jmin] = -INFINITY_INT;
          bi0[jmin] = -INFINITY_INT;
          bi1[jmin] = -INFINITY_INT;
          jmin++;
        }
      while (jmax > bandwidth
             && bm1[jmax] < -XDROP
             && bm0[jmax] < -XDROP)
        {
          bm0[jmax] = -INFINITY_INT;
          bm1[jmax] = -INFINITY_INT;
          bi0[jmax] = -INFINITY_INT;
          bi1[jmax] = -INFINITY_INT;
          jmax--;
        }
      if (jmax <= jmin)
        return -XDROP;

      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>jmin; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmin < bandwidth)
        bi1[jmin] = MAX(bi1[jmin+1] + gp_cost_gap_continued, bm1[jmin+1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<jmax; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmax > bandwidth)
        bi1[jmax] = MAX(bi1[jmax-1] + gp_cost_gap_continued, bm1[jmax-1] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */
  return bm1[bandwidth];
}





/* I.b) Alignment methods :
 */


#ifdef INLINE
inline
#endif
long int alignment_SG(char * data1, long int pos1, long int len1,
                      char * data2, long int pos2, long int len2,
                      Feature *feature)
{
  long int i,j;

  M2(0,0) = 0;
  I2(0,0) = gp_cost_gap_opened - gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M2(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I2(i,0) = I2(i-1,0) +  gp_cost_gap_continued;

  for (j = 1; j <= len2; j++) {

    M2(0,j) =  -INFINITY_INT;
    I2(0,j) =  I2(0,j-1) + gp_cost_gap_continued;

    for (i = 1; i <= len1;i++) {
      I2(i,j) = MAX (
                     MAX( M2(i,j-1) + gp_cost_gap_opened    , M2(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I2(i,j-1) + gp_cost_gap_continued ,  I2(i-1,j) + gp_cost_gap_continued )
                     );
      M2(i,j) = MAX ( M2(i-1,j-1) + S2(i,j), I2(i-1,j-1) + S2(i,j) );
    }
  }
  return MAX(M2(len1,len2),I2(len1,len2));
}


#ifdef INLINE
inline
#endif
long int alignment_SG_Strait(char * data1, long int pos1, char * data2, long int pos2, long int len,
                             Feature *feature)
{
  long int i,result = 0;
  for (i = 0; i < len; i++) {
    result += gp_substitution_matrix[(long int)(data2[pos2+i])][(long int)data1[pos1+i]];
#ifdef DEBUG_ALIGNMENT_STRAIT
    fprintf(stdout,"* score[%c/%c]=%ld\n",LOOKUP(data2[pos2+i]),LOOKUP(data1[pos1+i]),gp_substitution_matrix[(int)(data2[pos2+i])][(int)data1[pos1+i]]);
#endif
  }
  return result;
}


#ifdef INLINE
inline
#endif
long int alignment_SG_DROP(char * data1, long int pos1, long int len1,
                           char * data2, long int pos2, long int len2,
                           long int XDROP,
                           Feature *feature
                           )
{
  long int i,j;
  long int max;

  M2(0,0) = 0;
  I2(0,0) = gp_cost_gap_opened -  gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M2(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I2(i,0) = I2(i-1,0) +  gp_cost_gap_continued;

  for (j = 1; j <= len2; j++) {

    M2(0,j) =  -INFINITY_INT;
    I2(0,j) =  I2(0,j-1)  + gp_cost_gap_continued;
    max     =  M2(0,j);

    for (i = 1; i <= len1; i++) {
      I2(i,j) = MAX (
                     MAX( M2(i,j-1) + gp_cost_gap_opened    , M2(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I2(i,j-1) + gp_cost_gap_continued ,  I2(i-1,j) + gp_cost_gap_continued )
                     );
      M2(i,j) = MAX ( M2(i-1,j-1) + S2(i,j), I2(i-1,j-1) + S2(i,j) );
      max = MAX(max,M2(i,j));
    }
    if (max < -XDROP)
      return max;
  }
  return MAX(M2(len1,len2),I2(len1,len2));
}


#define M2opt_current(a)   (*(feature->buffer01 + d_current + (a)))
#define M2opt_previous(a)  (*(feature->buffer01 + d_previous + (a)))
#define I2opt_current(a)   (*(feature->buffer01 + I2_DELTA + d_current + (a)))
#define I2opt_previous(a)  (*(feature->buffer01 + I2_DELTA + d_previous + (a)))
#define S2opt(i,j) (gp_substitution_matrix[(long int)(data1[pos1+(i)-1])][(long int)(data2[pos2+(j)-1])])

#ifdef INLINE
inline
#endif
long int alignment_SG_DROP_opt(char * data1, long int pos1, long int len1,
                               char * data2, long int pos2, long int len2,
                               long int XDROP,
                               Feature *feature)
{
  long int i,j,s;
  long int max;
  long int I2_DELTA = (len1 + 1)<<1;
  long int DELTA    = (len1 + 1);
  long int d_previous = DELTA;
  long int d_current  = 0;

  M2opt_current(0) = 0;
  I2opt_current(0) = gp_cost_gap_opened -  gp_cost_gap_continued;

  for (i = 1; i <= len1; i++) {
    M2opt_current(i) = -INFINITY_INT;
    I2opt_current(i) = I2opt_current(i-1) +  gp_cost_gap_continued;
  }


  for (j = 1; j <= len2; j++) {
    if (j&1) {
      d_previous = 0;
      d_current  = DELTA;
    } else {
      d_previous = DELTA;
      d_current  = 0;
    }


    M2opt_current(0) =  -INFINITY_INT;
    I2opt_current(0) =  I2opt_previous(0)  + gp_cost_gap_continued;
    max              =  -INFINITY_INT;

    for (i = 1; i <= len1; i++) {
      d_previous ++;
      d_current  ++;
      I2opt_current(0) = MAX (
                              MAX( M2opt_previous(0) + gp_cost_gap_opened    , M2opt_current(-1) + gp_cost_gap_opened )
                              ,
                              MAX( I2opt_previous(0) + gp_cost_gap_continued , I2opt_current(-1) + gp_cost_gap_continued )
                              );
      s = S2opt(i,j);
      M2opt_current(0) = MAX ( M2opt_previous(-1) + s, I2opt_previous(-1) + s );
      max = MAX (max,M2opt_current(0));
    }
    if (max < -XDROP)
      return max;
  }
  return MAX(M2opt_current(0),I2opt_current(0));
}






#ifdef INLINE
inline
#endif
long int  alignment_SG_Lz(char * data1, long int pos1, long int len1,
                          char * data2, long int pos2, long int len2,
                          long int XDROP,
                          Feature *feature)
{

  long int last_score  = 0, score  = 0;
  long int last_indels = 0, indels = 0;
  long int last_i      = 0, i      = 0;

  if (len1 < len2) {

    for (i = 0; i < len1; i++) {
      score += D_S2(i,i+indels);
      if (score > last_score) {
        last_score  = score;
        last_i      = i;
        last_indels = indels;
      } else {
        if (score < -XDROP) {
          indels++;
          if (indels > len2 - len1)
            return -XDROP;
          score  = last_score + gp_cost_gap_opened + (indels-last_indels) * gp_cost_gap_continued;
          i      = last_i;
        }
      }
    }/* for i */
    return last_score;
  } else {
    for (i = 0; i < len2; i++) {
      score += D_S2(i+indels,i);
      if (score > last_score) {
        last_score  = score;
        last_i      = i;
        last_indels = indels;
      }
      if (score < -XDROP) {
        indels++;
        if (indels > len1 - len2)
          return -XDROP;
        score  = last_score + gp_cost_gap_opened + (indels-last_indels) * gp_cost_gap_continued;
        i      = last_i;
      }
    }/* for i */
    return last_score;
  }
}




/*
 * band alignment
 *
 *
 *          ------------------
 *          |^               |
 *          | bandwidth      |
 *          |v               |
 *          |----------------|
 *          |I dlen          |
 *          |----------------|
 *          |^               |
 *          | bandwidth      |
 *          |v               |
 *          ------------------
 *
 *
 *
 *  Classical DP matrix :
 *       (a) ->
 *     1  2  4  7  11 16 22
 * (b) 3  5  8  12 17 23
 *  |  6  9  13 18 24
 *  v  10 14 19 25
 *     15 20 26
 *     21 27
 *     28
 *
 *
 *
 *  Border DP matrix
 *
 *               (i)->
 *     22
 *     16----------        (a,b-1)
 *(j)  11 23                      \
 * |   7  17            (a-1,b-1)-[*]
 * V   4  12 24                    |                   (a,b-1)
 *     2  8  18                  (a-1,b)                    \
 *---- 1  5  13  25 ----------------------------- (a-1,b-1)-[*]
 *     3  9  19                  (a,b-1)                    /
 *     6  14 26                    |                   (a-1,b)
 *     10 20            (a-1,b-1)-[*]
 *     15 27                     /
 *     21----------        (a-1,b)
 *     28
 */

#ifdef INLINE
inline
#endif
long int alignment_SG_Border(char * data1, long int pos1, long int len1,
                             char * data2, long int pos2, long int len2,
                             long int XDROP,
                             Feature *feature) {

  /* band alignment length anw width*/
  long int dwidth       = ABS(len2-len1);
  long int bandwidth    = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + dwidth + 1;
  long int row2       = 2*row1;

  /* dp tables */
  long int * bm0 = (long int *)feature->buffer01;
  long int * bm1 = (long int *)feature->buffer01+row1;
  long int * bi0 = (long int *)feature->buffer01+row2;
  long int * bi1 = (long int *)feature->buffer01+row2+row1;

  /* maximal and minimal positions taken into account */
  long int jmin = 0, jmax = row1-1;
  {
    long int j;
    for (j = 0; j < row1; j++) {
      bm1[j] = -INFINITY_INT;
      bi1[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
      bm0[j] = -INFINITY_INT;
      bi0[j] = -INFINITY_INT;
    }
    bm1[bandwidth] = 0;
    bi1[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      /* switch tables */
      {
        long int * bmt = bm0; long int * bit = bi0;
        bm0 = bm1; bi0 = bi1;
        bm1 = bmt; bi1 = bit;
      }

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(jmin,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(jmax,len2+bandwidth-i-1); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }

      /* MinMax evaluations */
      while (jmin < bandwidth && bm1[jmin] < -XDROP && bm0[jmin] < -XDROP)
        {
          bm0[jmin] = -INFINITY_INT;
          bm1[jmin] = -INFINITY_INT;
          bi0[jmin] = -INFINITY_INT;
          bi1[jmin] = -INFINITY_INT;
          jmin++;
        }
      while (jmax > bandwidth && bm1[jmax] < -XDROP && bm0[jmax] < -XDROP)
        {
          bm0[jmax] = -INFINITY_INT;
          bm1[jmax] = -INFINITY_INT;
          bi0[jmax] = -INFINITY_INT;
          bi1[jmax] = -INFINITY_INT;
          jmax--;
        }
      if (jmax <= jmin)
        return -XDROP;

      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>jmin; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmin < bandwidth)
        bi1[jmin] = MAX(bi1[jmin+1] + gp_cost_gap_continued, bm1[jmin+1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<jmax; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmax > bandwidth)
        bi1[jmax] = MAX(bi1[jmax-1] + gp_cost_gap_continued, bm1[jmax-1] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      /* switch tables */
      {
        long int * bmt = bm0; long int * bit = bi0;
        bm0 = bm1; bi0 = bi1;
        bm1 = bmt; bi1 = bit;
      }

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(jmin,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(jmax,len1+bandwidth-i-1); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }

      /* MinMax evaluations */
      while (jmin < bandwidth && bm1[jmin] < -XDROP && bm0[jmin] < -XDROP)
        {
          bm0[jmin] = -INFINITY_INT;
          bm1[jmin] = -INFINITY_INT;
          bi0[jmin] = -INFINITY_INT;
          bi1[jmin] = -INFINITY_INT;
          jmin++;
        }
      while (jmax > bandwidth && bm1[jmax] < -XDROP && bm0[jmax] < -XDROP)
        {
          bm0[jmax] = -INFINITY_INT;
          bm1[jmax] = -INFINITY_INT;
          bi0[jmax] = -INFINITY_INT;
          bi1[jmax] = -INFINITY_INT;
          jmax--;
        }
      if (jmax <= jmin)
        return -XDROP;

      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>jmin; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmin < bandwidth)
        bi1[jmin] = MAX(bi1[jmin+1] + gp_cost_gap_continued, bm1[jmin+1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<jmax; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      if (jmax > bandwidth)
        bi1[jmax] = MAX(bi1[jmax-1] + gp_cost_gap_continued, bm1[jmax-1] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */
  return MAX(bm1[bandwidth + dwidth],bi1[bandwidth + dwidth]);
}



#ifdef INLINE
inline
#endif
long int alignment_SG_on_MA(char * data1,  char * data2,
                            MA * ma, long int left_correction,
                            Feature *feature)
{
  long int i,result = 0;
  long int left_size,right_size;
  tuple * t = ma->first_tuple, * t_prev = ma->first_tuple;

  /* [1] Compute score before first tuple */
  left_size  = TBL_POS(t) - ma->left_pos_begin;
  right_size = TBR_POS(t) - ma->right_pos_begin;

  if (left_size > 0 && right_size > 0)
    result += alignment_SG(data1, ma->left_pos_begin,  left_size,
                           data2, ma->right_pos_begin, right_size,
                           feature
                           );

  /* [2] Compute score for several tuples */
  while (t != NULL) {
    /* (1) compute score inside tuple */
    for (i = 0; i < TSIZE(t); i++)
      result += gp_substitution_matrix[(long int)(data1[i+TBL_POS(t)])][(long int)(data2[i+TBR_POS(t)])];

    if (t->next) {
      /* (2) compute score between tuple (t) and (t->next) */
      result += alignment_SG(data1, TEL_POS(t), TGAP_L(t,t->next),
                             data2, TER_POS(t), TGAP_R(t,t->next),
                             feature
                             );
    }/* if */
    t_prev = t;
    t = t->next;
  }/* while */

  /* [3] Compute score after last tuple */
  left_size  =  ma->left_pos_end  - TEL_POS(t_prev);
  right_size =  ma->right_pos_end - TER_POS(t_prev);
  if (left_size > 0 && right_size > 0)
    result += alignment_SG(data1, TEL_POS(t_prev), left_size,
                           data2, TER_POS(t_prev), right_size,
                           feature
                           );

  return result;
}



/*
 * II) Statistical functions:
 *  - (compute transition/transversion ratio)
 *  - (compute triplet mutation frequency)
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */



typedef enum backtrack_SG_t { _BSG_EMPTY_,
                              M_TO_M_SUBS,
                              I_TO_M_SUBS,
                              FIRST_DEL_LEFT,
                              FIRST_DEL_RIGHT,
                              CONT_DEL_LEFT,
                              CONT_DEL_RIGHT } backtrack_SG;


#define M3(i,j) (mt[(i)+ (j)*(len1+1)])
#define I3(i,j) (in[(i)+ (j)*(len1+1)])
#define S3(i,j) S2(i,j)
#define D3(k)   (bk[(k)])


#define WORDCOUNTANDSTAT(count,letter,word) {(word) <<= 2; (word) |= backcode[(long int)(letter)]; (word) &= 0xfc00003f;  if (((word) & 0x3f) == (word)) count[(word)]++;}
#define PAIROFWORDCOUNANDSTAT(count,word1,word2) { if ((word1) > 0 && (word1) < 64 && (word2) > 0 && (word2) < 64) ((count)[word1][word2])++;}



/* II.a) Extensions ("left" or "right") methods :
 */


long int left_alignment_SG_stats_Border(char * data1, long int pos1, long int len1,
                                        char * data2, long int pos2, long int len2,
                                        unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                        MA   * ma,
                                        Feature * feature) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;

 /* various switch elements to access information */
  long int width       = 2*bandwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < row1; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 1; i <= len1; i++) {

      long int * bm0 = mt[i-1];
      long int * bm1 = mt[i];
      long int * bi0 = in[i-1];
      long int * bi1 = in[i];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i-bandwidth+j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i-bandwidth+j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(row1-1,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 1; i <= len2; i++) {

      long int * bm0 = mt[i-1];
      long int * bm1 = mt[i];
      long int * bi0 = in[i-1];
      long int * bi1 = in[i];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(row1-1,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i+bandwidth-j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i+bandwidth-j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  /* [3] display alignment */

  if (len1  < len2) {

    long int i = len1, j = len2 - len1;
    long int i1 = -len1, i2 = -len2;
    backtrack_SG bk;

    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk = I_TO_M_SUBS;
    else
      bk = M_TO_M_SUBS;

    if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
      (ma->trinomial_count1[(pos1+i1)%3])++;
      (ma->trinomial_count2[(pos2+i2)%3])++;
      *nonmutatedword = 0xfc000000;
      if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
        (ma->transindels[0])++; /* transition */
      else
        (ma->transindels[1])++; /* transversion */
    } else {
      WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
    }
    WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
    WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
    PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

    i1++;
    i2++;
    i--;

    while (i > 0) {

      if (bk == M_TO_M_SUBS || bk == FIRST_DEL_LEFT || bk == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-j])]:
             gp_substitution_matrix[(long int)(data1[pos1-i+j])][(long int)(data2[pos2-i])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk = M_TO_M_SUBS;

          if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
            (ma->trinomial_count1[(pos1+i1)%3])++;
            (ma->trinomial_count2[(pos2+i2)%3])++;
            *nonmutatedword = 0xfc000000;
            if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
              (ma->transindels[0])++; /* transition */
            else
              (ma->transindels[1])++; /* transversion */
          } else {
            WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
          }
          WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
          WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
          PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

          i1++;
          i2++;

        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-j])]:
               gp_substitution_matrix[(long int)(data1[pos1-i+j])][(long int)(data2[pos2-i])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk = I_TO_M_SUBS;

            if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
              (ma->trinomial_count1[(pos1+i1)%3])++;
              (ma->trinomial_count2[(pos2+i2)%3])++;
              *nonmutatedword = 0xfc000000;
              if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
                (ma->transindels[0])++; /* transition */
              else
                (ma->transindels[1])++; /* transversion */
            } else {
              WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
            }
            WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
            WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
            PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

            i1++;
            i2++;

          } else {
            _WARNING("left_alignment_SG_stats_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk =  FIRST_DEL_RIGHT;

          ma->transindels[6]++; /* block of indels text */
          *nonmutatedword = 0xfc000000;
          ma->transindels[2]++; /* indels */
          ma->transindels[4]++; /* indels text */

          i2++;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk =  FIRST_DEL_LEFT;

            ma->transindels[5]++; /* block of indels query */
            *nonmutatedword = 0xfc000000;
            ma->transindels[2]++; /* indels */
            ma->transindels[3]++; /* indels query */

            i1++;

          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk =  CONT_DEL_RIGHT;

              *nonmutatedword = 0xfc000000;
              ma->transindels[2]++; /* indels */
              ma->transindels[4]++; /* indels text */

            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk =  CONT_DEL_LEFT;

                *nonmutatedword = 0xfc000000;
                ma->transindels[2]++; /* indels */
                ma->transindels[3]++; /* indels query */

                i1++;

              } else {
                _WARNING("left_alignment_SG_stats_Border() 2");
              }
            }
          }
        }
      }
    }/* while */


    if (j > 0) {

      while (j > 1) {
        j--;
        bk = CONT_DEL_RIGHT;

        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[4]++; /* indels text */

        i2++;
      }

      bk = FIRST_DEL_RIGHT;

      ma->transindels[6]++; /* block of indels text */
      *nonmutatedword = 0xfc000000;
      ma->transindels[2]++; /* indels */
      ma->transindels[4]++; /* indels text */

      i2++;
    }

    if (j < 0) {

      while (j < -1) {
        j++;
        bk = CONT_DEL_LEFT;

        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[3]++; /* indels query */

        i1++;
      }

      bk = FIRST_DEL_LEFT;

      ma->transindels[5]++; /* block of indels query */
      *nonmutatedword = 0xfc000000;
      ma->transindels[2]++; /* indels */
      ma->transindels[3]++; /* indels query */
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = len1 - len2;
    long int i1 = -len1, i2 = -len2;
    backtrack_SG bk;

    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk = I_TO_M_SUBS;
    else
      bk = M_TO_M_SUBS;

    if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
      (ma->trinomial_count1[(pos1+i1)%3])++;
      (ma->trinomial_count2[(pos2+i2)%3])++;
      *nonmutatedword = 0xfc000000;
      if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
        (ma->transindels[0])++; /* transition */
      else
        (ma->transindels[1])++; /* transversion */
    } else {
      WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
    }
    WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
    WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
    PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

    i1++;
    i2++;
    i--;

    while (i > 0) {

      if (bk == M_TO_M_SUBS || bk == FIRST_DEL_LEFT || bk == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1-i-j])][(long int)(data2[pos2-i])]:
             gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+j])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk = M_TO_M_SUBS;

          if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
            (ma->trinomial_count1[(pos1+i1)%3])++;
            (ma->trinomial_count2[(pos2+i2)%3])++;
            *nonmutatedword = 0xfc000000;
            if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
              (ma->transindels[0])++; /* transition */
            else
              (ma->transindels[1])++; /* transversion */
          } else {
            WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
          }
          WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
          WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
          PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

          i1++;
          i2++;

        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1-i-j])][(long int)(data2[pos2-i])]:
               gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+j])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk = I_TO_M_SUBS;

            if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
              (ma->trinomial_count1[(pos1+i1)%3])++;
              (ma->trinomial_count2[(pos2+i2)%3])++;
              *nonmutatedword = 0xfc000000;
              if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
                (ma->transindels[0])++; /* transition */
              else
                (ma->transindels[1])++; /* transversion */
            } else {
              WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
            }
            WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
            WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
            PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

            i1++;
            i2++;

          } else {
            _WARNING("left_alignment_SG_stats_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk =  FIRST_DEL_LEFT;

          ma->transindels[5]++; /* block of indels query */
          *nonmutatedword = 0xfc000000;
          ma->transindels[2]++; /* indels */
          ma->transindels[3]++; /* indels query */

          i1++;

        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk =  FIRST_DEL_RIGHT;

            ma->transindels[6]++; /* block of indels text */
            *nonmutatedword = 0xfc000000;
            ma->transindels[2]++; /* indels */
            ma->transindels[4]++; /* indels text */

            i2++;

          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk =  CONT_DEL_LEFT;

              *nonmutatedword = 0xfc000000;
              ma->transindels[2]++; /* indels */
              ma->transindels[3]++; /* indels query */

              i1++;

            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk =  CONT_DEL_RIGHT;

                *nonmutatedword = 0xfc000000;
                ma->transindels[2]++; /* indels */
                ma->transindels[4]++; /* indels text */

                i2++;

              } else {
                _WARNING("left_alignment_SG_stats_Border() 2");
              }
            }
          }
        }
      }
    }/* while */


    if (j < 0) {

      while (j < -1) {
        j++;
        bk = CONT_DEL_RIGHT;

        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[3]++; /* indels query */

        i2++;
      }

      bk = FIRST_DEL_RIGHT;

      ma->transindels[6]++; /* block of indels text */
      *nonmutatedword = 0xfc000000;
      ma->transindels[2]++; /* indels */
      ma->transindels[4]++; /* indels text */

      i2++;
    }

    if (j > 0) {

      while (j > 1) {
        j--;
        bk = CONT_DEL_LEFT;

        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[3]++; /* indels query */

        i1++;
      }

      bk = FIRST_DEL_LEFT;

      ma->transindels[5]++; /* block of indels query */
      *nonmutatedword = 0xfc000000;
      ma->transindels[2]++; /* indels */
      ma->transindels[3]++; /* indels query */

      i1++;
    }
  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);

  return 0;
}



long int right_alignment_SG_stats_Border(char * data1, long int pos1, long int len1,
                                         char * data2, long int pos2, long int len2,
                                         unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                         MA   * ma,
                                         Feature * feature) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);
  long int k = 0;
  backtrack_SG *  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,right_alignment_SG_stats_Border);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < row1; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(row1,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(row1,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  if (len1  < len2) {

    long int i = len1, j = len2 - len1;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("right_alignment_SG_stats_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_RIGHT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_LEFT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_RIGHT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_LEFT;
              } else {
                _WARNING("right_alignment_SG_stats_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = len1 - len2;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("right_alignment_SG_stats_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_LEFT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_RIGHT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_LEFT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_RIGHT;
              } else {
                _WARNING("right_alignment_SG_stats_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }
  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);




  /* [3] display alignment */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {

      switch(D3(k--)) {
      case M_TO_M_SUBS:
      case I_TO_M_SUBS:

        if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
          (ma->trinomial_count1[(pos1+i1)%3])++;
          (ma->trinomial_count2[(pos2+i2)%3])++;
          *nonmutatedword = 0xfc000000;
          if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
            (ma->transindels[0])++; /* transition */
          else
            (ma->transindels[1])++; /* transversion */
        } else {
          WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
        }
        WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
        WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
        PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
        ma->transindels[5]++; /* block of indels query */
      case CONT_DEL_LEFT:
        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[3]++; /* indels query */
        i1++;
        break;

      case FIRST_DEL_RIGHT:
        ma->transindels[6]++; /* block of indels text */
      case CONT_DEL_RIGHT:
        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[4]++; /* indels text */
        i2++;
        break;

      default:
        _WARNING("right_alignment_SG_stats_Border() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return 0;
}





/* II.b) Alignment methods :
 */

long int alignment_SG_stats_Border(char * data1, long int pos1, long int len1,
                                   char * data2, long int pos2, long int len2,
                                   unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                   MA   * ma,
                                   Feature * feature) {

  /* band alignment length anw width*/
  long int dwidth       = ABS(len2-len1);
  long int bandwidth    = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + dwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);
  long int k = 0;
  backtrack_SG *  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,alignment_SG_stats_Border);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < width; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(width,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<width-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[width-1] = MAX(bi1[width-2] + gp_cost_gap_continued, bm1[width-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(width,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<width-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[width-1] = MAX(bi1[width-2] + gp_cost_gap_continued, bm1[width-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  if (len1  < len2) {

    long int i = len1, j = dwidth;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("alignment_SG_stats_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_RIGHT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_LEFT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_RIGHT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_LEFT;
              } else {
                _WARNING("alignment_SG_stats_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = dwidth;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("alignment_SG_stats_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_LEFT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_RIGHT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_LEFT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_RIGHT;
              } else {
                _WARNING("alignment_SG_stats_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);




  /* [3] display alignment */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {

      switch(D3(k--)) {
      case M_TO_M_SUBS:
      case I_TO_M_SUBS:

        if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
          (ma->trinomial_count1[(pos1+i1)%3])++;
          (ma->trinomial_count2[(pos2+i2)%3])++;
          *nonmutatedword = 0xfc000000;
          if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
            (ma->transindels[0])++; /* transition */
          else
            (ma->transindels[1])++; /* transversion */
        } else {
          WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
        }
        WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
        WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
        PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
        ma->transindels[5]++; /* block of indels query */
      case CONT_DEL_LEFT:
        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[3]++; /* indels query */
        i1++;
        break;

      case FIRST_DEL_RIGHT:
        ma->transindels[6]++; /* block of indels text */
      case CONT_DEL_RIGHT:
        *nonmutatedword = 0xfc000000;
        ma->transindels[2]++; /* indels */
        ma->transindels[4]++; /* indels text */
        i2++;
        break;

      default:
        _WARNING("alignment_SG_stats_Border() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return 0;
}


long int alignment_SG_stats(char * data1, long int pos1, long int len1,
                            char * data2, long int pos2, long int len2,
                            unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                            MA   * ma,
                            Feature * feature) {

  long int i,j,k;
  long int i1,i2;
  long int result;
  long int *mt, *in;
  backtrack_SG * bk;

  mt = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(mt,alignment_SG_stats);
  in = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(in,alignment_SG_stats);
  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,alignment_SG_stats);

  /* [1] compute */

  M3(0,0) = 0;

  I3(0,0) = gp_cost_gap_opened - gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M3(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I3(i,0) = I3(i-1,0) + gp_cost_gap_continued;


  for (j = 1; j <= len2; j++) {

    M3(0,j) =  -INFINITY_INT;
    I3(0,j) =  gp_cost_gap_opened  + (j-1)*gp_cost_gap_continued;

    for (i = 1; i <= len1; i++) {
      I3(i,j) = MAX (
                     MAX( M3(i,j-1) + gp_cost_gap_opened    , M3(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I3(i,j-1) + gp_cost_gap_continued , I3(i-1,j) + gp_cost_gap_continued )
                     );
      M3(i,j) = MAX ( M3(i-1,j-1) + S3(i,j), I3(i-1,j-1) + S3(i,j) );
    }
  }

  result = MAX( M3(len1,len2) , I3(len1,len2) );

  /* [2] backtrace */

  i = len1;
  j = len2;
  k = 0;

  if (M3(len1,len2) < I3(len1,len2))
    D3(k++) = I_TO_M_SUBS;
  else
    D3(k++) = M_TO_M_SUBS;


  while (i > 0 && j > 0) {

    if (D3(k-1) == M_TO_M_SUBS || D3(k-1) == FIRST_DEL_LEFT || D3(k-1) == FIRST_DEL_RIGHT) {

      /* no gap to open  */
      if (M3(i,j) == M3(i-1,j-1) + S3(i,j)) {
        i--;
        j--;
        D3(k++) = M_TO_M_SUBS;
      } else {
        if (M3(i,j) == I3(i-1,j-1) + S3(i,j)) {
          i--;
          j--;
          D3(k++) = I_TO_M_SUBS;
        } else {
          _WARNING("alignment_SG_stats() 1");
        }
      }
    } else {
      /* need to open a gap */
      if (I3(i,j) == M3(i,j-1) + gp_cost_gap_opened) {
        j--;
        D3(k++) =  FIRST_DEL_RIGHT;
      } else {
        if (I3(i,j) == M3(i-1,j) + gp_cost_gap_opened) {
          i--;
          D3(k++) =  FIRST_DEL_LEFT;
        } else {
          if (I3(i,j) == I3(i,j-1) + gp_cost_gap_continued) {
            j--;
            D3(k++) =  CONT_DEL_RIGHT;
          } else {
            if (I3(i,j) == I3(i-1,j) + gp_cost_gap_continued) {
              i--;
              D3(k++) =  CONT_DEL_LEFT;
            } else {
              _WARNING("alignment_SG_stats() 2");
            }
          }
        }
      }
    }
  }/* while */


  if (j > 0) {
    while (j > 1) {
      j--;
      D3(k++) = CONT_DEL_RIGHT;
    }
    D3(k++) = FIRST_DEL_RIGHT;
  }

  if (i > 0) {
    while (i > 1) {
      i--;
      D3(k++) = CONT_DEL_LEFT;
    }
    D3(k++) = FIRST_DEL_LEFT;
  }

  FREE(mt,(len1+1) * (len2+1) * sizeof(long int));
  FREE(in,(len1+1) * (len2+1) * sizeof(long int));


  /* [3] backtrace to compute some statistics */

  i1 = 0; i2 = 0; k--;

  while (i1 < len1 || i2 < len2) {

    switch(D3(k--)) {
    case M_TO_M_SUBS:
    case I_TO_M_SUBS:

      if (TNORM(data1[pos1+i1]) != TNORM(data2[pos2+i2])) { /* substitution */
        (ma->trinomial_count1[(pos1+i1)%3])++;
        (ma->trinomial_count2[(pos2+i2)%3])++;
        *nonmutatedword = 0xfc000000;
        if ((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2])) == 2)
          (ma->transindels[0])++; /* transition */
        else
          (ma->transindels[1])++; /* transversion */
      } else {
        WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i1],*nonmutatedword);
      }
      WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i1],*mutatedword1);
      WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i2],*mutatedword2);
      PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);

      i1++;
      i2++;
      break;

    case FIRST_DEL_LEFT:
      ma->transindels[5]++; /* block of indels query */
    case CONT_DEL_LEFT:
      *nonmutatedword = 0xfc000000;
      ma->transindels[2]++; /* indels */
      ma->transindels[3]++; /* indels query */
      i1++;
      break;

    case FIRST_DEL_RIGHT:
      ma->transindels[6]++; /* block of indels text */
    case CONT_DEL_RIGHT:
      *nonmutatedword = 0xfc000000;
      ma->transindels[2]++; /* indels */
      ma->transindels[4]++; /* indels text */
      i2++;
      break;

    default:
      _WARNING("alignment_SG_stats() 3");
    }/* case */
  }/* while */
  FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  return result;
}






long int alignment_SG_stats_Strait(char * data1, long int pos1,
                                   char * data2, long int pos2, long int len,
                                   unsigned * mutatedword1, unsigned * mutatedword2, unsigned * nonmutatedword,
                                   MA * ma,
                                   Feature * feature) {

  long int i;
  for (i = 0; i < len; i++) {
    if (TNORM(data1[pos1+i]) != TNORM(data2[pos2+i])) { /* substitution */
      (ma->trinomial_count1[(pos1+i)%3])++;
      (ma->trinomial_count2[(pos2+i)%3])++;
      *nonmutatedword = 0xfc000000;
      if ((TNORM(data1[pos1+i])^TNORM(data2[pos2+i])) == 2)
        (ma->transindels[0])++; /* transition */
      else
        (ma->transindels[1])++; /* transverstion */
    } else {
      WORDCOUNTANDSTAT(feature->nb_non_mutated_triplet_count,data1[pos1+i],*nonmutatedword);
    }
    WORDCOUNTANDSTAT(feature->nb_triplet_count1,data1[pos1+i],*mutatedword1);
    WORDCOUNTANDSTAT(feature->nb_triplet_count2,data2[pos2+i],*mutatedword2);
    PAIROFWORDCOUNANDSTAT(feature->nb_pair_of_triplets,*mutatedword1,*mutatedword2);
  }
  return 0;
}



long int alignment_SG_stats_on_MA(char * data1, char * data2,
                                  long int left_correction,
                                  MA * ma,
                                  Feature * feature,
                                  long int full) {


  tuple * t = ma->first_tuple, * t_prev = ma->first_tuple;
  unsigned mutatedword1    = 0xfc000000;
  unsigned mutatedword2    = 0xfc000000;
  unsigned nonmutatedword  = 0xfc000000;

  /* [1] Compute alignment before first tuple */
  if (full) {
    long int left_size  = TBL_POS(t) - ma->left_pos_begin;
    long int right_size = TBR_POS(t) - ma->right_pos_begin;
    if (left_size > 0 && right_size > 0)
      left_alignment_SG_stats_Border(data1,  TBL_POS(t), left_size,
                                     data2,  TBR_POS(t), right_size,
                                     &mutatedword1,&mutatedword2,&nonmutatedword,
                                     ma, feature);

#ifdef DEBUG_STATS
    fprintf(stderr,"Stat MA=(%ld,%ld)-(%ld,%ld) -> trinomial1=(%ld,%ld,%ld) tr=%ld tv=%ld id=%ld trinomial2=(%ld,%ld,%ld)\n",
            ma->left_pos_begin,ma->left_pos_end,
            ma->right_pos_begin,ma->right_pos_end,
            ma->trinomial_count1[0],ma->trinomial_count1[1],ma->trinomial_count1[2],
            ma->transindels[0],ma->transindels[1],ma->transindels[2],
            ma->trinomial_count2[0],ma->trinomial_count2[1],ma->trinomial_count2[2]
            );
#endif
  }


  /* [2] Compute alignment for several tuples */
  while (t != NULL) {
    /* (1) compute inside a tuple */
    alignment_SG_stats_Strait(data1,TBL_POS(t),
                              data2,TBR_POS(t),TSIZE(t),
                              &mutatedword1,&mutatedword2,&nonmutatedword,
                              ma, feature);

#ifdef DEBUG_STATS
    fprintf(stderr,"Stat MA=(%ld,%ld)-(%ld,%ld) -> trinomial1=(%ld,%ld,%ld) tr=%ld tv=%ld id=%ld trinomial2=(%ld,%ld,%ld)\n",
            ma->left_pos_begin,ma->left_pos_end,
            ma->right_pos_begin,ma->right_pos_end,
            ma->trinomial_count1[0],ma->trinomial_count1[1],ma->trinomial_count1[2],
            ma->transindels[0],ma->transindels[1],ma->transindels[2],
            ma->trinomial_count2[0],ma->trinomial_count2[1],ma->trinomial_count2[2]
            );
#endif
    if (full) {
      if (t->next) {
        /* (2) alignment between (t) et (t->next) */
        alignment_SG_stats_Border(data1, TEL_POS(t), TGAP_L(t,t->next),
                                  data2, TER_POS(t), TGAP_R(t,t->next),
                                  &mutatedword1, &mutatedword2, &nonmutatedword,
                                  ma, feature);


#ifdef DEBUG_STATS
        fprintf(stderr,"Stat MA=(%ld,%ld)-(%ld,%ld) -> trinomial1=(%ld,%ld,%ld) tr=%ld tv=%ld id=%ld trinomial2=(%ld,%ld,%ld)\n",
                ma->left_pos_begin,ma->left_pos_end,
                ma->right_pos_begin,ma->right_pos_end,
                ma->trinomial_count1[0],ma->trinomial_count1[1],ma->trinomial_count1[2],
                ma->transindels[0],ma->transindels[1],ma->transindels[2],
                ma->trinomial_count2[0],ma->trinomial_count2[1],ma->trinomial_count2[2]
                );
#endif
      }
    } /* full */
    t_prev = t;
    t = t->next;
  }/* while */


  /* [3] Compute alignement after last tuple */
  if (full) {
    long int left_size  =  ma->left_pos_end  - TEL_POS(t_prev);
    long int right_size =  ma->right_pos_end - TER_POS(t_prev);
    if (left_size > 0 && right_size > 0)
      right_alignment_SG_stats_Border(data1, TEL_POS(t_prev), left_size,
                                      data2, TER_POS(t_prev), right_size,
                                      &mutatedword1,&mutatedword2,&nonmutatedword,
                                      ma, feature);
#ifdef DEBUG_STATS
    fprintf(stderr,"Stat MA=(%ld,%ld)-(%ld,%ld) -> trinomial1=(%ld,%ld,%ld) tr=%ld tv=%ld id=%ld trinomial2=(%ld,%ld,%ld)\n",
            ma->left_pos_begin,ma->left_pos_end,
            ma->right_pos_begin,ma->right_pos_end,
            ma->trinomial_count1[0],ma->trinomial_count1[1],ma->trinomial_count1[2],
            ma->transindels[0],ma->transindels[1],ma->transindels[2],
            ma->trinomial_count2[0],ma->trinomial_count2[1],ma->trinomial_count2[2]
            );
#endif
  }
  return 0;
}




/*-----------------------------------------------------------------*/
/*
 * III) PSL output functions:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */

#define CHOICE_LEFT_POS_BEGIN_LIST  1
#define CHOICE_RIGHT_POS_BEGIN_LIST 2
#define CHOICE_BLOCK_LENGTH_LIST    3
#define CHOICE_BLOCK_NB             4


/* III.a) Extensions ("left" or "right") methods :
 */


long int left_alignment_SG_PSL_Border(char * data1, long int pos1, long int len1,
                                      char * data2, long int pos2, long int len2,
                                      /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                      /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                      long int choice
                                      ) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;

 /* various switch elements to access information */
  long int width       = 2*bandwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < row1; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 1; i <= len1; i++) {

      long int * bm0 = mt[i-1];
      long int * bm1 = mt[i];
      long int * bi0 = in[i-1];
      long int * bi1 = in[i];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i-bandwidth+j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i-bandwidth+j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(row1-1,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 1; i <= len2; i++) {

      long int * bm0 = mt[i-1];
      long int * bm1 = mt[i];
      long int * bi0 = in[i-1];
      long int * bi1 = in[i];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(row1-1,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i+bandwidth-j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i+bandwidth-j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  /* [3] display alignment */

  if (len1  < len2) {

    long int i = len1, j = len2 - len1;
    long int i1 = -len1, i2 = -len2;
    backtrack_SG bk;

    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk = I_TO_M_SUBS;
    else
      bk = M_TO_M_SUBS;
    (*pPSL_length)++;
    i1++;
    i2++;
    i--;

    while (i > 0) {

      if (bk == M_TO_M_SUBS || bk == FIRST_DEL_LEFT || bk == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-j])]:
             gp_substitution_matrix[(long int)(data1[pos1-i+j])][(long int)(data2[pos2-i])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk = M_TO_M_SUBS;
          (*pPSL_length)++;
          i1++;
          i2++;

        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-j])]:
               gp_substitution_matrix[(long int)(data1[pos1-i+j])][(long int)(data2[pos2-i])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk = I_TO_M_SUBS;
            (*pPSL_length)++;
            i1++;
            i2++;

          } else {
            _WARNING("left_alignment_SG_PSL_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk =  FIRST_DEL_RIGHT;
          switch (choice) {
          case CHOICE_LEFT_POS_BEGIN_LIST:
            fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
            break;
          case CHOICE_RIGHT_POS_BEGIN_LIST:
            fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
            break;
          case CHOICE_BLOCK_LENGTH_LIST:
            fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
            break;
          }
          i2++;

          (*pPSL_left_pos_begin)  =  pos1 + i1;
          (*pPSL_right_pos_begin) =  pos2 + i2;
          (*pPSL_nb)++;
          (*pPSL_length) = 0;


        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk =  FIRST_DEL_LEFT;
            switch (choice) {
            case CHOICE_LEFT_POS_BEGIN_LIST:
              fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
              break;
            case CHOICE_RIGHT_POS_BEGIN_LIST:
              fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
              break;
            case CHOICE_BLOCK_LENGTH_LIST:
              fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
              break;
            }
            i1++;

            (*pPSL_left_pos_begin)  =  pos1 + i1;
            (*pPSL_right_pos_begin) =  pos2 + i2;
            (*pPSL_nb)++;
            (*pPSL_length) = 0;


          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk =  CONT_DEL_RIGHT;
              i2++;

            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk =  CONT_DEL_LEFT;
                i1++;

              } else {
                _WARNING("left_alignment_SG_PSL_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {

      while (j > 1) {
        j--;
        bk = CONT_DEL_RIGHT;
        i2++;
      }

      bk = FIRST_DEL_RIGHT;
      switch (choice) {
      case CHOICE_LEFT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
        break;
      case CHOICE_RIGHT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
        break;
      case CHOICE_BLOCK_LENGTH_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
        break;
      }
      i2++;

      (*pPSL_left_pos_begin)  =  pos1 + i1;
      (*pPSL_right_pos_begin) =  pos2 + i2;
      (*pPSL_nb)++;
      (*pPSL_length) = 0;

    }

    if (j < 0) {

      while (j < -1) {
        j++;
        bk = CONT_DEL_LEFT;
        i1++;
      }

      bk = FIRST_DEL_LEFT;
      switch (choice) {
      case CHOICE_LEFT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
        break;
      case CHOICE_RIGHT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
        break;
      case CHOICE_BLOCK_LENGTH_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
        break;
      }
      i1++;

      (*pPSL_left_pos_begin)  =  pos1 + i1;
      (*pPSL_right_pos_begin) =  pos2 + i2;
      (*pPSL_nb)++;
      (*pPSL_length) = 0;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = len1 - len2;
    long int i1 = -len1, i2 = -len2;
    backtrack_SG bk;

    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk = I_TO_M_SUBS;
    else
      bk = M_TO_M_SUBS;
    (*pPSL_length)++;
    i1++;
    i2++;
    i--;

    while (i > 0) {

      if (bk == M_TO_M_SUBS || bk == FIRST_DEL_LEFT || bk == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1-i-j])][(long int)(data2[pos2-i])]:
             gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+j])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk = M_TO_M_SUBS;
          (*pPSL_length)++;
          i1++;
          i2++;

        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1-i-j])][(long int)(data2[pos2-i])]:
               gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+j])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk = I_TO_M_SUBS;
            (*pPSL_length)++;
            i1++;
            i2++;

          } else {
            _WARNING("left_alignment_SG_PSL_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk =  FIRST_DEL_LEFT;
          switch (choice) {
          case CHOICE_LEFT_POS_BEGIN_LIST:
            fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
            break;
          case CHOICE_RIGHT_POS_BEGIN_LIST:
            fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
            break;
          case CHOICE_BLOCK_LENGTH_LIST:
            fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
            break;
          }
          i1++;

          (*pPSL_left_pos_begin)  =  pos1 + i1;
          (*pPSL_right_pos_begin) =  pos2 + i2;
          (*pPSL_nb)++;
          (*pPSL_length) = 0;

        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk =  FIRST_DEL_RIGHT;
            switch (choice) {
            case CHOICE_LEFT_POS_BEGIN_LIST:
              fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
              break;
            case CHOICE_RIGHT_POS_BEGIN_LIST:
              fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
              break;
            case CHOICE_BLOCK_LENGTH_LIST:
              fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
              break;
            }
            i2++;

            (*pPSL_left_pos_begin)  =  pos1 + i1;
            (*pPSL_right_pos_begin) =  pos2 + i2;
            (*pPSL_nb)++;
            (*pPSL_length) = 0;


          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk =  CONT_DEL_LEFT;
              i1++;

            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk =  CONT_DEL_RIGHT;
                i2++;

              } else {
                _WARNING("left_alignment_SG_PSL_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {

      while (j > 1) {
        j--;
        bk = CONT_DEL_LEFT;
        i1++;
      }

      bk = FIRST_DEL_LEFT;
      switch (choice) {
      case CHOICE_LEFT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
        break;
      case CHOICE_RIGHT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
        break;
      case CHOICE_BLOCK_LENGTH_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
        break;
      }
      i1++;

      (*pPSL_left_pos_begin)  =  pos1 + i1;
      (*pPSL_right_pos_begin) =  pos2 + i2;
      (*pPSL_nb)++;
      (*pPSL_length) = 0;

    }

    if (j < 0) {

      while (j < -1) {
        j++;
        bk = CONT_DEL_RIGHT;
        i2++;
      }

      bk = FIRST_DEL_RIGHT;
      switch (choice) {
      case CHOICE_LEFT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
        break;
      case CHOICE_RIGHT_POS_BEGIN_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
        break;
      case CHOICE_BLOCK_LENGTH_LIST:
        fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
        break;
      }
      i2++;

      (*pPSL_left_pos_begin)  =  pos1 + i1;
      (*pPSL_right_pos_begin) =  pos2 + i2;
      (*pPSL_nb)++;
      (*pPSL_length) = 0;

    }
  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);

  return 0;
}




long int right_alignment_SG_PSL_Border(char * data1, long int pos1, long int len1,
                                       char * data2, long int pos2, long int len2,
                                       /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                       /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                       long int choice
                                       ) {


  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);
  long int k = 0;
  backtrack_SG *  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,right_alignment_SG_PSL_Border_noflush);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < row1; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(row1,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(row1,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  if (len1  < len2) {

    long int i = len1, j = len2 - len1;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("right_alignment_SG_PSL_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_RIGHT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_LEFT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_RIGHT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_LEFT;
              } else {
                _WARNING("right_alignment_SG_PSL_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = len1 - len2;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("right_alignment_SG_PSL_Border() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_LEFT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_RIGHT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_LEFT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_RIGHT;
              } else {
                _WARNING("right_alignment_SG_PSL_Border() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);


  /* [3] display alignment */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {
      switch(D3(k--)) {

      case I_TO_M_SUBS:
        (*pPSL_left_pos_begin)  =  pos1 + i1;
        (*pPSL_right_pos_begin) =  pos2 + i2;
      case M_TO_M_SUBS:
        (*pPSL_length)++;
        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
        switch (choice) {
        case CHOICE_LEFT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
          break;
        case CHOICE_RIGHT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
          break;
        case CHOICE_BLOCK_LENGTH_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
          break;
        }
        (*pPSL_nb)++;
        (*pPSL_length) = 0;
      case CONT_DEL_LEFT:
        i1++;
        break;

      case FIRST_DEL_RIGHT:
        switch (choice) {
        case CHOICE_LEFT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
          break;
        case CHOICE_RIGHT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
          break;
        case CHOICE_BLOCK_LENGTH_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
          break;
        }
        (*pPSL_nb)++;
        (*pPSL_length) = 0;
      case CONT_DEL_RIGHT:
        i2++;
        break;

      default:
        _WARNING("right_alignment_SG_PSL_Border() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return 0;
}


/* III.b) Alignment methods :
 */


long int alignment_SG_PSL(char * data1, long int pos1, long int len1,
                          char * data2, long int pos2, long int len2,
                          /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                          /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                          long int choice
                          ) {

  long int i,j,k;
  long int result;
  long int *mt, *in;
  backtrack_SG * bk;

  mt = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(mt,alignment_SG_PSL);
  in = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(in,alignment_SG_PSL);
  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,alignment_SG_PSL);

  /* [1] compute */

  M3(0,0) = 0;

  I3(0,0) = gp_cost_gap_opened - gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M3(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I3(i,0) = I3(i-1,0) + gp_cost_gap_continued;


  for (j = 1; j <= len2; j++) {

    M3(0,j) =  -INFINITY_INT;
    I3(0,j) =  gp_cost_gap_opened  + (j-1)*gp_cost_gap_continued;

    for (i = 1; i <= len1; i++) {
      I3(i,j) = MAX (
                     MAX( M3(i,j-1) + gp_cost_gap_opened    , M3(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I3(i,j-1) + gp_cost_gap_continued , I3(i-1,j) + gp_cost_gap_continued )
                     );
      M3(i,j) = MAX ( M3(i-1,j-1) + S3(i,j), I3(i-1,j-1) + S3(i,j) );
    }
  }

  result = MAX( M3(len1,len2) , I3(len1,len2) );

  /* [2] backtrace */

  i = len1;
  j = len2;
  k = 0;

  if (M3(len1,len2) < I3(len1,len2))
    D3(k++) = I_TO_M_SUBS;
  else
    D3(k++) = M_TO_M_SUBS;


  while (i > 0 && j > 0) {

    if (D3(k-1) == M_TO_M_SUBS || D3(k-1) == FIRST_DEL_LEFT || D3(k-1) == FIRST_DEL_RIGHT) {

      /* no gap to open  */
      if (M3(i,j) == M3(i-1,j-1) + S3(i,j)) {
        i--;
        j--;
        D3(k++) = M_TO_M_SUBS;
      } else {
        if (M3(i,j) == I3(i-1,j-1) + S3(i,j)) {
          i--;
          j--;
          D3(k++) = I_TO_M_SUBS;
        } else {
          _WARNING("alignment_SG_PSL() 1");
        }
      }
    } else {
      /* need to open a gap */
      if (I3(i,j) == M3(i,j-1) + gp_cost_gap_opened) {
        j--;
        D3(k++) =  FIRST_DEL_RIGHT;
      } else {
        if (I3(i,j) == M3(i-1,j) + gp_cost_gap_opened) {
          i--;
          D3(k++) =  FIRST_DEL_LEFT;
        } else {
          if (I3(i,j) == I3(i,j-1) + gp_cost_gap_continued) {
            j--;
            D3(k++) =  CONT_DEL_RIGHT;
          } else {
            if (I3(i,j) == I3(i-1,j) + gp_cost_gap_continued) {
              i--;
              D3(k++) =  CONT_DEL_LEFT;
            } else {
              _WARNING("alignment_SG_PSL() 2");
            }
          }
        }
      }
    }
  }/* while */


  if (j > 0) {
    while (j > 1) {
      j--;
      D3(k++) = CONT_DEL_RIGHT;
    }
    D3(k++) = FIRST_DEL_RIGHT;
  }

  if (i > 0) {
    while (i > 1) {
      i--;
      D3(k++) = CONT_DEL_LEFT;
    }
    D3(k++) = FIRST_DEL_LEFT;
  }

  FREE(mt,(len1+1) * (len2+1) * sizeof(long int));
  FREE(in,(len1+1) * (len2+1) * sizeof(long int));


  /* [3] backtrace to compute some statistics */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {

      switch(D3(k--)) {
      case I_TO_M_SUBS:
        (*pPSL_left_pos_begin)  =  pos1 + i1;
        (*pPSL_right_pos_begin) =  pos2 + i2;
      case M_TO_M_SUBS:
        (*pPSL_length)++;
        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
        switch (choice) {
        case CHOICE_LEFT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
          break;
        case CHOICE_RIGHT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
          break;
        case CHOICE_BLOCK_LENGTH_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
          break;
        }
        (*pPSL_nb)++;
        (*pPSL_length) = 0;
      case CONT_DEL_LEFT:
        i1++;
        break;

      case FIRST_DEL_RIGHT:
        switch (choice) {
        case CHOICE_LEFT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
          break;
        case CHOICE_RIGHT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
          break;
        case CHOICE_BLOCK_LENGTH_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
          break;
        }
        (*pPSL_nb)++;
        (*pPSL_length) = 0;
      case CONT_DEL_RIGHT:
        i2++;
        break;

      default:
        _WARNING("alignment_SG_PSL() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return result;
}




long int alignment_SG_PSL_Border(char * data1, long int pos1, long int len1,
                                 char * data2, long int pos2, long int len2,
                                 /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                 /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                 long int choice) {

  /* band alignment length anw width*/
  long int dwidth       = ABS(len2-len1);
  long int bandwidth    = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + dwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);
  long int k = 0;
  backtrack_SG *  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,alignment_SG_PSL_Border);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < width; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(width,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<width-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[width-1] = MAX(bi1[width-2] + gp_cost_gap_continued, bm1[width-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(width,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-+bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<width-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[width-1] = MAX(bi1[width-2] + gp_cost_gap_continued, bm1[width-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */





  /* [2] backtrace */
  if (len1  < len2) {

    long int i = len1, j = dwidth;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("alignment_SG_PSL_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_RIGHT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_LEFT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_RIGHT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_LEFT;
              } else {
                _WARNING("alignment_SG_PSL_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = dwidth;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("alignment_SG_PSL_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_LEFT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_RIGHT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_LEFT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_RIGHT;
              } else {
                _WARNING("alignment_SG_PSL_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }
  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);


  /* [3] display alignment */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {
      switch(D3(k--)) {

      case I_TO_M_SUBS:
        (*pPSL_left_pos_begin)  =  pos1 + i1;
        (*pPSL_right_pos_begin) =  pos2 + i2;
      case M_TO_M_SUBS:
        (*pPSL_length)++;
        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
        switch (choice) {
        case CHOICE_LEFT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
          break;
        case CHOICE_RIGHT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
          break;
        case CHOICE_BLOCK_LENGTH_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
          break;
        }
        (*pPSL_nb)++;
        (*pPSL_length) = 0;
      case CONT_DEL_LEFT:
        i1++;
        break;

      case FIRST_DEL_RIGHT:
        switch (choice) {
        case CHOICE_LEFT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_left_pos_begin));
          break;
        case CHOICE_RIGHT_POS_BEGIN_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_right_pos_begin));
          break;
        case CHOICE_BLOCK_LENGTH_LIST:
          fprintf(OUTSTREAM,"%ld,",(*pPSL_length));
          break;
        }
        (*pPSL_nb)++;
        (*pPSL_length) = 0;
      case CONT_DEL_RIGHT:
        i2++;
        break;

      default:
        _WARNING("alignment_SG_PSL_Border_noflush() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return 0;
}






long int alignment_SG_PSL_Strait(char * data1, long int pos1,
                                 char * data2, long int pos2, long int len,
                                 /* out */ long int * pPSL_left_pos_begin, /* out */ long int * pPSL_right_pos_begin,
                                 /* out */ long int * pPSL_length,         /* out */ long int * pPSL_nb,
                                 long int choice
                                 ) {

  long int i;
  for (i = 0; i < len; i++) {
    (*pPSL_length)++;
  }
  return 0;
}




long int alignment_SG_PSL_on_MA(char * data1, char * data2,
                                long int left_correction,
                                MA * ma,
                                long int choice
                                ) {

  tuple * t = ma->first_tuple, * t_prev = ma->first_tuple;

  /* PSL properties */
  long int PSL_left_pos_begin  = ma->left_pos_begin;
  long int PSL_right_pos_begin = ma->right_pos_begin;
  long int PSL_length          = 0;
  long int PSL_nb              = 1;


  /* [1] Compute alignment before first tuple */
  long int left_size  = TBL_POS(t) - ma->left_pos_begin;
  long int right_size = TBR_POS(t) - ma->right_pos_begin;
  if (left_size > 0 && right_size > 0)
    left_alignment_SG_PSL_Border(data1, TBL_POS(t), left_size,
                                 data2, TBR_POS(t), right_size,
                                 &PSL_left_pos_begin, &PSL_right_pos_begin,
                                 &PSL_length,  &PSL_nb,
                                 choice
                                 );

  /* [2] Compute alignment for several tuples */
  while (t != NULL) {
    /* (1) compute inside a tuple */
    alignment_SG_PSL_Strait(data1, TBL_POS(t),
                            data2, TBR_POS(t), TSIZE(t),
                            &PSL_left_pos_begin, &PSL_right_pos_begin,
                            &PSL_length,  &PSL_nb,
                            choice
                            );
    if (t->next) {
      /* (2) alignment between (t) et (t->next) */
      alignment_SG_PSL_Border(data1, TEL_POS(t), TGAP_L(t,t->next),
                              data2, TER_POS(t), TGAP_R(t,t->next),
                              &PSL_left_pos_begin, &PSL_right_pos_begin,
                              &PSL_length,  &PSL_nb,
                              choice
                              );
    }
    t_prev = t;
    t = t->next;
  }/* while */


  /* [3] Compute alignement after last tuple */
  {
    long int left_size  =  ma->left_pos_end  - TEL_POS(t_prev);
    long int right_size =  ma->right_pos_end - TER_POS(t_prev);
    if (left_size > 0 && right_size > 0)
      right_alignment_SG_PSL_Border(data1, TEL_POS(t_prev), left_size,
                                    data2, TER_POS(t_prev), right_size,
                                    &PSL_left_pos_begin, &PSL_right_pos_begin,
                                    &PSL_length,  &PSL_nb,
                                    choice
                                    );
  }

  /* last block */
  switch (choice) {
  case CHOICE_LEFT_POS_BEGIN_LIST:
    fprintf(OUTSTREAM,"%ld,",(PSL_left_pos_begin));
    break;
  case CHOICE_RIGHT_POS_BEGIN_LIST:
    fprintf(OUTSTREAM,"%ld,",(PSL_right_pos_begin));
    break;
  case CHOICE_BLOCK_LENGTH_LIST:
    fprintf(OUTSTREAM,"%ld,",(PSL_length));
    break;
  case CHOICE_BLOCK_NB:
    fprintf(OUTSTREAM,"%ld",(PSL_nb));
    break;
  }

  return 0;
}









/*-----------------------------------------------------------------*/

/* IV) YASS output functions:
 *
 * WARNING : set in "data1" its reversed complementary "datarev" when needed.
 * coordonninate "pos1" must then be explicitely expressed.
 */



/*
 *  display a buffer composed of 5 lines with NBL chars
 */


long int line_pos = 0;
char line_dpos1[NBL+1];
char line_d1[NBL+1];
char line_op[NBL+1];
char line_d2[NBL+1];
char line_dpos2[NBL+1];



void flush()
{
  long int i;
  if (line_pos > 0) {
    /* close strings */
    i = line_pos;
    if (line_dpos1[i]!='\0') {
      while (line_dpos1[i]!=' ' && i>0)
        i--;
      line_dpos1[i]='\0';
    }

    i = line_pos;
    if (line_dpos2[i]!='\0') {
      while (line_dpos2[i]!=' ' && i>0)
        i--;
      line_dpos2[i]='\0';
    }

    line_dpos1[line_pos] = '\0';
    line_dpos2[line_pos] = '\0';
    line_d1[line_pos]    = '\0';
    line_d2[line_pos]    = '\0';
    line_op[line_pos]    = '\0';
    /* print the 5 lines */
    fprintf(OUTSTREAM,"%s",line_dpos1);fprintf(OUTSTREAM,"\n");
    fprintf(OUTSTREAM,"%s",line_d1   );fprintf(OUTSTREAM,"\n");
    fprintf(OUTSTREAM,"%s",line_op   );fprintf(OUTSTREAM,"\n");
    fprintf(OUTSTREAM,"%s",line_d2   );fprintf(OUTSTREAM,"\n");
    fprintf(OUTSTREAM,"%s",line_dpos2);fprintf(OUTSTREAM,"\n");
    fprintf(OUTSTREAM,"\n");
    /* reset buffers */
    line_pos = 0;
    memset(line_dpos1,0,NBL+1);
    memset(line_dpos2,0,NBL+1);
    memset(line_d1   ,0,NBL+1);
    memset(line_d2   ,0,NBL+1);
    memset(line_op   ,0,NBL+1);
  }
}

void push(char d1, long int dpos1,char op,char d2, long int dpos2)
{
  /* put positions in the buffer */
  if ((d1 != CHAR_DUM) && !(dpos1%10) && size_lint(dpos1) + line_pos < NBL) {
    sprintf(line_dpos1+line_pos,"|%ld",dpos1);
  } else {
    if (line_dpos1[line_pos]=='\0')
      line_dpos1[line_pos]=' ';
  }

  if ((d2 != CHAR_DUM) && !(dpos2%10) && size_lint(dpos2) + line_pos < NBL) {
    sprintf(line_dpos2+line_pos,"|%ld",dpos2);
  } else {
    if (line_dpos2[line_pos]=='\0')
      line_dpos2[line_pos]=' ';
  }

  /* memorize alignment */
  line_d1[line_pos]=d1;
  line_d2[line_pos]=d2;
  line_op[line_pos]=op;

  line_pos++;

  /* print if end of line */
  if (line_pos == NBL) {
    flush();
  }
}




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
                                                  char * data2, long int pos2, long int len2) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < row1; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 1; i <= len1; i++) {

      long int * bm0 = mt[i-1];
      long int * bm1 = mt[i];
      long int * bi0 = in[i-1];
      long int * bi1 = in[i];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i-bandwidth+j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i-bandwidth+j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(row1-1,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 1; i <= len2; i++) {

      long int * bm0 = mt[i-1];
      long int * bm1 = mt[i];
      long int * bi0 = in[i-1];
      long int * bi1 = in[i];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i == data2+pos2-i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<=MIN(row1-1,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1-i+bandwidth-j == data2+pos2-i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1-i+bandwidth-j])][(long int)(data2[pos2-i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  /* [3] display alignment */

  if (len1  < len2) {

    long int i = len1, j = len2 - len1;
    long int i1 = -len1, i2 = -len2;
    backtrack_SG bk;

    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk = I_TO_M_SUBS;
    else
      bk = M_TO_M_SUBS;

    push(LOOKUP(data1[pos1+i1]),
         BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
         (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
         LOOKUP(data2[pos2+i2]),
         BCOUNT + pos2+i2
         );

    i1++;
    i2++;
    i--;

    while (i > 0) {

      if (bk == M_TO_M_SUBS || bk == FIRST_DEL_LEFT || bk == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-j])]:
             gp_substitution_matrix[(long int)(data1[pos1-i+j])][(long int)(data2[pos2-i])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk = M_TO_M_SUBS;

          push(LOOKUP(data1[pos1+i1]),
               BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
               (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
               LOOKUP(data2[pos2+i2]),
               BCOUNT + pos2+i2
               );
          i1++;
          i2++;

        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i-j])]:
               gp_substitution_matrix[(long int)(data1[pos1-i+j])][(long int)(data2[pos2-i])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk = I_TO_M_SUBS;

            push(LOOKUP(data1[pos1+i1]),
                 BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                 (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
                 LOOKUP(data2[pos2+i2]),
                 BCOUNT + pos2+i2
                 );
            i1++;
            i2++;

          } else {
            _WARNING("display_left_alignment_SG_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk =  FIRST_DEL_RIGHT;

          push(CHAR_DUM,
               BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
               CHAR_INS,
               LOOKUP(data2[pos2+i2]),
               BCOUNT + pos2+i2
               );
          i2++;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk =  FIRST_DEL_LEFT;

            push(LOOKUP(data1[pos1+i1]),
                 BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                 CHAR_DEL,
                 CHAR_DUM,
                 BCOUNT + pos2+i2
                 );
            i1++;

          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk =  CONT_DEL_RIGHT;

              push(CHAR_DUM,
                   BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                   CHAR_INS,
                   LOOKUP(data2[pos2+i2]),
                   BCOUNT + pos2+i2
                   );
              i2++;

            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk =  CONT_DEL_LEFT;

                push(LOOKUP(data1[pos1+i1]),
                     BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                     CHAR_DEL,
                     CHAR_DUM,
                     BCOUNT + pos2+i2
                     );
                i1++;

              } else {
                _WARNING("display_left_alignment_SG_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */


    if (j > 0) {

      while (j > 1) {
        j--;
        bk = CONT_DEL_RIGHT;
        push(CHAR_DUM,
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_INS,
             LOOKUP(data2[pos2+i2]),
             BCOUNT + pos2+i2
             );
        i2++;
      }

      bk = FIRST_DEL_RIGHT;
      push(CHAR_DUM,
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           CHAR_INS,
           LOOKUP(data2[pos2+i2]),
           BCOUNT + pos2+i2
           );
      i2++;
    }

    if (j < 0) {

      while (j < -1) {
        j++;
        bk = CONT_DEL_LEFT;
        push(LOOKUP(data1[pos1+i1]),
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_DEL,
             CHAR_DUM,
             BCOUNT + pos2+i2
             );
        i1++;
      }

      bk = FIRST_DEL_LEFT;
      push(LOOKUP(data1[pos1+i1]),
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           CHAR_DEL,
           CHAR_DUM,
           BCOUNT + pos2+i2
           );
      i1++;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = len1 - len2;
    long int i1 = -len1, i2 = -len2;
    backtrack_SG bk;

    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk = I_TO_M_SUBS;
    else
      bk = M_TO_M_SUBS;

    push(LOOKUP(data1[pos1+i1]),
         BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
         (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
         LOOKUP(data2[pos2+i2]),
         BCOUNT + pos2+i2
         );

    i1++;
    i2++;
    i--;

    while (i > 0) {

      if (bk == M_TO_M_SUBS || bk == FIRST_DEL_LEFT || bk == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1-i-j])][(long int)(data2[pos2-i])]:
             gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+j])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk = M_TO_M_SUBS;

          push(LOOKUP(data1[pos1+i1]),
               BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
               (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
               LOOKUP(data2[pos2+i2]),
               BCOUNT + pos2+i2
               );
          i1++;
          i2++;

        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1-i-j])][(long int)(data2[pos2-i])]:
               gp_substitution_matrix[(long int)(data1[pos1-i])][(long int)(data2[pos2-i+j])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk = I_TO_M_SUBS;

            push(LOOKUP(data1[pos1+i1]),
                 BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                 (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
                 LOOKUP(data2[pos2+i2]),
                 BCOUNT + pos2+i2
                 );
            i1++;
            i2++;

          } else {
            _WARNING("display_left_alignment_SG_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk =  FIRST_DEL_LEFT;

          push(LOOKUP(data1[pos1+i1]),
               BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
               CHAR_DEL,
               CHAR_DUM,
               BCOUNT + pos2+i2
               );
          i1++;

        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk =  FIRST_DEL_RIGHT;

            push(CHAR_DUM,
                 BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                 CHAR_INS,
                 LOOKUP(data2[pos2+i2]),
                 BCOUNT + pos2+i2
                 );
            i2++;

          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk =  CONT_DEL_LEFT;

              push(LOOKUP(data1[pos1+i1]),
                   BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                   CHAR_DEL,
                   CHAR_DUM,
                   BCOUNT + pos2+i2
                   );
              i1++;

            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk =  CONT_DEL_RIGHT;

                push(CHAR_DUM,
                     BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
                     CHAR_INS,
                     LOOKUP(data2[pos2+i2]),
                     BCOUNT + pos2+i2
                     );
                i2++;

              } else {
                _WARNING("display_left_alignment_SG_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */


    if (j < 0) {

      while (j < -1) {
        j++;
        bk = CONT_DEL_RIGHT;
        push(CHAR_DUM,
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_INS,
             LOOKUP(data2[pos2+i2]),
             BCOUNT + pos2+i2
             );
        i2++;
      }

      bk = FIRST_DEL_RIGHT;
      push(CHAR_DUM,
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           CHAR_INS,
           LOOKUP(data2[pos2+i2]),
           BCOUNT + pos2+i2
           );
      i2++;
    }

    if (j > 0) {

      while (j > 1) {
        j--;
        bk = CONT_DEL_LEFT;
        push(LOOKUP(data1[pos1+i1]),
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_DEL,
             CHAR_DUM,
             BCOUNT + pos2+i2
             );
        i1++;
      }

      bk = FIRST_DEL_LEFT;
      push(LOOKUP(data1[pos1+i1]),
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           CHAR_DEL,
           CHAR_DUM,
           BCOUNT + pos2+i2
           );
      i1++;
    }
    /*<<*/
  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);

  return 0;
}




long int display_right_alignment_SG_Border_noflush(long int querychunk, long int reverse,
                                                   char * data1, long int pos1, long int len1,
                                                   char * data2, long int pos2, long int len2
                                                   ) {

  /* band alignment length anw width*/
  long int bandwidth  = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int row1       = 2*bandwidth + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);
  long int k = 0;
  backtrack_SG *  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,display_right_alignment_SG_Border_noflush);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < row1; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(row1,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i=0; i < len2; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(row1,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<row1-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[row1-1] = MAX(bi1[row1-2] + gp_cost_gap_continued, bm1[row1-2] + gp_cost_gap_opened);
    }/* i */
  } /* len1 <> len2 */



  /* [2] backtrace */
  if (len1  < len2) {

    long int i = len1, j = len2 - len1;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("display_right_alignment_SG_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_RIGHT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_LEFT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_RIGHT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_LEFT;
              } else {
                _WARNING("display_right_alignment_SG_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = len1 - len2;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("display_right_alignment_SG_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_LEFT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_RIGHT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_LEFT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_RIGHT;
              } else {
                _WARNING("display_right_alignment_SG_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }
  } /* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);




  /* [3] display alignment */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {

      switch(D3(k--)) {
      case M_TO_M_SUBS:
      case I_TO_M_SUBS:
        push(LOOKUP(data1[pos1+i1]),
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
             LOOKUP(data2[pos2+i2]),
             BCOUNT + pos2+i2
             );
        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
      case CONT_DEL_LEFT:
        push(LOOKUP(data1[pos1+i1]),
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_DEL,
             CHAR_DUM,
             BCOUNT + pos2+i2
             );
        i1++;
        break;

      case FIRST_DEL_RIGHT:
      case CONT_DEL_RIGHT:
        push(CHAR_DUM,
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_INS,
             LOOKUP(data2[pos2+i2]),
             BCOUNT + pos2+i2
             );
        i2++;
        break;

      default:
        _WARNING("display_right_alignment_SG_Border_noflush() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return 0;
}



/* IV.b) Alignment methods :
 */


long int display_alignment_SG_Border_noflush(long int querychunk, long int reverse,
                                             char * data1, long int pos1, long int len1,
                                             char * data2, long int pos2, long int len2) {
  /* band alignment length anw width*/
  long int dwidth       = ABS(len2-len1);
  long int bandwidth    = gp_delta_stat + 1;

  /* various switch elements to access information */
  long int width       = 2*bandwidth + dwidth + 1;
  long int height      = MIN(len1,len2) + 1;

  /* dp tables */
  long int ** mt = lint_directtable(height,width,-INFINITY_INT);
  long int ** in = lint_directtable(height,width,-INFINITY_INT);
  long int k = 0;
  backtrack_SG *  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,display_alignment_SG_Border_noflush);


  /* [1] compute */

  /* init table */
  {
    long int * bm0 = mt[0];
    long int * bi0 = in[0];
    long int j;
    for (j = 0; j < width; j++) {
      bm0[j] = -INFINITY_INT;
      bi0[j] = (gp_cost_gap_opened-gp_cost_gap_continued) + ABS(j-bandwidth)*gp_cost_gap_continued;
    }
    bm0[bandwidth] = 0;
    bi0[bandwidth] = -INFINITY_INT;
  }

  /* len1 < len2 */
  if (len1 < len2) {
    long int i;
    for (i = 0; i < len1; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j = MAX(0,i+bandwidth-len1+1); j < bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i+bandwidth-j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i+bandwidth-j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(width,len2+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i-bandwidth+j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i-bandwidth+j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<width-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[width-1] = MAX(bi1[width-2] + gp_cost_gap_continued, bm1[width-2] + gp_cost_gap_opened);
    }/* i */
  } else { /* len2 < len1 */
    long int i;
    for (i = 0; i < len2; i++) {

      long int * bm0 = mt[i];
      long int * bm1 = mt[i+1];
      long int * bi0 = in[i];
      long int * bi1 = in[i+1];

      /* M(j) 2 full vectorisable loops */
      {
        long int j;
        for (j=MAX(0,i+bandwidth-len2+1); j<bandwidth; j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i == data2+pos2+i+bandwidth-j)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i])][(long int)(data2[pos2+i+bandwidth-j])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }
      {
        long int j;
        for (j=bandwidth; j<MIN(width,len1+bandwidth-i); j++) {
          /* [SFS] : single file stop */
          long int d;
          if (data1 == data2 && data1+pos1+i-bandwidth+j == data2+pos2+i)
            d = -INFINITY_INT;
          else
            d = gp_substitution_matrix[(long int)(data1[pos1+i-bandwidth+j])][(long int)(data2[pos2+i])];
          bm1[j] = MAX( bm0[j] , bi0[j] ) + d;
        }
      }


      /* I(j) in two half loops */
      bi1[bandwidth] = MAX(
                           MAX(bi0[bandwidth-1],bi0[bandwidth+1]) + gp_cost_gap_continued,
                           MAX(bm0[bandwidth-1],bm0[bandwidth+1]) + gp_cost_gap_opened
                           );
      /* first upper half-loop */
      {
        long int j;
        for (j=bandwidth-1; j>0; j--) {
          bi1[j] = MAX(
                       MAX(bi0[j-1],bi1[j+1]) + gp_cost_gap_continued,
                       MAX(bm0[j-1],bm1[j+1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[0] = MAX(bi1[1] + gp_cost_gap_continued, bm1[1] + gp_cost_gap_opened);

      /* second lower half-loop */
      {
        long int j;
        for (j=bandwidth+1; j<width-1; j++) {
          bi1[j] = MAX(
                       MAX(bi0[j+1],bi1[j-1]) + gp_cost_gap_continued,
                       MAX(bm0[j+1],bm1[j-1]) + gp_cost_gap_opened
                       );
        }
      }
      bi1[width-1] = MAX(bi1[width-2] + gp_cost_gap_continued, bm1[width-2] + gp_cost_gap_opened);
    }/* i */
  }/* len1 <> len2 */



  /* [2] backtrace */
  if (len1  < len2) {

    long int i = len1, j = dwidth;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i+j-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-j-1])][(long int)(data2[pos2+i-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("display_alignment_SG_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_RIGHT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_LEFT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_RIGHT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_LEFT;
              } else {
                _WARNING("display_alignment_SG_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

  } else { /* len1 >= len2 */

    long int i = len2, j = dwidth;
    if (mt[i][bandwidth + j] < in[i][bandwidth + j])
      bk[k++] = I_TO_M_SUBS;
    else
      bk[k++] = M_TO_M_SUBS;


    while (i > 0) {

      if (bk[k-1] == M_TO_M_SUBS || bk[k-1] == FIRST_DEL_LEFT || bk[k-1] == FIRST_DEL_RIGHT) {

        /* no gap to open */

        if (
            mt[i][bandwidth + j] == mt[i-1][bandwidth + j] +
            (
             (j>=0)?
             gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
             gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
             )
            /* M3(i,j) == M3(i-1,j-1) + S3(i,j) */
            ) {
          i--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(MM) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] = M_TO_M_SUBS;
        } else {
          if (
              mt[i][bandwidth + j] == in[i-1][bandwidth + j] +
              (
               (j>=0)?
               gp_substitution_matrix[(long int)(data1[pos1+i+j-1])][(long int)(data2[pos2+i-1])]:
               gp_substitution_matrix[(long int)(data1[pos1+i-1])][(long int)(data2[pos2+i-j-1])]
               )
              /* M3(i,j) == I3(i-1,j-1) + S3(i,j) */
                ) {
            i--;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(MI) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] = I_TO_M_SUBS;
          } else {
            _WARNING("display_alignment_SG_Border_noflush() 1");
          }
        }
      } else {
        /* need to open a gap */
        if ( j + bandwidth > 0 &&
             in[i][bandwidth + j] == mt[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_opened
             /* I3(i,j) == M3(i,j-1) + gp_cost_gap_opened */) {
          if (j<=0)
            i--;
          j--;
#ifdef DEBUG_BORDER_BACKTRACK
          printf("(II1) i:%ld j:%ld\n",i,j);
#endif
          bk[k++] =  FIRST_DEL_LEFT;
        } else {
          if ( j + bandwidth < width - 1 &&
               in[i][bandwidth + j] == mt[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_opened
               /* I3(i,j) == M3(i-1,j) + gp_cost_gap_opened */) {
            if (j>=0)
              i--;
            j++;
#ifdef DEBUG_BORDER_BACKTRACK
            printf("(II2) i:%ld j:%ld\n",i,j);
#endif
            bk[k++] =  FIRST_DEL_RIGHT;
          } else {
            if ( j + bandwidth > 0 &&
                 in[i][bandwidth + j] == in[i+(j<=0?-1:0)][bandwidth + j - 1] + gp_cost_gap_continued
                 /* I3(i,j) == I3(i,j-1) + gp_cost_gap_continued */) {
              if (j<=0)
                i--;
              j--;
#ifdef DEBUG_BORDER_BACKTRACK
              printf("(II3) i:%ld j:%ld\n",i,j);
#endif
              bk[k++] =  CONT_DEL_LEFT;
            } else {
              if ( j + bandwidth < width - 1 &&
                   in[i][bandwidth + j] == in[i+(j>=0?-1:0)][bandwidth + j + 1] + gp_cost_gap_continued
                   /* I3(i,j) == I3(i-1,j) + gp_cost_gap_continued */) {
                if (j>=0)
                  i--;
                j++;
#ifdef DEBUG_BORDER_BACKTRACK
                printf("(II4) i:%ld j:%ld\n",i,j);
#endif
                bk[k++] =  CONT_DEL_RIGHT;
              } else {
                _WARNING("display_alignment_SG_Border_noflush() 2");
              }
            }
          }
        }
      }
    }/* while */

    if (j > 0) {
      while (j > 1) {
        j--;
        bk[k++] = CONT_DEL_LEFT;
      }
      bk[k++] = FIRST_DEL_LEFT;
    }

    if (j < 0) {
      while (j < -1) {
        j++;
        bk[k++] = CONT_DEL_RIGHT;
      }
      bk[k++] = FIRST_DEL_RIGHT;
    }
  }/* len1 <> len2 */

  lint_free_directtable(mt,height,width);
  lint_free_directtable(in,height,width);




  /* [3] display alignment */
  {
    long int i1 = 0, i2 = 0;
    k--;

    while (i1 < len1 || i2 < len2) {

      switch(D3(k--)) {
      case M_TO_M_SUBS:
      case I_TO_M_SUBS:
        push(LOOKUP(data1[pos1+i1]),
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
             LOOKUP(data2[pos2+i2]),
             BCOUNT + pos2+i2
             );
        i1++;
        i2++;
        break;

      case FIRST_DEL_LEFT:
      case CONT_DEL_LEFT:
        push(LOOKUP(data1[pos1+i1]),
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_DEL,
             CHAR_DUM,
             BCOUNT + pos2+i2
             );
        i1++;
        break;

      case FIRST_DEL_RIGHT:
      case CONT_DEL_RIGHT:
        push(CHAR_DUM,
             BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
             CHAR_INS,
             LOOKUP(data2[pos2+i2]),
             BCOUNT + pos2+i2
             );
        i2++;
        break;

      default:
        _WARNING("display_alignment_SG_Border_noflush() 3");
      }/* case */
    }/* while */
    FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  }
  return 0;
}





long int display_alignment_SG_noflush(long int querychunk, long int reverse,
                                      char * data1, long int pos1, long int len1,
                                      char * data2, long int pos2, long int len2)
{

  long int i,j;
  long int i1,i2;
  long int k;
  long int * mt, * in;
  backtrack_SG * bk;

  mt = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(mt,display_alignment_SG_noflush);
  in = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(in,display_alignment_SG_noflush);
  bk = (backtrack_SG*) MALLOC( (len1+len2+1) * sizeof(backtrack_SG));
  ASSERT(bk,display_alignment_SG_noflush);

  /* [1] compute */

  M3(0,0) = 0;

  I3(0,0) = gp_cost_gap_opened - gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M3(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I3(i,0) = I3(i-1,0) +  gp_cost_gap_continued;

  for (j = 1; j <= len2; j++) {

    M3(0,j) =  -INFINITY_INT;
    I3(0,j) =  I3(0,j-1)  + gp_cost_gap_continued;

    for (i = 1; i <= len1; i++) {
      I3(i,j) = MAX (
                     MAX( M3(i,j-1) + gp_cost_gap_opened    , M3(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I3(i,j-1) + gp_cost_gap_continued , I3(i-1,j) + gp_cost_gap_continued )
                     );
      if  (data1 + pos1 + i !=  data2 + pos2 + j)
        M3(i,j) = MAX ( M3(i-1,j-1) + S3(i,j), I3(i-1,j-1) + S3(i,j) );
      else
        M3(i,j) = -INFINITY_INT;
    }
  }


#ifdef DEBUGDISPLAYALIGN
  {
    long int i,j;
    fprintf(stdout,"   ");
    for (i = 0; i <= len1; i++)
      fprintf(stdout,"|  %c  ",LOOKUP((long int)data1[pos1+i]));
    fprintf(stdout,"|\n");
    for (j = 0; j <= len2; j++) {
      fprintf(stdout," %c |",LOOKUP((long int)data2[pos2+j]));
      for (i = 0; i <= len1; i++) {
        fprintf(stdout,"%2.1g|",(float)(M3(i,j)));
      }
      fprintf(stdout,"\n");
    }
  }
#endif

  /* [2] backtrace */
  i = len1;
  j = len2;
  k = 0;

  if (M3(len1,len2) < I3(len1,len2))
    D3(k++) = I_TO_M_SUBS;
  else
    D3(k++) = M_TO_M_SUBS;


  while (i > 0 && j > 0) {

    if (D3(k-1) == M_TO_M_SUBS || D3(k-1) == FIRST_DEL_LEFT || D3(k-1) == FIRST_DEL_RIGHT) {

      /* no gap to open */

      if (M3(i,j) == M3(i-1,j-1) + S3(i,j)) {
        i--;
        j--;
        D3(k++) = M_TO_M_SUBS;
      } else {
        if (M3(i,j) == I3(i-1,j-1) + S3(i,j)) {
          i--;
          j--;
          D3(k++) = I_TO_M_SUBS;
        } else {
          _WARNING("display_alignment_SG_noflush() 1");
        }
      }
    } else {

      /* need to open a gap */
      if (I3(i,j) == M3(i,j-1) + gp_cost_gap_opened) {
        j--;
        D3(k++) =  FIRST_DEL_RIGHT;
      } else {
        if (I3(i,j) == M3(i-1,j) + gp_cost_gap_opened) {
          i--;
          D3(k++) =  FIRST_DEL_LEFT;
        } else {
          if (I3(i,j) == I3(i,j-1) + gp_cost_gap_continued) {
            j--;
            D3(k++) =  CONT_DEL_RIGHT;
          } else {
            if (I3(i,j) == I3(i-1,j) + gp_cost_gap_continued) {
              i--;
              D3(k++) =  CONT_DEL_LEFT;
            } else {
              _WARNING("display_alignment_SG_noflush() 2");
            }
          }
        }
      }
    }
  }/* while */

  if (j > 0) {
    while (j > 1) {
      j--;
      D3(k++) = CONT_DEL_RIGHT;
    }
    D3(k++) = FIRST_DEL_RIGHT;
  }

  if (i > 0) {
    while (i > 1) {
      i--;
      D3(k++) = CONT_DEL_LEFT;
    }
    D3(k++) = FIRST_DEL_LEFT;
  }

  FREE(mt, (len1+1) * (len2+1) * sizeof(long int));
  FREE(in, (len1+1) * (len2+1) * sizeof(long int));


#ifdef DEBUGDISPLAYALIGN
  fprintf(stdout,"backtrack : ");
  for (i = k; i >= 0; i--) {
    fprintf(stdout,"[%ld]",D3(k));
  }
  fprintf(stdout,"\n");
#endif


  /* [3] display alignment */

  i1 = 0; i2 = 0; k--;

  while (i1 < len1 || i2 < len2) {

    switch(D3(k--)) {
    case M_TO_M_SUBS:
    case I_TO_M_SUBS:
      push(LOOKUP(data1[pos1+i1]),
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           (TNORM(data1[pos1+i1])==TNORM(data2[pos2+i2]))?CHAR_EQU:(((TNORM(data1[pos1+i1])^TNORM(data2[pos2+i2]))==2)?CHAR_SS:CHAR_SV),
           LOOKUP(data2[pos2+i2]),
           BCOUNT + pos2+i2
           );
      i1++;
      i2++;
      break;

    case FIRST_DEL_LEFT:
    case CONT_DEL_LEFT:
      push(LOOKUP(data1[pos1+i1]),
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           CHAR_DEL,
           CHAR_DUM,
           BCOUNT + pos2+i2
           );
      i1++;
      break;

    case FIRST_DEL_RIGHT:
    case CONT_DEL_RIGHT:
      push(CHAR_DUM,
           BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i1) : pos1+i1),
           CHAR_INS,
           LOOKUP(data2[pos2+i2]),
           BCOUNT + pos2+i2
           );
      i2++;
      break;

    default:
      _WARNING("display_alignment_SG_noflush() 3");
    }/* case */
  }/* while */
  FREE(bk, (len1+len2+1) * sizeof(backtrack_SG));
  return 0;
}



long int display_alignment_SG_Strait_noflush(long int querychunk, long int reverse,
                                             char *data1, long int pos1,
                                             char *data2, long int pos2, long int len) {
  long int i;
  for (i = 0; i < len; i++) {
    push(
         LOOKUP(data1[pos1+i]),
         BCOUNT + ((reverse) ? (gp_chunksize_query[querychunk] - 1) - (pos1 + i) : pos1+i),
         (TNORM(data1[pos1+i]) == TNORM(data2[pos2+i]))?CHAR_EQU:((TNORM(data1[pos1+i])^TNORM(data2[pos2+i]))==2)?CHAR_SS:CHAR_SV,
         LOOKUP(data2[pos2+i]),
         BCOUNT + pos2 + i
         );
  }
  return 0;
}




long int display_alignment_SG(long int querychunk, long int reverse,
                              char * data1, long int pos1, long int len1,
                              char * data2, long int pos2, long int len2)
{
  display_alignment_SG_noflush(querychunk, reverse,
                               data1,pos1,len1,
                               data2,pos2,len2);
  flush();
  return 0;
}


long int display_alignment_SG_on_MA(char * data1,char * data2,
                                    MA * ma)

{
  long int left_size,right_size;
  long int left_correction = 0;
  tuple * t = ma->first_tuple, * t_prev = ma->first_tuple;

  /* [1] Display alignment before first tuple */
  left_size  = TBL_POS(t) - ma->left_pos_begin;
  right_size = TBR_POS(t) - ma->right_pos_begin;

#ifdef DEBUG_DISPLAY
#define DEBUG_DISPLAY_PADDING 40
  display_alignment_SG_noflush(ma->j_chunk, ma->reverse,
                               data1, MAX(0, ma->left_pos_begin  - DEBUG_DISPLAY_PADDING), MIN(DEBUG_DISPLAY_PADDING,ma->left_pos_begin),
                               data2, MAX(0, ma->right_pos_begin - DEBUG_DISPLAY_PADDING), MIN(DEBUG_DISPLAY_PADDING,ma->right_pos_begin));
  push('{',1,'{','{',1);
#endif

  if (left_size > 0 && right_size > 0)
    display_left_alignment_SG_Border_noflush(ma->j_chunk, ma->reverse,
                                             data1, TBL_POS(t), left_size,
                                             data2, TBR_POS(t), right_size
                                             );

  /* [2] Display alignment for several tuples */
  while (t != NULL) {

    /* (1) display alignment of tuple "t" */
#ifdef DEBUG_DISPLAY
  push('<',1,'<','<',1);
#endif
    display_alignment_SG_Strait_noflush(ma->j_chunk, ma->reverse,
                                        data1, TBL_POS(t),
                                        data2, TBR_POS(t),
                                        TSIZE(t));
#ifdef DEBUG_DISPLAY
  push('>',1,'>','>',1);
#endif

    if (t->next) {
      /* (2) display alignment between (t) and (t->next) */
      display_alignment_SG_Border_noflush(ma->j_chunk, ma->reverse,
                                          data1, TEL_POS(t), TGAP_L(t,t->next),
                                          data2, TER_POS(t), TGAP_R(t,t->next));
    }
    t_prev = t;
    t = t->next;
  }


  /* [3] Display alignment after first tuple */
  left_size  =  ma->left_pos_end  - TEL_POS(t_prev);
  right_size =  ma->right_pos_end - TER_POS(t_prev);

  if (left_size > 0 && right_size > 0)
    display_right_alignment_SG_Border_noflush(ma->j_chunk, ma->reverse,
                                              data1, TEL_POS(t_prev), left_size,
                                              data2, TER_POS(t_prev), right_size);

#ifdef DEBUG_DISPLAY
  push('}',1,'}','}',1);
  display_alignment_SG_noflush(ma->j_chunk, ma->reverse,
                               data1, ma->left_pos_end,  MIN(DEBUG_DISPLAY_PADDING, gp_chunksize_query[ma->j_chunk] - ma->left_pos_end),
                               data2, ma->right_pos_end, MIN(DEBUG_DISPLAY_PADDING, gp_chunksize_text[ma->i_chunk]  - ma->right_pos_end));
#endif
  flush();
  return 0;
}




/*******************************
 *
 *  Compute the real score
 *
 ********************************/



long int alignment_score_on_MAs(char * data1, char * data2,MA * firstMa,long int left_correction) {
  MA * ma=firstMa;

  while (ma != NULL) {
    alignment_SG_score_on_MA(data1,data2,ma,left_correction);
#ifdef DEBUG_SCORE_MA
    fprintf(stderr,"(%ld,%ld)(%ld,%ld)  score=%ld\n",ma->left_pos_begin,ma->left_pos_end,ma->right_pos_begin,ma->right_pos_end,ma->blastscore);
#endif
    ma = ma->next;
  }
  return 1;
}




long int alignment_SG_score_on_MA(char * data1, char * data2,
                                  MA * ma,long int left_correction) {
  long int left_size,right_size;
  tuple * t = ma->first_tuple, * t_prev = ma->first_tuple;
  ma->blastscore = 0;

  /* [1] Compute alignment before first tuple */
  left_size  = TBL_POS(t) - ma->left_pos_begin;
  right_size = TBR_POS(t) - ma->right_pos_begin;

  if (left_size > 0 && right_size > 0)
    alignment_SG_score(data1, ma->left_pos_begin,  left_size,
                       data2, ma->right_pos_begin, right_size,
                       &(ma->blastscore));
#ifdef DEBUG_SCORE_MA
  fprintf(stderr,"left score : %ld {len1=%ld len2=%ld}\n",ma->blastscore,left_size,right_size);
#endif

  /* [2] Compute alignment for several tuples */
  while (t != NULL) {
    /* (1) compute inside a tuple */
    alignment_SG_score_Strait(data1,TBL_POS(t),
                              data2,TBR_POS(t),TSIZE(t),
                              &(ma->blastscore));
#ifdef DEBUG_SCORE_MA
    fprintf(stderr,"after strait score : %ld  {length=%ld,(%ld-%ld)(%ld-%ld)}\n",ma->blastscore,TSIZE(t),
            TBL_POS(t), TEL_POS(t),
            TBR_POS(t), TER_POS(t));
#endif
    if (t->next) {
      /* (2) alignment between (t) et (t->next) */
      alignment_SG_score(data1, TEL_POS(t), TGAP_L(t,t->next),
                         data2, TER_POS(t), TGAP_R(t,t->next),
                         &(ma->blastscore));
#ifdef DEBUG_SCORE_MA
      fprintf(stderr,"after between tuples : %ld\n",ma->blastscore);
#endif
    }/* if */
    t_prev = t;
    t = t->next;
  }/* while */

  /* [3] Compute alignement after last tuple */
  left_size  =  ma->left_pos_end  - TEL_POS(t_prev);
  right_size =  ma->right_pos_end - TER_POS(t_prev);
#ifdef DEBUG_SCORE_MA
  fprintf(stderr,"before last right score : %ld {len1=%ld len2=%ld}\n",ma->blastscore,left_size,right_size);
#endif

  if (left_size > 0 && right_size > 0)
    alignment_SG_score(data1, TEL_POS(t_prev), left_size,
                       data2, TER_POS(t_prev), right_size,
                       &(ma->blastscore));
#ifdef DEBUG_SCORE_MA
  fprintf(stderr,"after last right score : %ld\n",ma->blastscore);
#endif
  return 0;
}


long int alignment_SG_score_Strait(char * data1,  long int pos1, char * data2,  long int pos2, long int len,
                                   /*out*/long int    *p_score) {
  long int i;

  for (i = 0; i < len; i++) {
    *p_score += gp_substitution_matrix[((long int)(data1 + pos1)[i])][(long int)(data2 + pos2)[i]];
#ifdef DEBUG_SCORE_MA
    fprintf(stderr,"STRAIT: %c<->%c score=%ld\n",
            LOOKUP(data1[pos1+i]),LOOKUP(data2[pos2+i]),
            gp_substitution_matrix[((int)(data1 + pos1)[i])][(int)(data2 + pos2)[i]]);
#endif
  }
  return 0;
}




long int alignment_SG_score(char * data1, long int pos1, long int len1,
                            char * data2, long int pos2, long int len2,
                            /*out*/long int    *p_score
                            ) {
  long int i,j;
  long int gain_score;
  long int * mt, * in;

  mt = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(mt,alignment_SG_score);
  in = (long int*) MALLOC( (len1+1) * (len2+1) * sizeof(long int));
  ASSERT(in,alignment_SG_score);


  M3(0,0) = 0;

  I3(0,0) = gp_cost_gap_opened - gp_cost_gap_continued;

  for (i = 1; i <= len1; i++)
    M3(i,0) = -INFINITY_INT;

  for (i = 1; i <= len1; i++)
    I3(i,0) = I3(i-1,0) + gp_cost_gap_continued;


  for (j = 1; j <= len2; j++) {

    M3(0,j) =  -INFINITY_INT;
    I3(0,j) =  I3(0,j-1)  + gp_cost_gap_continued;

    for (i = 1; i <= len1; i++) {
      I3(i,j) = MAX (
                     MAX( M3(i,j-1) + gp_cost_gap_opened    , M3(i-1,j) + gp_cost_gap_opened )
                     ,
                     MAX( I3(i,j-1) + gp_cost_gap_continued , I3(i-1,j) + gp_cost_gap_continued )
                     );
      M3(i,j) = MAX ( M3(i-1,j-1) + S3(i,j), I3(i-1,j-1) + S3(i,j) );
    }
  }

  gain_score = MAX( M3(len1,len2) , I3(len1,len2) );
  *p_score += gain_score;
  return gain_score;
}

