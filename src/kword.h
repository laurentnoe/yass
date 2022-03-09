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

#ifndef __KWORD_H_
#define __KWORD_H_

#include "threads.h"

/* compute the keyvalue,and the pend, given the pos,data and seed number  */
#define KEY(__keyvalue__,__pend__,__data__,__pos__,__seed__) {              \
    long int _u_;                                                           \
    long int _j_;                                                           \
    __pend__     = __pos__;                                                 \
    __keyvalue__ = 0;                                                       \
    for ( _j_ = 0 ; _j_ < gp_seeds_span[__seed__] ; __pend__++,_j_++ ) {    \
        long unsigned _data_element_ = (long unsigned)(__data__ [__pend__]);\
        char          _seed_element_ = gp_motifs[__seed__][_j_];            \
        if  (unindexable[_data_element_]) {                                 \
            __keyvalue__ =  -1;                                             \
            goto key_end;                                                   \
        }                                                                   \
        switch (_seed_element_) {                                           \
          case 'N' :                                                        \
          case '#' :                                                        \
              __keyvalue__ <<= 2;                                           \
              __keyvalue__  |= (backcode[_data_element_]);                  \
           break;                                                           \
           case '@' :                                                       \
              __keyvalue__ <<= 1;                                           \
              __keyvalue__  |= (backcode[_data_element_] & 0x1);            \
           break;                                                           \
           case 'R' :                                                       \
             if ((backcode[_data_element_] & 0x1) == 0x0) { /* purine */    \
               __keyvalue__ <<= 1;                                          \
               __keyvalue__  |= ((backcode[_data_element_] >> 1) & 0x1);    \
             } else {                                                       \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'Y' :                                                       \
             if ((backcode[_data_element_] & 0x1) == 0x1) { /* pyrimid. */  \
               __keyvalue__ <<= 1;                                          \
               __keyvalue__  |= ((backcode[_data_element_] >> 1) & 0x1);    \
             } else {                                                       \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'S' :                                                       \
             _u_  = backcode[_data_element_] >> 1;                          \
             _u_ ^= backcode[_data_element_] & 0x1;                         \
             if (_u_ == 0x1) { /* strong */                                 \
               __keyvalue__ <<= 1;                                          \
               __keyvalue__  |= (backcode[_data_element_] & 0x1);           \
             } else {                                                       \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'W' :                                                       \
             _u_  = backcode[_data_element_] >> 1;                          \
             _u_ ^= backcode[_data_element_] & 0x1;                         \
             if (_u_ == 0x0) { /* weak */                                   \
               __keyvalue__ <<= 1;                                          \
               __keyvalue__  |= (backcode[_data_element_] & 0x1);           \
             } else {                                                       \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'M' :                                                       \
             if ((backcode[_data_element_] & 0x2) == 0x0) { /* amino */     \
               __keyvalue__ <<= 1;                                          \
               __keyvalue__  |= (backcode[_data_element_] & 0x1);           \
             } else {                                                       \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'K' :                                                       \
             if ((backcode[_data_element_] & 0x2) == 0x2) { /* keto */      \
               __keyvalue__ <<= 1;                                          \
               __keyvalue__  |= (backcode[_data_element_] & 0x1);           \
             } else {                                                       \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'r' :                                                       \
             if ((backcode[_data_element_] & 0x1) == 0x1) { /* not purine */\
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'y' :                                                       \
             if ((backcode[_data_element_] & 0x1) == 0x0) { /* not pyrimi.*/\
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 's' :                                                       \
             _u_  = backcode[_data_element_] >> 1;                          \
             _u_ ^= backcode[_data_element_] & 0x1;                         \
             if (_u_ == 0x0) { /* not strong */                             \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'w' :                                                       \
             _u_  = backcode[_data_element_] >> 1;                          \
             _u_ ^= backcode[_data_element_] & 0x1;                         \
             if (_u_ == 0x1) { /* not weak */                               \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'm' :                                                       \
             if ((backcode[_data_element_] & 0x2) == 0x2) { /* not amino */ \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'k' :                                                       \
             if ((backcode[_data_element_] & 0x2) == 0x0) { /* not keto */  \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'a' :                                                       \
           case 'A' :                                                       \
             if (backcode[_data_element_] != 0x0) { /* A */                 \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'c' :                                                       \
           case 'C' :                                                       \
             if (backcode[_data_element_] != 0x1) { /* C */                 \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 'g' :                                                       \
           case 'G' :                                                       \
             if (backcode[_data_element_] != 0x2) { /* G */                 \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
           case 't' :                                                       \
           case 'T' :                                                       \
             if (backcode[_data_element_] != 0x3) { /* T */                 \
               __keyvalue__   = -1;                                         \
               goto key_end;                                                \
             }                                                              \
           break;                                                           \
       }                                                                    \
    }                                                                       \
  key_end:;                                                                 \
  }


/* Read multifasta sequence */
long int CreateData(/* in */  char *filename,
                    /* out */ char **p_data,          /* out */ long int *p_datasize,
                    /* out */ long int *p_nbchunks,   /* out */ char ***p_chunkname,
                    /* out */ long int **p_chunksize, /* out */ long int **p_chunkstart,
                    /* out */ long int * nbletters  /* [4]  */ ,
                    /* out */ long int * nbtriplets /* [64] */ );

/* Sequence reverse complement */
long int CreateReverseComplement(
                                 /* in  */ char *   data,
                                 /* in  */ long int datasize,
                                 /* in  */ long int nbchunks,
                                 /* in  */ long int * chunksize,
                                 /* in  */ long int * chunkstart,
                                 /* out */ char ** p_data_rev_out);

/* Compute the number of dots per query chunk */
long int * ComputeDotsTable(long int nbchunks, long int * chunksize);

/* Build the Index : Some part of it is now also in feature ... */
long int CreateKeyList( /* in */  char *data,         /* in  */ long int datasize,      /* in */ long int seed,
                        /* out */ long int **p_list,  /* out */ long int * p_list_size, /* out */ Feature * feature);

long int ComputeLengthAndSortSeeds();

long int lowercaseCount(char *data, long int pos, long int len);
long int nCount(char *data, long int pos, long int len);

double Entropy3mer(char *data, long int pos, long int len);

/* Debugging */
long int DisplaySequence(char *data, long int pos, long int len);
long int DisplayQuery();
long int DisplayText();

long int DisplayData(char * filename, char * data, long int size, long int nbchunks, char ** chunkname, long int * chunksize, long int * chunkstrt);

long int DisplaySeeds();

#endif
