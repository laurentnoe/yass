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

#ifndef __KWORD_H_
#define __KWORD_H_

#include "threads.h"

/* compute the keyvalue,and the pend, given the pos,data and seed number  */
#define KEY(__keyvalue__,__pend__,__data__,__pos__,__seed__) {              \
    long int _j_;                                                           \
    __pend__     = __pos__;                                                 \
    __keyvalue__ = 0;                                                       \
    for ( _j_ = 0 ; _j_ < gp_seeds_span[__seed__] ; __pend__++,_j_++ ) {    \
        long unsigned _data_element_ = (long unsigned)(__data__ [__pend__]);\
        char          _seed_element_ = gp_motifs[__seed__][_j_];            \
        if  (unindexable[_data_element_]) {                                 \
            __keyvalue__ =  -1;                                             \
            break;                                                          \
        }                                                                   \
        switch (_seed_element_) {                                           \
           case '#' :                                                       \
              __keyvalue__ <<= 2;                                           \
              __keyvalue__  |= (backcode[_data_element_]);                  \
           break;                                                           \
           case '@' :                                                       \
              __keyvalue__ <<= 1;                                           \
              __keyvalue__  |= (backcode[_data_element_]&0x1);              \
           break;                                                           \
        }                                                                   \
    }                                                                       \
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
