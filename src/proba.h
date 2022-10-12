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

#ifndef __PROBA_H_
#define __PROBA_H_

long int statistical_bound_of_waiting_time1(double p, long int k, double alpha);
long int statistical_bound_of_waiting_time2(double p, long int k, double alpha);
double *randomwalk_probability_of_pos1(double pI, long int L);
double *randomwalk_probability_of_pos2(double pI, long int L);
double *randomwalk_probability_of_pos3(double pI, long int L);
long int statistical_bound_of_randomwalk1(double pI, long int L, double alpha);
long int statistical_bound_of_randomwalk2(double pI, long int L, double alpha);
double Evalue(double K, double Lambda, long int query_size, long int text_size,
            long int score);
long int MinScore(double K, double Lambda, long int query_size, long int text_size,
           double MaxEvalue);
double BitScore(double K, double Lambda, long int score);
double P_mutation_bias(long int freq[]);

void computeLettersFrequency(long int nb_letters[2][4], /* out */ double letters_frequency[2][4]);
void computeBackgroundFrequency(double letters_frequency[2][4], /* out */ double background_frequency[4][4]);
void computeBackgroundTripletFrequency(long int nb_triplets[2][64], /* out */ double probabilities[64][64]);

double computeLambda(double freq_background[4][4]);
double computeK(double freq_background[4][4], double lambda);

double entropyTriplet(long int * count /* [4^3]*/);
double mutualInformationTriplet(long int ** count /* [4^3] x [4^3] */);

#endif
