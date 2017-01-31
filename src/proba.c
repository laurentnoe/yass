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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "global_var.h"
#include "proba.h"



/*****************************
 *
 * Waiting Time Distribution
 *
 *****************************/


/*
 * bound of waiting time
 *
 */

#ifdef INLINE
inline
#endif
long int statistical_bound_of_waiting_time1(double p, long int k, double alpha)
{
    double a1 = dpow(p, k);
    double a2 = (1 - p) * a1;
    long int x = 0;
    double Sum = 0.00;

    double *last_k_prob =
        (double *) MALLOC((k + 1) * sizeof(double));
    ASSERT(last_k_prob, waiting_time_distrib_1);


    while (Sum < (1 - alpha)) {

        if (x < k) {
            last_k_prob[x % (k + 1)] = 0.00;
        } else if (x == k) {
            last_k_prob[x % (k + 1)] = a1;
        } else if (x <= 2 * k) {
            last_k_prob[x % (k + 1)] = a2;
        } else {
            last_k_prob[x % (k + 1)] =
                last_k_prob[(x + k) % (k + 1)] -
                a2 * last_k_prob[x % (k + 1)];
        }
        Sum += last_k_prob[x % (k + 1)];
        x++;
    }
    FREE(last_k_prob, (k + 1) * sizeof(double));

    return x;
}


#ifdef INLINE
inline
#endif
long int statistical_bound_of_waiting_time2(double p, long int k, double alpha)
{
    double p_k = dpow(p, k);
    double qp_k = (1 - p) * p_k;
    long int x = 0;
    double Sum_0_xk1 = 0.00;
    double Sum = 0;

    double *last_k_prob =
        (double *) MALLOC((k + 1) * sizeof(double));
    ASSERT(last_k_prob, waiting_time_distrib_2);

    while (Sum < (1 - alpha)) {
        if (x < k) {
            last_k_prob[x % (k + 1)] = 0.00;
        } else if (x == k) {
            last_k_prob[x % (k + 1)] = p_k;
        } else {
            Sum_0_xk1 += last_k_prob[x % (k + 1)];
            last_k_prob[x % (k + 1)] = qp_k * (1 - Sum_0_xk1);
        }

        Sum += last_k_prob[x % (k + 1)];
        x++;
    }
    FREE(last_k_prob, (k + 1) * sizeof(double));

    return x;
}




/*****************************
 *
 * Random Walk distribution:
 *
 *****************************/


/*
 *
 * table of probabilities to reach
 * position "i" at step L
 * with the rdwalk
 * -> pr pI for left shift
 * -> pr pI for right shift
 * (pr 1-2pI dont move)
 */


#ifdef INLINE
inline
#endif
double *randomwalk_probability_of_pos1(double pI, long int L)
{

    double *C_N_K;              /* Cnk table */
    double *P_for_score;        /* random walk scores */
    double Sum_of_C_N_K = 1;
    double P_dep;

    long int i, nb_dep;

    /* (1) allocation and initialisation */

    C_N_K = (double *) MALLOC((2 * L + 3) * sizeof(double));
    ASSERT(C_N_K, randomwalk_probability_of_pos1);
    for (i = 0; i < 2 * L + 3; i++)
        C_N_K[i] = 0;


    P_for_score =
        (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(P_for_score, randomwalk_probability_of_pos1);
    for (i = 0; i < 2 * L + 1; i++)
        P_for_score[i] = 0;



    /* (2) poly computed
     */

    C_N_K[L] = 1;


    for (nb_dep = 0; nb_dep <= L; nb_dep++) {

      P_dep =
        (C(L, nb_dep) * dpow(pI, nb_dep)) * dpow((1 - pI), (L - nb_dep));


    /***************************************/
        if (P_dep > 1 || 1) {   /* infinity problems ??? */
            goto end;
        }


        /* (2.1) proba computed here */
        for (i = -nb_dep; i <= nb_dep; i++) {
            P_for_score[i + L] += ((C_N_K[i + L]) * P_dep / Sum_of_C_N_K);
        }

        /* (2.2) C_N_K update */
        if (nb_dep < L) {
            for (i = -nb_dep - 1; i <= nb_dep + 1; i += 2) {
                C_N_K[i + L] = C_N_K[i - 1 + L] + C_N_K[i + 1 + L];
            }
            for (i = -nb_dep; i <= nb_dep; i += 2) {
                C_N_K[i + L] = 0;
            }
            Sum_of_C_N_K *= 2;
        }
    }

  end:
    /* (3) free unused stuff */
    FREE(C_N_K, (2 * L + 3) * sizeof(double));

    return P_for_score;
}


#ifdef INLINE
inline
#endif
double *randomwalk_probability_of_pos2(double pI, long int L)
{

    long int i, j;
    double a = pI / 2;
    double b = 1 - pI;
    double *poly, *polytmp;

    poly    = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(poly, randomwalk_probability_of_pos2);
    polytmp = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(polytmp, randomwalk_probability_of_pos2);


    /* (1)  table inits */

    for (i = 0; i < 2 * L + 1; i++) {
        poly[i] = 0;
    }

    /* (2) polynom computed
     * [ p*X^2 + b*X + p ] ^ L
     *
     */

    poly[0] = 1;
    for (i = 0; i < L; i++) {

        /* clear tmp buffer */
        for (j = 0; j < (2 * L + 1); j++) {
            polytmp[j] = 0;
        }

        /* badly have [ p*X^2 + b*X + p ] ^ i+1 from [ p*X^2 + b*X + p ] ^ i */
        for (j = 0; j < (2 * L + 1); j++)
            polytmp[j] += poly[j] * a;
        for (j = 0; j < (2 * L); j++)
            polytmp[j + 1] += poly[j] * b;
        for (j = 0; j < (2 * L - 1); j++)
            polytmp[j + 2] += poly[j] * a;

        /* copy [ p*X^2 + b*X + p ] ^ i+1 into poly[] */
        for (j = 0; j < (2 * L + 1); j++)
            poly[j] = polytmp[j];

    }                           /* i */
    FREE(polytmp, (2 * L + 1) * sizeof(double));

    return poly;

}


#ifdef INLINE
inline
#endif
double *randomwalk_probability_of_pos3(double pI, long int L)
{

    long int i, j;
    long int P = L;
    double a = pI * 0.50;
    double b = 1 - pI;
    double *u, *f, *t, *s;

    u = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(u, randomwalk_probability_of_pos3);
    f = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(f, randomwalk_probability_of_pos3);
    t = (double *) MALLOC((2 * L + 1) * sizeof(double));
    ASSERT(t, randomwalk_probability_of_pos3);

    /* (1) tables inits */

    for (i = 0; i < 2 * L + 1; i++)
        u[i] = 0;
    for (i = 0; i < 2 * L + 1; i++)
        f[i] = 0;

    f[0] = 1.00;
    u[0] = a;
    u[1] = b;
    u[2] = a;

    /* (2) qpow */
    while (P > 0) {
        if (P & 1) {
            /* f = f*u */
            for (i = 0; i < 2 * L + 1; i++) {
                t[i] = 0;
                for (j = 0; j <= i; j++)
                    t[i] += u[i - j] * f[j];
            }
            s = t;
            t = f;
            f = s;
        }
        /* u = u * u; */
        for (i = 0; i < 2 * L + 1; i++) {
            t[i] = 0;
            for (j = 0; j <= i; j++)
                t[i] += u[i - j] * u[j];
        }
        s = t;
        t = u;
        u = s;
        P >>= 1;
    }

    FREE(u, (2 * L + 1) * sizeof(double));
    FREE(t, (2 * L + 1) * sizeof(double));
    return f;

}







/*
 *
 * give rd walk bounds in  (1-alpha)%
 *
 */



#ifdef INLINE
inline
#endif
long int statistical_bound_of_randomwalk1(double pI, long int L, double alpha)
{
    double *RDW_Bound = randomwalk_probability_of_pos1(pI, L);
    double Sum = RDW_Bound[L];
    long int bound = 0;
    while (Sum < (1 - alpha) && bound < L) {
        bound++;
        Sum += RDW_Bound[L - bound] + RDW_Bound[L + bound];
    }

    /* overflow in table */
    if (bound == L)
        bound = (long int) (pow((double) pI * L, (double) (1 / 2.3))) + 1;      /* asympt value */

    FREE(RDW_Bound, (2 * L + 1) * sizeof(double));
    return bound;
}




#ifdef INLINE
inline
#endif
long int statistical_bound_of_randomwalk2(double pI, long int L, double alpha)
{
    double *RDW_Bound = randomwalk_probability_of_pos3(pI, L);
    double Sum = RDW_Bound[L];
    long int bound = 0;
    while (Sum < (1 - alpha) && bound < L) {
        bound++;
        Sum += RDW_Bound[L - bound] + RDW_Bound[L + bound];
    }

    FREE(RDW_Bound, (2 * L + 1) * sizeof(double));
    return bound;
}



/*
 * Scores,Evalue,and Bitscore ...
 *
 */

double
Evalue(double K, double Lambda, long int query_size, long int text_size, long int score)
{
    return K * query_size * text_size * exp(-Lambda * score);
}

double BitScore(double K, double Lambda, long int score)
{
  return ((Lambda * ((double) score)) - log(K)) / log((double)2.0);
}

long int
MinScore(double K, double Lambda, long int query_size, long int text_size,
         double MaxEvalue)
{
  return (long int) floor((1 / Lambda) * (log(K * query_size * text_size) - log(MaxEvalue)));
}



/*
 *  trinomial prob
 */

#ifdef INLINE
inline
#endif
double P_mutation_bias(long int freq[])
{
    long int i;
    long int n = freq[0] + freq[1] + freq[2];
    double result = 1.00;

    if (n == 0)
        return 1.00;

    for (i = 1; i <= n; i++) {
        result *= ((double) i) / 3;
    }
    for (i = 1; i <= freq[0]; i++)
        result /= (double) i;
    for (i = 1; i <= freq[1]; i++)
        result /= (double) i;
    for (i = 1; i <= freq[2]; i++)
        result /= (double) i;

    return result;
}







/*
 *  Compute letters frequency
 */



void computeLettersFrequency(long int nb_letters[2][4], /* out */ double letters_frequency[2][4])
{
  unsigned long int nbL[2];

  nbL[0] =  nb_letters[0][0] + nb_letters[0][1] + nb_letters[0][2] + nb_letters[0][3];
  nbL[1] =  nb_letters[1][0] + nb_letters[1][1] + nb_letters[1][2] + nb_letters[1][3];

  if (nbL[0] == 0) {
    nb_letters[0][0] = nb_letters[0][1] = nb_letters[0][2] = nb_letters[0][3] = 1;
    nbL[0] = 4;
  }

  if (nbL[1] == 0) {
    nb_letters[1][0] = nb_letters[1][1] = nb_letters[1][2] = nb_letters[1][3] = 1;
    nbL[1] = 4;
  }

  /* freq computed */
  {
    long int i,j;
    for(i=0;i<2;i++)
      for(j=0;j<4;j++)
        letters_frequency[i][j] = (double) (nb_letters[i][j]) / (nbL[i]);
  }
}


/*
 *  Compute mutation probability on single nucleotides
 */

void computeBackgroundFrequency(double letters_frequency[2][4], /* out */ double background_frequency[4][4])
{

  /* table allocation */
  {
    long int i,j;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        background_frequency[i][j] =  letters_frequency[0][i] * letters_frequency[1][j];
  }
}





/*
 *  Compute "mutation" probabilty on triple nucleotides words
 */


void computeBackgroundTripletFrequency(long int nb_triplets[2][64], /* out */ double probabilities[64][64])
{

  /* Number of words in each sequence */
  unsigned long int nbL[2] = {0,0};
  double   sizeproduct;

  {
    long int i;
    for(i=0 ; i<64; i++){
      nbL[0] +=  nb_triplets[0][i];
      nbL[1] +=  nb_triplets[1][i];
    }
  }

  sizeproduct = ((double)nbL[0]) * ((double)nbL[1]);

  /* Compute probability to have this triplet conserved on both sequences */
  {
    long int i;
    for(i=0 ; i<64; i++){
      long int j;
      for(j=0 ; j<64; j++){
        (probabilities[i])[j] = (((double)nb_triplets[0][i]) * ((double)nb_triplets[1][j])) / (sizeproduct);
      }
    }
  }
}





/*
 *
 *  Methods for assessing the statistical signifiance of molecular
 *  sequence features by using general scoring schemes, writting by Samuel Karlin
 *  and F. Atschul in 1989-90.
 *
 */




#define COMPUTELAMBDA(lambda)                                                              \
  {                                                                                        \
    long int i,j;                                                                          \
    S = 0.0;                                                                               \
    for (i=0; i<4; i++)                                                                    \
      for (j=0; j<4; j++)                                                                  \
        S += (freq_background[i][j]) * exp(lambda*(double)(gp_substitution_matrix[i][j])); \
  }


/*
 *
 * Compute Lambda
 *
 */


double computeLambda(double freq_background[4][4]) {

  double lambda_lower = 0.0;
  double lambda_upper = 1.0;
  double lambda = 0.0;
  double S = 0.0;

  /*1) check feasability  */
  {
    long int i,j;
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
        S += (freq_background[i][j]) * (double) (gp_substitution_matrix[i][j]);
  }

  if (S >= 0.0) {
    _ERROR("Lambda value cannot be computed : a common reason is a strong AT/GC sequence bias, so please fix it with the \'-L\' option, or change the scoring system with the \'-C\' parameter");
  }

  /*2) looking for lambda upper */
  COMPUTELAMBDA(lambda_upper);
  while (S < 1) {
    lambda_lower = lambda_upper;
    lambda_upper *= 2.0;
    COMPUTELAMBDA(lambda_upper);
  }


  /*3) looking for lambda exactly */
  while (lambda_upper - lambda_lower > SENSITIVITY_LAMBDA) {
    lambda =  ( lambda_upper + lambda_lower ) * 0.5;
    COMPUTELAMBDA(lambda);
    if (S > 1.0) {
      lambda_upper = lambda;
    } else {
      lambda_lower = lambda;
    }
  }
  return lambda;
}


/*
 *
 *   compute K
 *
 */

 double computeK(double freq_background[4][4], double lambda) {

    double *pr = NULL, *pr_new = NULL, *pr_tmp = NULL;
    double denominator = 0.0;
    double numerator = 0.0;
    double Pb = 0.0 , Ek = 0.0 , C = 0.0;
    long int i=0, k=0, scoreMin = 0, scoreMax = 0;
    long int lengthProb, zeroPosition;


    /*
     * Tests values from BLAST:
     */

    /*
     * Blast values that must be checked
     *  -C +1/-3 (and 50 %GC)
     *  lambda = 1.371;
     *  K      = 0.711;
     */


    /* Find minScore,MaxScore to compute the tab length */
    {
      long int i,j;
      for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
          if ( gp_substitution_matrix[i][j] > scoreMax )
            scoreMax = gp_substitution_matrix[i][j];
          if ( gp_substitution_matrix[i][j] < scoreMin)
            scoreMin =  gp_substitution_matrix[i][j];
        }
      }
    }

    lengthProb   = (scoreMax - scoreMin + 1) * NBLOOPS_K + 1;
    zeroPosition = ABS(scoreMin) * NBLOOPS_K;


    /*1) initialise vector probability tables */
    if (lengthProb > 4 * MEGA) {
      _ERROR("compute K is not possible, because the scores given (-C parameter) are too large");
    }

    pr     = (double *)MALLOC(lengthProb*sizeof(double));
    ASSERT(pr,computeK);
    pr_new = (double *)MALLOC(lengthProb*sizeof(double));
    ASSERT(pr_new,computeK);

    for (i = 0; i < lengthProb; i++) {
        pr[i] = 0.0;
        pr_new[i] = 0.0;
    }

    /* Initial step for K = 0 */
    pr[zeroPosition] = 1.0;

    /*
     * 2) Compute denominator lambda*E[ S_1 * e^(lambda*S_1)]
     */
    denominator = 0.0;
    {
      long int i,j;
      for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
          denominator += freq_background[i][j] * gp_substitution_matrix[i][j] * exp(lambda * gp_substitution_matrix[i][j]);
        }
      }
    }
    denominator *= lambda;


    /*
     * 3) numerator 1/k*SUM (E[e^lambda*S_k;S_k<0]+Prob(S_k>=0),k,1,SENSIBILITY_K)
     */

    for (k = 1; k < NBLOOPS_K; k++) {
      long int x;

      /*
       * 3.1) compute a new probability
       */

      for (x = ABS(scoreMin); x < lengthProb - ABS(scoreMax); x++){
         long int i,j;
         for (i = 0; i < 4; i++) {
           for (j = 0; j < 4; j++) {
             pr_new[x +  gp_substitution_matrix[i][j]] += pr[x] * freq_background[i][j];
           }
         }
      }


      /*
       * 3.2) permutation pr and pr_new
       */
      pr_tmp = pr;
      pr = pr_new;
      pr_new = pr_tmp;

      for (i = 0; i < lengthProb; i++)
        pr_new[i] = 0.0;


      Pb = 0.0;         /* Probabilité Prob(S_k) */
      Ek = 0.0;         /* Esperance   E[e^(lambda*S_k);S_k<0 */


      /*
       * 3.4) esperance on the negative side (<0)
       */
      for (i = 0; i < zeroPosition; i++) {
        long int Sk = i - zeroPosition;

        if (exp(lambda * (double)Sk) > DBL_MAX) {
          _ERROR("compute K,  exponent is too large");
        }
        Ek += pr[i] * exp(lambda * (double)Sk);
      }

      /*
       * 3.5) probabilities on the postif sidee (>=0)
       */
      for (i = zeroPosition ; i < lengthProb; i++)
        Pb += pr[i];

      numerator += (Pb + Ek) / (double) k;
    }

    C = exp(-2 * numerator) / (denominator);

    FREE(pr,     lengthProb*sizeof(double));
    FREE(pr_new, lengthProb*sizeof(double));

    return C * (lambda * DELTA_K) / (1 - exp(-lambda * DELTA_K));
}




/*
 * Entropy
 */

double entropyTriplet(long int * count /* [4^3]*/)  {
  double result = 0;
  long int    i,sum  = 0;
  for (i = 0; i < 64; i++) {
    sum += count[i];
  }

  if (sum != 0) {
    for (i = 0; i < 64; i++) {
      if (count[i] != 0) {
        double p = (double) (count[i]) / (double) (sum);
        result -= p * log(p);
      }
    }
    result /= log(2.0);
  }
  return result;
}


/*
 * Mutual information
 */

double mutualInformationTriplet(long int ** count /* [4^3] x [4^3] */)  {
  double result = 0;
  long int    sum  = 0;
  long int * rcount  = count[0];                        /* indirection du to "directtable" */
  double *rprob = gp_freq_tripletbackground[0];    /* indirection du to "directtable" */
  double rprobsum = 0.0;

  {
    long int i;
    for (i = 0; i < 64*64; i++){
      if (rcount[i]){
        sum += rcount[i];
        rprobsum += rprob[i];
      }
    }

  }


  if (sum != 0) {
    long int i;
    for (i = 0; i < 64*64; i++) {
      double q = (double) rcount[i] / sum;

      if (rcount[i] != 0) {
        if ((i>>6) == (i&0x3f))
          result +=  q * log(q);
        else
          result +=  rprob[i] * log(rprob[i]);

      }
    }
  }

  result /= log(2.0);

  return result;
}

