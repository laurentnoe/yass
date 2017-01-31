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

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "util.h"
#include "global_var.h"
#include "kword.h"
#include "tuple.h"
#include "assemble.h"
#include "prdyn.h"
#include "proba.h"
#include "display.h"
#include "regroup.h"
#include "threads.h"


static void usage()
{
  SETCOLOR(stderr,GREEN);
  fprintf(stderr, "* Usage :\n");
  RESET(stderr);
  fprintf(stderr,
        "  " PACKAGE_NAME " [options] { file.mfas | file1.mfas file2.mfas }\n");

  fprintf(stderr,
        "      -h       display this Help screen\n");
  fprintf(stderr,
        "      -d <int>    0 : Display alignment positions (kept for compatibility)\n"
        "                  1 : Display alignment positions + alignments + stats (default)\n"
        "                  2 : Display blast-like tabular output\n"
        "                  3 : Display light tabular output (better for post-processing)\n"
        "                  4 : Display BED file output\n"
        "                  5 : Display PSL file output\n");
  fprintf(stderr,
        "      -r <int>    0 : process forward (query) strand\n"
        "                  1 : process Reverse complement strand\n"
        "                  2 : process both forward and Reverse complement strands (default)\n");
  fprintf(stderr,
        "      -o <str> Output file\n");
  fprintf(stderr,
        "      -l       mask Lowercase regions (seed algorithm only)\n");
  fprintf(stderr,
        "      -s <int> Sort according to\n"
        "                  0 : alignment scores\n"
        "                  1 : entropy\n"
        "                  2 : mutual information (experimental)\n"
        "                  3 : both entropy and score\n"
        "                  4 : positions on the 1st file\n"
        "                  5 : positions on the 2nd file\n"
        "                  6 : alignment %% id\n"
        "                  7 : 1st file sequence %% id\n"
        "                  8 : 2nd file sequence %% id\n"
        );
  fprintf(stderr,
        "                  10-18 : (0-8) + sort by \"first fasta-chunk order\"\n"
        "                  20-28 : (0-8) + sort by \"second fasta-chunk order\"\n"
        "                  30-38 : (0-8) + sort by \"both first/second chunks order\"\n");
  fprintf(stderr,
        "                  40-48 : equivalent to (10-18) where \"best score for fst fasta-chunk\" replaces \"...\"\n"
        "                  50-58 : equivalent to (20-28) where \"best score for snd fasta-chunk\" replaces \"...\"\n");
  fprintf(stderr,
        "                  60-68 : equivalent to (70-78), where \"sort by best score of fst fasta-chunk\" replaces \"...\"\n"
        "                  70-78 : BLAST-like behavior : \"keep fst fasta-chunk order\", then trigger all hits for all good snd fasta-chunk\n"
        "                  80-88 : equivalent to (30-38) where \"best score for fst/snd fasta-chunk\"\" replaces \"...\"\n"
        );
  fprintf(stderr,
        "      -v       display the current Version\n");
  fprintf(stderr,
        "                                                                               \n");
  fprintf(stderr,
        "      -M <int> select a scoring Matrix (default 3):\n"
        "               [Match,Transversion,Transition],(Gopen,Gext)\n");
  {
    long int i;
    for (i = 0; i < NBMATRICES; i += 2) {
      if (i == NBMATRICES-1) {
      fprintf(stderr,
            "               %2ld : [%3ld,%3ld,%3ld],(%3ld,%3ld)\n",
            i,
            SUBMATRIXTABLE[i][0],
            SUBMATRIXTABLE[i][1],
            SUBMATRIXTABLE[i][2],
            INDELSTABLE[i][0],
            INDELSTABLE[i][1]
            );
      } else {
      fprintf(stderr,
            "               %2ld : [%3ld,%3ld,%3ld],(%3ld,%3ld)  %2ld : [%3ld,%3ld,%3ld],(%3ld,%3ld)\n",
            i,
            SUBMATRIXTABLE[i][0],
            SUBMATRIXTABLE[i][1],
            SUBMATRIXTABLE[i][2],
            INDELSTABLE[i][0],
            INDELSTABLE[i][1],
            i+1,
            SUBMATRIXTABLE[i+1][0],
            SUBMATRIXTABLE[i+1][1],
            SUBMATRIXTABLE[i+1][2],
            INDELSTABLE[i+1][0],
            INDELSTABLE[i+1][1]
            );
      }
    }
  }
  fprintf(stderr,
        "      -C <int>[,<int>[,<int>[,<int>]]]\n"
        "               reset match/mismatch/transistion/other Costs (penalties)\n"
        "               you can also give the 16 values of matrix (ACGT order)\n");
  fprintf(stderr,
        "      -G <int>,<int> reset Gap opening/extension penalties\n");
  fprintf(stderr,
        "      -L <real>,<real> reset Lambda and K parameters of Gumbel law\n");
  fprintf(stderr,
        "      -X <int>  Xdrop threshold score (default %ld)\n",
        gp_xdrop);
  fprintf(stderr,
        "      -E <int>  E-value threshold (default %1.2g)\n",
        gp_expectation_value);
  fprintf(stderr,
        "      -e <real> low complexity filter :\n"
        "                minimal allowed Entropy of trinucleotide distribution\n"
        "                ranging between 0 (no filter) and 6 (default %1.2f)\n",
        gp_entropy_min);
  fprintf(stderr,
        "                                                                               \n");
  fprintf(stderr,
        "      -O <int> limit number of Output alignments (default %ld)\n",
        gp_nbmaxlines);
  fprintf(stderr,
        "      -S <int> Select sequence from the first multi-fasta file (default %ld)\n"
          "                 * use 0 to select the full first multi-fasta file\n",
        gp_selection_fasta);
  fprintf(stderr,
        "      -T <int> forbid aligning too close regions (e.g. Tandem repeats)\n"
        "               valid for single sequence comparison only (default %ld bp)\n",
        gp_distdiag);
  fprintf(stderr,
        "                                                                               \n");
  fprintf(stderr,
        "      -p <str> seed Pattern(s)\n"
        "                 * use \'#\' for match\n"
        "                 * use \'@\' for match or transition\n"
        "                 * use \'-\' or \'_\' for joker\n"
        "                 * use \',\' for seed separator (max: %ld seeds)\n"
        "                 - example with one seed :\n"
        "                    " PACKAGE_NAME " file.fas -p  \"#@#--##--#-##@#\"\n"
        "                 - example with two complementary seeds :\n"
        "                    " PACKAGE_NAME " file.fas -p \"##-#-#@#-##@-##,##@#--@--##-#--###\"\n",
        (long int) MAX_SEED);
  fprintf(stderr,
        "                 (default  \"");
  {
    long int i;
    for (i = 0; i < gp_nb_seeds; i++)
      fprintf(stderr,"%s%c",gp_motifs[i],(i == (gp_nb_seeds-1))?'"':',');
  }
  fprintf(stderr,")\n");
  fprintf(stderr,
        "      -c <int> seed hit Criterion : 1 or 2 seeds to consider a hit (default %ld)\n",
        gp_hitcriterion);
  fprintf(stderr,
        "                                                                               \n");
  fprintf(stderr,
        "      -t <real> Trim out over-represented seeds codes\n"
        "                ranging between 0.0 (no trim) and +inf (default %1.2g)\n",
        gp_t);

  fprintf(stderr,
        "      -a <int> statistical tolerance Alpha (%%) (default %ld%%)\n",
        gp_alpha_percent);
  fprintf(stderr,
        "      -i <int> Indel rate (%%)                  (default %ld%%)\n",
        gp_indels_percent);
  fprintf(stderr,
        "      -m <int> Mutation rate (%%)               (default %ld%%)\n",
        gp_mutations_percent);
  fprintf(stderr,
        "                                                                               \n");
  fprintf(stderr,
        "      -W <int>,<int> Window <min,max> range for post-processing and grouping\n"
        "                     alignments (default <%ld,%ld>)\n",
        gp_win_min, gp_win_max);
  fprintf(stderr,
        "      -w <real> Window size coefficient for post-processing and grouping\n"
        "                alignments (default %g)\n"
        "                NOTE : -w 0 disables post-processing\n",
        gp_win_mul);

#ifdef MEM_ALLOCATED
  fprintf(stderr,"\n"
        "      -Alloc <int> maximal total Allocation allowed (default %lu)\n",
        gp_max_mem_allocated);
#endif

  RESET(stderr);
  RESET(stdout);
  exit(0);
}


static void vers()
{
    fprintf(stderr, "version " PACKAGE_VERSION "\n");
    exit(0);
}



static void error(char *E1, char *E2, char *E3)
{
    SETCOLOR(stderr, RED);
    fprintf(stderr, "* Error : ");
    RESET(stderr);
    fprintf(stderr,"\" ");
    fprintf(stderr,"%s", E1);
    fprintf(stderr,"%s", E2);
    fprintf(stderr,"%s", E3);
    fprintf(stderr," \"\n");
    usage();
}

/*
 * read several seed patterns separated with ','
 */
static void parsepattern(int argc, char ** argv, long int * i) {
  char * s_begin = NULL;
  char * s_end   = NULL;
  (*i)++;
  if (*i >= argc)
    error("\"",argv[*i-1],"\" found without motif");
  if (!(strchr(argv[*i],'#') || strchr(argv[*i],'@') || strchr(argv[*i],'1') || strchr(argv[*i],'T') || strchr(argv[*i],'X') || strchr(argv[*i],'x') ))
    error("pattern \"",argv[*i],"\"  doesn't have any \"#,1\" / \"@,T,X,x\" symbol");
  /* free previous seeds and set gp_nb_seeds to zero */
  while (gp_nb_seeds > 0) {
    gp_nb_seeds--;
    FREE(gp_motifs[gp_nb_seeds],sizeof(strlen()+1));
  }
  s_begin = argv[*i];
  while (1) {
    if (gp_nb_seeds == MAX_SEED)
      error("pattern \"",argv[*i-1],"\" has too many seeds");
    s_end = strchr(s_begin,',');
    if (s_end) {
      *s_end = '\0';
      gp_motifs[gp_nb_seeds++] = (char *) strdup(s_begin);
      s_begin = s_end + 1;
    } else {
      gp_motifs[gp_nb_seeds++] = (char *) strdup(s_begin);
      break;
    }
  }
  ComputeLengthAndSortSeeds();
}

/*
 * read the output filename
 */
static void parseoutfile(int argc, char ** argv, long int * i) {
  (*i)++;
  if ((*i) >= argc)
    error("\"",argv[(*i)-1],"\" found without filename");
  gv_outstream = fopen(argv[*i],"w");
  if (!gv_outstream)
    error("Can't create output file \"",argv[*i],"\"");
}

/*
 * read one flag
 */
static void parseflag(int argc, char ** argv, long int * i, long int * var, long int val) {
  *var = val;
}

/*
 * read one integer -<Option> <int>
 */
static void parselint(int argc, char ** argv, long int * i, long int * var, long int min, long int max) {
  long int i_tmp;
  (*i)++;
  if ((*i)>=argc)
    error("\"",argv[(*i)-1],"\" found without argument");
  i_tmp = strtol(argv[(*i)],NULL,10);
  if (errno != 0)
    error("\"",argv[(*i)-1],"\" is not followed by an integer or by an integer outside the range");
  if (i_tmp < min || i_tmp  > max)
    error("\"",argv[(*i)-1],"\" is followed by an integer outside the range");
  *var = i_tmp;
}


#ifdef MEM_ALLOCATED

/*
 * read one unsigned long integer -<Option> <int>
 */

static void parseulint(int argc, char ** argv, long int * i, unsigned long int * var, unsigned long int min, unsigned long int max) {
  unsigned long int uli_tmp;
  (*i)++;
  if ((*i)>=argc)
    error("\"",argv[(*i)-1],"\" found without argument");
  uli_tmp = strtoul(argv[(*i)],NULL,10);
  if (errno != 0)
    error("\"",argv[(*i)-1],"\" is not followed by an integer or by an integer outside the range");
  if (uli_tmp < min || uli_tmp > max)
    error("\"",argv[(*i)-1],"\" is followed by an integer outside the range");
  *var = uli_tmp;
}

#endif

/*
 * read two integers -<Option> <int>,<int>
 */
static void parsedlint(int argc, char ** argv, long int * i, long int * var1, long int min1, long int max1, long int * var2, long int min2, long int max2) {
  long int i_tmp1, i_tmp2, err;
  (*i)++;
  if ((*i)>=argc)
    error("\"",argv[(*i)-1],"\" found without argument");
  err = sscanf(argv[(*i)],"%ld,%ld",&i_tmp1,&i_tmp2);
  if (err != 2)
    error("\"",argv[(*i)-1],"\" is not followed by an <int>,<int>");
  if (i_tmp1 < min1 || i_tmp1 > max1 || i_tmp2 < min2 || i_tmp2 > max2)
    error("\"",argv[(*i)-1],"\" is followed by integers  outside the range");
  *var1 = i_tmp1;
  *var2 = i_tmp2;
}


static void parsedouble(int argc, char ** argv, long int * i, double * var, double min, double max) {
  double d_tmp;
  (*i)++;
  if ((*i)>=argc)
    error("\"",argv[(*i)-1],"\" found without argument");
  d_tmp = strtod(argv[(*i)],NULL);
  if (errno != 0)
    error("\"",argv[(*i)-1],"\" is not followed by a double") ;
  if (d_tmp>max || d_tmp <min)
    error("\"",argv[(*i)-1],"\" is followed by a double outside the range");
  *var = d_tmp;
}

static void parseddouble(int argc, char ** argv, long int * i, double * var1, double min1, double max1, double * var2, double min2, double max2) {
  double d_tmp1, d_tmp2;
  int    err;
  (*i)++;
  if (*i>=argc)
    error("\"",argv[(*i)-1],"\" found without argument");
  err = sscanf(argv[*i],"%lf,%lf",&d_tmp1,&d_tmp2);
  if (err != 2)
    error("\"",argv[(*i)-1],"\" is not followed by an <real>,<real>");
  if (d_tmp1 < min1 || d_tmp1 > max1 || d_tmp2 < min2 || d_tmp2 > max2)
    error("\"",argv[*i-1],"\" is followed by integers outside the range");
  *var1 = d_tmp1;
  *var2 = d_tmp2;
}


/*
 * fillmatrixmttrtv
 */

void fillmatrixmttrtv() {
  long int x,y;
  for (x = 0; x < 32; x++) {
    for (y = 0; y < 32; y++) {
      if (((x >= 0 && x < 6) || (x >= 16 && x< 22)) && ((y >= 0 && y < 6) || (y >= 16 && y < 22)))
        if ((x % 16) == (y % 16))
          gp_substitution_matrix[x][y] = (int) gp_costs[0]; /* match */
        else
          if ((x % 2) == (y % 2))
            gp_substitution_matrix[x][y]= (int) gp_costs[2]; /* transition */
          else
            gp_substitution_matrix[x][y]= (int) gp_costs[1]; /* transversion */
      else
        gp_substitution_matrix[x][y] = (int) gp_costs[3]; /* other */

      gp_cost_max_substitution_matrix =  MAX(
                                             ABS(gp_substitution_matrix[x][y]),
                                             gp_cost_max_substitution_matrix
                                             );
      /* symbols U (7) ~ T (3) */
      if (y == 7)
        gp_substitution_matrix[x][y]=gp_substitution_matrix[x][3];
      if (x == 7)
        gp_substitution_matrix[x][y]=gp_substitution_matrix[3][y];
      if (y == 7 && x == 7)
        gp_substitution_matrix[x][y]=gp_substitution_matrix[3][3];
      if (y == 23)
        gp_substitution_matrix[x][y]=gp_substitution_matrix[x][3];
      if (x == 23)
        gp_substitution_matrix[x][y]=gp_substitution_matrix[3][y];
      if (y == 23 && x == 23)
        gp_substitution_matrix[x][y]=gp_substitution_matrix[3][3];
    }
  }
}

int correctmatrix() {
  long int modified = 0;
  long int di,dj;
  for (di = 0 ; di < 32 ; di += 16) {
    for (dj = 0 ; dj < 32 ; dj += 16) {
      long int i,j;
      for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) {
          long int score = (long int)(

            gp_substitution_matrix[i + di][j + dj] * ABS( i == j
                                                ?
                                                log(0.01 + gp_freq_background[i][j])/log(0.0625)
                                                :
                                                log(0.0625)/log(0.01 + gp_freq_background[i][j])
                                          ) + (i == j?0.5:-0.5));

          if (score != gp_substitution_matrix[i + di][j + dj]) {
            modified = 1;
          }
          gp_substitution_matrix[i + di][j + dj] = score ;

          if (i == 3)
            gp_substitution_matrix[7 + di][j + dj] = gp_substitution_matrix[i + di][j + dj];
          if (j == 3)
            gp_substitution_matrix[i + di][7 + dj] = gp_substitution_matrix[i + di][j + dj];
          if (i == 3 && j == 3)
            gp_substitution_matrix[7 + di][7 + dj] = gp_substitution_matrix[i + di][j + dj];
        }
      }
    }
  }
  return modified;
}


/*
 * read matrix substitution costs  -<Option> <int> ( ,<int>  ( ,<int>  ( ,<int> )))
 */

static void parsecost(int argc, char ** argv, long int * i) {
  char * s_begin,* s_end;
  long int  cost_granularity = 0;
  long int  cost_table[16];

  (*i)++;
  if ((*i) >= argc)
    error("\"",argv[(*i)-1],"\" found without argument");

  /* 1) read parameters */
  s_begin = argv[(*i)];
  cost_granularity = 0;

  while (1) {
    if (cost_granularity == 16)
      error(" -\"",argv[*i-1],"\" has too many arguments");

    s_end = strchr(s_begin,',');

    cost_table[cost_granularity] = strtol(s_begin, NULL, 10);
    if (errno != 0)
      error("\"",argv[(*i)-1],"\" is not followed by a list of integers");

    /* update tables */
    if (cost_granularity < 4)
      gp_costs[cost_granularity] = cost_table[cost_granularity];

    cost_granularity++;


    if (s_end)
      s_begin = s_end + 1;
    else
      break;
  }

  if (cost_granularity != 16 && cost_granularity > 4)
    error(" -\"",argv[*i-1],"\" has an incorrect number of arguments");


  /* 2) select two cases :
   *    either the full 16 matrix
   *    or the transition/transversion one
   */

  if (cost_granularity <= 4) {

    /* 2.1) warnings on positives/negatives values */
    if (gp_costs[0] < 0) {
      _WARNING("match costs are converted into positive value");
      gp_costs[0] =   ABS(gp_costs[0]);
    }

    if (gp_costs[1] > 0 || gp_costs[2] > 0 || gp_costs[3] > 0) {
      _WARNING("mismatch costs (transition/transversion/other) are converted into negative values");
      gp_costs[1] = - ABS(gp_costs[1]);
      gp_costs[2] = - ABS(gp_costs[2]);
      gp_costs[3] = - ABS(gp_costs[3]);
    }

    gp_adhoc_matrix = 0;


    /* 2.2) affect costs depending on parameters granularity */
    switch (cost_granularity) {
    case 1: /* only match cost given : other costs are -match */
      gp_costs[1] =  gp_costs[2] = gp_costs[3] = - gp_costs[0];
      break;
    case 2: /* match/mismatch costs : other costs are mismatches */
      gp_costs[2] = gp_costs[3] = gp_costs[1];
      break;
    case 3: /* match/mismatch/transition : other is ~ to mismatch */
      gp_costs[3] = gp_costs[1];
      break;
    }
    /* 2.3  fill matrix */
    fillmatrixmttrtv();

  } else {
    /* cost granularity is 16 */
     long int x,y;
     for (x = 0; x < 32; x++) {
       for (y = 0; y < 32; y++) {
         if (((x >= 0 && x < 6) || (x >= 16 && x < 22)) && ((y >=0 && y < 6) || (y >= 16 && y< 22)))
           gp_substitution_matrix[x][y] = (long int) cost_table[(x%4)*4 + (y%4)]; /* match */
         else
           gp_substitution_matrix[x][y] = (long int) - gp_cost_max_substitution_matrix;

         gp_cost_max_substitution_matrix =  MAX(
                                                ABS(gp_substitution_matrix[x][y]),
                                                gp_cost_max_substitution_matrix
                                                );
         /* symbols U (7) ~ T (3) */
         if (y == 7)
           gp_substitution_matrix[x][y]=gp_substitution_matrix[x][3];
         if (x == 7)
           gp_substitution_matrix[x][y]=gp_substitution_matrix[3][y];
         if (y == 7 && x == 7)
           gp_substitution_matrix[x][y]=gp_substitution_matrix[3][3];
         if (y == 23)
           gp_substitution_matrix[x][y]=gp_substitution_matrix[3][3];
         if (x == 23)
           gp_substitution_matrix[x][y]=gp_substitution_matrix[3][y];
         if (y == 23 && x == 23)
           gp_substitution_matrix[x][y]=gp_substitution_matrix[3][3];
       }
     }
     gp_adhoc_matrix = 1;
  }
}


/*
 * read complete parameters (unless for files)
 */

static void ParseParams(int argc, char ** argv, long int * i) {
  if (!strcmp(argv[*i],"-h") || !strcmp(argv[*i],"--help"))
    usage();
  else if (!strcmp(argv[*i],"-d") || !strcmp(argv[*i],"--display"))
    parselint(argc,argv,i,
             &gp_display,0,5);
  else if (!strcmp(argv[*i],"-r") || !strcmp(argv[*i],"--reverse"))
    parselint(argc,argv,i,
             &gp_reverse,0,2);
  else if (!strcmp(argv[*i],"-l") || !strcmp(argv[*i],"--lowercase"))
    parseflag(argc,argv,i,
              &gp_lowercase,1);
  else if (!strcmp(argv[*i],"-s") || !strcmp(argv[*i],"--sort")) {
    parselint(argc,argv,i,
             &gp_sortcriterion,0,(NBPOTENTIALCRITERIA*NBSORTBLOCKSCRITERIA-1));
    gp_sortblockscriterion = gp_sortcriterion/NBPOTENTIALCRITERIA;
    if (gp_sortcriterion%NBPOTENTIALCRITERIA >= NBSORTCRITERIA)
      error("\"",argv[*i-1],"\" is followed by integers outside the range");
    gp_sortcriterion_func       = sortcriteria[gp_sortcriterion%NBPOTENTIALCRITERIA];
    gp_sortblockscriterion_func = sortblockscriteria[gp_sortcriterion/NBPOTENTIALCRITERIA];
  }
  else if (!strcmp(argv[*i],"-o") || !strcmp(argv[*i],"--output"))
    parseoutfile(argc,argv,i);
  else if (!strcmp(argv[*i],"-T") || !strcmp(argv[*i],"--Tandems"))
    parselint(argc,argv,i,
             &gp_distdiag,0,0x7fffffff);
  else if (!strcmp(argv[*i],"-O") || !strcmp(argv[*i],"--nbO"))
    parselint(argc,argv,i,
             &gp_nbmaxlines,1,0x7fffffff);
  else if (!strcmp(argv[*i],"-S") || !strcmp(argv[*i],"--nbs"))
    parselint(argc,argv,i,
             &gp_selection_fasta,0,0x7fffffff);
  else if (!strcmp(argv[*i],"-M") || !strcmp(argv[*i],"--Matrix")) {
    parselint(argc,argv,i,
             &gp_matrix,0,(NBMATRICES-1));
    gp_costs[0]            = SUBMATRIXTABLE[gp_matrix][0];
    gp_costs[1]            = SUBMATRIXTABLE[gp_matrix][1];
    gp_costs[2]            = SUBMATRIXTABLE[gp_matrix][2];
    gp_costs[3]            = SUBMATRIXTABLE[gp_matrix][3];
    gp_cost_gap_opened     =    INDELSTABLE[gp_matrix][0];
    gp_cost_gap_continued  =    INDELSTABLE[gp_matrix][1];
    gp_adhoc_matrix        = 0;
    fillmatrixmttrtv();
  }
  else if (!strcmp(argv[*i],"-C") || !strcmp(argv[*i],"--Costs")) {
    parsecost(argc,argv,i);
  }
  else if (!strcmp(argv[*i],"-G") || !strcmp(argv[*i],"--Gaps")) {
    parsedlint(argc,argv,i,
              &gp_cost_gap_opened   ,-100000,+100000,
              &gp_cost_gap_continued,-100000,+100000);
    if ( gp_cost_gap_opened > 0 ||  gp_cost_gap_continued > 0) {
      _WARNING("gaps costs are converted into negative values");
      gp_cost_gap_opened     = -ABS(gp_cost_gap_opened);
      gp_cost_gap_continued  = -ABS(gp_cost_gap_continued);
    }

  }

  else if (!strcmp(argv[*i],"-X") || !strcmp(argv[*i],"--Xdrop")) {
    parselint(argc,argv,i,
             &gp_xdrop,0,100);
  }
  else if (!strcmp(argv[*i],"-E") || !strcmp(argv[*i],"--Expectation"))
    parsedouble(argc,argv,i,
                &gp_expectation_value,1e-100,1e+10);
  else if (!strcmp(argv[*i],"-L") || !strcmp(argv[*i],"--Lambda"))
    parseddouble(argc,argv,i,
                 &gp_lambda_blast,1e-9,1e+9,
                 &gp_k_blast,     1e-9,1e+9);
  else if (!strcmp(argv[*i],"-W") || !strcmp(argv[*i],"--Windows"))
    parsedlint(argc,argv,i,
              &gp_win_min,0,0x7fffffff,
              &gp_win_max,0,0x7fffffff);
  else if (!strcmp(argv[*i],"-w") || !strcmp(argv[*i],"--window"))
    parsedouble(argc,argv,i,
                &gp_win_mul,0,1e+16);
  else if (!strcmp(argv[*i],"-p") || !strcmp(argv[*i],"--pattern"))
    parsepattern(argc,argv,i);
  else if (!strcmp(argv[*i],"-c") || !strcmp(argv[*i],"--hitcriterion"))
    parselint(argc,argv,i,
             &gp_hitcriterion,1,2);
  else if (!strcmp(argv[*i],"-e") || !strcmp(argv[*i],"--entropy"))
    parsedouble(argc,argv,i,
                &gp_entropy_min,0.0,6.0);
  else if (!strcmp(argv[*i],"-t") || !strcmp(argv[*i],"--trim"))
    parsedouble(argc,argv,i,
             &gp_t,0.0,1e+100);
  else if (!strcmp(argv[*i],"-v") || !strcmp(argv[*i],"--version"))
    vers();
  else if (!strcmp(argv[*i],"-m") || !strcmp(argv[*i],"--mutations"))
    parselint(argc,argv,i,
             &gp_mutations_percent,0,99);
  else if (!strcmp(argv[*i],"-a") || !strcmp(argv[*i],"--alpha"))
    parselint(argc,argv,i,
              &gp_alpha_percent,1,99);
  else if (!strcmp(argv[*i],"-i") || !strcmp(argv[*i],"--indels"))
    parselint(argc,argv,i,
             &gp_indels_percent,0,100);

#ifdef MEM_ALLOCATED
  else if (!strcmp(argv[*i],"-Alloc"))
    parseulint(argc,argv,i,
               &gp_max_mem_allocated,0,~0);
#endif
  else
    error("\"",argv[*i],"\" is not a valid argument");
}


/*
 * read filenames
 */
static void ParseName(int argc, char ** argv, long int * i) {
  FILE * f = NULL;
  if (gp_nbfiles == 2)
    error("Too many files : \"",argv[*i],"\" is the third one ...");
  strncpy(gp_files[gp_nbfiles],argv[*i],16384);
  if (!(f=fopen(gp_files[gp_nbfiles],"r")))
    error("Can't open file \"",gp_files[gp_nbfiles],"\"");
  fclose(f);
  gp_nbfiles++;
}

/*
 * scan arguments function
 */

static void ScanArg(int argc, char ** argv) {
  long int i;
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      ParseParams(argc,argv,&i);
    } else {
      ParseName(argc,argv,&i);
    }
  }

  /* 1) parameters post-checking */
  if (gp_nbfiles == 0)
    error("No file found"," ... ","give at least one file ...");

  if (gp_lowercase) { /* mask lowercase codes */
    long int i = 0;
    for (i = 16 ; i < 32; i++)
      unindexable[i] = 1;
  }
}



/*
 *   Function to be lauched by threads
 */

static
#ifdef THREAD_FORWARD_REVERSE
#if defined(WIN32) || defined(WIN64)
DWORD WINAPI thread_work_assemble(PVOID fvoid)
#else
     void * thread_work_assemble(void * fvoid)
#endif
#else
     void * thread_work_assemble(void * fvoid)
#endif
{

  /*
   * [5] Assemble work
   */

  /* global variable for one thread */
  Feature * f = (Feature *) fvoid;

  if (gp_nbfiles == 2) {
    MultiAssemble_Double(f->chunk_query, f->chunk_query_size,
                         gp_text, gp_textsize, gp_nbchunks_text,
                         gp_chunkname_text, gp_chunksize_text,
                         gp_chunkstrt_text, f);

  } else {
    if (gp_nbchunks_query > 1) {
      /* FIXME : need  MultiAssemble_SingleRev which does no report comparison on the same diagonal neither the half top matrix */
      MultiAssemble_Double(f->chunk_query, f->chunk_query_size,
                           gp_text, gp_textsize, gp_nbchunks_text,
                           gp_chunkname_text, gp_chunksize_text,
                           gp_chunkstrt_text, f);
    } else {
      if (f->reverse)
      Assemble_SingleRev(f->chunk_query,
                         gp_query + gp_chunkstrt_query[f->j_chunk],
                         f->chunk_query_size, f);
      else
        Assemble_Single   (f->chunk_query,
                           f->chunk_query_size, f);
    }
  }


#ifdef DEBUG_MEMORY
  fprintf (stderr,"Memorized MA size before post-processing :\n");SizeMA (f->first_MA);
#endif

  /*
   * [6]
   */

  /*
   * [6.2] Regrouping MA list
   */
  DISPLAY_BEGIN(regrouping_mutex_begin,regrouping_begin ,"\n\n grouping   ....  ");
  if (f->first_MA)
    Regroup_MAList(f);
  DISPLAY_END(regrouping_mutex_end,regrouping_end);

  /*
   * [6.3] Filtering MA list
   */

  DISPLAY_BEGIN(filtering_mutex_begin, filtering_begin ," filtering   ...  ");
  if (f->first_MA)
    EntropyAndScoreFilter_MAList(f,
                         f->chunk_query,
                         gp_text,gp_chunkstrt_text);
  DISPLAY_END(filtering_mutex_end, filtering_end);

  /*
   * [6.4] Sorting MA list according to score, entropy
   */

  DISPLAY_BEGIN(sorting_mutex_begin,sorting_begin ," sorting   .....  ");
  if (f->first_MA)
    Sort_MAList(f);
  DISPLAY_END(sorting_mutex_end,sorting_end);

  STATS_ADD_CLOCK(f,clock_post);
#ifdef THREAD_FORWARD_REVERSE
  END_THREAD();
#else
  return 0;
#endif
}




/*
 *   Function to be lauched at first by Query threads
 */

static
#ifdef THREAD_QUERY_CHUNK
#if defined(WIN32) || defined(WIN64)
DWORD WINAPI thread_query_chunk(PVOID num)
#else
void * thread_query_chunk(void *num)
#endif
#else
void * thread_query_chunk(void *num)
#endif
{
  long int test = 0;
  long int current_chunk_nb = 0;
  long int q_size, minscore;
  char * q_chunk     = NULL;
  char * q_chunk_rev = NULL;
  long int i = (*((long int *)num));

  do {
    LOCK(query_chunk_mutex);
    if (gv_chunk_nb >= gv_chunk_nb_end) {
      test = 1;
      UNLOCK(query_chunk_mutex);
    } else {
      current_chunk_nb = gv_chunk_nb++;
      UNLOCK(query_chunk_mutex);

      q_size      = gp_chunksize_query[current_chunk_nb];
      q_chunk     = gp_query     + gp_chunkstrt_query[current_chunk_nb];
      q_chunk_rev = gp_query_rev + gp_chunkstrt_query[current_chunk_nb];

      if (q_size < gp_seeds_span_min) {
        _ERROR("selected query chunk of the first file is too small to be indexed, please change the seed with the \'-p \"<seed pattern>\"\' parameter");
      }

      minscore =  MinScore(gp_k_blast, gp_lambda_blast,
                     gp_selection_fasta?q_size:gp_querysize,
                     gp_textsize, gp_expectation_value);

#ifdef DEBUG_DATA
      fprintf(stderr, "querychunk_size: %ld, text_size: %ld\n", q_size, gp_textsize);
      fprintf(stderr, "minscore : %ld\n", minscore);
#endif

      if (gp_reverse != 1) {
        RESETFEATURE(gv_feature[i][0], minscore, q_chunk, q_size, 0, current_chunk_nb);
#ifdef THREAD_FORWARD_REVERSE
        CREATE_THREAD(gv_feature[i][0]->thread_assemble,
                      thread_work_assemble,gv_feature[i][0]);
#else
        thread_work_assemble((void *) gv_feature[i][0]);
#endif
      }


      if (gp_reverse) {
        RESETFEATURE(gv_feature[i][1], minscore, q_chunk_rev, q_size, 1, current_chunk_nb);
#ifdef THREAD_FORWARD_REVERSE
        CREATE_THREAD(gv_feature[i][1]->thread_assemble,
                      thread_work_assemble,gv_feature[i][1]);
#else
        thread_work_assemble((void *) gv_feature[i][1]);
#endif
      }

      /*
       * [7.0] Waiting for threads to end ...
       */

#ifdef THREAD_FORWARD_REVERSE
      if (gp_reverse != 1) {
        WAIT_THREAD(gv_feature[i][0]->thread_assemble);
      }

      if (gp_reverse) {
        WAIT_THREAD(gv_feature[i][1]->thread_assemble);
      }
#endif

      /*
       * [8.0] Merge results
       */

      {
        MA * first_MA = NULL, * last_MA = NULL;
        /* the function merges the forward and reverse part */
        MergeSort_forward_reverse(gv_feature[i][0],gv_feature[i][1],&first_MA,&last_MA);

        /* [FIXME : should be modified for 6x, so this part is subject to some changes] */
        if ((gp_sortblockscriterion >= 5))
          ListSort_MAList(&first_MA,&last_MA,0x40000000,TRUE);

        /* this function locks the main "gv_first/last_MA" list */
        MergeSort_MAList(first_MA,last_MA);
      }
    } /* else ... */
  } while (test == 0); /* while ... */

#ifdef THREAD_QUERY_CHUNK
  END_THREAD();
#else
  return 0;
#endif

}



/************************
   START OF HOSTILITIES
 ************************/

int main(int argc, char *argv[]) {

    int k;

    /*[0] Default parameters */

    /* Init seeds */
    gp_motifs[0] = (char *) strdup("###-#@-##@##");
    gp_motifs[1] = (char *) strdup("###--#-#--#-###");
    gp_seeds_bitweight[0] = gp_seeds_bitweight[1] = 18;
    gp_seeds_span_min = gp_seeds_span[0] = 12;
    gp_seeds_span_max = gp_seeds_span[1] = 15;
    gp_nb_seeds = 2;


    /* Init scoring matrix */
    gp_substitution_matrix = lint_directtable(32,32,0);
    fillmatrixmttrtv();

    /* Init semaphores */
    INIT_MUTEX();
    INIT_QUERY_MUTEX();

    /* [1] command line read */
    ScanArg(argc,argv);

    /* [2] initialize parameters and statistics  */
    gp_mutations = (double) gp_mutations_percent / 100.00;
    gp_indels    = (double) gp_indels_percent    / 100.00;
    gp_alpha     = (double) gp_alpha_percent     / 100.00;
    gp_rho_stat  = statistical_bound_of_waiting_time2(
                                          (1 - gp_mutations),
                                          gp_seeds_bitweight[0] / 2,
                                          gp_alpha
                                          );
    gp_delta_stat = statistical_bound_of_randomwalk2(gp_indels,
                                         gp_rho_stat,
                                         gp_alpha);

    gp_border     = gp_rho_stat + gp_delta_stat + 1;


    /* [3] chaining algorithm memory required */
    initialise_deltashift();


    /* [4] read files and set their attributes */
    if (gp_nbfiles == 2) {
      /* two files */
      CreateData(gp_files[0], &gp_query, &gp_querysize, &gp_nbchunks_query,
             &gp_chunkname_query, &gp_chunksize_query, &gp_chunkstrt_query, gp_nb_letters[0], gp_nb_triplets[0]);

      CreateData(gp_files[1], &gp_text, &gp_textsize, &gp_nbchunks_text,
             &gp_chunkname_text, &gp_chunksize_text, &gp_chunkstrt_text, gp_nb_letters[1], gp_nb_triplets[1]);


    } else {
      /* one file */
      CreateData(gp_files[0], &gp_text, &gp_textsize, &gp_nbchunks_text,
             &gp_chunkname_text, &gp_chunksize_text, &gp_chunkstrt_text, gp_nb_letters[0], gp_nb_triplets[0]);

      gp_query           = gp_text;
      gp_querysize       = gp_textsize;
      gp_nbchunks_query  = gp_nbchunks_text;
      gp_chunkname_query = gp_chunkname_text;
      gp_chunksize_query = gp_chunksize_text;
      gp_chunkstrt_query = gp_chunkstrt_text;

      { /* copy stats to the second file */
      int i;
      for (i = 0; i < 4; i++)
        gp_nb_letters[1][i]  = gp_nb_letters[0][i];
      for (i = 0; i < 64; i++)
        gp_nb_triplets[1][i] = gp_nb_triplets[0][i];
      }
    }

    /* check sizes of the sequences compared to */
    /* first sequence */
    if (gp_nbfiles > 0) {
      if (gp_nbchunks_query == 1) {
        if (gp_chunksize_query[0] < gp_seeds_span_min) {
          _ERROR("first file is too small to be indexed with the current seed. Please change the seed with the \'-p \"<seed pattern>\"\' parameter");
        }
      } else {
        if (gp_nbchunks_query > 1) {
          int i;
          for (i = 0; i < gp_nbchunks_query; i++) {
            if (gp_chunksize_query[i] < gp_seeds_span_min) {
              _ERROR("first file contains too small chunks to be indexed with the current seed. Please change the seed with the \'-p \"<seed pattern>\"\' parameter");
            }
          }
        } else {
          _ERROR("first file is empty or incorrect fasta/multifasta file");
        }
      }
    }

    /* second sequence */
    if (gp_nbfiles > 1) {
      if (gp_nbchunks_text == 1) {
        if (gp_chunksize_text[0] < gp_seeds_span_min) {
          _ERROR("second file is too small to be indexed with the current seed. Please change the seed with the \'-p \"<seed pattern>\"\' parameter");
        }
      } else {
        if (gp_nbchunks_text > 1) {
          int i;
          for (i = 0; i < gp_nbchunks_text; i++) {
            if (gp_chunksize_text[i] < gp_seeds_span_min) {
              _WARNING("second file contains small chunks that will be ignored using the current seed. Please change the seed with the \'-p \"<seed pattern>\"\' parameter\n");
              break;
            }
          }
        } else {
          _ERROR("second file is empty or incorrect fasta/multifasta file");
        }
      }
    }

#ifdef TRACE
    gp_dots = ComputeDotsTable(gp_nbchunks_query,gp_chunksize_query);
#endif

#ifdef DEBUG_COST
    fprintf(stderr, "costs: {match=%d, mis=%d, trans=%d, other=%d} - {gap-open:%d, gap-continue:%d}\n",
          gp_costs[0],gp_costs[1],
          gp_costs[2],gp_costs[3],
          gp_cost_gap_opened,gp_cost_gap_continued
          );
#endif

    /* compute probabilities of single nucleotides
     * and triple nucleotide words in the background
     * ("all genome") model
     */

    /* a) triple letters */
    computeBackgroundTripletFrequency(gp_nb_triplets, gp_freq_tripletbackground);

    /* b) single letter */
    computeLettersFrequency(gp_nb_letters, gp_freq_letters);
    computeBackgroundFrequency(gp_freq_letters, gp_freq_background);

    /* c) correct matrix according to background single letter probabilities */
    if (gp_adhoc_matrix == 0) {
      if (correctmatrix()) {
        _WARNING("substitution matrix was modified due to sequences composition bias");
      }
    }

    /*
     * compute lambda and K
     */

    if (gp_lambda_blast < 0 && gp_k_blast < 0) {
      gp_lambda_blast      = computeLambda(gp_freq_background);
      gp_k_blast           = computeK     (gp_freq_background,gp_lambda_blast);
    }

#ifdef DEBUG_PROBA
    fprintf(stderr, "seq 1: {#A:%d, #T:%d, #G:%d, #C:%d}\n",
          gp_nb_letters[0][0], gp_nb_letters[0][1],
          gp_nb_letters[0][2], gp_nb_letters[0][3]);
    fprintf(stderr, "seq 2: {#A:%d, #T:%d, #G:%d, #C:%d}\n",
          gp_nb_letters[1][0], gp_nb_letters[1][1],
          gp_nb_letters[1][2], gp_nb_letters[1][3]);

    fprintf(stderr, "word  : { ");
    {
      long int i,j;
      for (i = 0; i < 64; i++)
      for (j = 0; j < 64; j++)
        fprintf(stderr,"pr(%c%c%c,%c%c%c):%e\n",
              lookup[(i>>4)&0x3],lookup[(i>>2)&0x3],lookup[(i)&0x3],
              lookup[(j>>4)&0x3],lookup[(j>>2)&0x3],lookup[(j)&0x3],
              gp_freq_tripletbackground[i][j]);
    }
    fprintf(stderr, "Lambda: %f, K: %f \n",gp_lambda_blast,gp_k_blast);
#endif


#ifdef DEBUG_SEED
    DisplaySeeds();
#endif


    if (gp_reverse) {
      CreateReverseComplement(gp_query,
                        gp_querysize,
                        gp_nbchunks_query,
                        gp_chunksize_query,
                        gp_chunkstrt_query,
                        &gp_query_rev);
    }

#ifdef STATS
    gv_time_spent = time(NULL);
#endif

    /* check which query chunk to use */
    if (gp_selection_fasta == 0) {
      gv_chunk_nb = 0;
      gv_chunk_nb_end = gp_nbchunks_query;
      /* if more than 10 chunks ... "-S 0" warning */
      if (gp_nbchunks_query > 10) {
        _WARNING("the -S 0 parameter enumerates all the first file chunks, indexing one at a time : it can be slow !!");
        if (gp_nbchunks_text < gp_nbchunks_query) {
          fprintf(stderr,"  reversing the order of the two files would help...\n");
        }
      }
    } else {
      if (gp_nbchunks_query < gp_selection_fasta) {
        _ERROR("selected fasta chunk does not exist in the first file, please check your -S parameter");
      }
      gv_chunk_nb     = gp_selection_fasta-1;
      gv_chunk_nb_end = gp_selection_fasta;
    }


    for (k = 0; k < MAX_QUERY_CHUNK_THREADS; k++) {
      AllocInitFeature(&gv_feature[k][0]);
      AllocInitFeature(&gv_feature[k][1]);
    }

    /* (a) main part of the program : create k threads on the "thread_query_chunk" function */
#ifdef THREAD_QUERY_CHUNK
    for (k = 0; k < MAX_QUERY_CHUNK_THREADS; k++) {
      gv_thread_num[k] = k;
      CREATE_THREAD(gv_threads[k], thread_query_chunk, (void *)(&(gv_thread_num[k])));
    }
    for (k = 0; k < MAX_QUERY_CHUNK_THREADS; k++) {
      WAIT_THREAD(gv_threads[k]);
    }
#else
    /* (b) otherwise, call the function of one thread by simple call */
    gv_thread_num[0] = 0;
    thread_query_chunk((void *)(&(gv_thread_num[0])));
#endif

    /* [FIXME : should be modified for 6x, so this part is subject to some changes] */
    if (gp_sortblockscriterion == 6)
      gp_sortblockscriterion_func = SortBlocksCriterionQueryNumber;

    if ((gp_sortblockscriterion == 4) || (gp_sortblockscriterion == 6) || (gp_sortblockscriterion == 8))
      ListSort_MAList(&gv_first_MA,&gv_last_MA,0x40000000,TRUE);



#ifdef DEBUG_LISTMA
    DisplayListMA(gv_first_MA);
#endif

    Display_Alignements(gv_first_MA);


#ifdef STATS
    gv_time_spent = time(NULL) - gv_time_spent;
#endif


    Display_Params();
    Display_Stats();

    return 0;
}
