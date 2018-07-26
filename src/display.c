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
#include <math.h>
#include <time.h>


/* 1) include utils macro */
#include "util.h"
/* 2) include global variables */
#include "global_var.h"
/* 3) the current file defs */
#include "display.h"
/* 4) other files */
#include "align.h"
#include "prdyn.h"
#include "tuple.h"
#include "proba.h"
#include "kword.h"

void Display_Alignements(MA * first_MA) {

  long int nb_displayed_alignments = 0;
  long int previous_query_chunk    = -1;
  MA *ma = NULL;

  /* [0] header if needed */
  if (gp_display == 2) {
    /* This one is for BLAST tabular files */
    fprintf(OUTSTREAM,
            "# Fields: Query id, Subject id, %% identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n");

  } else if (gp_display == 3) {
    /* This one is mine */
    fprintf(OUTSTREAM,
            "# q. start,\tq. end,\ts. start,\ts. end,\tq. length,\ts.length,\tf/r\tbit score\te-value\n");

  } else if (gp_display == 4) {
    /* This one is for BED files */
    fprintf(OUTSTREAM,"# s. id,\ts. start,\ts. end,\tstrand,\tq. start,\tq. end,\tbit score\te-value\t%%identity\n");
  } else if (gp_display == 5) {
    /* This one is for PSL files */
    /* FIXME : add a header for PSL output (but needed ??) */
  }

  /* [1] Effective display of each Memorized Alignment */
  for (ma = first_MA;
       ma != NULL && nb_displayed_alignments < gp_nbmaxlines;
       ma = ma->next) {

    /* count the output size in number of alignments */
    nb_displayed_alignments++;


#ifdef DEBUG
      DisplayMA(ma);
#endif


      /* output */
      switch (gp_display) {
         case 0:
         case 1:
            /*  (1.1) alignment position repeat */
            if (!(ma->reverse))
                fprintf(OUTSTREAM, "*(%ld-%ld)(%ld-%ld)"  ,
                        BCOUNT + ma->left_pos_begin   ,
                        BCOUNT + ma->left_pos_end  - 1,
                        BCOUNT + ma->right_pos_begin  ,
                        BCOUNT + ma->right_pos_end - 1);

            else {
                fprintf(OUTSTREAM, "*(%ld-%ld)(%ld-%ld)",
                        BCOUNT + gp_chunksize_query[ma->j_chunk] - ma->left_pos_begin - 1,
                        BCOUNT + gp_chunksize_query[ma->j_chunk] - ma->left_pos_end      ,
                        BCOUNT + ma->right_pos_begin               ,
                        BCOUNT + ma->right_pos_end              - 1);
            }

            /* (1.2) print details (size, evalue, etc ...) */
            fprintf(OUTSTREAM, " Ev: %g s: %ld/%ld %c\n",
                    Evalue(gp_k_blast, gp_lambda_blast,
                           gp_selection_fasta?gp_chunksize_query[ma->j_chunk]:gp_querysize,
                           gp_textsize, ma->blastscore),
                    ma->left_pos_end - ma->left_pos_begin   ,
                    ma->right_pos_end - ma->right_pos_begin ,
                    (ma->reverse) ? 'r' : 'f');

            if (gp_display == 0)
              break;


            /* (2.1) print chunk names */
            fprintf(OUTSTREAM, "* \"%s\" (%ld bp) / \"%s\" (%ld bp)\n",
                    gp_chunkname_query[ma->j_chunk], gp_chunksize_query[ma->j_chunk],
                    gp_chunkname_text[ma->i_chunk],  gp_chunksize_text[ma->i_chunk]
                    );
            /* (2.2) score bitscore and statistics */
            fprintf(OUTSTREAM, "* score = %ld : bitscore = %.2f\n",
                    (ma->blastscore), (float) BitScore(gp_k_blast,
                                                       gp_lambda_blast,
                                                       (ma->blastscore)));


            fprintf(OUTSTREAM,
                    "* mutations per triplet %ld, %ld, %ld (%.2e) | ts : %ld tv : %ld | entropy : %g\n\n",
                    ma->trinomial_count1[0], ma->trinomial_count1[1], ma->trinomial_count1[2],
                    (float) (P_mutation_bias(ma->trinomial_count1)), ma->transindels[0],
                    ma->transindels[1],
                    ma->entropy
                    /* , ma->mutual */
                    );


            /* (2.3) display alignment */
            display_alignment_SG_on_MA((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],
                                        gp_text + gp_chunkstrt_text[ma->i_chunk], ma);

#ifdef DEBUG
            DisplayMA(ma);
#endif

            break;              /* case 1: */

        case 2:
            /* (3) blast like tabular format */
            /* (3.1) fasta name */
            fprintf(OUTSTREAM, "%s\t%s\t", gp_chunkname_query[ma->j_chunk],
                    gp_chunkname_text[ma->i_chunk]);


            /* (3.2) alignment percent identity */
            fprintf(OUTSTREAM, "%3.2f\t",
                    100.00 -
                    100.00 *
                    (
                     (double)
                     2.00 * (ma->transindels[0] + ma->transindels[1] + ma->transindels[2])
                     /
                     (
                      (ma->left_pos_end - ma->left_pos_begin) +
                      (ma->right_pos_end - ma->right_pos_begin) +
                      (ma->transindels[2])
                      )
                     )
                    );
            /* (3.3) alignment length */
            fprintf(OUTSTREAM, "%ld\t",
                    (
                     (ma->left_pos_end - ma->left_pos_begin) +
                     (ma->right_pos_end - ma->right_pos_begin) +
                     (ma->transindels[2])
                     )/2
                    );
            /* (3.5) mutations indels */
            fprintf(OUTSTREAM, "%ld\t%ld\t", ma->transindels[0] + ma->transindels[1],
                    ma->transindels[2]);
            /* (3.6) alignment position repeat */
            if (!(ma->reverse))
                fprintf(OUTSTREAM, "%ld\t%ld\t%ld\t%ld\t",
                        BCOUNT + ma->left_pos_begin,
                        BCOUNT + ma->left_pos_end - 1,
                        BCOUNT + ma->right_pos_begin,
                        BCOUNT + ma->right_pos_end - 1);
            else
                fprintf(OUTSTREAM, "%ld\t%ld\t%ld\t%ld\t",
                        BCOUNT + gp_chunksize_query[ma->j_chunk] - ma->left_pos_end      ,  /* inverted pair */
                        BCOUNT + gp_chunksize_query[ma->j_chunk] - ma->left_pos_begin - 1,
                        BCOUNT + ma->right_pos_end - 1,               /* inverted pair : blast database reversed */
                        BCOUNT + ma->right_pos_begin );
            fprintf(OUTSTREAM, "%.2g\t%.3g\n",
                    Evalue(gp_k_blast, gp_lambda_blast,
                           gp_selection_fasta?gp_chunksize_query[ma->j_chunk]:gp_querysize,
                           gp_textsize, ma->blastscore),
                    BitScore(gp_k_blast, gp_lambda_blast,
                             ma->blastscore));
            break;              /* case 2: */
        case 3:
          /*  (4.1) alignment position repeat */
          if (!(ma->reverse))
            fprintf(OUTSTREAM, "%ld\t%ld\t%ld\t%ld",
                    BCOUNT + ma->left_pos_begin,
                    BCOUNT + ma->left_pos_end - 1,
                    BCOUNT + ma->right_pos_begin,
                    BCOUNT + ma->right_pos_end - 1);
          else
            fprintf(OUTSTREAM, "%ld\t%ld\t%ld\t%ld",
                    BCOUNT + gp_chunksize_query[ma->j_chunk] - ma->left_pos_end,       /* inverted pair */
                    BCOUNT + gp_chunksize_query[ma->j_chunk] - ma->left_pos_begin - 1,
                    BCOUNT + ma->right_pos_begin,                /* normal pair */
                    BCOUNT + ma->right_pos_end              - 1);

            /* (4.2) print details (size, f/b , bitscore , evalue ) */
          fprintf(OUTSTREAM, "\t%ld\t%ld\t%c\t%2g\t%3g\n",
                  ma->left_pos_end  - ma->left_pos_begin,
                  ma->right_pos_end - ma->right_pos_begin,
                  (ma->reverse) ? 'r' : 'f',
                  BitScore(gp_k_blast, gp_lambda_blast, ma->blastscore),
                  Evalue(gp_k_blast, gp_lambda_blast,
                         gp_selection_fasta?gp_chunksize_query[ma->j_chunk]:gp_querysize,
                         gp_textsize, ma->blastscore)
                  );
          break;

        case 4:
          /* (5) potential header */
          if (previous_query_chunk != ma->j_chunk) {
            previous_query_chunk = ma->j_chunk;
            fprintf(OUTSTREAM,"track name=\"%s\" description=\"%s\"\n",gp_chunkname_query[ma->j_chunk]/* queryname */,"empty");
          }

          /*  (5.0) bed file name */
          fprintf(OUTSTREAM, "%s\t", gp_chunkname_query[ma->j_chunk] );

          /*  (5.1) alignment position repeat */
          if (!(ma->reverse))
            fprintf(OUTSTREAM, "%ld\t%ld\t%c\t%ld\t%ld\t",
                    ma->left_pos_begin,
                    ma->left_pos_end,
                    '+',
                    ma->right_pos_begin,
                    ma->right_pos_end);
          else
            fprintf(OUTSTREAM, "%ld\t%ld\t%c\t%ld\t%ld\t",
                    gp_chunksize_query[ma->j_chunk] - ma->left_pos_end,       /* inverted pair */
                    gp_chunksize_query[ma->j_chunk] - ma->left_pos_begin,
                    '-',
                    ma->right_pos_begin,                /* normal pair */
                    ma->right_pos_end);

          /* (5.2) bitscore / evalue */
          fprintf(OUTSTREAM, "%.3g\t%.2g\t",
                  BitScore(gp_k_blast, gp_lambda_blast,
                           ma->blastscore),
                  Evalue(gp_k_blast, gp_lambda_blast,
                         gp_selection_fasta?gp_chunksize_query[ma->j_chunk]:gp_querysize,
                         gp_textsize, ma->blastscore)
                  );


          /* (5.3) percent identity */
          fprintf(OUTSTREAM, "%3.2f\n",
                  100.00 -
                  100.00 *
                  (
                   (double)
                   2.00 * (ma->transindels[0] + ma->transindels[1] + ma->transindels[2])
                   /
                   (
                    (ma->left_pos_end - ma->left_pos_begin) +
                    (ma->right_pos_end - ma->right_pos_begin) +
                    (ma->transindels[2])
                    )
                   )
                  );
          break;
      case 5:
        /* (6) potential header */
        if (previous_query_chunk != ma->j_chunk) {
          previous_query_chunk = ma->j_chunk;
          fprintf(OUTSTREAM,"track name=\"%s\" description=\"%s\"\n",gp_chunkname_query[ma->j_chunk]/* queryname */,"empty");
        }
        /* (6.1) psl : base matches */
        fprintf(OUTSTREAM,"%ld\t",
                (
                 (ma->left_pos_end - ma->left_pos_begin) +
                 (ma->right_pos_end - ma->right_pos_begin) -
                 (ma->transindels[2])
                 )/2 - (
                 (ma->transindels[0] +  ma->transindels[1])
                 )
                );
        /* (6.1) base mismatches */
        fprintf(OUTSTREAM,"%ld\t",(ma->transindels[0] +  ma->transindels[1]));
        /* (6.2) repmatches, N count */
        fprintf(OUTSTREAM,"%ld\t%ld\t",
                lowercaseCount((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],  ma->left_pos_begin, ma->left_pos_end - ma->left_pos_begin + 1),
                nCount((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],  ma->left_pos_begin,  ma->left_pos_end - ma->left_pos_begin + 1));
        /* (6.3) nb indels blocks / indels in the query */
        fprintf(OUTSTREAM,"%ld\t%ld\t",ma->transindels[5],ma->transindels[3]);
        /* (6.4) nb indels blocks / indels in the text */
        fprintf(OUTSTREAM,"%ld\t%ld\t",ma->transindels[6],ma->transindels[4]);
        /* (6.5) strands (-+) of (++) */
        fprintf(OUTSTREAM,"%c\t",ma->reverse?'-':'+');
        /* (6.6) qName - qSize - qAlignStart - qAlignEnd */
        if (!(ma->reverse))
          fprintf(OUTSTREAM, "%s\t%ld\t%ld\t%ld\t",
                  gp_chunkname_query[ma->j_chunk],
                  gp_selection_fasta?gp_chunksize_query[ma->j_chunk]:gp_querysize,
                  ma->left_pos_begin,
                  ma->left_pos_end);
        else
          fprintf(OUTSTREAM, "%s\t%ld\t%ld\t%ld\t",
                  gp_chunkname_query[ma->j_chunk],
                  gp_selection_fasta?gp_chunksize_query[ma->j_chunk]:gp_querysize,
                  gp_chunksize_query[ma->j_chunk] - ma->left_pos_end,
                  gp_chunksize_query[ma->j_chunk] - ma->left_pos_begin
                  );

        /* (6.7) tName - tSize - tAlignStart - tAlignEnd */
        fprintf(OUTSTREAM, "%s\t%ld\t%ld\t%ld\t",
                gp_chunkname_text[ma->i_chunk],
                gp_chunksize_text[ma->i_chunk],
                ma->right_pos_begin,
                ma->right_pos_end);

        /* (6.8) block count */
        alignment_SG_PSL_on_MA((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],
                               gp_text + gp_chunkstrt_text[ma->i_chunk],
                               0, ma, 4);
        fprintf(OUTSTREAM, "\t");

        /* (6.9) block sizes (list with separator ',' even at the end) */
        alignment_SG_PSL_on_MA((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],
                               gp_text + gp_chunkstrt_text[ma->i_chunk],
                               0, ma, 3);
        fprintf(OUTSTREAM, "\t");

        /* (6.10) qStarts (list with separator ',' even at the end) */

        /* SAME "BUG" AS IN PSL : dont ask me why ...
         * bugs are programs specific features in some softs :
         * they are documented but not solved
         *
         * see http://genome.ucsc.edu/FAQ/FAQformat#format2 ,
         * the "Be aware that" paragraph ...
         *
         * :o)
         */

        alignment_SG_PSL_on_MA((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],
                               gp_text + gp_chunkstrt_text[ma->i_chunk],
                               0, ma, 1);
        fprintf(OUTSTREAM, "\t");




        /* (6.10) tStarts (list with separator ',' even at the end) */
        alignment_SG_PSL_on_MA((ma->reverse?gp_query_rev:gp_query) + gp_chunkstrt_query[ma->j_chunk],
                               gp_text + gp_chunkstrt_text[ma->i_chunk],
                               0, ma, 2);
        fprintf(OUTSTREAM, "\n");

        break;


        } /* switch */
    } /* for */
}




void Display_Stats() {
    fprintf(stderr, "* Statistical parameters\n");
    fprintf(stderr, "  lambda = %g\n", gp_lambda_blast);
    fprintf(stderr, "  K      = %g\n", gp_k_blast);
    fprintf(stderr, "  rho    = %ld\n", gp_rho_stat);
    fprintf(stderr, "  delta  = %ld\n", gp_delta_stat);

#ifdef STATS
    {
      int k;
      for (k = 0; k < ((gp_selection_fasta == 0)?MAX_QUERY_CHUNK_THREADS:1); k++) {

        Feature ** feature = gv_feature[k];

        if (feature[0]->nb_seeds + feature[1]->nb_seeds) {

          fprintf(stderr, "  ----\n");
          /* preprocess */
          fprintf(stderr, "  # keys removed     =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_keys_removed,
                  feature[1]->nb_keys_removed,
                  feature[0]->nb_keys_removed + feature[1]->nb_keys_removed);
          /* hits */
          fprintf(stderr, "  -\n");
          fprintf(stderr, "  # single hit       =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_seeds,
                  feature[1]->nb_seeds,
                  feature[0]->nb_seeds + feature[1]->nb_seeds);
          if (gp_hitcriterion == 1 ) {
            fprintf(stderr, "  # single hit tests =   f:%10ld   r:%10ld   total:%10ld \n",
                    feature[0]->nb_single_tests,
                    feature[1]->nb_single_tests,
                    feature[0]->nb_single_tests + feature[1]->nb_single_tests);
            fprintf(stderr, "  # single hit ext.s =   f:%10ld   r:%10ld   total:%10ld \n",
                    feature[0]->nb_single_hits,
                    feature[1]->nb_single_hits,
                    feature[0]->nb_single_hits + feature[1]->nb_single_hits);
          }
          fprintf(stderr, "  # chains built     =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_chains_built ,
                  feature[1]->nb_chains_built ,
                  feature[0]->nb_chains_built + feature[1]->nb_chains_built
                  );
          fprintf(stderr, "  # align. detected  =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_ma ,
                  feature[1]->nb_ma ,
                  feature[0]->nb_ma + feature[1]->nb_ma
                  );
          /* regroup */
          fprintf(stderr, "  -\n");
          fprintf(stderr, "  # grouping tests   =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_postprocessed_grouping_tests,
                  feature[1]->nb_postprocessed_grouping_tests,
                  feature[0]->nb_postprocessed_grouping_tests + feature[1]->nb_postprocessed_grouping_tests
                  );
          fprintf(stderr, "  # align. linked    =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_postprocessed_grouping_links,
                  feature[1]->nb_postprocessed_grouping_links,
                  feature[0]->nb_postprocessed_grouping_links + feature[1]->nb_postprocessed_grouping_links
                  );
          fprintf(stderr, "  # align. postpros  =   f:%10ld   r:%10ld   total:%10ld \n",
                  feature[0]->nb_postprocessed_ma,
                  feature[1]->nb_postprocessed_ma,
                  feature[0]->nb_postprocessed_ma + feature[1]->nb_postprocessed_ma
                  );

          /* cpu-times */
          fprintf(stderr, "  ----\n");
          fprintf(stderr, "  pre-proc. cpu-time =   f:%9.2fs   r:%9.2fs   total:%9.2fs [%ld%%]\n",
                  (double) feature[0]->clock_pre  /CLOCKS_PER_SEC,
                  (double) feature[1]->clock_pre  /CLOCKS_PER_SEC,
                  (double)(feature[0]->clock_pre   + feature[1]->clock_pre  )/CLOCKS_PER_SEC,
                  (long int)(
                             (feature[0]->clock_pre   + feature[1]->clock_pre  )/
                             (feature[0]->clock_pre   + feature[1]->clock_pre   +
                              feature[0]->clock_chain + feature[1]->clock_chain +
                              feature[0]->clock_align + feature[1]->clock_align +
                              feature[0]->clock_post  + feature[1]->clock_post  + 1e-6)*100)
                  );
          fprintf(stderr, "  %s cpu-time =   f:%9.2fs   r:%9.2fs   total:%9.2fs [%ld%%]\n",


#ifdef THREAD_ASSEMBLE_ALIGN
                  (gp_hitcriterion == 1)?"ch./h.(+)":"chain.(+)",
#else
                  (gp_hitcriterion == 1)?"ch./hit. ":"chain.   ",
#endif
                  (double)feature[0]->clock_chain/CLOCKS_PER_SEC,
                  (double)feature[1]->clock_chain/CLOCKS_PER_SEC,
                  (double)(feature[0]->clock_chain + feature[1]->clock_chain)/CLOCKS_PER_SEC,
                  (long int)(
                             (feature[0]->clock_chain + feature[1]->clock_chain)/
                             (feature[0]->clock_pre   + feature[1]->clock_pre   +
                              feature[0]->clock_chain + feature[1]->clock_chain +
                              feature[0]->clock_align + feature[1]->clock_align +
                              feature[0]->clock_post  + feature[1]->clock_post  + 1e-6)*100)
                  );
          fprintf(stderr,
#ifdef THREAD_ASSEMBLE_ALIGN
                  "  align.(-) cpu-time =   f:%9.2fs   r:%9.2fs   total:%9.2fs [%ld%%]\n",
#else
                  "  align.    cpu-time =   f:%9.2fs   r:%9.2fs   total:%9.2fs [%ld%%]\n",
#endif
                  (double)feature[0]->clock_align/CLOCKS_PER_SEC,
                  (double)feature[1]->clock_align/CLOCKS_PER_SEC,
                  (double)(feature[0]->clock_align + feature[1]->clock_align)/CLOCKS_PER_SEC,
                  (long int)(
                             (feature[0]->clock_align + feature[1]->clock_align)/
                             (feature[0]->clock_pre   + feature[1]->clock_pre   +
                              feature[0]->clock_chain + feature[1]->clock_chain +
                              feature[0]->clock_align + feature[1]->clock_align +
                              feature[0]->clock_post  + feature[1]->clock_post   + 1e-6)*100)
                  );
          fprintf(stderr, "  post-pro. cpu-time =   f:%9.2fs   r:%9.2fs   total:%9.2fs [%ld%%]\n",
                  (double)feature[0]->clock_post/CLOCKS_PER_SEC,
                  (double)feature[1]->clock_post/CLOCKS_PER_SEC,
                  (double)(feature[0]->clock_post + feature[1]->clock_post)/CLOCKS_PER_SEC,
                  (long int)(
                             (feature[0]->clock_post  + feature[1]->clock_post )/
                             (feature[0]->clock_pre   + feature[1]->clock_pre   +
                              feature[0]->clock_chain + feature[1]->clock_chain +
                              feature[0]->clock_align + feature[1]->clock_align +
                              feature[0]->clock_post  + feature[1]->clock_post   + 1e-6)*100)
                  );
#ifdef DEBUG
          DisplayHistoScore(feature);
#endif
        } /* if (working thread) */
      } /* for (k < MAX_THREADS) */
    } /* int k; */
    fprintf(stderr, "  ----\n");
    fprintf(stderr, "     real-time spent = %ds\n", (unsigned int)gv_time_spent);
#endif
}




void Display_Params() {
  long int i;
  long int querysize_chunk = 0;
  fprintf(stderr, "\n");
  fprintf(stderr, RULE);

  if (gp_selection_fasta > 0)
    querysize_chunk = gp_chunksize_query[gp_selection_fasta-1];

  /* filenames */
  if (gp_nbfiles <= 1) {
    fprintf(stderr, "  DNA file : ");
    if (gp_selection_fasta == 0) {
      fprintf(stderr, "\"%s\", %ld chunk(s), (overall size : %ld bp, GC: %2.2f%%)\n", gp_files[0], gp_nbchunks_query, gp_querysize, 100*(gp_freq_letters[0][1/*C*/] +gp_freq_letters[0][2/*G*/]));
    } else {
      fprintf(stderr, "\"%s\", %ld chunk(s), selected chunk: %ld (size : %ld bp, GC: %2.2f%%)\n", gp_files[0], gp_nbchunks_query, gp_selection_fasta, querysize_chunk,  100*(gp_freq_letters[0][1/*C*/] +gp_freq_letters[0][2/*G*/]));
    }
  } else {
    fprintf(stderr, "  first DNA file : ");
    if (gp_selection_fasta == 0) {
      fprintf(stderr, "\"%s\", %ld chunk(s), (overall size : %ld bp, GC: %2.2f%%)\n", gp_files[0], gp_nbchunks_query, gp_querysize,  100*(gp_freq_letters[0][1/*C*/] +gp_freq_letters[0][2/*G*/]));
    } else {
      fprintf(stderr, "\"%s\", %ld chunk(s), selected chunk: %ld (size : %ld bp, GC: %2.2f%%)\n", gp_files[0], gp_nbchunks_query, gp_selection_fasta, querysize_chunk,  100*(gp_freq_letters[0][1/*C*/] +gp_freq_letters[0][2/*G*/]));
    }

    fprintf(stderr, "  second DNA file : ");
    fprintf(stderr, "\"%s\"  %ld chunk(s), (overall size : %ld bp,  GC: %2.2f%%)\n\n", gp_files[1], gp_nbchunks_text, gp_textsize, 100*(gp_freq_letters[1][1/*C*/] +gp_freq_letters[1][2/*G*/]));
  }

  /* parameters */
  fprintf(stderr, "* Parameters\n");
  for (i = 0; i < gp_nb_seeds; i++)
    fprintf(stderr,
            "  seed[%ld]: pattern = \"%s\", weight = %ld bits, span = %ld bases\n",
            i, gp_motifs[i], gp_seeds_bitweight[i], gp_seeds_span[i]
            );
  fprintf(stderr,
          " original substitution costs (match,mismatch,transition,other) = %ld,%ld,%ld,%ld\n",
          gp_costs[0], gp_costs[1], gp_costs[2], gp_costs[3]
            );
  DisplaySubstitutionMatrix();
  fprintf(stderr,
          "  gap costs (open,extend) = %ld,%ld\n",
          gp_cost_gap_opened, gp_cost_gap_continued
          );
  fprintf(stderr,
          "  x-drop cost = %ld\n",
          gp_xdrop
          );
  fprintf(stderr,
          "  mutation rate = %ld%%\n",
          gp_mutations_percent
          );
  fprintf(stderr,
          "  indels rate = %ld%%\n",
          gp_indels_percent
          );
  fprintf(stderr,
          "  tolerance  = %ld%%\n",
          gp_alpha_percent
          );
}



/* Dot progress bar :
 * "launch and lets take a coffee"(TM)  progress bar
 */

void Display_Progress(long int pos, long int maxpos, Feature  *feature)
{

  if (gp_selection_fasta) {
    /* one query chunk */
    if (gp_reverse == 2) {
      while ( (double) feature->last_point * maxpos <(double) pos * (NBL / 2) && (feature->last_point < (NBL / 2)) ) {
        fprintf(stderr, ".");
        fflush(NULL);
        (feature->last_point)++;
      }
    } else {
      while ( (double) feature->last_point * maxpos <(double) pos * (NBL) && (feature->last_point < (NBL)) ) {
        fprintf(stderr, ".");
        fflush(NULL);
        (feature->last_point)++;
      }
    }
  } else {
    /* multiple query chunks */
    while ( (double) feature->last_point * maxpos <(double) pos * (gp_dots[feature->j_chunk]) && feature->last_point < gp_dots[feature->j_chunk]) {
      fprintf(stderr, ".");
      fflush(NULL);
      (feature->last_point)++;
    }
  }
}


void DisplayHistoScore(Feature ** feature /* feature[2] */) {
  long int i;

  fprintf(stderr," Scoring Histogram : \n");
  for (i = 0; i < sizeof(long int)*8-2; i++) {



    long int lines_f = 20;
    long int lines_r = 20;

    fprintf(stderr,"\t score 2^%2ld to 2^%2ld \t",i,i+1);

    if (feature[0]->MAcount) {
      long int count_f = feature[0]->MAcount[i];
      while (count_f > 0) {
        count_f >>= 1;
        lines_f--;
        fprintf(stderr,".");
      }
      while (lines_f > 0) {
        lines_f--;
        fprintf(stderr," ");
      }
      fprintf(stderr,"(%8ld)",feature[0]->MAcount[i]);
    }else {
      fprintf(stderr,
              "                                  "
              );
    }

    if (feature[1]->MAcount) {
      long int count_r = feature[1]->MAcount[i];
      while (count_r > 0) {
        count_r >>= 1;
        lines_r--;
        fprintf(stderr,".");
      }
      while (lines_r > 0) {
        lines_r--;
        fprintf(stderr," ");
      }
      fprintf(stderr,"(%8ld)\n",feature[1]->MAcount[i]);
    }else {
      fprintf(stderr,
              "                                  "
              );
    }
  }
}



void DisplaySubstitutionMatrix() {
  long int i,j;
  for (i = 0; i < 4; i++) {
    fprintf(stderr,"\t");
    for (j = 0; j < 4; j++) {
      fprintf(stderr,"%2ld",gp_substitution_matrix[i][j]);
      if (j == 3)
        fprintf(stderr,"\n");
      else
        fprintf(stderr,",");
    }
  }
}
