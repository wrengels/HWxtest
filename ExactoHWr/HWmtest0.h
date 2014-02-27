//
//  HWmtest.h
//  ExactoHW
//
//  Created by Bill Engels on 2/23/14.
//
//

#ifndef ExactoHW_HWmtest_h
#define ExactoHW_HWmtest_h



#endif


void mtest0 (int * m,
            int * k,
            double * observedVals, // observed stats: LLR, Prob, U, X2
            double * Pvals, // computed P values: LLR, Prob, U, X2
            int * rstatID, // which statistic to use for histogram (1-4)
            int * histobins, // number of bins for histogram. (no histogram if 0)
            double * histobounds, // Two values indicating the range for histogram
            double * histoData, // histogram data. length = histobounds.
            int * safeSecs, // abort calculation after this many seconds
            double * nTrials // the number of tables to bo examined
);