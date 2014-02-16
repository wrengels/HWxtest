//
//  HWxtest.c
//  ExactoHW
//
//  Created by Bill Engels on 2/11/14.
//
//
/************************************************************
 *                  For R Package `HWxtest`
 *
 * DESCRIPTION
 *
 *Function to perform an exact test for Hardy-
 *Weinberg proportions. The method uses recursion as described
 *by
 *
 *           Engels, 2009, Genetics 183, pp1431-1441
 *
 *************************************************************/

#define COUNTTYPE unsigned short

#ifdef NOT_READY_FOR_R    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "HWxcount.h"

    #define R_pow_di pow  // The following defines convert R calls to normal C
    #define Calloc(x,y) calloc(x,sizeof(y))
    #define Free free
    #define Rprintf printf
    #define lgammafn lgamma
#else
    #include <R.h>
    #include <Rmath.h>
#endif

#include <time.h>

// Globals
COUNTTYPE * Rarray;
unsigned nAlleles, Rbytes;
time_t start;
int timeLimit;
int histobins;

void xtest (int * rm,
            int * rk,
            double * robservedVals, // observed stats: LLR, Prob, U, X2
            double * rPvals, // computed P values: LLR, Prob, U, X2
            int * rstatID, // which statistic to use for histogram (1-4)
            int * rhistobins, // number of bins for histogram. (no histogram if 0)
            double * rhistobounds, // Two values indicating the range for histogram
            double * rhistoData, // histogram data. length = histobounds.
            int * rsafeSecs // abort calculation after this many seconds
            )
{
    // Set up global variables used during recursion
    nAlleles = *rk;
    Rbytes = *rk * sizeof(COUNTTYPE);
    Rarray = Calloc(*rk * *rk * (*rk-1)/2, COUNTTYPE);
    for (int i = 0; i < nAlleles; i++) Rarray[i] = rm[i];
    
    rPvals[0] = 42.4242;
    
}