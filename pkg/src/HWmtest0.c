//
//  HWmtest.c
//  ExactoHW
//
//  Created by Bill Engels on 2/23/14.
//
//

//
//  HWxtest.c
//  ExactoHW
//
//  Created by Bill Engels on 2/11/14.
//
/*
    Perform Monte Carlo test for Hardy-Weinberg
 
 Copyright 2009-2014  William R. Engels
 
 This file used in computing an exact test for Hardy-Weinberg
 frequencies. It is part of the 'HWxtest' package for R
 and is made available under the terms of the GNU General Public
 License, version 2, or at your option, any later version,
 incorporated herein by reference.
 
 This program is distributed in the hope that it will be
 useful, but WITHOUT ANY WARRANTY; without even the implied
 warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the GNU General Public License for more
 details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, a copy is available at
 http://www.r-project.org/Licenses/
 
 */


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

/*     NOTES
 
 permutation via MWC takes about half as long as via unif_rand
 when called from R
 (about 3 seconds versus 7 seconds for 100,000,000 swaps)
 However, this might be a smaller fraction of the time when all
 other calculations are taken into account.
 
*/

#define COUNTTYPE unsigned short

#ifdef NOT_READY_FOR_R    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "HWmtest0.h"

    #define R_pow_di pow  // The following defines convert R calls to normal C
    #define Calloc(x,y) calloc((x),sizeof(y))
    #define Free free
    #define Rprintf printf
    #define lgammafn lgamma
    #define unif_rand drand48
    #define GetRNGstate()
    #define PutRNGstate()
#else
    #include <R.h>
    #include <Rmath.h>
#endif

#include <time.h>

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
) {
double * xlnx, *lnFact, *exa; // Lookup tables
    int hdex, HN = *histobins;
    int nAlleles = *k;
    int statID = *rstatID;
    double pLLR, pU, pPr, pX2;  // P values
    double maxLLR, maxlPr, minmaxU, minX2; // cutoff values
    double problT, statlT, uT, x2T, prob; // temp values for each trial
    double statSpan, constProbTerm, constLLRterm, probSum, leftStat;
    unsigned as, d; // for figuring stats
    double x;
    pU = pLLR = pPr = pX2 =0;
    GetRNGstate();
    
    // Make lookup tables
    unsigned * aij = Calloc(nAlleles * nAlleles, unsigned); // to hold the allele counts;
    xlnx = Calloc(m[0] + 1, double);
    lnFact = Calloc(m[0] + 1, double);
    xlnx[0] = 0;
    lnFact[0] = 0;
    double lni;
    for (int i = 1; i <= m[0]; i++) {
        lni = log(i);
        xlnx[i] = lni * i;
        lnFact[i] = lnFact[i-1] + lni;
    }
    size_t nsq = fmax(2, nAlleles * nAlleles);
    exa = Calloc(nsq, double); // Expected numbers. Array uses extra space but saves time
    int  nGenes = 0;
    for(int i = 0; i < nAlleles; i++) nGenes += m[i];
    int ntotal = nGenes/2;
    for(int i = 0; i < nAlleles; i++) {
        exa[i * nAlleles + i] = (double)(m[i] * m[i])/(2.0 * nGenes);
        for (int j = 0; j < i; j++) {
            exa[i * nAlleles + j] = (double)(m[i] * m[j])/nGenes;
        }
    }
    
    // Get constant terms for LLR and Prob
    constProbTerm = constLLRterm = 0;
    for (int i = 0; i < nAlleles; i++) {
        constProbTerm +=  lgammafn(m[i] + 1); //lnFact[m[i]];
        constLLRterm += xlnx[m[i]];
    }
    constProbTerm += log(2)*ntotal + lgammafn(ntotal+1) - lgammafn(nGenes +1);
    constLLRterm += -log(2)*ntotal - log(ntotal) * ntotal;
    
    // Get cutoffs for the four test statistics
    double oneMinus = 0.9999999; // Guards against floating-point-equality-test errors
    if(observedVals[0] > 0.000000000001) observedVals[0] = 0; // positive values are rounding errors
    maxLLR = observedVals[0] * oneMinus;
    maxlPr = log(observedVals[1]) * oneMinus;
    minmaxU = observedVals[2] * oneMinus;
    minX2 = observedVals[3] * oneMinus;

    // Set up histogram
    if (HN) {
        for(int i = 0; i < HN; i++) histoData[i] = 0;
        switch (*rstatID) {
            case 0: // LLR -- histobounds gives bounds for -2LLR
                leftStat = histobounds[0]/(-2.0);
                statSpan = -2.0 * HN/(histobounds[1] - histobounds[0]);
                break;
            case 1: // Prob -- histobounds gives bounds for -2ln(pr)
                leftStat = histobounds[0]/(-2.0);
                statSpan = -2.0 * HN/(histobounds[1] - histobounds[0]);
                break;
            case 2: // U score  -- histobounds is actual bounds
                leftStat = histobounds[0];
                statSpan = (double)HN/(histobounds[1] - histobounds[0]);
                break;
            case 3: // X2 -- histobounds is actual bounds
                leftStat = histobounds[0];
                statSpan = (double)HN/(histobounds[1] - histobounds[0]);
                break;
            default:
                break;
        }
    }
    
    // Set up genes array;
    unsigned char * genes = Calloc(nGenes, unsigned char);
    unsigned char * g = genes, *gt, tchar;
    for (int j = 0; j < nAlleles; j++) {
        for (int i = 0; i < m[j]; i++) {
            *(g++) = j;
        }
    }
    
    // Scramble the genes a few times for good measure
    for (int scramble = 0; scramble < 10; scramble++) {
        g = genes;
        for (int i = nGenes; i > 0; i--) {
            gt = g + (int)(unif_rand() * i);
            tchar = *g;
            *g = *gt;
            *gt = tchar;
            g++;
        }
    }
    
    // Perform the actual trials
    unsigned char * lastgene;
    for (int trial = 0; trial < *nTrials; trial++) {
        
        // Permute one half of the genes
        g = genes;
        for (int i = nGenes; i > ntotal; i--) {
            gt = g + (int)(unif_rand() * i);
            tchar = *g;
            *g = *gt;
            *gt = tchar;
            g++;
        } // permute
        
        // Build aij array
        memset(aij, 0, nAlleles * nAlleles * sizeof(unsigned));
        g = genes;
        gt = genes + ntotal;
        lastgene = genes + nGenes;
        while (gt < lastgene) {
            aij[*g * nAlleles + *gt]++;
            g++; gt++;
        }
        
        // Make aij triangular
        for (int i = 0; i < nAlleles; i++) {
            for (int j = 0; j < i; j++) {
                aij[i * nAlleles + j] += aij[j * nAlleles + i];
            }
        }
        
        // Find statistics
        uT = 0; problT = 0; statlT = 0; x2T = 0;
        d = 0;
        for (int i = 0; i < nAlleles; i++) {
            as = aij[i * nAlleles + i];
            d += as;  // total homozygotes
            uT += (double)as/m[i];
            for (int j = 0; j <= i; j++) {
                as = aij[i * nAlleles + j];
                problT += lnFact[as];
                statlT += xlnx[as];
                x = exa[i * nAlleles + j];
                x2T += R_pow_di(as - x,2)/x;
            } // for j
        } // for i
        uT = nGenes * uT - ntotal;
        problT = constProbTerm - problT - d * M_LN2;
        prob = exp(problT);
        statlT = constLLRterm - statlT - d * M_LN2;
  
//        Rprintf("\n\nx2T = %.5f\n", x2T);
//        for (int i = 0; i < nAlleles; i++) {
//            Rprintf("\n");
//            for (int j = 0; j < nAlleles; j++) {
//                Rprintf("%d\t", aij[i * nAlleles + j]);
//            }
//        }
        
        //Now process the new values of prob and stat
		if(statlT <= maxLLR) pLLR += 1;
        if(problT <= maxlPr) pPr += 1;
        if (minmaxU < 0) {
            if(uT <= minmaxU) pU += 1;
        } else {
            if(uT >= minmaxU) pU += 1;
        }
        if(x2T >= minX2) pX2 += 1;
        
        // Update histogram if needed
        if (HN) {
            switch (statID) {
                case 0:
                    x = statlT;
                    break;
                case 1:
                    x = problT;
                    break;
                case 2:
                    x = uT;
                    break;
                case 3:
                    x = x2T;
                default:
                    break;
            }
            hdex = statSpan * (x - leftStat);
            if ((hdex >= 0) && (hdex < HN)) {
                histoData[hdex] += 1;
            }
        }
} // for trial
    
    Pvals[0] = pLLR/ *nTrials;
    Pvals[1] = pPr/ *nTrials;
    Pvals[2] = pU/ *nTrials;
    Pvals[3] = pX2/ *nTrials;

    PutRNGstate();
    Free(xlnx);Free(lnFact);Free(genes);
    Free(exa);
    Free(aij);
}

