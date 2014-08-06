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
 *Function to perform Monte Carlo 'exact' test for Hardy-
 *Weinberg proportions. The method uses permutation as described
 *by
 *
 *           Engels, 2009, Genetics 183, pp1431-1441
 *
 *         NOTE ABOUT PERMUTATION METHOD
 * This algorithm performs the permutations on
 * 'long long' types, each containing 8 'unsigned char'
 * variables, holding the ID of one allele. Therefore, it is
 * doing 8 permutations at a time. Each of the 8 positions
 * had been previously scrambled independently. So every
 * trial is a combination of one of the original 8 permutations
 * and a new permutation of half the genes. That is, it is
 * the product of two independent permutations even though
 * both of those two permutations are used in other trials.
 * In one sense, that means the trials are not truly 
 * independent, but for practical purposes they are. (The
 * same can be said for pseudorandom numbers in
 * general.) This technique results in a significant
 * increase in the permutation speed which is most
 * beneficial in large samples. When the sample size is
 * much larger than the number of genotypes, the speed
 * increase is about 5-fold. For small samples, it's about
 * 1.5 fold.
 *************************************************************/

/*     NOTES
 
 permutation via MWC takes about half as long as via unif_rand
 when called from R
 (about 3 seconds versus 7 seconds for 100,000,000 swaps)
 However, this might be a smaller fraction of the time when all
 other calculations are taken into account.
 
*/

#define COUNTTYPE unsigned short

// Check to see if long long are 64-bit
#define LONG64 LONG_LONG_MAX > 4294967296

#ifdef NOT_READY_FOR_R    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "HWmtest.h"

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


//This call (MWCranx) generates a random unsigned long named ranx
#define MWCranx  t=a*Q[++ii]+c; c=(t>>32);  ranx=(unsigned)(t+c); if(ranx<c){ranx++;c++;} Q[ii]= ranx;
// Turn it into a random integer in the range 0 to k-1 with (unsigned)((unsigned short)ranx * (k)) >> 16;
#define MWCrank(k) (unsigned)(((unsigned long long)ranx * (unsigned long long)(k)) >> 32)


void mtest (int * m,
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
double * xlnx, *lnFact, *exa, *uTerm1, *uTerm2, *x211, *x221, *x222; // Lookup tables
    int hdex, HN = *histobins;
    int nAlleles = *k;
    int statID = *rstatID;
    double pLLR, pU, pPr, pX2;  // P values
    double maxLLR, maxlPr, minmaxU, minX2; // cutoff values
    double problT, statlT, uT, x2T; // temp values for each trial
    double statSpan, constProbTerm, constLLRterm, leftStat;
    unsigned as, d, aij0, aij1, aij2; // for figuring stats
    double x;
    pU = pLLR = pPr = pX2 =0;
    GetRNGstate();
    
    // Set up MWC
    unsigned Q[256];
    unsigned long long c;
    unsigned char ii;
    c = 362436;
    ii = 255;
    Q[0] = 0;
    unsigned long long t, a = 1540315826LL;
    unsigned ranx;
    // Get "random" Q
    for (int i = 0; i < 256; i++) {
        Q[i] = (unsigned)((double)4294967295 * unif_rand());
    }
    
    //Run off some random MWCs just to mix things up
#if LONG64
        for (int i=0; i<1000; i++) {
            MWCranx
        }
#endif
    
    // Make lookup tables
    unsigned * aij = Calloc(nAlleles * nAlleles, unsigned); // to hold the allele counts;
    xlnx = Calloc(m[0] + 1, double);
    lnFact = Calloc(m[0] + 1, double);
    uTerm1 = Calloc(m[0]/2 + 1, double);
    uTerm2 = Calloc(m[1]/2 + 1, double);
    int biggesta11 = m[0]/2;
    int biggesta22 = m[1]/2;
    int biggesta21 = m[1];
    x211 = Calloc((biggesta11+1), double);
    x222 = Calloc((biggesta22 + 1), double);
    x221 = Calloc((biggesta21+1), double);
    xlnx[0] = 0;
    lnFact[0] = 0;
    double lni;
    for (int i = 1; i <= m[0]; i++) {
        lni = log(i);
        xlnx[i] = lni * i;
        lnFact[i] = lnFact[i-1] + lni;
    }
    for(int i = 0; i <= m[0]/2; i++) uTerm1[i] = (double)i/m[0];
    for(int i = 0; i <= m[1]/2; i++) uTerm2[i] = (double)i/m[1];
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
    for(int i = 0; i <= biggesta11; i++) x211[i] = R_pow_di(exa[0] - i, 2)/exa[0];
    for(int i = 0; i <= biggesta21; i++) x221[i] = R_pow_di(exa[nAlleles] - i, 2)/exa[nAlleles];
    for(int i = 0; i <= biggesta22; i++) x222[i] = R_pow_di(exa[nAlleles + 1] - i, 2)/exa[nAlleles + 1];
    
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
    
    // Set up genes8 array
//    long long *genes8 = Calloc(nGenes, long long);
    int mplex = sizeof(long long)/sizeof(char);  // So we can permute mplex at a time
    // Set up genes array;
    long long * genes = Calloc(nGenes,long long);
    unsigned char * g = (unsigned char *)genes, *gt, tchar;
    for (int j = 0; j < nAlleles; j++) {
        for (int i = 0; i < m[j]; i++) {
            for (int k = 0; k < mplex; k++) {
                g[k] = j;
            }
            g += mplex;
        }
    }
    
    // Scramble the genes a few times for good measure
    long long * gplex, *gtplex, tplex;
    for (int scramble = 0; scramble < 10; scramble++) {
        for (int k = 0; k < mplex; k++) {
            // scramble the kth position
            gplex = genes;
            for (int i = nGenes; i > 0; i--) {
                gtplex = gplex + (int)(unif_rand() * i);
                tchar = *(((unsigned char*)gplex) + k);
                *(((unsigned char*)gplex) + k) = *(((unsigned char*)gtplex) + k);
                *(((unsigned char*)gtplex) + k) = tchar;
                gplex++;
            }
        }
    }
    
    // Perform the actual trials
    time_t start = time(NULL);
    unsigned char * lastgene;
    int trialPlex = (*nTrials)/mplex;
    int otherGene = mplex * ntotal;
    for (int trial = 0; trial < trialPlex; trial++) {
        if (time(NULL) - start >= *safeSecs) {
            *nTrials = trial * mplex; // to signify that the process was stopped
            trial = trialPlex; // To stop the process
            break;
        }
        // Permute one half of the genes in mplex chunks
        gplex = genes;
        for (int i = nGenes; i > ntotal; i--) {
//            gt = g + (int)(unif_rand() * i);
#if LONG64
            MWCranx
            gtplex = gplex + MWCrank(i);
#else
            // Use R built-in RNG as suggested by B. Ripley if long long is 32-bit
            gtplex = gplex+ (int)(unif_rand() * i);
#endif
            tplex = *gplex;
            *gplex = *gtplex;
            *gtplex = tplex;
            gplex++;
        } // permute
        
        for (int k = 0; k < mplex; k++) {
            
            // Build aij array
            memset(aij, 0, nAlleles * nAlleles * sizeof(unsigned));
            g = ((unsigned char *)genes) + k;
            gt = g + otherGene;
            lastgene = (unsigned char *)(genes + nGenes) + k;
            while (gt < lastgene) {
                aij[*g * nAlleles + *gt]++;
                g += mplex;
                gt += mplex;
            }
            // Make aij triangular
            for (int i = 0; i < nAlleles; i++) {
                for (int j = 0; j < i; j++) {
                    aij[i * nAlleles + j] += aij[j * nAlleles + i];
                }
            }
            // Find statistics
            // Get the first 2 rows from lookup tables, etc. This cuts time by maybe 5%
            aij0 = aij[0];
            aij1 = aij[nAlleles];
            aij2 = aij[nAlleles + 1];
            d = aij0 + aij2;
            problT = lnFact[aij0] + lnFact[aij1] + lnFact[aij2];
            statlT = xlnx[aij0] + xlnx[aij1] + xlnx[aij2];
            uT = uTerm1[aij0] + uTerm2[aij2];
            x2T = x211[aij0] + x221[aij1] + x222[aij2];
            // Get remaining rows the slower way
            for (int i = 2; i < nAlleles; i++) {
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
            statlT = constLLRterm - statlT - d * M_LN2;
            
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
            } // if HN
        } // for k
    } // for trial
    
    Pvals[0] = pLLR/ *nTrials;
    Pvals[1] = pPr/ *nTrials;
    Pvals[2] = pU/ *nTrials;
    Pvals[3] = pX2/ *nTrials;

    PutRNGstate();
    Free(xlnx);Free(lnFact);Free(genes);
    Free(exa); Free(uTerm1); Free(uTerm2);
    Free(x211); Free(x221); Free(x222);
    Free(aij);
}

