//
//  HWxtest.c
//  ExactoHW
//
//  Created by Bill Engels on 2/11/14.
//
 /*
   Full Enumeration Test for Hardy-Weinberg Proportions
  
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

#define COUNTTYPE unsigned short

#ifdef NOT_READY_FOR_R    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "HWxtest.h"

    #define R_pow_di pow  // The following defines convert R calls to normal C
    #define Calloc(x,y) calloc((x),sizeof(y))
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
unsigned nAlleles, Rbytes, ntotal;
time_t start;
int timeLimit;
int histobins, HN;
int statID;
int * mi;
double * xlnx, *lnFact, *exa, *uTerm1, *uTerm2, *x211, *x221, *x222; // Lookup tables
double tableCount;
double pLLR, pU, pPr, pX2;  // P values
double maxLLR, maxlPr, minmaxU, minX2; // cutoff values
double statSpan, constProbTerm, constLLRterm, probSum, leftStat;
double * hProb;
double umean, uvariance;

static void heterozygote (unsigned r, unsigned c, double probl, double statl, double u, double x2, COUNTTYPE * R);

static void homozygote (unsigned r, double probl, double statl, double u, double x2, COUNTTYPE * R)
{
    // If the process takes longer than `timeLimit` seconds, set
    // `tableCount` negative to signify that the job is aborted
    if(tableCount < 0) return;
    if(time(NULL) - start >= timeLimit) tableCount = -tableCount;
    
	COUNTTYPE * res, *resn;
	int lower, upper, exindix;
	unsigned i, arr;
    double arrln2;
	COUNTTYPE * Rnew = R + nAlleles;
	memcpy(Rnew, R, Rbytes);
	//Find upper and lower limits for arr.
    res = R-1;  // So res is a 1-based version of R
    resn = Rnew-1; // resn is 1 based for Rnew
	lower = res[r];
	for (i = 1; i <= r-1; i++) lower -= res[i];
    lower = lower < 2 ? 0 : lower/2;
    upper = res[r]/2;
    //For each possible value of arr, examine the heterozygote at r, r-1
    for(arr = lower; arr <= upper; arr++) {
        resn[r] = res[r] - 2*arr;
        arrln2 = arr * M_LN2;
        exindix = (r-1)*nAlleles + r - 1;  // index of homozygote
        heterozygote(r,
                     r-1,
                     probl + lnFact[arr] + arrln2,
                     statl + xlnx[arr] + arrln2,
                     u + (double)arr/mi[r],
                     x2 + R_pow_di(arr - exa[exindix],2)/exa[exindix],
                     Rnew);
    }
}

static void heterozygote (unsigned r, unsigned c, double probl, double statl, double u, double x2, COUNTTYPE * R)
{
    if(tableCount < 0) return;
    COUNTTYPE *res, *resn;
	int lower, upper, exindex;
	unsigned i, arc, ar1, ar2, a32, a31, a11, a21, a22;
	unsigned res1, res2, resTemp, dT;
	int hdex;
	double probl3, statl3, x23, problT, statlT, uT, x2T, prob, x=0;
    COUNTTYPE * Rnew = R + nAlleles;
    
    res = R-1; // to make res a 1-based version of R
    resn = Rnew-1; // so resn is 1-based for Rnew
	lower = res[r];
	for (i = 1; i < c; i++) lower -= res[i];
    lower = fmax(0, lower);
    upper = fmin(res[r], res[c]);
    if(c > 2) for (arc = lower; arc <= upper; arc++) {
        memcpy(Rnew, R, Rbytes); // Put a fresh set of residuals from R into Rnew
        
        // decrement residuals for the current value of arc.
        resn[r] -= arc;
        resn[c] -= arc;
        exindex = (r-1)*nAlleles + c - 1;
        heterozygote(r, c-1,
                     probl+lnFact[arc],
                     statl + xlnx[arc],
                     u,
                     x2 + R_pow_di(arc - exa[exindex], 2)/exa[exindex],
                     Rnew);
    } // for arc
    if(c==2){
		if(r > 3) for (ar2= lower; ar2 <= upper; ar2++) {
			 memcpy(Rnew, R, Rbytes); // Put a fresh set of residuals from R into Rnew
			// decrement residuals for the current value of arc.
			resn[r] -= ar2;
			resn[c] -= ar2;
			// The value of ar1 is now fixed, so no need for any more calls to heterozygote in this row
			ar1 = fmin(resn[r], resn[1]);
			resn[1] -= ar1;
			resn[r] -= ar1;
            exindex = (r-1)*nAlleles;
            homozygote(r-1,
                       probl + lnFact[ar2] + lnFact[ar1],
                       statl + xlnx[ar2] + xlnx[ar1],
                       u,
                        x2 + R_pow_di(ar1 - exa[exindex], 2)/exa[exindex]+ R_pow_di(ar2 - exa[exindex+1], 2)/exa[exindex+1] ,
                       Rnew);
        } // if r > 3
		if(r==3) // and c = 2, then we can handle a series of two-allele cases with no deeper recursion
		{
			double * uT1, *uT2, *x11, *x22;
			for(a32 = lower; a32 <= upper; a32++) {
				a31 = fmin(res[1], res[3]-a32); //Value of a31 is now fixed for each a32
				probl3 = probl + lnFact[a32] + lnFact[a31];
				statl3 = statl + xlnx[a32] + xlnx[a31];
                exindex = 2*nAlleles;
                x23 = x2 + R_pow_di(a31 - exa[exindex], 2)/exa[exindex]+ R_pow_di(a32 - exa[exindex+1], 2)/exa[exindex+1] ;
				// get residual allele counts for two-allele case
				res1 = res[1] - a31;
				res2 = res[2] - a32;
                // make pointers to lookups in case they need to be swapped
				uT1 = uTerm1;
				uT2 = uTerm2;
                x11 = x211;
                x22 = x222;
                
				if(res1 > res2) {            // make sure res1 <= res2. If they need swapping, then swap lookups too
					resTemp = res2;
					res2 = res1;
					res1 = resTemp;
					uT1 = uTerm2;
					uT2 = uTerm1;
                    x11 = x222;
                    x22 = x211;
				}
				
				// Now process two-allele case with allele counts res1 and res2
                tableCount += res1/2 + 1;
                    for(a11 = 0; a11 <= res1/2; a11++) {
					a21 = res1-a11*2; // integer arithmetic rounds down
					a22 = (res2-a21)/2;
					problT = probl3 + lnFact[a11] + lnFact[a21] + lnFact[a22];
					statlT = statl3 + xlnx[a11] + xlnx[a21] + xlnx[a22];
					dT = a11 + a22;
					
					// Here come the actual probability and LLR and X2 and U values
					problT = constProbTerm - problT -dT * M_LN2;
					prob = exp(problT);
					statlT = constLLRterm - statlT - dT * M_LN2;
					uT = 2 * ntotal * (u + uT1[a11] + uT2[a22]) - ntotal;
                    x2T = x23 + x221[a21] + x11[a11] + x22[a22];
                        
//                    umean += prob * uT;
//                    uvariance += prob * uT * uT;
                    
                    //Now process the new values of prob and stat
                    probSum += prob;
                    if(statlT <= maxLLR) pLLR += prob;
                    if(problT <= maxlPr) pPr += prob;
                    if (minmaxU < 0) {
                        if(uT <= minmaxU) pU += prob;
                    } else {
                        if(uT >= minmaxU) pU += prob;
                    }
                    if(x2T >= minX2) pX2 += prob;
                    
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
                            hProb[hdex] += prob;
                        }
                    }
                } // for a11
			} // for a32
		} // if r == 3
	} // if c == 2
}

static void twoAlleleSpecialCase() {
    unsigned a11, a21, a22, res1, res2;
	double problT, statlT, prob, uT, x2T, x;
	unsigned dT;
	int hdex;
	res1 = mi[2]; // because they come ordered largest to smallest
	res2 = mi[1];
	for(a11 = 0; a11 <= res1/2; a11++) {
		a21 = res1-a11*2; // integer arithmetic rounds down
		a22 = (res2-a21)/2;
		problT =  lnFact[a11] + lnFact[a21] + lnFact[a22];
		statlT = xlnx[a11] + xlnx[a21] + xlnx[a22];
		dT =  a11 + a22;
		
		// Here come the actual probability and LLR values
		problT = constProbTerm - problT -dT * M_LN2;
		prob = exp(problT);
		statlT = constLLRterm - statlT - dT * M_LN2;
        x2T = x211[a22] + x221[a21] + x222[a11];
		uT = ntotal *2 * (uTerm1[a22] + uTerm2[a11]) - ntotal;
        
		//Now process the new values of prob and stat
		probSum += prob;
		if(statlT <= maxLLR) pLLR += prob;
        if(problT <= maxlPr) pPr += prob;
        if (minmaxU < 0) {
            if(uT <= minmaxU) pU += prob;
        } else {
            if(uT >= minmaxU) pU += prob;
        }
        if(x2T >= minX2) pX2 += prob;
        
        // Update histogram if needed
        x = statlT;
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
                hProb[hdex] += prob;
            }
        }
        
    } // for a11
}


void xtest (int * rm,
            int * rk,
            double * robservedVals, // observed stats: LLR, Prob, U, X2
            double * rPvals, // computed P values: LLR, Prob, U, X2
            int * rstatID, // which statistic to use for histogram (1-4)
            int * rhistobins, // number of bins for histogram. (no histogram if 0)
            double * rhistobounds, // Two values indicating the range for histogram
            double * rhistoData, // histogram data. length = histobounds.
            int * rsafeSecs, // abort calculation after this many seconds
            double * tables // the number of tables examined
            )
{
    // Set up global variables used during recursion
    nAlleles = *rk;
    pU = pLLR = pPr = pX2 =probSum = 0;
    hProb = rhistoData;
    Rbytes = *rk * sizeof(COUNTTYPE);
    statID = *rstatID;
    timeLimit = *rsafeSecs;
    HN = *rhistobins;
    start = time(NULL);
    Rarray = Calloc(*rk * *rk * (*rk-1)/2, COUNTTYPE);
    for (int i = 0; i < nAlleles; i++) Rarray[i] = rm[i];
    mi = rm-1; // 1-based list of allele counts
    tableCount = 0;
    umean = 0; uvariance = 0;
    
    // Make lookup tables
    xlnx = Calloc(rm[0] + 1, double);
    lnFact = Calloc(rm[0] + 1, double);
    uTerm1 = Calloc(rm[0]/2 + 1, double);
    uTerm2 = Calloc(rm[1]/2 + 1, double);
    int biggesta11 = rm[0]/2;
    int biggesta22 = rm[1]/2;
    int biggesta21 = rm[1];
    x211 = Calloc((biggesta11+1), double);
    x222 = Calloc((biggesta22 + 1), double);
    x221 = Calloc((biggesta21+1), double);
    xlnx[0] = 0;
    lnFact[0] = 0;
    double lni;
    for (int i = 1; i <= rm[0]; i++) {
        lni = log(i);
        xlnx[i] = lni * i;
        lnFact[i] = lnFact[i-1] + lni;
    }
    for(int i = 0; i <= rm[0]/2; i++) uTerm1[i] = (double)i/rm[0];
    for(int i = 0; i <= rm[1]/2; i++) uTerm2[i] = (double)i/rm[1];
    size_t nsq = fmax(2, nAlleles * nAlleles);
    exa = Calloc(nsq, double); // Expected numbers. Array uses extra space but saves time
    int  nGenes = 0;
    for(int i = 0; i < nAlleles; i++) nGenes += rm[i];
    ntotal = nGenes/2;
    for(int i = 0; i < nAlleles; i++) {
        exa[i * nAlleles + i] = (double)(rm[i] * rm[i])/(2.0 * nGenes);
        for (int j = 0; j < i; j++) {
            exa[i * nAlleles + j] = (double)(rm[i] * rm[j])/nGenes;
        }
    }
    for(int i = 0; i <= biggesta11; i++) x211[i] = R_pow_di(exa[0] - i, 2)/exa[0];
    for(int i = 0; i <= biggesta21; i++) x221[i] = R_pow_di(exa[nAlleles] - i, 2)/exa[nAlleles];
    for(int i = 0; i <= biggesta22; i++) x222[i] = R_pow_di(exa[nAlleles + 1] - i, 2)/exa[nAlleles + 1];
    
    // Get constant terms for LLR and Prob
    constProbTerm = constLLRterm = 0;
    for (int i = 0; i < nAlleles; i++) {
        constProbTerm +=  lgammafn(rm[i] + 1); //lnFact[rm[i]];
        constLLRterm += xlnx[rm[i]];
    }
    constProbTerm += log(2)*ntotal + lgammafn(ntotal+1) - lgammafn(nGenes +1);
    constLLRterm += -log(2)*ntotal - log(ntotal) * ntotal;
    
    // Get cutoffs for the four test statistics
    double oneMinus = 0.9999999; // Guards against floating-point-equality-test errors
    if(robservedVals[0] > 0.000000000001) robservedVals[0] = 0; // positive values are rounding errors
    maxLLR = robservedVals[0] * oneMinus;
    maxlPr = log(robservedVals[1]) * oneMinus;
    minmaxU = robservedVals[2] * oneMinus;
    minX2 = robservedVals[3] * oneMinus;
    
    // Set up histogram
    if (HN) {
        switch (*rstatID) {
            case 0: // LLR -- histobounds gives bounds for -2LLR
                leftStat = rhistobounds[0]/(-2.0);
                statSpan = -2.0 * HN/(rhistobounds[1] - rhistobounds[0]);
                break;
            case 1: // Prob -- histobounds gives bounds for -2ln(pr)
                leftStat = rhistobounds[0]/(-2.0);
                statSpan = -2.0 * HN/(rhistobounds[1] - rhistobounds[0]);
                break;
            case 2: // U score  -- histobounds is actual bounds
                leftStat = rhistobounds[0];
                statSpan = (double)HN/(rhistobounds[1] - rhistobounds[0]);
                break;
            case 3: // X2 -- histobounds is actual bounds
                leftStat = rhistobounds[0];
                statSpan = (double)HN/(rhistobounds[1] - rhistobounds[0]);
                break;
            default:
                break;
        }
        hProb = rhistoData;
        for(int i = 0; i < HN; i++) hProb[i] = 0;
    }
    start = time(NULL);
    if (nAlleles == 2) {
        twoAlleleSpecialCase();
    } else {
        homozygote(nAlleles, 0, 0, 0, 0, Rarray);
    }
    
    *tables = tableCount;
    rPvals[0] = pLLR;
    rPvals[1] = pPr;
    rPvals[2] = pU;
    rPvals[3] = pX2;
    if (tableCount < 0) for(int i = 0; i < 4; i++) rPvals[i] = -1; // Process timed out and p values are meaningless
    
//    printf("\nU mean = %.8f", umean);
//    printf("\nU variance = %.8f\n", uvariance - umean * umean);
    
    Free(xlnx);Free(lnFact);Free(Rarray);
    Free(exa); Free(uTerm1); Free(uTerm2);
    Free(x211); Free(x221); Free(x222);
    
}