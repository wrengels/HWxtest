//
//  HWxcount.c
//  ExactoHW
//
//  Created by William R. Engels 2014.
//
/*
 Calculate the number of genotype tables for a given set of allele counts
 
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
 *Function to count the number of sets of genotype numbers for
 *a given set of allele numbers. That is, it determines how many
 *tables must be examined to perform an exact test for Hardy-
 *Weinberg proportions. The method uses recursion as described
 *by
 *
 *           Engels, 2009, Genetics 183, pp1431-1441
 *
 *The function `xcount` is designed to be called from R using
 *the .C method. Functions `homozygote` and `heterozygote`
 *perform the recursive traversion of the tree of tables.
 *If there are only two alleles, then the number can be
 *found without recursion by calling `twoAlleleSpecialCase`
 *The `makeHash` function assigned a unique number to each node
 *so that those tables downstream from it do not have to be 
 *re-counted. Only up to `MAXNODE` of them are stored.
 *Recursion depth reaches nAlleles*(nAlleles-1)/2 but the number
 *of nodes can be much larger.
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
#define MAXNODE 2000

// Globals
COUNTTYPE * Rarray;
unsigned nAlleles, Rbytes;
double tableCount, countLimit;
struct node {
    double count;
    unsigned long long hash;
} nodez[MAXNODE];
int nextNode;
time_t start;
int timeLimit;
unsigned long long * hashCoefs;

// Compute a unique "hash" to identify a node, defined by it's level (r) and the residuals (R)
unsigned long long makeHash (unsigned r, COUNTTYPE * R) {
    unsigned long long h = R[0];
    for (int i = 1; i < nAlleles; i++) {
        h += R[i] * hashCoefs[i-1];
    }
    h += r * hashCoefs[nAlleles -1];
    return h;
}


void heterozygote (unsigned r, unsigned c, COUNTTYPE * R);

void homozygote (unsigned r, COUNTTYPE * R)
{
    // If the process takes longer than `timeLimit` seconds, set
    // `tableCount` negative to signify that the job is aborted
    if(tableCount < 0) return;
    if(time(NULL) - start >= timeLimit) tableCount = -tableCount;
    if (tableCount > countLimit) tableCount = -tableCount;
	COUNTTYPE * res, *resn;
	int lower, upper;
	unsigned i, arr;
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
            heterozygote(r, r-1, Rnew);
        }
}

void heterozygote (unsigned r, unsigned c, COUNTTYPE * R)
{
    if(tableCount < 0) return;  // aborted because of time limit
	COUNTTYPE *res, *resn;
	int lower, upper;
	unsigned ntables;
	unsigned i, arc, ar1, ar2, a32, a31;
	COUNTTYPE * Rnew = R + nAlleles;
	double countsSoFar;
    unsigned long long hash;
    //	NSNumber * n;
	
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
        heterozygote(r, c-1, Rnew);
    }
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
            // Before calling homozygote, see if we have visited this node before by comparing its hash tag.
            hash = makeHash(r-1, Rnew);
            i = 0;
            // Search list of old nodes
            while (hash != nodez[i].hash && i < nextNode) i++;
            if(i < nextNode) {
				// old node was found, no need to go any further.
				tableCount += nodez[i].count;
			} else {
				// new node
				countsSoFar =  tableCount;
                homozygote(r-1, Rnew);
                if (nextNode < MAXNODE) {
                    // Make a new node
                    nodez[i].hash = hash;
                    nodez[i].count = tableCount - countsSoFar;
                    nextNode++;
                }
  			} // new node
        }
		if(r==3) // and c = 2, then we can handle a series of two-allele cases with no deeper recursion
		{
			for(a32 = lower; a32 <= upper; a32++) {
				a31 = fmin(res[1], res[3]-a32); //Value of a31 is now fixed for each a32
				ntables = (fmin(res[1] - a31, res[2]-a32))/2 + 1;
				tableCount += ntables;
			}
		} // if r == 3
	} // if c == 2
} // heterozygote

double twoAlleleSpecialCase(int * m)  {
	unsigned ntables = fmin(m[0], m[1])/2 + 1;
	return (double)ntables;
}


void xcount (int * m, int * k, double * count, int * safeSecs) {
    countLimit = *count;
    if (*count < 1) countLimit = 1e100; // Set limit to a huge number unless user asked for a limit
    tableCount = 0;
    // If there are only 2 alleles, that's a special case.
    if (*k==2) {
        *count = twoAlleleSpecialCase(m);
        return;
    }
    
    // Set up global variables used during recursion
    nAlleles = *k;
    Rbytes = *k * sizeof(COUNTTYPE);
    Rarray = Calloc(*k * *k * (*k-1)/2, COUNTTYPE);
    for (int i = 0; i < nAlleles; i++) {
        Rarray[i] = m[i];
    }
    hashCoefs = Calloc(nAlleles, unsigned long long);
    hashCoefs[0] = m[0] + 1;
    for (int i = 1; i < nAlleles; i++) {
        hashCoefs[i] = hashCoefs[i-1] * (m[i] + 1);
    }
    nextNode = 0;
    
    
    timeLimit = *safeSecs;
    start = time(NULL);
    
//    *****************      This is the call to do all the work!
    homozygote(nAlleles, Rarray);
//
    *count = tableCount;
    Free(hashCoefs);
    Free(Rarray);
}
