<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{HWxtest}
-->







<STYLE type="text/css">
  h1,h2,h3,h4,h5 { 
    font-family: palatino, georgia, serif;
    color: Maroon;
  }
    h1, h6{
    text-align: center;
  }
  h1{line-height: 50px}
  body{
    font-size: 0.9em;
    line-height: 23px;
  }
  h3{
  font-weight: normal;
  font-size: 1.8em;
  }
  h6{
        font-size: 0.9em;
        font-weight: normal;
        line-height: 5px;      
   }
   hr{
     border-top-style: solid;
     border-top-width: medium;
   }
  code {
    font-size: 80%;
    line-height: 140%;
    border: 1px solid #ccc;
  }
   @media print{
  hr { 
      visibility: inherit;
      page-break-before: auto;
    }
   }
 </STYLE>


* * *
HWxtest -- Exact Tests for Hardy-Weinberg Proportions With Multiple Alleles
========================================================

###### William R. Engels  --  <wrengels@wisc.edu>  
###### University of Wisconsin, Madison -- Genetics Department
* * *
<center>
* [Main features](#intro)
* [Quick example](#qex)
* [Plot of test statistic](#plot1)
* [Extended example: Bowhead whales](#wex)
* [Inaccuracy in asymptotic tests](#asym)
* [Bibliography](#bib)
</center>

### <a name="intro">Main features</a>
The __HWxtest__ package tests whether a set of genotype counts fits Hardy-Weinberg proportions using the methods described by <a href="http://dx.doi.org/10.1534/genetics.109.108977">Engels (2009)</a>. The package's main features are:
* Performs fast exact test with any number of alleles
* The P value is determined using a choice of test statistics (see below). One option is to test specifically for an excess of homozygotes or heterozygotes.
* Can handle data sets with multiple loci or populations. It accepts data from other population genetics software: *GenePop* (<a href="http://dx.doi.org/10.1111/j.1471-8286.2007.01931.x">Rousset, 2008</a>), *adegenet* (<a href="http://dx.doi.org/10.1093/bioinformatics/btn129">Jombart, 2008</a>), *pegas* (<a href="http://dx.doi.org/10.1093/bioinformatics/btp696">Paradis, 2010</a>)
* Can perform either full enumeration test or Monte Carlo, depending on which is optimal for a given data set.
* Includes functions for determining the number of tables for a given set of allele counts

### <a name="qex">Quick example</a>
Suppose you sample 279 individuals from a population and determine their genotypes at a particular 3-allele locus. If the allele names are $a_1, a_2$ and $a_3$, the genotype counts might be:
$$latex
{\bf{a}} = \left[ {\begin{array}{*{20}{c}}
{{a_{11}}}&{}&{}\\
{{a_{21}}}&{{a_{22}}}&{}\\
{{a_{31}}}&{{a_{32}}}&{{a_{33}}}
\end{array}} \right] = \left[ {\begin{array}{*{20}{c}}
{83}&{}&{}\\
{49}&{18}&{}\\
{74}&{34}&{21}
\end{array}} \right]
$$


That is, we observed 83 homozygotes for allele $a_1$, 49 heterozygotes of genotype $a_1/a_2$, and so on. The genotype counts are in the lower half of a $3 \times 3$ matrix. (These values come from <a href="http://dx.doi.org/10.3354/esr00459">Morin et al. (2012)</a>) We'll enter the numbers by rows in an *R* vector:

```r
obs <- c(83, 49, 18, 74, 34, 21)
```

To test whether these numbers fit the Hardy-Weinberg proportions, call the function `hw.test`

```r
result <- hw.test(obs)
result
```

```
## 
## *****    Sample of 279 diploids with 3 alleles
## Full enumeration of 204350 tables to test for HW
## 
## P value (LLR)   = 0.116908
## P value (Prob)  = 0.098767
## P value (U)     = 0.030499 (test for homozygote excess)
## P value (Chisq) = 0.109703
```

The output tells us that there are 204350 possible samples of 279 diploids with the same allele counts as the observed data. The P values were computed by looking at each of those possibilities and adding up the total probability of all tables which deviate from HW by at least as much as the observed. 

But why does it report *four* different P values? The reason is that there are several ways one can measure just how deviant a particular outcome is compared to the HW expectation. These measures are explained in more detail below. Briefly, if you had no prior expectation that there might be too many homozygotes, then the first (*LLR*) P value is recommended. However, if you had a prior suspicion that homozygotes would be in excess, such as from inbreeding, population admixture, null alleles, *etc.*, then the third choice (*U*) would be prefereable.

### <a name="plot1">Plot a histogram of the test statistic</a>

The `hw.test` function can also show a plot of how the test statistic is actually distributed, as opposed to its theoretical distribution for large samples. To do this, specify the number of bins to use in the histogram:

```r
hw.test(obs, histobins = T)
```

![plot of chunk plot1](figure/plot1.png) 

```
## 
## *****    Sample of 279 diploids with 3 alleles
## Full enumeration of 204350 tables to test for HW
## 
## P value (LLR)   = 0.116908
## P value (Prob)  = 0.098767
## P value (U)     = 0.030499 (test for homozygote excess)
## P value (Chisq) = 0.109703
```

The red area represents those outcomes found to be more extreme than the observed data by the *LLR* criterion. In other words, the red area represents the P value. The blue curve shows the theoretical $\chi^2$ distribution for comparison. In this case the fit is pretty good.

We can do the same thing for the *U* score by specifying the `statName`:

```r
hw.test(obs, histobins = T, detail = 0, statName = "U")
```

![plot of chunk plot2](figure/plot2.png) 

Note that setting the parameter `detail` to zero suppresses the re-printing of the numbers.

### <a name="wex">Example with real data: bowhead whales</a>
A set of data from P. <a href="http://dx.doi.org/10.3354/esr00459">Morin et al. (2012)</a> includes genotypes of 280 whales from one population and 49 from another. Each individual was classified at 51 genetic loci. Some of these loci had more than 20 alleles. We can test for Hardy-Weinberg fit for each locus and each population:

```r
data(whales.genind)
wtest <- hw.test(whales.genind)
```

Note that the format of the data is of class `genind` from the package `adegenet` which must be installed to perform the test. The result includes 102 individual tests, so let's not print them all out. Instead, we'll use the function `hwdf` to put the results in the form of a data frame, and look at just the first 10 in the list:

```r
dfwhales <- hwdf(wtest)
dfwhales[1:10, ]
```

```
##                         P-val(LLR)        ± SE  obs-LLR Method  Tables Trials   N k
## P1.BmC5R700_Y9101        6.276e-01          NA  -0.9044  exact    6160     NA 280 3
## P1.BmCATR205_R212        4.808e-01          NA  -1.5072  exact   15200     NA 280 3
## P1.BmCHYY286_R417        3.695e-01          NA  -2.9063  exact 1292854     NA 278 4
## P1.BmMPOR184_R284        1.169e-01          NA  -2.9735  exact  204350     NA 279 3
## P1.Bmys28Y154_R162       7.816e-01          NA  -1.3733  exact  176208     NA 277 4
## P1.Bmys387R245_R361      7.593e-01          NA  -0.6841  exact   11410     NA 279 3
## P1.Bmys404Y286_K316      3.593e-02          NA  -4.3185  exact  212551     NA 272 3
## P1.Bmys412R79_R463       2.347e-06          NA -14.0552  exact   12901     NA 280 3
## P1.Bmys42aK46_R225_K232  2.000e-05 ± 1.414e-05 -30.0616  monte      NA  1e+05 263 7
## P1.Bmys43Y237_Y377       8.540e-03 ± 2.910e-04  -8.6796  monte      NA  1e+05 277 4
```

The first column shows the population, such as "P1," and locus, such as "BmC5R700_Y9101" of each sample. The second column is the P value from the *LLR* criterion by default. *N* is the number of diploid individuals in the sample and *k* is the number of alleles.

Note that there is a standard error associated with the last two tests, but not the others. This is because those two tests were done by the Monte Carlo method rather than full enumeration. The `hw.test` function decides which method to use based on how many tables would need to be examined in order to do the full enumeration. In these two samples, the number of tables would be 3.8215 &times; 10<sup>22</sup> and 152168317  respectively, and that would take too long. Instead, `hw.test` performed 100,000 random trials to estimate the P value in each case.

If you need a more precise P value than the one provided by 100,000 trials of the Monte Carlo test, there are two ways to do it. One way is to simply increase the Monte Carlo sample size. Suppose we want to do that for the 7-allele case, P1.Bmys42aK46_R225_K232, but not for the whole data set. We can start by pulling out the genotype counts for that specific test. These genotype counts are contained in the `wtest` object and can be picked out as follows:

```r
counts1 <- wtest$P1$Bmys42aK46_R225_K232$genotypes
counts1
```

```
##     001 002 003 004 005 006 007
## 001   1  NA  NA  NA  NA  NA  NA
## 002   1  19  NA  NA  NA  NA  NA
## 003   1  46  18  NA  NA  NA  NA
## 004   1  24  33   3  NA  NA  NA
## 005   1   7  17  10   5  NA  NA
## 006   0  14  19   7  14   4  NA
## 007   0   0   6   0   8   1   3
```

Now we can call `hw.test` again, but increase the sample size to one million by changing the parameter `B`:

```r
hw.test(counts1, detail = 1, B = 1e+06)
```

```
## P value (LLR) = 6e-06 ± 2.4495e-06
```

The standard error of the P value is now considerably smaller. For the other sample, P1.Bmys43Y237_Y377, the number of tables is not as huge, and we can obtain an exact P value by increasing the cutoff number. The following tells `hw.test` to perform the full enumeration test instead of the Monte Carlo unless the number of tables exceeds 2 billion:

```r
counts2 <- wtest$P1$Bmys43Y237_Y377$genotypes
hw.test(counts2, detail = 1, cutoff = 2e+08)
```

```
## P value (LLR) = 0.008886
```

This time, the P value is reported without a standard error, indicating that an exact test with full enumeration was performed.

The results contained in `wtest` includes P values for all four criteria. So, for example, if we were interested in the *U* test rather than the *LLR* we could have specified that in the data frame call:

```r
hwdf(wtest, statName = "U")[1:10, ]
```

```
##                          P-val(U)        ± SE    obs-U Method  Tables Trials   N k
## P1.BmC5R700_Y9101       3.043e-01          NA -14.7089  exact    6160     NA 280 3
## P1.BmCATR205_R212       1.132e-01          NA -29.1880  exact   15200     NA 280 3
## P1.BmCHYY286_R417       3.130e-01          NA   9.2964  exact 1292854     NA 278 4
## P1.BmMPOR184_R284       3.050e-02          NA  43.7794  exact  204350     NA 279 3
## P1.Bmys28Y154_R162      3.124e-01          NA -13.6940  exact  176208     NA 277 4
## P1.Bmys387R245_R361     5.305e-01          NA  -0.8355  exact   11410     NA 279 3
## P1.Bmys404Y286_K316     2.630e-02          NA -45.3983  exact  212551     NA 272 3
## P1.Bmys412R79_R463      1.404e-07          NA -90.6985  exact   12901     NA 280 3
## P1.Bmys42aK46_R225_K232 4.040e-03 ± 0.0002006 128.7425  monte      NA  1e+05 263 7
## P1.Bmys43Y237_Y377      3.198e-01 ± 0.0014749  10.5665  monte      NA  1e+05 277 4
```

When using the *U* test, it is important to pay attention to the sign of *U*. If it is negative, as in the first two samples, the test is specifically for heterozygote excess. If it is positive, then we are testing for homozygote excess. This point will be discussed further below. 

### <a name="asym">Standard asymptotic tests are often inaccurate</a>
Most introductory genetics classes teach that the way to test for Hardy-Weinberg proportions is the standard $\chi^2$ test with $k(k-1)/2$ degrees of freedom. (I am guilty of this myself in teaching Genetics 466 at the University of Wisconsin!) The problem is that with real data and multiple alleles, there are often rare alleles which make the asymptotic test unreliable. 


### <a name="bib">Bibliography</a>

- W. R. Engels,   (2009) Exact Tests For Hardy-Weinberg Proportions.  <em>Genetics</em>  <strong>183</strong>  1431-1441  <a href="http://dx.doi.org/10.1534/genetics.109.108977">10.1534/genetics.109.108977</a>
- F. Rousset,   (2008) genepop'007: a complete re-implementation of the genepop software for Windows and Linux.  <em>Molecular Ecology Resources</em>  <strong>8</strong>  (1)   103-106-NA  <a href="http://dx.doi.org/10.1111/j.1471-8286.2007.01931.x">http://dx.doi.org/10.1111/j.1471-8286.2007.01931.x</a>
- T. Jombart,   (2008) Adegenet: A R Package For The Multivariate Analysis of Genetic Markers.  <em>Bioinformatics</em>  <strong>24</strong>  1403-1405  <a href="http://dx.doi.org/10.1093/bioinformatics/btn129">10.1093/bioinformatics/btn129</a>
- PA Morin, FI Archer, VL Pease, BL Hancock-Hanser, KM Robertson, RM Huebinger, KK Martien, JW Bickham, JC George, LD Postma, BL Taylor,   (2012) Empirical comparison of single nucleotide polymorphisms and microsatellites for population and demographic analyses of bowhead whales.  <em>Endangered Species Research</em>  <strong>19</strong>  (2)   129-147-NA  <a href="http://dx.doi.org/10.3354/esr00459">10.3354/esr00459</a>
- E. Paradis,   (2010) Pegas: an R Package For Population Genetics With an Integrated-Modular Approach.  <em>Bioinformatics</em>  <strong>26</strong>  419-420  <a href="http://dx.doi.org/10.1093/bioinformatics/btp696">10.1093/bioinformatics/btp696</a>

