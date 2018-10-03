---
title: "Randomisation Inference - worksheet"
author: Sebastián Martínez
date: "2018-10-03"
showtoc: false
---

Randomisation Inference - Worksheet
===========================

This document looks to understand the examples in the Randomisation Inference package by [Samii and Aronow](https://cran.r-project.org/web/packages/ri/ri.pdf). 

## Initial examples

We are first looking at the example in the `dispdist` help section. The description of the function is for "estimated ATE distribution display, summary and significance testing"

We define the observed outcomes, *y*<sub>*i*</sub>, and the randomisation from the intervention, *Z*<sub>*i*</sub>.

``` r
# Estimated Results
y <- c(8,6,2,0,3,1,1,1,2,2,0,1,0,2,2,4,1,1)
Z <- c(1,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0)
df <- data.frame(y, Z)
```

An initial exploration of the data shows that individuals who received treatment have a larger distribution of outcomes than those that did not. Some jitter is added to be able to differentiate some points with the same *y* value.

``` r
g <- ggplot(data = df, aes(x = Z, y = y)) + geom_jitter(width = 0.1, height = 0) + scale_x_discrete(limits = c(0,1), labels = c("Control", "Treatment"), name ="Treatment type")
```
{{< figure src="/resources/treat-scatter.png" >}}


We now what to divide the observations further into clusters and blocks.

``` r
cl <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)
bl <- c(rep(1,4),rep(2,6),rep(3,8))
cluster <- as.character(cl)
block <- as.character(bl)

df <- cbind(df, cl, bl, cluster, block)
```

Let us see what happens when we divide the observations by clusters (color), and by block (shape).

``` r
g <- ggplot(data = df, aes(x = Z, y = y)) + geom_jitter(aes(color = cluster, shape = block), width = 0.1, height = 0) + scale_x_discrete(limits = c(0,1), labels = c("Control", "Treatment"), name ="Treatment type")
g
```

![](/resources/unnamed-chunk-4-1.png)

We are going to generate all the different permutations available considering the clusters and blocks the observations belong to.

Comparing the permutation sets when considering only the blocks, only the clusters, or both variations we can see how the sample gets reduced when we impose more restrictions.

| Blocks | Clusters | Both |
|--------|----------|------|
| 6300   | 126      | 36   |

Now we generate the probabilities of treatment considering the clusters and the blocks for observation.

``` r
probs_cl <- genprobexact(Z, clustvar = cluster)
probs_bl <- genprobexact(Z, blockvar = block)
probs <- genprobexact(Z, blockvar = block, clustvar = cluster)
```

The variable `probs` contains the probability of being part of the treatment group for each cluster/block group. Internally, the function makes use of a calculation of how many different permutations for each cluster-block combination are available given the observed treatment *Z*. If there is no cluster variable, each observation is assumed to be in its own cluster, similarly with the blocks.

The function `estate` estimate the average treatment effect given an outcome variable (in our case, *y*), and a treatment variable (in our case *Z*), considering the probability of treatment assignment

``` r
# ATE assuming equal probability
ate_noprobs <- estate(y, Z)
```

    ## Warning in estate(y, Z): Probabilities not specified. Assuming equal
    ## probabilities.

``` r
# ATE assuming only blocks
ate_bl <- estate(y, Z, prob = probs_bl)
# ATE assuming only cluster
ate_cl <- estate(y, Z, prob = probs_cl)
# ATE assuming both clusters and blocks
ate <- estate(y, Z, prob = probs)
```

| Equal Prob | Blocks | Clusters | Both |
|------------|--------|----------|------|
| 1.9        | 2      | 1.9      | 2    |

We are now interested in comparing the observed outcome variable *y* with a simulated outcome variable *y*<sub>*s*</sub> which assumes there is no effect after being treated (the 'sharp null' case). \`\`

``` r
# Calculating the assumed outcome under the given ATE, in this case, a sharp null is assumed.
# ys$Y0 stores the values that answer: what would have happened to the TREATED      observations had they NOT  been treated
# ys$Y1 stores the values that answer: what would have happened to the NON-TREATED  observations had they      been treated
# So, two scenarios: treated units do not improve, or non-treated units do improve. 
sharp_null_ate <- 0
ys <- genouts(y, Z, ate = sharp_null_ate) 
```

Now that we have the simulated outcome scenarios, we use the `gendist` function to calculate the average treatment effect for each one of the permutations calculated before. However, we are interested in estimating the difference between the distributions using the Horvitz-Thompson estimators (`HT = TRUE`) and the inverse-probability weighted regression estimator (`HT = FALSE`)

``` r
## Blocks
# Inverse-probability weighted regression estimators
distout_bl <- gendist(ys, perms_bl, prob = probs_bl, HT = FALSE) 
# Horvits-Thompson estimators
distout_bl_HT <- gendist(ys, perms_bl, prob = probs_bl, HT = TRUE) 

## Clusters
# Inverse-probability weighted regression estimators
distout_cl <- gendist(ys, perms_cl, prob = probs_cl, HT = FALSE) 
# Horvits-Thompson estimators
distout_cl_HT <- gendist(ys, perms_cl, prob = probs_cl, HT = TRUE) 


## Blocks and Clusters
# Inverse-probability weighted regression estimators
distout <- gendist(ys, perms, prob = probs, HT = FALSE) 
# Horvits-Thompson estimators
distout_HT <- gendist(ys, perms, prob = probs, HT = TRUE) 

ate_dist_bl <- data.frame(ipw = distout_bl, ht = distout_bl_HT)
ate_dist_cl <- data.frame(ipw = distout_cl, ht = distout_cl_HT)
ate_dist <- data.frame(ipw = distout, ht = distout_HT)
```

We use the output from the previous sequence to feed into the `dispdist` function, which generates a set of results that include the *p* values and the standard deviations of the statistical tests from the randomisation. Additionally, it prints out a histogram with the density of results from the different permutations assuming the null hypothesis, and compares this with the observed Average Treatment Effect statistic. 
``` r
dispdist(distout, ate) # display characteristics of sampling dist. for inference
```

![](/resources/ate-null-dist.png)

    ## $two.tailed.p.value
    ## [1] 0.1666667
    ## 
    ## $two.tailed.p.value.abs
    ## [1] 0.1944444
    ## 
    ## $greater.p.value
    ## [1] 0.08333333
    ## 
    ## $lesser.p.value
    ## [1] 0.9444444
    ## 
    ## $quantile
    ##      2.5%     97.5% 
    ## -2.055556  2.222222 
    ## 
    ## $sd
    ## [1] 1.440879
    ## 
    ## $exp.val
    ## [1] 1.048454e-16

Using these results we cannot completely reject the null hypothesis that the treatment had a causal effect on the sampled population. 
