# KMD
Kernel measure of multi-sample dissimilarity (KMD)
measures the dissimilarity between multiple samples, based on the observations from them. It converges to the population quantity (depending on the kernel) which is between 0 and 1. A small value indicates the multiple samples are from the same distribution, and a large value indicates the corresponding distributions are different. The population quantity is 0 if and only if all distributions are the same, and 1 if and only if all distributions are mutually singular.

This package implements the computation of the sample KMD between several distributions based on independent observations from them, using K-nearest neighbor graphs and minimum spanning trees. It also implements the tests based on KMD for H0: the M distributions are equal against H1: not all the distributions are equal. Both permutation test and asymptotic test are available. The tests are consistent against all alternatives where at least two samples have different distributions.

## Installation

This package depends on R (>= 4.0.0). You can install the package KMD by (the package `devtools` needs to be installed first):

``` r
devtools::install_github("zh2395/KMD")
```

Alternatively, download the entire folder, and execute the following command in R:
``` r
install.packages("~/Downloads/KMD-main", repos=NULL, type="source")
```

You can uninstall the package by:
```r
remove.packages("KMD")
```

## Usage of the functions
Here we briefly introduce the functions in this package.
See the documentation (help page) of the R package for more details.

`KMD` implements the KMD estimator based on geometric graphs.
The inputs are:
`X`: the data matrix (n by dx) or the distance/similarity matrix (n by n);
`Y`: a vector of length n, indicating the labels (from 1 to M) of the data;
`M`: the number of possible labels;
`Knn`: the number of nearest neighbors to use, or "MST".
The recommended default value for `Knn` is 1;
`Kernel`: an M by M kernel matrix with row i and column j being the kernel value k(i, j); or "discrete" which indicates using the discrete kernel.
``` r
library(KMD)
n = 60
d = 2
set.seed(1)
X1 = matrix(runif(n*d/2),ncol = d)
X2 = matrix(runif(n*d/2),ncol = d)
X2[,1] = X2[,1] + 1
X = rbind(X1,X2)
Y = c(rep(1,n/2),rep(2,n/2))
print(KMD(X, Y, M = 2, Knn = 1, Kernel = "discrete"))
# 0.9344444. X1 and X2 are mutually singular, so the theoretical KMD is 1.
print(KMD(X, Y, M = 2, Knn = 1, Kernel = base::diag(c(1,1))))
# 0.9344444. This is essentially the same as specifying the discrete kernel above.
print(KMD(X, Y, M = 2, Knn = 2, Kernel = "discrete"))
print(KMD(X, Y, M = 2, Knn = "MST", Kernel = "discrete"))
# 0.9508333, 0.9399074. One can also use other geometric graphs (2-NN graph and MST here) to estimate the same theoretical quantity.
```

`KMD_test` implements the tests based on KMD.
Both permutation test and asymptotic test are available.
The tests are consistent against all alternatives where at least two samples have different distributions.
A small KMD value indicates the multiple samples are from the same distribution, and a large KMD value indicates the corresponding distributions are different.
The null hypothesis that all samples are from the same distribution is rejected for large KMD value. The permutation test returns the p-value given by (sum(KMD_i >= KMD_0) + 1)/(B + 1), where KMD_i is the KMD computed after a random permutation on the Y labels, and B is the total number of permutations that have been performed.
The asymptotic test first normalizes the KMD by the square root of the permutation variance, and then returns the p-value given by: P(N(0,1) > normalized KMD).

The inputs of `KMD_test` are:
`X`: the data matrix (n by dx) or the distance/similarity matrix (n by n);
`Y`: a vector of length n, indicating the labels (from 1 to M) of the data;
`M`: the number of possible labels;
`Knn`: the number of nearest neighbors to use, or "MST".
The recommended default value for `Knn` is 0.1n;
`Kernel`: an M by M kernel matrix with row i and column j being the kernel value k(i, j); or "discrete" which indicates using the discrete kernel;
`Permutation`: TRUE or FALSE; whether to perform permutation test or the asymptotic test;
`B`: the number of permutations to perform, only used for permutation test.

``` r
d = 2
set.seed(1)
X1 = matrix(rnorm(100*d), nrow = 100, ncol = d)
X2 = matrix(rnorm(100*d,sd=sqrt(1.5)), nrow = 100, ncol = d)
X3 = matrix(rnorm(100*d,sd=sqrt(2)), nrow = 100, ncol = d)
X = rbind(X1,X2,X3)
Y = c(rep(1,100),rep(2,100),rep(3,100))

print(KMD_test(X, Y, M = 3, Knn = 1, Kernel = "discrete"))
# A small p-value since the three distributions are not the same.
print(KMD_test(X, Y, M = 3, Knn = 1, Kernel = "discrete", Permutation = FALSE))
# p-value of the asymptotic test is similar to that of the permutation test
print(KMD_test(X, Y, M = 3, Knn = 1, Kernel = diag(c(10,1,1))))
# p-value is improved by using a different kernel
print(KMD_test(X, Y, M = 3, Knn = 30, Kernel = "discrete"))
# The suggested choice Knn = 0.1n yields a very small p-value.
print(KMD_test(X, Y, M = 3, Knn = "MST", Kernel = "discrete"))
# One can also use the MST.
print(KMD_test(X, Y, M = 3, Knn = 2, Kernel = "discrete"))
# MST has similar performance as 2-NN, which is between 1-NN and 30-NN

# Check null distribution of the z values
ni = 100
n = 3*ni
d = 2
Null_KMD = function(id){
   set.seed(id)
   X = matrix(rnorm(n*d), nrow = n, ncol = d)
   Y = c(rep(1,ni),rep(2,ni),rep(3,ni))
   return(KMD_test(X, Y, M = 3, Knn = "MST", Kernel = "discrete", Permutation = FALSE)[1,1])
}
hist(sapply(1:1000, Null_KMD), breaks = c(-Inf,seq(-5,5,length=50),Inf), freq = FALSE,
   xlim = c(-4,4), ylim = c(0,0.5), main = expression(paste(n[i]," = 100")),
   xlab = expression(paste("normalized ",hat(eta))))
 lines(seq(-5,5,length=1000),dnorm(seq(-5,5,length=1000)),col="red")
# The histogram of the normalized KMD is close to that of a standard normal distribution.
```





