# Obtain KNN with ties broken at random
#
# @param X Matrix (n by dx)
# @param Knn number of nearest neighbors
# @return an n by Knn matrix showing the indices of KNN
get_neighbors = function(X,Knn) {
  if (!is.matrix(X)) X = as.matrix(X)
  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = Knn + 2)
  nn_index_X = nn_X$nn.idx[, 2:(Knn+1), drop=F]

  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)
  if (length(repeat_data) > 0) {
    gp_size = id = NULL
    df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
    df_X[, gp_size := length(id), by = "group"]

    for (i in 1:length(repeat_data)) {
      if (df_X$gp_size[i] > Knn) {
        # The number of repeated data is more than Knn
        group_indices = df_X$id[df_X$group==df_X$group[i]]
        if (Knn == 1 & length(group_indices) == 2) {
          nn_index_X[df_X$id[i],] = setdiff(group_indices, df_X$id[i]) # avoid length 1 vector in sample function
        }
        else {
          nn_index_X[df_X$id[i],] = sample(setdiff(group_indices, df_X$id[i]),Knn)
        }
      }
      else {
        if (nn_X$nn.dists[df_X$id[i], Knn+1] < nn_X$nn.dists[df_X$id[i], Knn+2]) {
          # The number of repeated data is less than Knn
          # but there is no tie at the KNN
          nn_index_X[df_X$id[i],] = setdiff(nn_X$nn.idx[df_X$id[i], 1:(Knn+1)], df_X$id[i])
        }
        else {
          # The number of repeated data is less than Knn
          # There are ties at the kNN
          distances <- proxy::dist(matrix(X[df_X$id[i], ], ncol = ncol(X)), matrix(X[-df_X$id[i], ], ncol = ncol(X)))
          tie_dist <- sort(distances, partial = Knn)[Knn]
          id_small <- which(distances < tie_dist)
          id_small = id_small + (id_small >= df_X$id[i])
          nn_index_X[df_X$id[i],1:length(id_small)] = id_small
          id_equal = sample(which(distances == tie_dist),Knn-length(id_small))
          id_equal = id_equal + (id_equal >= df_X$id[i])
          nn_index_X[df_X$id[i],(1+length(id_small)):Knn] = id_equal
        }
      }
    }
  }
  ties = which(nn_X$nn.dists[, Knn+1] == nn_X$nn.dists[, Knn+2])
  ties = setdiff(ties, repeat_data)
  if (length(ties) > 0) {
    for (i in ties) {
      distances <- proxy::dist(matrix(X[i, ], ncol = ncol(X)), matrix(X[-i, ], ncol = ncol(X)))
      tie_dist <- sort(distances, partial = Knn)[Knn]
      id_small <- which(distances < tie_dist)
      if (length(id_small) > 0) {
        id_small = id_small + (id_small >= i)
        nn_index_X[i,1:length(id_small)] = id_small
      }
      id_equal = sample(which(distances == tie_dist),Knn-length(id_small))
      id_equal = id_equal + (id_equal >= i)
      nn_index_X[i,(1+length(id_small)):Knn] = id_equal
    }
  }
  return(nn_index_X)
}

#' Kernel Measure of Multi-sample Dissimilarity
#'
#' Compute the kernel measure of multi-sample dissimilarity (KMD) with directed K-nearest neighbor (K-NN) graph or minimum spanning tree (MST).
#'
#' The kernel measure of multi-sample dissimilarity (KMD) measures the dissimilarity between
#' multiple samples, based on the observations from them.
#' It converges to the population quantity (depending on the kernel) which is between 0 and 1.
#' A small value indicates the multiple samples are from the same distribution,
#' and a large value indicates the corresponding distributions are different.
#' The population quantity is 0 if and only if all distributions are the same,
#' and 1 if and only if all distributions are mutually singular.
#'
#' If X is an n by n matrix, it will be interpreted as a distance/similarity matrix.
#' In such case, MST requires it to be symmetric (an undirected graph).
#' K-NN graph does not require it to be symmetric, with the nearest neighbors of point i computed based on the i-th row, and ties broken at random.
#' The diagonal terms (self-distances) will be ignored.
#' If X is an n by dx data matrix, Euclidean distance will be used for computing the K-NN graph (ties broken at random) and the MST.
#'
#' @param X the data matrix (n by dx) or the distance/similarity matrix (n by n)
#' @param Y a vector of length n, indicating the labels (from 1 to M) of the data
#' @param M the number of possible labels
#' @param Knn the number of nearest neighbors to use, or "MST"
#' @param Kernel an M by M kernel matrix with row i and column j being the kernel value \eqn{k(i, j)}; or "discrete" which indicates using the discrete kernel.
#'
#' @import data.table
#' @export
#' @return The algorithm returns a real number which is the sample KMD and is asymptotically between 0 and 1.
#' @seealso \code{\link{KMD_test}}
#' @examples
#' n = 60
#' d = 2
#' set.seed(1)
#' X1 = matrix(runif(n*d/2),ncol = d)
#' X2 = matrix(runif(n*d/2),ncol = d)
#' X2[,1] = X2[,1] + 1
#' X = rbind(X1,X2)
#' Y = c(rep(1,n/2),rep(2,n/2))
#' print(KMD(X, Y, M = 2, Knn = 1, Kernel = "discrete"))
#' # 0.9344444. X1 and X2 are mutually singular, so the theoretical KMD is 1.
#' print(KMD(X, Y, M = 2, Knn = 1, Kernel = base::diag(c(1,1))))
#' # 0.9344444. This is essentially the same as specifying the discrete kernel above.
#' print(KMD(X, Y, M = 2, Knn = 2, Kernel = "discrete"))
#' print(KMD(X, Y, M = 2, Knn = "MST", Kernel = "discrete"))
#' # 0.9508333, 0.9399074. One can also use other geometric graphs (2-NN graph and MST here)
#' # to estimate the same theoretical quantity.
KMD = function(X, Y, M = length(unique(Y)), Knn = 1, Kernel = "discrete") {
  if (!is.numeric(Y)) stop("Input Y is not numeric.")
  if ((sum(Y < 1) >= 1) || (sum(Y > M) >= 1)) stop("Label Y should be in 1,...,M.")
  if (is.matrix(Kernel)) discrete_kernel = FALSE
  else if (Kernel == "discrete") {
    discrete_kernel = TRUE
    Kernel = diag(rep(1,M))
  }
  else stop("The Kernel argument should be a matrix.")

  if (!is.matrix(X)) X = as.matrix(X)
  n = dim(X)[1]
  n_i = rep(0,M)
  for (i in 1:M) n_i[i] = sum(Y == i)

  if (Knn != "MST") {
    # First case: Knn graph
    if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
    if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
    if (nrow(X) == ncol(X)) {
      if (sum(diag(X)) != 0) warning("The distance of some data point to itself is non-zero. Self-distances are ignored when computing the nearest neighbors.")
      node_neighbors = function(i) {
        tie_dist <- sort(X[i,-i], partial = Knn)[Knn]
        id_small <- which(X[i,-i] < tie_dist)
        id_small = id_small + (id_small >= i)
        id_equal = which(X[i,-i] == tie_dist)
        if (length(id_equal) == 1) return(c(id_small,id_equal + (id_equal >= i)))
        if (length(id_equal) > 1) return(c(id_small, sample(id_equal + (id_equal >= i), Knn-length(id_small))))
      }
      nn_index_X = do.call(rbind, lapply(1:nrow(X),node_neighbors))
    }
    else {
      nn_index_X = get_neighbors(X,Knn)
    }
  }
  else {
    # Second case: MST
    return(KMD_MST(X, Y, M, discrete_kernel, Kernel, n_i))
  }

  # Compute U_stats
  if (discrete_kernel) U_stats = sum(n_i*(n_i-1))/n/(n-1)
  else U_stats = (t(n_i)%*%Kernel%*%n_i - sum(diag(Kernel)*n_i))/n/(n-1)

  # Compute 1/n sum K(\Delta_i,\Delta_i)
  if (discrete_kernel) mean_Kii = 1
  else mean_Kii = sum(diag(Kernel)*n_i)/n

  node_calculator = function(j) sum(Kernel[Y[j],Y[nn_index_X[j,]]])
  return((mean(sapply(1:n, node_calculator))/Knn - U_stats)/(mean_Kii - U_stats))
}



# Calculate KMD using minimum spanning tree (MST).
KMD_MST = function(X, Y, M, discrete_kernel, Kernel, n_i) {
  n = dim(X)[1]
  if (discrete_kernel) {
    U_stats = sum(n_i*(n_i-1))/n/(n-1)
    Kernel = diag(rep(1,M))
    mean_Kii = 1
  }
  else {
    U_stats = (t(n_i)%*%Kernel%*%n_i - sum(diag(Kernel)*n_i))/n/(n-1)
    mean_Kii = sum(diag(Kernel)*n_i)/n
  }
  # One-dimensional MST is simply a sorting
  if (dim(X)[2] == 1) {
    Y = Y[order(X)]
    node_calculator = function(j) return(Kernel[Y[j],Y[j-1]] + Kernel[Y[j],Y[j+1]])
    res = Kernel[Y[1],Y[2]] + Kernel[Y[n-1],Y[n]]
    return(((sum(sapply(2:(n-1), node_calculator))/2 + res)/n - U_stats)/(mean_Kii - U_stats))
  }

  # MST in general dimension
  if (nrow(X) == ncol(X)) {
    # MST from distance matrix
    X = igraph::graph_from_adjacency_matrix(X,
                  mode = "undirected", weighted = TRUE)
    mst = igraph::mst(X)
    out  = igraph::get.edges(mst,1:(igraph::ecount(mst)))
    # (n-1) by 2 matrix; the first and second columns correspond to "from" and "to" indices
  }
  else {
    # Euclidean MST
    out = mlpack::emst(as.data.frame(X))$output
    # (n-1) by 3 matrix; the first and second columns correspond to "from" and "to" indices
    # the index starts from 0
    out[,1] = out[,1] + 1
    out[,2] = out[,2] + 1
  }

  tmp = matrix(0,n,2)
  # the first column is the degree of node i
  # the second column is the sum of k(xi,x_{N(i)})
  for (i in 1:(n-1)) {
    tmp[out[i,1],1] = tmp[out[i,1],1] + 1
    tmp[out[i,2],1] = tmp[out[i,2],1] + 1
    tmp[out[i,1],2] = tmp[out[i,1],2] + Kernel[Y[out[i,1]],Y[out[i,2]]]
    tmp[out[i,2],2] = tmp[out[i,2],2] + Kernel[Y[out[i,2]],Y[out[i,1]]]
  }
  return((mean(tmp[,2]/tmp[,1]) - U_stats)/(mean_Kii - U_stats))
}


#' Testing based on KMD
#'
#' Testing based on the kernel measure of multi-sample dissimilarity (KMD).
#' Both permutation test and asymptotic test are available.
#' The tests are consistent against all alternatives where at least
#' two samples have different distributions.
#'
#' The kernel measure of multi-sample dissimilarity (KMD) measures the dissimilarity between
#' multiple samples using geometric graphs such as K-nearest neighbor (K-NN) graph and minimum spanning tree (MST),
#' based on the observations from them.
#' A small value indicates the multiple samples are from the same distribution,
#' and a large value indicates the corresponding distributions are different.
#' The test rejects the null hypothesis that all samples are from the same distribution for
#' large value of sample KMD. The permutation test returns the p-value given by
#' (sum(KMD_i >= KMD_0) + 1)/(B + 1), where KMD_i is the KMD computed after a random permutation on the Y labels,
#' and B is the total number of permutations that have been performed.
#' The asymptotic test first normalizes the KMD by the square root of the permutation variance,
#' and then returns the p-value given by: P(N(0,1) > normalized KMD).
#'
#' If X is an n by n matrix, it will be interpreted as a distance/similarity matrix.
#' In such case, MST requires it to be symmetric (an undirected graph).
#' K-NN graph does not require it to be symmetric, with the nearest neighbors of point i computed based on the i-th row, and ties broken at random.
#' The diagonal terms (self-distances) will be ignored.
#' If X is an n by dx data matrix, Euclidean distance will be used for computing the K-NN graph (ties broken at random) and the MST.
#'
#' @param X the data matrix (n by dx) or the distance/similarity matrix (n by n)
#' @param Y a vector of length n, indicating the labels (from 1 to M) of the data
#' @param M the number of possible labels
#' @param Knn the number of nearest neighbors to use, or "MST"
#' @param Kernel an M by M kernel matrix with row i and column j being the kernel value \eqn{k(i, j)}; or "discrete" which indicates using the discrete kernel.
#' @param Permutation TRUE or FALSE; whether to perform permutation test or the asymptotic test.
#' @param B the number of permutations to perform, only used for permutation test.
#'
#' @import data.table
#' @export
#' @return If Permutation == TRUE, permutation test is performed and the algorithm returns a p-value for testing H0: the M distributions are equal against H1: not all the distributions are equal.
#' If Permutation == FALSE, asymptotic test is performed and a 1 by 2 matrix: (z value, p-value) is returned.
#' @seealso \code{\link{KMD}}
#' @examples
#' d = 2
#' set.seed(1)
#' X1 = matrix(rnorm(100*d), nrow = 100, ncol = d)
#' X2 = matrix(rnorm(100*d,sd=sqrt(1.5)), nrow = 100, ncol = d)
#' X3 = matrix(rnorm(100*d,sd=sqrt(2)), nrow = 100, ncol = d)
#' X = rbind(X1,X2,X3)
#' Y = c(rep(1,100),rep(2,100),rep(3,100))
#' print(KMD_test(X, Y, M = 3, Knn = 1, Kernel = "discrete"))
#' # A small p-value since the three distributions are not the same.
#' print(KMD_test(X, Y, M = 3, Knn = 1, Kernel = "discrete", Permutation = FALSE))
#' # p-value of the asymptotic test is similar to that of the permutation test
#' print(KMD_test(X, Y, M = 3, Knn = 1, Kernel = diag(c(10,1,1))))
#' # p-value is improved by using a different kernel
#' print(KMD_test(X, Y, M = 3, Knn = 30, Kernel = "discrete"))
#' # The suggested choice Knn = 0.1n yields a very small p-value.
#' print(KMD_test(X, Y, M = 3, Knn = "MST", Kernel = "discrete"))
#' # One can also use the MST.
#' print(KMD_test(X, Y, M = 3, Knn = 2, Kernel = "discrete"))
#' # MST has similar performance as 2-NN, which is between 1-NN and 30-NN
#'
#' # Check null distribution of the z values
#' ni = 100
#' n = 3*ni
#' d = 2
#' Null_KMD = function(id){
#'   set.seed(id)
#'   X = matrix(rnorm(n*d), nrow = n, ncol = d)
#'   Y = c(rep(1,ni),rep(2,ni),rep(3,ni))
#'   return(KMD_test(X, Y, M = 3, Knn = "MST", Kernel = "discrete", Permutation = FALSE)[1,1])
#' }
#' hist(sapply(1:500, Null_KMD), breaks = c(-Inf,seq(-5,5,length=50),Inf), freq = FALSE,
#'   xlim = c(-4,4), ylim = c(0,0.5), main = expression(paste(n[i]," = 100")),
#'   xlab = expression(paste("normalized ",hat(eta))))
#' lines(seq(-5,5,length=1000),dnorm(seq(-5,5,length=1000)),col="red")
#' # The histogram of the normalized KMD is close to that of a standard normal distribution.
KMD_test = function(X, Y, M = length(unique(Y)), Knn = ceiling(length(Y)/10), Kernel = "discrete", Permutation = TRUE, B = 500) {
  if (!is.numeric(Y)) stop("Input Y is not numeric.")
  if ((sum(Y < 1)>=1) || (sum(Y > M)>=1)) stop("Label Y should be in 1,...,M.")
  if (is.matrix(Kernel)) discrete_kernel = FALSE
  else if (Kernel == "discrete") {
    discrete_kernel = TRUE
    Kernel = diag(rep(1,M))
  }
  else stop("The Kernel argument should be a matrix.")

  if (!is.matrix(X)) X = as.matrix(X)
  n = dim(X)[1]
  n_i = rep(0,M)
  for (i in 1:M) n_i[i] = sum(Y == i)

  ####################################################
  # If Permutation == TRUE, perform permutation test
  if (Permutation) {
    # Under the permutation distribution, only the first term in the numerator has randomness
    # Case 1: Knn graph
    if (Knn != "MST") {
      if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
      if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
      if (nrow(X) == ncol(X)) {
        if (sum(diag(X)) != 0) warning("The distance of some data point to itself is non-zero. Self-distances are ignored when computing the nearest neighbors.")
        node_neighbors = function(i) {
          tie_dist <- sort(X[i,-i], partial = Knn)[Knn]
          id_small <- which(X[i,-i] < tie_dist)
          id_small = id_small + (id_small >= i)
          id_equal = which(X[i,-i] == tie_dist)
          if (length(id_equal) == 1) return(c(id_small,id_equal + (id_equal >= i)))
          if (length(id_equal) > 1) return(c(id_small, sample(id_equal + (id_equal >= i), Knn-length(id_small))))
        }
        nn_index_X = do.call(rbind, lapply(1:nrow(X),node_neighbors))
      }
      else {
        nn_index_X = get_neighbors(X,Knn)
      }

      Perm_stat = function(Y,resample_vector) {
        Y = Y[resample_vector]
        node_calculator = function(j) return(sum(Kernel[Y[j],Y[nn_index_X[j,]]]))
        return(sum(sapply(1:n, node_calculator)))
      }
    }
    else {
      # Case 2: MST
      # Case 2.1: One-dimensional MST
      if (dim(X)[2] == 1) {
        order_of_X = order(X)
        Perm_stat = function(Y,resample_vector) {
          Y = Y[resample_vector]
          node_calculator = function(j) return(Kernel[Y[order_of_X[j-1]],Y[order_of_X[j]]] + Kernel[Y[order_of_X[j]], Y[order_of_X[j+1]]])
          res = Kernel[Y[order_of_X[1]],Y[order_of_X[2]]] + Kernel[Y[order_of_X[n-1]],Y[order_of_X[n]]]
          return(sum(sapply(2:(n-1), node_calculator))/2 + res)
        }
      }
      else {
        # Case 2.2: MST in general dimension
        if (nrow(X) == ncol(X)) {
          # MST from distance matrix
          X = igraph::graph_from_adjacency_matrix(X,
                                                  mode = "undirected", weighted = TRUE)
          mst = igraph::mst(X)
          out  = igraph::get.edges(mst,1:(igraph::ecount(mst)))
          # (n-1) by 2 matrix; the first and second columns correspond to "from" and "to" indices
        }
        else {
          # Euclidean MST
          out = mlpack::emst(as.data.frame(X))$output
          # (n-1) by 3 matrix; the first and second columns correspond to "from" and "to" indices
          # the index starts from 0
          out[,1] = out[,1] + 1
          out[,2] = out[,2] + 1
        }


        Perm_stat = function(Y,resample_vector) {
          Y = Y[resample_vector]
          tmp = matrix(0,n,2)
          # the first column is the degree of node i
          # the second column is the sum of k(xi,x_{N(i)})
          for (i in 1:(n-1)) {
            tmp[out[i,1],1] = tmp[out[i,1],1] + 1
            tmp[out[i,2],1] = tmp[out[i,2],1] + 1
            tmp[out[i,1],2] = tmp[out[i,1],2] + Kernel[Y[out[i,1]],Y[out[i,2]]]
            tmp[out[i,2],2] = tmp[out[i,2],2] + Kernel[Y[out[i,2]],Y[out[i,1]]]
          }
          return(sum(tmp[,2]/tmp[,1]))
        }
      }
    }
    b <- boot::boot(data = Y, statistic = Perm_stat, R = B, sim = "permutation")
    return((sum(b$t >= b$t0) + 1)/(B + 1))
  }

  ######################################################
  # If Permutation == FALSE, perform asymptotic test
  # n > 3 is needed to compute tilde_c
  if (n < 4) stop("At least 4 observations are needed for the asymptotic test.")

  # Compute tilde_g_i
  if (Knn != "MST") {
    # First case: Knn graph
    if ((floor(Knn) != Knn) || (Knn <= 0)) stop("Knn should be a positive integer or the string MST.")
    if (Knn + 2 > nrow(X)) stop("n should be greater than Knn + 1")
    n = dim(X)[1]
    if (nrow(X) == ncol(X)) {
      #node_neighbors = function(i) FastKNN::k.nearest.neighbors(i,X,k=Knn)
      if (sum(diag(X)) != 0) warning("The distance of some data point to itself is non-zero. Self-distances are ignored when computing the nearest neighbors.")
      node_neighbors = function(i) {
        tie_dist <- sort(X[i,-i], partial = Knn)[Knn]
        id_small <- which(X[i,-i] < tie_dist)
        id_small = id_small + (id_small >= i)
        id_equal = which(X[i,-i] == tie_dist)
        if (length(id_equal) == 1) return(c(id_small,id_equal + (id_equal >= i)))
        if (length(id_equal) > 1) return(c(id_small, sample(id_equal + (id_equal >= i), Knn-length(id_small))))
      }
      nn_index_X = do.call(rbind, lapply(1:nrow(X),node_neighbors))
    }
    else {
      nn_index_X = get_neighbors(X,Knn)
    }
    node_calculator = function(j) sum(Kernel[Y[j],Y[nn_index_X[j,]]])
    First_term_in_numerator = mean(sapply(1:n, node_calculator))/Knn

    num_in_neighbors = rep(0,n)
    g3 = 0
    for (i in 1:n) {
      for (j in nn_index_X[i,]) {
        num_in_neighbors[j] = num_in_neighbors[j] + 1
        if (i %in% nn_index_X[j,]) g3 = g3 + 1
      }
    }
    g3 = g3/n/Knn^2
    g2_g1 = mean(num_in_neighbors*(num_in_neighbors - 1))/Knn^2
    g1 = 1/Knn
  }
  else {
    # Second case: MST
    # Case 2.1: One-dimensional MST
    if (dim(X)[2] == 1) {
      Y = Y[order(X)]
      node_calculator = function(j) return(Kernel[Y[j-1],Y[j]] + Kernel[Y[j],Y[j+1]])
      First_term_in_numerator = (sum(sapply(2:(n-1), node_calculator))/2 + Kernel[Y[1],Y[2]] + Kernel[Y[n-1],Y[n]])/n

      g1 = 0.5 + 1/n
      g2_g1 = 0.5
      g3 = 0.5 + 0.5/n
    }
    else {
      # Case 2.2: MST in general dimension
      if (nrow(X) == ncol(X)) {
        # MST from distance matrix
        X = igraph::graph_from_adjacency_matrix(X,
                  mode = "undirected", weighted = TRUE)
        mst = igraph::mst(X)
        out  = igraph::get.edges(mst,1:(igraph::ecount(mst)))
        # (n-1) by 2 matrix; the first and second columns correspond to "from" and "to" indices
      }
      else {
        # Euclidean MST
        out = mlpack::emst(as.data.frame(X))$output
        # (n-1) by 3 matrix; the first and second columns correspond to "from" and "to" indices
        # the index starts from 0
        out[,1] = out[,1] + 1
        out[,2] = out[,2] + 1
      }

      tmp = matrix(0,n,2)
      # the first column is the degree of node i
      # the second column is the sum of k(xi,x_{N(i)})
      in_neighbor_indices <- vector("list", n)
      for (i in 1:(n-1)) {
        tmp[out[i,1],1] = tmp[out[i,1],1] + 1
        tmp[out[i,2],1] = tmp[out[i,2],1] + 1
        tmp[out[i,1],2] = tmp[out[i,1],2] + Kernel[Y[out[i,1]],Y[out[i,2]]]
        tmp[out[i,2],2] = tmp[out[i,2],2] + Kernel[Y[out[i,2]],Y[out[i,1]]]
        in_neighbor_indices[[out[i,1]]] = c(in_neighbor_indices[[out[i,1]]], out[i,2])
        in_neighbor_indices[[out[i,2]]] = c(in_neighbor_indices[[out[i,2]]], out[i,1])
      }
      First_term_in_numerator = mean(tmp[,2]/tmp[,1])

      g1 = mean(1/tmp[,1])
      node_calculator = function(j) return(sum(1/tmp[in_neighbor_indices[[j]],1])^2 - sum(1/tmp[in_neighbor_indices[[j]],1]^2))
      g2_g1 = mean(sapply(1:n, node_calculator))
      node_calculator = function(j) return(sum(1/tmp[in_neighbor_indices[[j]],1])/tmp[j,1])
      g3 = mean(sapply(1:n, node_calculator))
    }
  }

  # Compute tilde_a, tilde_b, tilde_c
  if (discrete_kernel) {
    tilde_a = sum(n_i*(n_i-1))/n/(n-1)
    tilde_b = sum(n_i*(n_i-1)*(n_i-2))/n/(n-1)/(n-2)
    tilde_c = (sum(n_i*(n_i-1))^2 - sum(n_i^2*(n_i-1)^2) + sum(n_i*(n_i-1)*(n_i-2)*(n_i-3)))/n/(n-1)/(n-2)/(n-3)
  }
  else {
    tilde_a = t(n_i)%*%(Kernel^2)%*%n_i - sum(diag(Kernel)^2*n_i)
    # sum_{i,j,k} K_{ij}K_{jk}n_in_jn_k = t(n)Kdiag(n)Kn
    tilde_b = t(n_i)%*%Kernel%*%(n_i * (Kernel%*%n_i)) - 2*t(n_i)%*%Kernel%*%(diag(Kernel)*n_i) - tilde_a + sum(diag(Kernel)^2*n_i)
    tilde_c = (t(n_i)%*%Kernel%*%n_i - sum(diag(Kernel)*n_i))^2 - 4*tilde_b - 2*tilde_a

    tilde_a = tilde_a/n/(n-1)
    tilde_b = tilde_b/n/(n-1)/(n-2)
    tilde_c = tilde_c/n/(n-1)/(n-2)/(n-3)
  }


  Sn = tilde_a*(g1+g3-2/(n-1)) + tilde_b*(g2_g1 - g1 - 2*g3 - 1 + 4/(n-1)) + tilde_c*(g3-g2_g1+1-2/(n-1))

  # Compute U_stats
  if (discrete_kernel) U_stats = sum(n_i*(n_i-1))/n/(n-1)
  else U_stats = (t(n_i)%*%Kernel%*%n_i - sum(diag(Kernel)*n_i))/n/(n-1)

  Output = (First_term_in_numerator - U_stats)*sqrt(n/Sn)
  Output = t(c(Output, stats::pnorm(Output, lower.tail = FALSE)))
  colnames(Output) = c("z value", "p value")
  return(Output)
}















