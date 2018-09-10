+++
draft = false
title = "Dimensionality Reduction using Spectral Analysis"
date = 2018-08-01
+++


Dimensionality Reduction using Spectral Analysis
================================================

This short document presents an interesting alternative to Principal
Component Analysis by using spectral decomposition. The original content
of this document is based on a
[presentation](https://juanitorduz.github.io/laplacian_eigenmaps_dim_red.html)
given by Dr. Juan Orduz at the 2018 PyData Berlin conference. The linked
post is a practical application using Python to perform non-linear
dimensionality reduction based on spectral techniques. This post aims at
replicating the same results, but using R instead.

This is a list of all the packages you are going to need to carry out
this exercise in R.

-   <code>loe</code>: Carries out the main calculations of the
    spectral analysis. Functions used: So and so
-   <code>igraph</code>: Graph analysis. Functions used: so and so
-   <code>ggplot2</code>: Plotting package. Carries out the main
    calculations of the spectral analysis
-   <code>scatterploted</code>: 3D plotting package.
-   <code>reshape</code>: Tydiverse package used to change the shape
    of dataframes. Functions used: so and s

*NOTE:* As this is only a "translation" post, we are going to follow the
same examples from the original.

The Shortest Algebra Section Ever
---------------------------------

Consider *M*<sub>*n**x**n*</sub>(ℝ) the space of all *n**x**n* matrices
with real elements. Let $A\\inM\_{nxn}(\\mathbb{R})$ be a symmetric
matrix. *λ* ∈ ℂ is said to be an *eigenvalue* of A, with associated
*eigenvector* *f* ∈ ℝ<sup>*n*</sup>,  *f* ≠ 0, if
*A**f* = *λ**f*
 Additionally, we say that a set of vectors
ℬ = {*f*<sub>1</sub>, *f*<sub>2</sub>, ⋯, *f*<sub>*n*</sub>} is a basis
for ℝ<sup>*n*</sup> if 1. They are linearly independent; and 2. They
generate ℝ<sup>*n*</sup> And we can say that ℬ is an orthonormal basis
if for two elements *f*<sub>*i*</sub>, *f*<sub>*j*</sub> ∈ ℬ,
$&lt;f\_i,fj} = \\delta\_{ij}$.

The Spectral Theorem, at the center of this entire post, suggests that
there exists an orthonormal basis for ℝ<sup>*n*</sup> consisting of the
eigenvalues of *A*, where each eigenvalue is real.

Toy Example
-----------

Consider the adjacency matrix from the toy example:

    adj_mat <- matrix(c(0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0), nrow = 4, ncol = 4)
    adj_mat

    ##      [,1] [,2] [,3] [,4]
    ## [1,]    0    1    1    1
    ## [2,]    1    0    0    0
    ## [3,]    1    0    0    0
    ## [4,]    1    0    0    0

We can create the graph associated to this matrix by using the
<code>graph\_from\_adjacency\_matrix</code> function from the
<code>igraph</code> package. We also make the degree matrix, which
places the degree of each one of the nodes in the graph in its diagonal,
and the Laplacian matrix which calculates the difference between the
degree matrix and the adjacency matrix.

    # Setting the graph from the adjacency matrix
    g <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_mat)

    # Giving it a name
    g <- set_graph_attr(g, "name", "Toy Example")

    # Plotting the graph with smaller arrow sizes, and a title
    igraph::plot.igraph(x = g, edge.arrow.size = 0.5, main = "Toy Example")

![](SpectralDegreeReduction_files/figure-markdown_strict/unnamed-chunk-2-1.png)

    # Degree matrix and Laplacian Matrix
    deg_mat <- Diagonal(rowsum(adj_mat,group = c(1,1,1,1)), n = 4)
    laplacian <- as.matrix(deg_mat-adj_mat)

We are now going to check what this projection of this two dimensional
object onto ℝ looks like using the <code>spec.emb</code> function from
the <code>loe</code> package. We need to include the number of
dimensions of the desired space in which we want to project the graph.
In this case, *y* contains the

    y <- loe::spec.emb(adj_mat,1)
    y

    ## [1]  0.0000000 -0.5773503  0.7886751 -0.2113249

Using what we learned from the slides and the presentation available in
the original post, we would like to check if, analitically, the result
given by <code>loe::spec.emb</code> is what it should be.

Recall that the first eigenvalue for our matrix was 1. Using *y*

    lambda_1 <- 1

    a1 <- laplacian %*% y
    a2 <- lambda_1*(deg_mat %*% y)
    result <- a1 - a2

    #make sphere in R

    number_points <- 1000
    r <- runif(number_points, min = 0, max = 2*pi)
    rho <- runif(number_points, min = 0, max = 2*pi)

    x <- sin(rho)*cos(r)
    y <- sin(rho)*sin(r)
    z <- cos(rho)
    data <- as.data.frame(matrix(c(x, y, z), ncol =3))
    colnames(data) <- c("x_ax","y_ax", "z_ax")


    x <- data$x_ax[!(abs(data$z_ax)>0.99)]
    y <- data$y_ax[!(abs(data$z_ax)>0.99)]
    z <- data$z_ax[!(abs(data$z_ax)>0.99)]

    data_nopo <- as.data.frame(matrix(c(x, y, z), ncol =3))
    colnames(data_nopo) <- c("x_ax","y_ax", "z_ax")
    x1 <- rnorm(length(data_nopo$x_ax))

    s_colors <- rgb(abs(data_nopo$x_ax), abs(data_nopo$y_ax), abs(data_nopo$z_ax), maxColorValue = 1)
    scatterplot3d(data_nopo, pch = 16, main="Sphere",
                  xlab = "X",
                  ylab = "Y",
                  zlab = "Z",
                  box = FALSE,
                  color = s_colors)

![](SpectralDegreeReduction_files/figure-markdown_strict/unnamed-chunk-5-1.png)

    ggplot2::ggplot(data = data_nopo, aes(x = x_ax, y = y_ax)) + geom_point( color = s_colors)

![](SpectralDegreeReduction_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    n_components <- 2
    prin <- prcomp(x = data_nopo, scale. = T, rank. = n_components)
    #prin$center
    #biplot(prin)
    prinx <- prin$x[,1]
    priny <- prin$x[,2]
    prin_comp <- as.data.frame(matrix(c(prinx, priny), ncol = 2))
    colnames(prin_comp) <- c("first", "second")

    ggplot2::ggplot(data = prin_comp, aes(x = first, y = second)) + geom_point( color = s_colors)

![](SpectralDegreeReduction_files/figure-markdown_strict/unnamed-chunk-6-2.png)

    n_neighbors <- 50

    # Fit the object.
    nn <- get.knn(data = data_nopo, k = n_neighbors)$nn.index
    nn <- as.matrix(nn)
    nn <- as.data.frame(nn)
    nn$id <- rownames(nn)


    nn_melt <- melt(nn, id.vars = c("id"))
    nn_melt <- melt(nn)

    ## Using id as id variables

    nn_melt <- nn_melt[c("id", "value")]
    nn_melt$id <- as.integer(nn_melt$id)


    g <- graph.edgelist(as.matrix(nn_melt), directed=TRUE)
    #plot(g)

    g <- graph.edgelist(as.matrix(nn_melt), directed=FALSE)
    #plot(g)

    adj_mat <- as.data.frame(as.matrix(get.adjacency(graph = g)))
    adjmat <- ifelse(adj_mat == 2, 1, ifelse(adj_mat == 1, 1, 0))


    spec_emb <- spec.emb(A = as.matrix(adjmat), p = n_components)
    spec1 <- spec_emb[,1]
    spec2 <- spec_emb[,2]
    spec <- as.data.frame(matrix(c(spec1, spec2), ncol = 2))
    colnames(spec) <- c("first", "second")


    ggplot2::ggplot(data = spec, aes(x = first, y = second)) + geom_point( color = s_colors)

![](SpectralDegreeReduction_files/figure-markdown_strict/unnamed-chunk-6-3.png)
