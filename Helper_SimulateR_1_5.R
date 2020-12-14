#### Created by Sreeskandarajan Sutharzan (University of Michigan)##
#### on 09/24/19 ###################################################

#### Modified by Tian Zhang (University of Michigan) ##
#### on 10/31/19 ######################################

## Contains the helper functions

get_means <- function(data){
  # calculate means based on a real data set "data"
  return(apply(data,1,mean))
}

get_vars <- function(data){
  # calculate variances based on a real data set "data"
  return(apply(data,1,var))
}

calc_p <- function(means,vars){
  # point-estimate parameter p of a binomial distribution based on given means and variances 
  # requirement: length of means is the same as that of vars
  p = 1 - means / vars 
  return(p)
} 

calc_r <- function(means,vars){
  # point-estimate parameter r of a binomial distribution based on given means and variances 
  # requirement: length of means is the same as that of vars
  r = means^2 / (vars - means)
  r[(vars<means) | ((vars - means)<1e-8)] = -1 # mark as -1 if it follows poisson distribution
  return(r)
}

adjust_mu <- function(mu){
  # restrict mu from being Inf
  mu[mu == Inf] <- 1.797693e+100
  return(mu)
}

runSimulation <- function(x, beta, size, s = 1){
  # This function runs the simulation
  # Inputs
    ## x: design matrix
    ## beta: matrix containing GLM cofficients. Rows: genes columns samples 
    ## size: vector containing size parameter for each gene
    ## s: matrix containing scale values. Rows:genes and columns:samples
  # Output
    ## k: The counts matrix
  
  # Generating the scaled mean expression based on the GLM
  # log2(q) = xs %*% beta
  q <- t(2^(x %*% beta)) # Unscaled mean expression
  mu <- s * q # scaled mean expression
  mu<-adjust_mu(mu)
  # Generating the counts based on negative binomial distribution
  nrows <- nrow(mu)
  ncols <- ncol(mu)
  k <- matrix(nrow = nrows, ncol = ncols) # The counts matrix
  for (r in 1:nrows){
    for (c in 1:ncols){
      if(size[r] > 0){
        k[r,c] <- rnbinom(n = 1, 
                          size = size[r],
                          mu = mu[r,c])
      } else {
        k[r,c] <- rpois(n = 1,
                        lambda = mu[r,c])
      }
    }
  }
  # Returning the counts
  return(list(k,mu))
}

genDesignMat <- function(n.indv, 
                         n.cluster, 
                         perc.cluster.gene = 0.05, 
                         dist.cluster = rep(1,n.cluster), 
                         n.biop = 1, 
                         plot.heatmap = FALSE){
  # This function generates the design matrix
  # Inputs
    ## n.indv: Number of indiviuals
    ## n.cluster: Number of desired clusters
    ## perc.cluster.gene: Percentage of genes would contribute to clustering
    ## dist.cluster: Desired distribution of numbers of genes in each cluster. 
      ## It should be a vector.The default setting is to make genes
      ## evenly distributed in the clusters.
    ## n.biop: Number of biopsies per indiviuals
    ## plot.heatmap: If TRUE plots the matrix as a heatmap
  # Output
    ## x: Design matrix
  
  x <- matrix(nrow = n.indv * n.biop, ncol = n.cluster)
  x[,] <- 0
  x[,1] <- 1 # intercept indication
  ## x[((n.biop*n.indv/2)+1):(n.biop*n.indv),2] <- 1 # Treatment indication
  
  # Filing the indiviual indications
  # check invariant
  if(length(dist.cluster) != n.cluster){ 
    stop("The length of dist.cluster doesn't match n.cluster!")
  }
  
  for (i in 2:n.cluster){ 
    startIdx <- n.indv * sum(dist.cluster[1:(i-1)]) / sum(dist.cluster) * n.biop + 1
    stopIdx <- startIdx + n.indv * dist.cluster[i] / sum(dist.cluster) * n.biop - 1
    x[startIdx:stopIdx,i] <- 1
  }
  
  # OPTIONAL: Plotting the design matrix as a heatmap
  if (plot.heatmap){
    heatmap(x,Rowv = NA, Colv = NA, revC = T, scale = "none",
            xlab = "Beta", ylab = "Indiviuals")
  }
  
  # generate labels for testing classification
  dist.num = numeric(n.cluster)
  for(i in 1:(n.cluster-1)){
    dist.num[i] = dist.cluster[i] / sum(dist.cluster) * n.indv
  }
  dist.num[n.cluster] = n.indv - sum(dist.num[1:(n.cluster - 1)])
  
  cluster.labels = numeric(n.indv)
  start = 1
  for(i in 1:n.cluster){
    cluster.labels[start:sum(dist.num[1:i])] = i
    start = sum(dist.num[1:i]) + 1
  }
  
  return(list(x,cluster.labels))
}

# to generate multivariate normal data, we will use "mvtnorm"
library(mvtnorm)

genBetaMat <- function(n.indv,
                       n.gene,
                       n.cluster, 
                       means,
                       covar,
                       perc.cluster.gene,
                       scale.var.beta = 1,
                       scale.covar.beta = 1,
                       plot.heatmap){
  # This function generates the mean functions
  # Inputs
    ## n.indv: Number of indiviuals
    ## n.gene: Number of genes
    ## n.cluster: Number of desired clusters
    ## means: Means of the genes. They are obtained from a real data set
    ## covar: covariance matrix learned before
    ## perc.cluster.gene: Percentage of genes would contribute to clustering
    ## plot.heatmap: If TRUE plots the matrix as a heatmap
  # Output
    ## beta: A matrix representing the mean functions
  
  beta = matrix(nrow = n.cluster, ncol = n.gene)
  n.gene.cluster = round(n.gene * perc.cluster.gene)
  
  # beta_0s are the means from the real data set 
  beta[1,] = log(means)
  
  # scaling covar/var of beta
  covar[row(covar) == col(covar)] = scale.var.beta * covar[row(covar) == col(covar)]
  covar[!(row(covar) == col(covar))] = scale.covar.beta * covar[!(row(covar) == col(covar))]
  
  # beta_c for clustering genes are generated from N(0,some sigma)
  # here, some sigma is equal to a adjustment parameter * the real standard deviation of the gene 
  for(i in 2:n.cluster){
    beta[i,1:n.gene.cluster] = rmvnorm(1, mean = rep(0,n.gene.cluster), sigma = covar)
  }
  
  # beta_c for non-clustering genes are just 0s
  for(i in (n.gene.cluster + 1):n.gene){
    beta[2:n.cluster,i] = 0
  }
  
  # heat plot
  if (plot.heatmap){
    heatmap(beta,Rowv = NA, Colv = NA, revC = T, scale = "none",
            xlab = "Genes", ylab = "Beta")
  }
  
  return(beta);
}
  
gene.simulate <- function(real.data,
                          n.indv, 
                          n.gene, 
                          n.cluster,
                          perc.cluster.gene = 0.05, 
                          dist.cluster = rep(1,n.cluster), 
                          scale.var.beta = 1,
                          scale.covar.beta = 1,
                          n.biop = 1,
                          s = 1, 
                          seed = 0,
                          designMat.heatmap = F,
                          betaMat.heatmap = F) {
  # The function generates a design matrix as users specify and 
    ## produce a matrix of genes' expression levels across different cases and genes.
  # Inputs
    ## real.data: A real data set of genes' expression levels. The reference group of the simulated
      ## data would be generated according to the data set.
    ## n.indv: Number of indiviuals
    ## n.genes: Number of genes
    ## n.cluster: Number of clusters
    ## perc.cluster.gene: The percentage of genes that will 
      ## contribute to clustering (except the mass cluster)
    ## dist.cluster: Desired distribution of numbers of genes in each cluster. 
      ## It should be a vector.The default setting is to make genes
      ## evenly distributed in the clusters.
    ## n.biop: Number of biopsies
    ## s: Scale
    ## seed: random seed
    ## designMat.heatmap: whether plot the heatmap for the design matrix or not
    ## betaMat.heatmap: whether plot the heatmap for the beta matrix or not
  # Output
    ## k: The simulated data matrix that has expression levels as its entries 
    ## for cases(labels for each column) and genes(labels for each row)
  
  print("Simulation Program Start...")
  # Setting the seed
  if(seed != 0){
    set.seed(seed)
  }
 
  
  print("Generating Design Matrix...")
  # Generating the design matrix
  x_labels <- genDesignMat(n.indv = n.indv,
                    n.biop = n.biop,
                    n.cluster = n.cluster,
                    dist.cluster = dist.cluster,
                    plot.heatmap = designMat.heatmap)
  # design matrix
  x = x_labels[[1]]
  # labels
  cluster.labels = x_labels[[2]]
  
  print("Design Matrix Generated!")
  
  # extract information from the real data
  means = get_means(real.data)
  vars = get_vars(real.data)
  size = calc_r(means,vars)
  n.cluster.gene = round(n.gene * perc.cluster.gene)
  #n.cluster.gene = round(n.gene * perc.cluster.gene)
  covar = cov(t(real.data)[,1:n.cluster.gene])
  
  print("Generating Beta Matrix...")
  # generate beta matrix
  beta = genBetaMat(n.indv,
                    n.gene,
                    n.cluster, 
                    means,
                    covar,
                    perc.cluster.gene, 
                    scale.var.beta,
                    scale.covar.beta,
                    betaMat.heatmap)
  print("Beta Matrix Generated!")
  
  print("Simulating based on Learned Information...")
  # Running the simulation
  k_Mu <- runSimulation(x = x, beta = beta, size = size, s = s)
  # return results: 
    ## K_Mu[[1]] is k; 
    ## cluster.label is labels indicating clusters for patients;
    ## beta is beta matrix
    ## x is design matrix
    ## size is size parameter for each gene of negative binomial distribution
    ## K_Mu[[2]] is mean parameter ofr each gene of negative binomial distribution
  result <- list(k_Mu[[1]],cluster.labels,beta,x,size,k_Mu[[2]])
  print("Simulation Done!")
  return(result)
}