# Clustering Test on Simulated data
rm(list = ls())
# include 17 data set with different parameters
setwd("C:/Users/Cory/Desktop/UROP2019/Simulation/Simulation_data_sets2")
library(stats)
library(aricode)
set.seed(2020)

################### Parameter: Number of clusters #####################

# NMI vector for the two algorithm: K-Means and Hierarchy-Clustering 
B = 5 # repeated times

# sequence of parameters
number_of_cluster = c(2,3,5,7,10,15)

nmi_kmeans_num = matrix(nrow = B, ncol = length(number_of_cluster))
nmi_hclust_num = numeric(length(number_of_cluster))

# evaluate performance; repeat B times for each parameter value
for(i in 1:length(number_of_cluster)){
        # read in data
        file_name = paste0("num_clusters_",number_of_cluster[i],".csv")
        data = read.csv(file_name, header = T, row.names = 1)
        # shuffle the data
        data = data[,sample(1:ncol(data),ncol(data))]
  
        # currentnumber of clusters
        k_clust = number_of_cluster[i]
  
        # percentage of cluster-associate genes
        perc = 0.003
  
        # rearrange the data
        data.label = as.factor(data[nrow(data),])
        data = data[-nrow(data),]
        data = data[1:(nrow(data)*perc),]
        data = apply(data,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))})
  
        for(b in 1:B){
                # fit K-Means and compute NMI
                kmeans_clustering = kmeans(data, k_clust, iter.max = 1000)
                nmi_kmeans_num[b,i] = NMI(kmeans_clustering$cluster, data.label)
        }  
        # fit Hierarchy and compute NMI
        dist_mat = dist(data)
        hier_clustring = hclust(dist_mat)
        hier_cut = cutree(hier_clustring, k = k_clust)
        nmi_hclust_num[i] = NMI(hier_cut, data.label)
}

# average performance for different parameter values across B tests
nmi_kmeans_num_avg = apply(nmi_kmeans_num, 2, mean)

# plot NMI ~ Number of clusters
plot(number_of_cluster, nmi_kmeans_num_avg, type = "b", col = 'red',
     main = "NMI ~ Number of clusters",
     ylab = "NMI", xlab = "Number of Clusters",
     xaxt = "n")
lines(number_of_cluster, nmi_hclust_num, type = "b", col = 'blue')
axis(side = 1, at = number_of_cluster)
legend('bottomright',legend = c("kmeans","hclust"),col = c("red","blue"),cex = 0.6,lty = c(1,1))



########### Parameter: Percentage of cluster-associate genes ############

# parameter sequence for perentage of cluster-associate genes
perc_of_cluster_genes = c(0.0001,0.0002,0.0005,0.001,0.002,
                              0.005,0.01,0.02,0.05,0.1,0.2,0.5)

# NMI vector for the two algorithm: K-Means and Hierarchy-Clustering 
nmi_kmeans_perc = matrix(nrow = B, ncol = length(perc_of_cluster_genes))
nmi_hclust_perc = numeric(length(perc_of_cluster_genes))

# evaluate performance; repeat B times for each parameter value
for(i in 1:length(perc_of_cluster_genes)){
  # read in data
  file_name = paste0("perc_clusters_",perc_of_cluster_genes[i],".csv")
  data = read.csv(file_name, header = T, row.names = 1)
  # shuffle the data
  data = data[,sample(1:ncol(data),ncol(data))]
  
  # current number of clusters
  k_clust = 5
  
  # percentage of cluster-associate genes
  perc = perc_of_cluster_genes[i]
  
  # rearrange the data
  data.label = as.factor(data[nrow(data),])
  data = data[-nrow(data),]
  data = data[1:(nrow(data)*perc),]
  data = apply(data,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))})
  
  for(b in 1:B){
                
      # fit K-Means and compute NMI
      kmeans_clustering = kmeans(data, k_clust, iter.max = 1000)
      nmi_kmeans_perc[b,i] = NMI(kmeans_clustering$cluster, data.label)
                
  }   
  # fit Hierarchy and compute NMI
  dist_mat = dist(data)
  hier_clustring = hclust(dist_mat)
  hier_cut = cutree(hier_clustring, k = k_clust)
      nmi_hclust_perc[i] = NMI(hier_cut, data.label)
}

# average performance for different parameter values across B tests
nmi_kmeans_perc_avg = apply(nmi_kmeans_perc, 2, mean)

axis_labels = as.character(perc_of_cluster_genes)
axis_labels = paste0('log(',axis_labels,')')

# plot: NMI ~ Percentage of cluster-associate Genes
plot(log(perc_of_cluster_genes), nmi_kmeans_perc_avg, type = "b", col = "red",
     main = "NMI ~ Perc of Cluster-associate Genes",
     ylab = "NMI", xlab = "log(Perc of cluster-associate Genes)",
     ylim = c(min(c(nmi_kmeans_perc_avg,nmi_hclust_perc)),
              max(c(nmi_kmeans_perc_avg,nmi_hclust_perc))),
     xaxt = "n")
lines(log(perc_of_cluster_genes), nmi_hclust_perc, type = "b", col = 'blue')
axis(1, at = log(perc_of_cluster_genes),
     labels = axis_labels)
legend('bottomright',legend = c("kmeans","hclust"),col = c("red","blue"),cex = 0.6,lty = c(1,1))


#################### Parameter: scale of variances #####################

# NMI vector for the two algorithm: K-Means and Hierarchy-Clustering 

# parameter sequence for scales of gene expresion variances
scales_of_vars = c(0.1,0.2,0.5,1,1.5,2)

nmi_kmeans_var = matrix(nrow = B, ncol = length(scales_of_vars))
nmi_hclust_var = numeric(length(scales_of_vars))

# evaluate performance; repeat B times for each parameter value
for(i in 1:length(scales_of_vars)){
  # read in data
  file_name = paste0("scale_var_beta_",scales_of_vars[i],".csv")
  data = read.csv(file_name, header = T, row.names = 1)
  # shuffle the data
  data = data[,sample(1:ncol(data),ncol(data))]
  
  # current number of clusters
  k_clust = 5
  
  # percentage of cluster-associate genes
  perc = 0.003
  
  # rearrange the data
  data.label = as.factor(data[nrow(data),])
  data = data[-nrow(data),]
  data = data[1:(nrow(data)*perc),]
  data = apply(data,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))})
  
        for(b in 1:B){
          
                
                # fit K-Means and compute NMI
                kmeans_clustering = kmeans(data, k_clust, iter.max = 1000)
                nmi_kmeans_var[b,i] = NMI(kmeans_clustering$cluster, data.label)
                
               
        }   
  # fit Hierarchy and compute NMI
  dist_mat = dist(data)
  hier_clustring = hclust(dist_mat)
  hier_cut = cutree(hier_clustring, k = k_clust)
  nmi_hclust_var[i] = NMI(hier_cut, data.label)
}

# average performance for different parameter values across B tests
nmi_kmeans_var_avg = apply(nmi_kmeans_var, 2, mean)

# plot: NMI ~ scales of variances of Beta
plot(scales_of_vars, nmi_kmeans_var_avg, type = "b", col = "red",
     main = "NMI ~ scales of vars of Beta",
     ylab = "NMI", xlab = "scales of vars of Beta",
     ylim = c(min(c(nmi_kmeans_var_avg,nmi_hclust_var_avg)),
              max(c(nmi_kmeans_var,nmi_hclust_var_avg))),
     xaxt = "n")
lines(scales_of_vars, nmi_hclust_var_avg, type = "b", col = 'blue')
axis(1, at = scales_of_vars)
legend('bottomright',legend = c("kmeans","hclust"),col = c("red","blue"),cex = 0.6,lty = c(1,1))




#################### Parameter: scale of covariances #####################

# parameter sequence for scales of gene expresion covariances
scales_of_covars = c(0.1,0.2,0.5,1,1.5,2)

# NMI vector for the two algorithm: K-Means and Hierarchy-Clustering 
nmi_kmeans_covar = matrix(nrow = B, ncol = length(scales_of_covars))
nmi_hclust_covar = numeric(length(scales_of_covars))

# evaluate performance; repeat B times for each parameter value
for(i in 1:length(scales_of_covars)){
  # read in data
  file_name = paste0("scale_covar_beta_",scales_of_covars[i],".csv")
  data = read.csv(file_name, header = T, row.names = 1)
  # shuffle the data
  data = data[,sample(1:ncol(data),ncol(data))]
  
  # current number of clusters
  k_clust = 5
  
  # percentage of cluster-associate genes
  perc = 0.003
  
  # rearrange the data
  data.label = as.factor(data[nrow(data),])
  data = data[-nrow(data),]
  data = data[1:(nrow(data)*perc),]
  data = apply(data,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1))})
  
        for(b in 1:B){
                
                
                # fit K-Means and compute NMI
                kmeans_clustering = kmeans(data, k_clust, iter.max = 1000)
                nmi_kmeans_covar[b,i] = NMI(kmeans_clustering$cluster, data.label)
                
                
        }   
  # fit Hierarchy and compute NMI
  dist_mat = dist(data)
  hier_clustring = hclust(dist_mat)
  hier_cut = cutree(hier_clustring, k = k_clust)
  nmi_hclust_covar[i] = NMI(hier_cut, data.label)
}

# average performance for different parameter values across B tests
nmi_kmeans_covar_avg = apply(nmi_kmeans_covar, 2, mean)

# plot: NMI ~ Percentage of cluster-associate Genes
plot(scales_of_covars, nmi_kmeans_covar_avg, type = "b", col = "red",
     main = "NMI ~ scales of covars",
     ylab = "NMI", xlab = "scales of covars",
     ylim = c(min(c(nmi_kmeans_covar_avg,nmi_hclust_covar_avg)),
              max(c(nmi_kmeans_covar_avg,nmi_hclust_covar_avg))),
     xaxt = "n")
lines(scales_of_covars, nmi_hclust_covar, type = "b", col = 'blue')
axis(1, at = scales_of_covars)
legend('bottomright',legend = c("kmeans","hclust"),col = c("red","blue"),cex = 0.6,lty = c(1,1))

