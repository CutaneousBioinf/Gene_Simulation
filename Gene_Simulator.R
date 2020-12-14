gene_simulator = function(object_parameter, par_sequence, 
                          level = 'any', removed = T, normalize = T, seed = 1234,
                          cnum_default = 5, cperc_default = 0.0025,
                          vscale_default = 1, covscale_default = 1,
                          n.indv = 300){
  # remove non-expression genes according to argument level
  if(level == "all"){ # use genes with all non-zero expression obs 
    select = !apply(PPnormalizedcount, MARGIN = 1, function(x){any(x == 0)})
  } else if(level == "any"){ # use genes with at least one non-zero obs
    select = apply(PPnormalizedcount, MARGIN = 1, function(x){sum(x)>0})
  } else{ # TODO: most informative way
  }
  
  PPnormalizedcount = PPnormalizedcount[select,]
  
  if(removed)
  {
    # number of genes removed
    print(paste("number of removed genes:", sum(!select)))
  }
    
  # set up default values of parameters
  n.gene <- nrow(PPnormalizedcount)
  n.cluster <- cnum_default # fix number of clusters to 5
  perc.cluster.gene = cperc_default # where, according to emperical test, 
    # K-Means can achieve about 0.5 NMI and Hierarchy can achieve about 0.1 NMI 
  scale.var.beta = vscale_default
  scale.covar.beta = covscale_default
  file_name = ""
  # generate data
  for(i in 1:length(par_sequence)){
    if(object_parameter == 'cnum'){
      n.cluster = par_sequence[i]
      file_name = "num_clusters_"
    }else if(object_parameter == 'cperc'){
      perc.cluster.gene = par_sequence[i]
      file_name = "perc_clusters_"
    }else if(object_parameter == 'vscale'){
      scale.var.beta = par_sequence[i]
      file_name = "scale_var_beta_"
    }else if(object_parameter == 'covscale'){
      scale.covar.beta = par_sequence[i]
      file_name = "scale_covar_beta_"
    }else{
      stop('Invalid object_parameter! It can only be one of the follows:
           \n cnum: number of clusters
           \n cperc: percentage of signal genes
           \n vscale: scale of gene variances
           \n covscale: scale of gene covariances')
    }
    result <- gene.simulate(real.data = PPnormalizedcount,
                            n.indv = n.indv,
                            n.gene = n.gene,
                            n.cluster = n.cluster,
                            perc.cluster.gene = perc.cluster.gene,
                            dist.cluster = rep(1,n.cluster),
                            scale.var.beta = scale.var.beta,
                            scale.covar.beta = scale.covar.beta,
                            n.biop = 1,
                            s = 1,
                            seed = seed,
                            designMat.heatmap = F,
                            betaMat.heatmap = F)
    k = result[[1]]
    
    if(normalize){
      # normalization
      k.inv <- t(apply(k,1,function(x){ qnorm((rank(x)-(3/8))/(length(x)-2*(3/8)+1)) }))
      k.inv[nrow(k.inv),] = k[nrow(k.inv),]
    }
    
    
    # row names
    gene.names <- sapply(1:n.gene, function(x) paste0("gene",x))
    rownames(k.inv)[-nrow(k.inv)] <- gene.names
    rownames(k.inv)[nrow(k.inv)] <- "cluster_label"
    
    # col names
    indv.names <- sapply(1:n.indv, function(x) paste0("individual",x))
    colnames(k.inv) = indv.names
    
    # write data
    print("Writing data...")
    data_file_name = paste0(file_name, i, ".csv")
    write.csv(k.inv, file = data_file_name,
              row.names = T, col.names = T)
    print("Writing done!")
  }
}

