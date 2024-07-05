##========================
## Imports
##========================
rm(list=ls())
library(dplyr)
library(SRTsim)
library(SingleCellExperiment)
library(zellkonverter)
library(readr)

##========================
## Overwrite simulation functions to use GP data (with multiple source and target genes) as input
##========================
srtsim_cci_ref = function(
    zero_prop_in = 0,
    disper_in = 0.01,
    disper_fc = 1/1000,
    mu_in = 1,
    ref_p = 0.4,
    free_p = 0.6,
    EstParam = NULL,
    numGene = 1000,
    location_in,
    region_cell_map,
    GP_in,
    sim_seed = 1,
    numKNN = 6,
    numSingleCellType = 2000){	
  
  set.seed(sim_seed)
  if(is.null(EstParam)){stop("EstParam is NULL, consider srtsim_cci_free for reference free simulation!")}
  if(!is.list(EstParam)){stop("EstParam has to be a list!!")}
  if(length(EstParam)< ncol(region_cell_map)){
    warning("The length of EstParam is less than targeted celltype number.")
    idx1 				<- sample(1:length(EstParam),ncol(region_cell_map),replace=TRUE)
    newEstParam <- EstParam[idx1]
    names(newEstParam) <- paste0("Celltype",1:ncol(region_cell_map))
  }else if(length(EstParam)> ncol(region_cell_map)){
    idx1 <- sample(1:length(EstParam),ncol(region_cell_map),replace=FALSE)
    newEstParam <- EstParam[idx1]
    names(newEstParam) <- paste0("Celltype",1:ncol(region_cell_map))
    print(paste0(paste(names(EstParam)[idx1],collapse=",")," are selected for data generation"))
  }else if(length(EstParam)== ncol(region_cell_map)){
    newEstParam <- EstParam
    names(newEstParam) <- paste0("Celltype",1:ncol(region_cell_map))
    print(paste0(names(EstParam)," are renamed to ",paste0("Celltype",1:ncol(region_cell_map))))
  }
  
  
  numLoc 	= nrow(location_in)
  region_label = location_in$region_label
  
  # numEntry = numSingleCellType*numGene # 5000 gene x 2000 cells to be sampled for each cell type
  
  # assign cell types to regions
  print("Generate Gene Expression Data...")
  celltype_count_in = list()
  free_parm <- c(zero_prop_in,disper_in,mu_in)
  for(celltype in colnames(region_cell_map)){
    celltype_param 	<- newEstParam[[celltype]]$marginal_param1
    idx_gene 				<- sample(1:nrow(celltype_param),numGene,replace=TRUE)
    tmp_parm 				<- celltype_param[idx_gene,1:3] ## ignore the model
    
    count_simu_ref 			<- t(apply(tmp_parm,1,function(param_in){rbinom(numSingleCellType, 1, 1-param_in[1]) * rnbinom(numSingleCellType, size = param_in[2]*disper_fc, mu = param_in[3])}))
    count_simu_free     <- matrix(rbinom(numSingleCellType*numGene, 1, 1-free_parm[1]) * rnbinom(numSingleCellType*numGene, size = free_parm[2]*disper_fc, mu = free_parm[3]), numGene, numSingleCellType)
    
    # Randomly pick from reference and free counts based on provided probabilities
    num_rows <- nrow(count_simu_ref)
    num_cols <- ncol(count_simu_ref)
    mask <- matrix(sample(c(TRUE, FALSE), num_rows * num_cols, replace = TRUE, prob = c(free_p, ref_p)), nrow = num_rows, ncol = num_cols)
    count_simu <- count_simu_ref
    count_simu[mask] <- count_simu_free[mask]
    
    celltype_count_in[[celltype]] = count_simu
    rm(celltype,count_simu)
  }
  
  # assign cell types to regions
  print("Generate Cell Types...")
  ztrue 	<- region_label
  c_simu 	<- rep(NA, length(ztrue))
  
  for(z in unique(region_label)){
    zi_idx <- ztrue == z # zi_idx is index for region z
    c_simu[zi_idx] <- sample(colnames(region_cell_map), sum(zi_idx), prob = region_cell_map[z,], replace = T)
  }
  
  # assign count data
  sim_cnt <- array(NA, c(numGene, numLoc))
  for(ct in unique(c_simu)){
    c_size <- sum(c_simu == ct)  # true number of cells in cell type c
    if(dim(celltype_count_in[[ct]])[2]<c_size){
      stop("Targeted cell number of ", ct, " is greater than prespecified cell number in background pool. Increase the numSingleCellType value")
    }else{
      cells_select <- sample(1:dim(celltype_count_in[[ct]])[2], c_size, replace = F)
    }
    
    # sample same number of group c cells in real data from generated cells
    sim_cnt[, c_simu == ct] <- as.matrix(celltype_count_in[[ct]][, cells_select])
    # for positions of original cell type c cells, assign generated counts of group c
  }
  colnames(sim_cnt) <- paste0("Cell" ,1:numLoc)
  rownames(sim_cnt) <- paste0("Gene" ,1:numGene)
  
  knn.list = FNN::get.knn(as.matrix(location_in[,c("x","y")]), k = numKNN)[[1]]
  gene_names = unique(unlist(GP_in[,1:2]))
  
  simu_count = sim_cnt
  # gene_names = unique(unlist(LR_dat))
  
  # sample genes to keep from GP genes
  if(dim(simu_count)[1] < length(gene_names)){
    gene_names <- sample(gene_names, dim(simu_count)[1], replace = FALSE)
  }
  
  rownames(simu_count)[1:length(gene_names)] = gene_names
  # mat_foldchange_record = matrix(0, length(gene_names), dim(simu_count)[2])
  
  # assign cell types to regions
  print("Construct Communications...")
  for(celltype_A in unique(GP_in$celltypeA)){
    for(region_A in unique(GP_in$regionA)){
      cellid_A = which(c_simu==celltype_A & ztrue==region_A)
      for(cellid in cellid_A){
        nearest_cells = knn.list[cellid,]
        for(nearid in nearest_cells){
          nearest_celltypeB = c_simu[nearid]
          for(increment_param in unique(GP_in$increment_param[GP_in$increment_param != 1])){
            tmp_GP = GP_in[which(GP_in$celltypeA == celltype_A &
                                 GP_in$regionA == region_A &
                                 GP_in$increment_param == increment_param &
                                 GP_in$celltypeB %in% nearest_celltypeB),]
            source_gene_ind=na.omit(match(c(unlist(tmp_GP$sources)),rownames(simu_count)))
            target_gene_ind=na.omit(match(c(unlist(tmp_GP$targets)),rownames(simu_count)))
            if(length(target_gene_ind)>0){
              if(length(target_gene_ind)>1){
                increment_mu <- rowMeans(sim_cnt[target_gene_ind,]) * increment_param
              } else {
                increment_mu <- mean(sim_cnt[target_gene_ind,]) * increment_param
              }
              increment_prob <- min(1, increment_param/10)
              increment <- rbinom(length(target_gene_ind), 1, increment_prob) * rnbinom(length(target_gene_ind), size = Inf, mu = increment_mu)
              simu_count[target_gene_ind,nearid] = sim_cnt[target_gene_ind,nearid] + increment
            }
            if(length(source_gene_ind)>0){
              if(length(source_gene_ind)>1){
                increment_mu <- rowMeans(sim_cnt[source_gene_ind,]) * increment_param
              } else {
                increment_mu <- mean(sim_cnt[source_gene_ind,]) * increment_param
              }
              increment_prob <- min(1, increment_param/10)
              increment <- rbinom(length(source_gene_ind), 1, increment_prob) * rnbinom(length(source_gene_ind), size = Inf, mu = increment_mu)
              simu_count[source_gene_ind,nearid] = sim_cnt[source_gene_ind,nearid] + increment
            }
          }
        }
      }
    }
  }
  
  simInfo <- cbind.data.frame(location_in,celltype=c_simu)
  simsrt 	<- new(
    Class = "simSRT",
    simCounts = as(as.matrix(simu_count),"sparseMatrix"),
    simcolData = as(simInfo,"DataFrame"),
    metaParam = list(
      simSeed = sim_seed,
      simCCIKNN= numKNN,
      simParam = list(background_param = newEstParam),
      simType = "CCI_REF")
  )
  return(simsrt)
}

unlockBinding("srtsim_cci_ref", asNamespace("SRTsim"))
assign("srtsim_cci_ref", srtsim_cci_ref, envir = asNamespace("SRTsim"))
lockBinding("srtsim_cci_ref", asNamespace("SRTsim"))

##========================
## Run all simulations
##========================
st_data_bronze_folder_path = "./" # "./", "../../datasets/srt_data/bronze"
sim_folder_path = "./" # "./", "../../datasets/gp_data/data_simulation"

setwd("/Users/sebastian.birk/Downloads")

isid = 1
set.seed(isid)

for(increment_mode in c("medium")){
  n_genes = c(1090)
  gp_file_paths = c(paste(sim_folder_path, "sim_gps_filtered_", increment_mode, "increments.csv", sep = ""),
                    paste(sim_folder_path, "sim_gps_", increment_mode, "increments.csv", sep = ""))
  for(nLoc in c(10000)){
    for (i in seq_along(n_genes)){
      nGene = n_genes[i]  
      load(paste0(st_data_bronze_folder_path, "starmap_mouse_visual_cortex_simulation_reference/starmap_1020_0410_seurat_filter_layer.rds"))
      info2   <- info %>% select(c(x,y,label)) %>% 
        filter(label!="Smc") %>% as.data.frame()
      
      region_celltype_df <- matrix(0.1,8,8)
      diag(region_celltype_df) <- 0.7
      rownames(region_celltype_df) <- paste0("Region",1:8)
      colnames(region_celltype_df) <- paste0("Celltype",1:8)
      region_celltype_df[1, 5] <- 0.
      region_celltype_df[1, 6] <- 0.
      region_celltype_df[1, 7] <- 0.
      region_celltype_df[1, 8] <- 0.
      region_celltype_df[2, 5] <- 0.
      region_celltype_df[2, 6] <- 0.
      region_celltype_df[2, 7] <- 0.
      region_celltype_df[2, 8] <- 0.
      region_celltype_df[3, 3] <- 0.5
      region_celltype_df[3, 7] <- 0.
      region_celltype_df[3, 8] <- 0.
      region_celltype_df[4, 4] <- 0.5
      region_celltype_df[4, 7] <- 0.
      region_celltype_df[4, 8] <- 0.
      region_celltype_df[5, 1] <- 0.
      region_celltype_df[5, 2] <- 0.
      region_celltype_df[5, 3] <- 0.
      region_celltype_df[5, 4] <- 0.
      region_celltype_df[5, 5] <- 1
      region_celltype_df[5, 6] <- 0.
      region_celltype_df[5, 7] <- 0.
      region_celltype_df[5, 8] <- 0.
      region_celltype_df[6, 1] <- 0.
      region_celltype_df[6, 2] <- 0.
      region_celltype_df[6, 3] <- 0.
      region_celltype_df[6, 4] <- 0.
      region_celltype_df[6, 5] <- 0.
      region_celltype_df[6, 6] <- 1
      region_celltype_df[6, 7] <- 0.
      region_celltype_df[6, 8] <- 0.
      region_celltype_df[7, 1] <- 0.125
      region_celltype_df[7, 2] <- 0.125
      region_celltype_df[7, 3] <- 0.125
      region_celltype_df[7, 4] <- 0.125
      region_celltype_df[7, 5] <- 0.125
      region_celltype_df[7, 6] <- 0.125
      region_celltype_df[7, 7] <- 0.125
      region_celltype_df[7, 8] <- 0.125
      region_celltype_df[8, 1] <- 0.
      region_celltype_df[8, 2] <- 0.
      region_celltype_df[8, 3] <- 0.
      region_celltype_df[8, 4] <- 0.
      region_celltype_df[8, 5] <- 0.25
      region_celltype_df[8, 6] <- 0.25
      region_celltype_df[8, 7] <- 0.25
      region_celltype_df[8, 8] <- 0.25
      
      simLoc <- simNewLocs(newN= nLoc,lay_out="random",preLoc = info2)
      
      simLoc %<>% mutate(region_label = case_when(
        x <= 0.5*median(x) & y <= median(y) ~ "Region1",
        x <= 0.5*median(x) & y > median(y) ~ "Region2",
        x > 0.5*median(x) & x <= median(x) & y <= median(y) ~ "Region3",
        x > 0.5*median(x) & x <= median(x) & y > median(y) ~ "Region4",
        x > median(x) & x <= 1.5*median(x) & y <= median(y) ~ "Region5",
        x > median(x) & x <= 1.5*median(x) & y > median(y) ~ "Region6",
        x > 1.5*median(x) & y <= median(y) ~ "Region7",
        TRUE ~ "Region8"
      )) %>% as.data.frame()
      
      simSRT  <- createSRT(count_in=sp_count[,rownames(info2)],loc_in =info2)
      simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")
      
      gp_df <- read_csv(gp_file_paths[i])
      gp_df$sources <- lapply(gp_df$sources, function(x) unlist(strsplit(x, ",")))
      gp_df$targets <- lapply(gp_df$targets, function(x) unlist(strsplit(x, ","))) 
      
      example_CCI_ref = srtsim_cci_ref(
        EstParam = simSRT1@EstParam,
        numGene = nGene,
        location_in  = simLoc[,c("x","y","region_label")],
        # region_label = region_label,
        region_cell_map = region_celltype_df,
        sim_seed = isid,
        GP_in = gp_df,
        numSingleCellType = nLoc/5*2
      )
      
      gp_df$sources <- sapply(gp_df$sources, function(x) paste(x, collapse = ","))
      gp_df$targets <- sapply(gp_df$targets, function(x) paste(x, collapse = ","))
      
      sce_ref <- SingleCellExperiment(
        list(counts=example_CCI_ref@simCounts),
        colData=DataFrame(region_label=example_CCI_ref@simcolData$region_label,
                          cell_type=example_CCI_ref@simcolData$celltype,
                          x=example_CCI_ref@simcolData$x,
                          y=example_CCI_ref@simcolData$y))
      
      rownames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[1]]
      colnames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[2]]
      metadata(sce_ref)$gp_data <- gp_df
      
      dir.create(paste(st_data_bronze_folder_path, "sim_ref_", nGene, "genes_", nLoc, "locs_", increment_mode, "increments", sep = ""))
      save_path <- paste(st_data_bronze_folder_path, "sim_ref_", nGene, "genes_", nLoc, "locs_", increment_mode, "increments", "/sim_ref_", nGene, "genes_", nLoc, "locs_", increment_mode, "increments", ".h5ad", sep = "")
      writeH5AD(sce_ref, file = save_path)
    }
  }
}

