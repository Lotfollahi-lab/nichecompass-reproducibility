##========================
## Installation (MacOS)
##========================
# SRTsim requirements
install.packages("rgl") # on MacOS this requires brew install xquartz, brew install mesa, brew install libpng
install.packages("units") # on MacOS, this requires brew install udunits
Sys.setenv(PATH = paste("/opt/homebrew/bin", Sys.getenv("PATH"), sep = ":"))
install.packages("sf") # on MacOS, this requires brew install gdal
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("S4Vectors")
install.packages("nloptr") # on MacOS, this requires brew install cmake
install.packages("lme4")
install.packages("ggpubr")

# SRTsim
install.packages('devtools')
devtools::install_github('xzhoulab/SRTsim')

# SingleCellExperiment
Sys.setenv(PATH = paste("/opt/gcc/bin", Sys.getenv("PATH"), sep = ":"))
Sys.setenv(PATH = paste("/opt/llvm/bin", Sys.getenv("PATH"), sep = ":"))
BiocManager::install(version = "3.19")
BiocManager::install("zellkonverter") # on MacOS, use gcc compiler
# brew install gcc
# mkdir -p ~/.R
# nano ~/.R/Makevars
# CC = /usr/local/bin/gcc-10
# CXX = /usr/local/bin/g++-10
install.packages("readr")

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
## Overwrite simulation function to use ligand-receptor-target GP data as input
##========================
srtsim_cci_ref = function(
    EstParam = NULL,
    numGene = 1000,
    location_in,
    region_cell_map,
    GP_in,
    sim_seed = 1,
    numKNN = 4,
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
  for(celltype in colnames(region_cell_map)){
    celltype_param 	<- newEstParam[[celltype]]$marginal_param1
    idx_gene 				<- sample(1:nrow(celltype_param),numGene,replace=TRUE)
    tmp_parm 				<- celltype_param[idx_gene,1:3] ## ignore the model
    
    count_simu 			<- t(apply(tmp_parm,1,function(param_in){rbinom(numSingleCellType, 1, 1-param_in[1]) * rnbinom(numSingleCellType, size = param_in[2], mu = param_in[3])}))
    
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
      cellid_A = which(c_simu==celltype_A & c_simu==region_A)
      for(cellid in cellid_A){
        nearest_cells = knn.list[cellid,]
        for(nearid in nearest_cells){
          nearest_celltypeB = c_simu[nearid]
          for(fc in unique(GP_in$fold_change)){
            tmp_GP = GP_in[which(GP_in$celltypeA == celltype_A &
                                   GP_in$regionA == region_A &
                                   GP_in$fold_change == fc &
                                   GP_in$celltypeB %in% nearest_celltypeB),]
            source_gene_ind=na.omit(match(c(unlist(tmp_GP$sources)),rownames(simu_count)))
            target_gene_ind=na.omit(match(c(unlist(tmp_GP$targets)),rownames(simu_count)))
            if(length(source_gene_ind)>0)
              simu_count[target_gene_ind,cellid] = fc*sim_cnt[target_gene_ind,cellid]
            if(length(target_gene_ind)>0)
              simu_count[source_gene_ind,nearid] = fc*sim_cnt[source_gene_ind,nearid]
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
## Configuration
##========================
#setwd("/Users/sebastian.birk/Downloads")
isid = 1
#foldchange = 3 # foldchange, 1, 5, 10
set.seed(isid)

##========================
## Load reference data
##========================
load(paste0("./starmap_1020_0410_seurat_filter_layer.rds"))
info2   <- info %>% select(c(x,y,label)) %>% 
				filter(label!="Smc") %>% as.data.frame()

##========================
## Define regions and cell type proportions
##========================
region_celltype_df <- matrix(0.1,4,4) # dim (regions, cell types)
diag(region_celltype_df) <- 0.7
rownames(region_celltype_df) <- paste0("Region",1:4)
colnames(region_celltype_df) <- paste0("Celltype",1:4)

##========================
## Set simulation params
##========================
nGene = 20000 # unique sampled genes in prior GPs
nLoc = 50000
fc = "weak"

##========================
## Simulate spatial locations
##========================
simLoc <- simNewLocs(newN= nLoc,lay_out="random",preLoc = info2)

##========================
## Assign region labels
##========================
simLoc %<>% mutate(region_label = case_when(
		x <= 0.5*median(x) ~ "Region1",
		x > 0.5*median(x) & x <= median(x)~ "Region2",
		x > median(x) & x <= 1.5*median(x)~ "Region3",
		TRUE ~ "Region4"
	)) %>% as.data.frame()

##========================
## Load ligand-receptor data
##========================
#load("./LR_dat.RData")

##========================
## Generate data w/o reference
##========================
#example_CCI_free = srtsim_cci_free(
#					zero_prop_in = 0,
#					disper_in = Inf,
#					mu_in = 1, 
#					numGene = nGene,
#					location_in  = simLoc[,c("x","y","region_label")],
#					# region_label = region_label,
#					region_cell_map = region_celltype_df,
#					sim_seed = isid,
#					fc = foldchange,
#					LR_in = LR_dat)

##========================
## Learn domain-specific (here domains are cell types) count distributions of reference data
##========================
simSRT  <- createSRT(count_in=sp_count[,rownames(info2)],loc_in =info2)
simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")

##========================
## Generate data with reference
##========================
#example_CCI_ref = srtsim_cci_ref(
#					EstParam = simSRT1@EstParam,
#					numGene = nGene,
#					location_in  = simLoc[,c("x","y","region_label")],
#					# region_label = region_label,
#					region_cell_map = region_celltype_df,
#					sim_seed = isid,
#					fc = foldchange,
#					LR_in = LR_dat)

##========================
## Read ligand-receptor-target data (including de-novo GPs)
##========================
file_path <- paste("./combined_denovo_gps_simulations_", fc, "_fc.csv", sep = "")
gp_df <- read_csv(file_path)
head(gp_df)

gp_df$sources <- lapply(gp_df$sources, function(x) unlist(strsplit(x, ",")))
gp_df$targets <- lapply(gp_df$targets, function(x) unlist(strsplit(x, ",")))

##========================
## Run reference-based simulation
##========================
example_CCI_ref = srtsim_cci_ref(
  EstParam = simSRT1@EstParam,
  numGene = nGene,
  location_in  = simLoc[,c("x","y","region_label")],
  # region_label = region_label,
  region_cell_map = region_celltype_df,
  sim_seed = isid,
  GP_in = gp_df,
  numSingleCellType = nLoc/5*2)

gp_df$sources <- sapply(gp_df$sources, function(x) paste(x, collapse = ","))
gp_df$targets <- sapply(gp_df$targets, function(x) paste(x, collapse = ","))

##========================
## Convert to h5ad object
##========================
sce_ref <- SingleCellExperiment(
  list(counts=example_CCI_ref@simCounts),
  colData=DataFrame(region_label=example_CCI_ref@simcolData$region_label,
                    cell_type=example_CCI_ref@simcolData$celltype,
                    x=example_CCI_ref@simcolData$x,
                    y=example_CCI_ref@simcolData$y))

rownames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[1]]
colnames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[2]]
metadata(sce_ref)$gp_data <- gp_df
save_path <- paste("./simulated_ref_", nGene, "genes_", nLoc, "locs_", fc, "fc.h5ad", sep = "")

##========================
## Run all reference-based simulations
##========================
isid = 1
#foldchange = 3 # foldchange, 1, 5, 10
set.seed(isid)

for(fc in c("weak", "medium", "strong")){
  for(nLoc in c(5000, 50000)){
    for(nGene in c(500, 20000)){
      load(paste0("./starmap_1020_0410_seurat_filter_layer.rds"))
      info2   <- info %>% select(c(x,y,label)) %>% 
        filter(label!="Smc") %>% as.data.frame()
      
      region_celltype_df <- matrix(0.1,4,4) # dim (regions, cell types)
      diag(region_celltype_df) <- 0.7
      rownames(region_celltype_df) <- paste0("Region",1:4)
      colnames(region_celltype_df) <- paste0("Celltype",1:4)
      
      simLoc <- simNewLocs(newN= nLoc,lay_out="random",preLoc = info2)
      simLoc %<>% mutate(region_label = case_when(
        x <= 0.5*median(x) ~ "Region1",
        x > 0.5*median(x) & x <= median(x)~ "Region2",
        x > median(x) & x <= 1.5*median(x)~ "Region3",
        TRUE ~ "Region4"
      )) %>% as.data.frame()
      
      simSRT  <- createSRT(count_in=sp_count[,rownames(info2)],loc_in =info2)
      simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")
      
      file_path <- paste("./combined_denovo_gps_simulations_", fc, "_fc.csv", sep = "")
      gp_df <- read_csv(file_path)
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
      
      ##========================
      ## Convert to h5ad object
      ##========================
      sce_ref <- SingleCellExperiment(
        list(counts=example_CCI_ref@simCounts),
        colData=DataFrame(region_label=example_CCI_ref@simcolData$region_label,
                          cell_type=example_CCI_ref@simcolData$celltype,
                          x=example_CCI_ref@simcolData$x,
                          y=example_CCI_ref@simcolData$y))
      
      rownames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[1]]
      colnames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[2]]
      metadata(sce_ref)$gp_data <- gp_df
      
      save_path <- paste("./simulated_ref_", nGene, "genes_", nLoc, "locs_", fc, "fc.h5ad", sep = "")
      writeH5AD(sce_ref, file = save_path)
    }
  }
}