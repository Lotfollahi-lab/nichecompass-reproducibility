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
## Installation (MacOS)
##========================
rm(list=ls())
library(dplyr)
library(SRTsim)
library(SingleCellExperiment)
library(zellkonverter)
library(readr)

isid = 1
foldchange = 3 # foldchange, 1, 3, 5, 10
set.seed(isid)

load(paste0("./starmap_1020_0410_seurat_filter_layer.rds"))
info2   <- info %>% select(c(x,y,label)) %>% 
				filter(label!="Smc") %>% as.data.frame()

region_celltype_df <- matrix(0.1,4,4)
diag(region_celltype_df) <- 0.7
rownames(region_celltype_df) <- paste0("Region",1:4)
colnames(region_celltype_df) <- paste0("Celltype",1:4)

nGene = 2000
nLoc = 5000

##========================
## simulate spatial spot
##========================
simLoc <- simNewLocs(newN= nLoc,lay_out="random",preLoc = info2)

## assign the region labels
simLoc %<>% mutate(region_label = case_when(
		x <= 0.5*median(x) ~ "Region1",
		x > 0.5*median(x) & x <= median(x)~ "Region2",
		x > median(x) & x <= 1.5*median(x)~ "Region3",
		TRUE ~ "Region4"
	)) %>% as.data.frame()


## load the Ligand-Receptor data
load("./LR_dat.RData")

##========================
## Data Generation
##========================

## reference free
example_CCI_free = srtsim_cci_free(
					zero_prop_in = 0,
					disper_in = Inf,
					mu_in = 1, 
					numGene = nGene,
					location_in  = simLoc[,c("x","y","region_label")],
					# region_label = region_label,
					region_cell_map = region_celltype_df,
					sim_seed = isid,
					fc = foldchange,
					LR_in = LR_dat)

## reference based
simSRT  <- createSRT(count_in=sp_count[,rownames(info2)],loc_in =info2)
simSRT1 <- srtsim_fit(simSRT,sim_schem="domain")

example_CCI_ref = srtsim_cci_ref(
					EstParam = simSRT1@EstParam,
					numGene = nGene,
					location_in  = simLoc[,c("x","y","region_label")],
					# region_label = region_label,
					region_cell_map = region_celltype_df,
					sim_seed = isid,
					fc = foldchange,
					LR_in = LR_dat)

sce_free <- SingleCellExperiment(
  list(counts=example_CCI_free@simCounts),
  colData=DataFrame(region_label=example_CCI_free@simcolData$region_label,
                    cell_type=example_CCI_free@simcolData$celltype,
                    x=example_CCI_free@simcolData$x,
                    y=example_CCI_free@simcolData$y))

rownames(sce_free) <- example_CCI_free@simCounts@Dimnames[[1]]
colnames(sce_free) <- example_CCI_free@simCounts@Dimnames[[2]]
metadata(sce_free)$lr_data <- LR_dat
writeH5AD(sce_free, file = "./simulated_free_1.h5ad")

sce_ref <- SingleCellExperiment(
  list(counts=example_CCI_ref@simCounts),
  colData=DataFrame(region_label=example_CCI_ref@simcolData$region_label,
                    cell_type=example_CCI_ref@simcolData$celltype,
                    x=example_CCI_ref@simcolData$x,
                    y=example_CCI_ref@simcolData$y))

rownames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[1]]
colnames(sce_ref) <- example_CCI_ref@simCounts@Dimnames[[2]]
metadata(sce_ref)$lr_data <- LR_dat
writeH5AD(sce_ref, file = "./simulated_ref_1.h5ad")

original_function <- getAnywhere("srtsim_cci_ref")
print(original_function)

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
  
  print(GP_in)
  
  knn.list = FNN::get.knn(as.matrix(location_in[,c("x","y")]), k = numKNN)[[1]]
  gene_names = unique(unlist(GP_in[,1:2]))

  print(gene_names)

  simu_count = sim_cnt
  # gene_names = unique(unlist(LR_dat))
  rownames(simu_count)[1:length(gene_names)] = gene_names
  # mat_foldchange_record = matrix(0, length(gene_names), dim(simu_count)[2])
  
  # assign cell types to regions
  print("Construct Communications...")
  for(celltype_A in unique(GP_in$celltypeA)){
    print(celltype_A)
    cellid_A = which(c_simu==celltype_A)
    for(cellid in cellid_A){
      nearest_cells = knn.list[cellid,]
        for(nearid in nearest_cells){
          nearest_celltypeB = c_simu[nearid]
          for(fc in unique(GP_in$fold_change)){
              tmp_GP = GP_in[which(GP_in$celltypeA == celltype_A & GP_in$fold_change == fc & GP_in$celltypeB %in% nearest_celltypeB),]
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
  
  simInfo <- cbind.data.frame(location_in,celltype=c_simu)
  simsrt 	<- new(
    Class = "simSRT",
    simCounts = as(as.matrix(simu_count),"sparseMatrix"),
    simcolData = as(simInfo,"DataFrame"),
    metaParam = list(
      simSeed = sim_seed,
      simCCIKNN= numKNN,
      simParam = list(background_param = newEstParam, effect_size = fc),
      simType = "CCI_REF")
  )
  return(simsrt)
}

unlockBinding("srtsim_cci_ref", asNamespace("SRTsim"))
assign("srtsim_cci_ref", srtsim_cci_ref, envir = asNamespace("SRTsim"))
lockBinding("srtsim_cci_ref", asNamespace("SRTsim"))

gp_df <- read_csv("./omnipath_lr_simulations.csv")
head(gp_df)

data$sources <- lapply(data$sources, function(x) unlist(strsplit(x, ",")))
data$targets <- lapply(data$targets, function(x) unlist(strsplit(x, ",")))

example_CCI_ref = srtsim_cci_ref(
  EstParam = simSRT1@EstParam,
  numGene = nGene,
  location_in  = simLoc[,c("x","y","region_label")],
  # region_label = region_label,
  region_cell_map = region_celltype_df,
  sim_seed = isid,
  GP_in = gp_df)

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
writeH5AD(sce_ref, file = "./simulated_ref_omnipath.h5ad")