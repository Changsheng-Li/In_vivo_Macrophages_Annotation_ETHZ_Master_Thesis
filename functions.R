library(GEOquery)
library(Seurat)
library(data.table)
library(dplyr)
library(readxl)
library(stringr)

load_seurat <- function(accession = NA,
                        study.code = NA,
                        rawdata.dir = NA,
                        saveddata.dir,
                        force.redownload = F,
                        after_UMAP = F){
  # get/check the accession ID
  entry <- accession_code_matching(accession = accession, study.code = study.code)
  ID <- entry$accession
  # check if previously stored Seurat object exists
  seuratfile.dir <- check_seurat_existance(accession = ID, saveddata.dir = saveddata.dir, after_UMAP = after_UMAP)
  
  if (!after_UMAP) {
    if (!is.na(seuratfile.dir) && !force.redownload){               # load if exists
      seurat_object <- readRDS(file = seuratfile.dir)
    } else {                                          # download it if not exists
      if(askYesNo(msg = "No local stored Seurat Object found that matches the accession code.\n Do you want to download it now?",
                  default = T,
                  prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))){
        download_funcion <- paste0("download_", ID)
        f <- get(download_funcion)
        seurat_object <- f(rawdata.dir = rawdata.dir, saveddata.dir = saveddata.dir)
      } else {
        stop()
      }
    }
  } else {
    seurat_object <- readRDS(file = seuratfile.dir)
  }
  return(seurat_object) 
}


####### data downloading functions for all studies #############################

##### Lung ###############

### 1.1 GSE123902

download_GSE123902 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE123902"
  location = "Lung"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # meta data
  gse <- getGEO("GSE123902")[[1]]
  names(pData(gse))
  
  #download raw data
  link = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/GSE123902_RAW.tar"
  file.dir = paste0(study.dir, "/", basename(link))
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
  }
  untar(tarfile = file.dir, exdir = study.dir)
  
  
  file_list = untar(tarfile = file.dir, list = T)
  data_sets_info <- data.frame(t(vapply(strsplit(file_list,split = "_",), '[', c(1,4), FUN.VALUE=character(2))))
  colnames(data_sets_info) <- c("ID", "Type")
  data_sets_info$name <- file_list
  
  
  for (file in file_list){
    file_ID = data_sets_info$ID[data_sets_info$name == file]
    
    #expression data
    counts <- fread(paste0(study.dir, "/", file), data.table = FALSE)
    rownames(counts) <- counts$V1
    counts <- counts[,-1]
    counts <- counts[,which(!duplicated(colnames(counts)))]
    counts <- t(counts)
    
    #meta data
    df<- do.call("rbind", replicate(pData(phenoData(gse))[file_ID,c(2,42:45)],n=ncol(counts),simplify = FALSE))
    colnames(df) <- vapply(strsplit(colnames(df),":"), '[', 1, FUN.VALUE=character(1))
    rownames(df) <- colnames(counts)
    
    #create/merge seurat object
    if (!exists("seurat_lung_1")) {
      seurat_lung_1 <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5, project = file_ID)
    } else {
      object <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5, project = file_ID)
      seurat_lung_1 <- merge(seurat_lung_1, object)
    }
  }
  
  #rename layers
  seurat_lung_1 <- JoinLayers(seurat_lung_1)
  seurat_lung_1 <- split(seurat_lung_1, f = seurat_lung_1$geo_accession)
  
  #save seurat object
  if (save_object){
    saveRDS(seurat_lung_1, file = paste0(saveddata.dir, "/seurat_lung_1.RDS"))
  }
  
  return(seurat_lung_1)
}

### 1.2 GSE131907

download_GSE131907 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE131907"
  location = "Lung"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # expression data
  gse <- getGEO("GSE131907")[[1]]
  names(pData(gse))
  
  link = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds.gz"
  file.dir = paste0(study.dir, "/", basename(link))
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
    gunzip(file.dir, remove=FALSE)
  }
  counts <- readRDS(paste0(study.dir, "/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds"))
  
  # meta data
  link = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_cell_annotation.txt.gz"
  file.dir = paste0(study.dir, "/", basename(link))
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
  }
  meta.data <- read.csv(file.dir, sep = "\t")
  
  seurat_lung_2 <- CreateSeuratObject(counts = counts, meta.data = meta.data)
  
  #save seurat object
  if (save_object){
    saveRDS(seurat_lung_2, file = paste0(saveddata.dir, "/seurat_lung_2.RDS"))
  }
  
  return(seurat_lung_2)
}



##### Brain ###############

### 2.1 GSM4972211

download_GSM4972211 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSM4972211"
  location = "Brain"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)

  # meta data
  gsm <- getGEO(ID)

  
  #download raw data
  link = gsm@header[["supplementary_file_1"]]
  file.dir = paste0(study.dir, "/", basename(link))
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
  }
  
  counts <- read.table(file = file.dir, header = T, row.names = 1, sep = ",")
  colnames(counts) <- gsub(pattern = ".", replacement = "-",x = colnames(counts), fixed = TRUE)
  
  #download raw data
  link = gsm@header[["supplementary_file_2"]]
  file.dir = paste0(study.dir, "/", basename(link))
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
  }
  
  meta.data <- read.table(file = file.dir, header = T, row.names = 5, sep = ",")
  
  seurat_brain_1 <- CreateSeuratObject(counts = counts, meta.data = meta.data, min.cells = 5)
  
  #save seurat object
  if (save_object){
    saveRDS(seurat_brain_1, file = paste0(saveddata.dir, "/seurat_brain_1.RDS"))
  }

  return(seurat_brain_1)
}


# download_GSE131928 <- function(rawdata.dir, saveddata.dir, save_object = T){
#   ID = "GSE131928"
#   location = "Brain"
#   study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
#   
#   # meta data
#   gse <- getGEO("GSE131928")[[1]]
#   head(pData(gse))
#   
#   #download raw data
#   link = gse@phenoData@data[["supplementary_file_1"]][gse@phenoData@data[["title"]] == "10X_GBM_IDHwt"]
#   file.dir = paste0(study.dir, "/", basename(link))
#   if (!file.exists(file.dir)) {
#     download.file(url = link, destfile = file.dir)
#   }
#   counts <- read.table(file = file.dir, header = T, row.names = 1)
#   
#   #download meta data
#   link = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131928/suppl/GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx"
#   file.dir = paste0(study.dir, "/", basename(link))
#   if (!file.exists(file.dir)) {
#     download.file(url = link, destfile = file.dir)
#   }
#   
#   meta.data <- read_excel(file.dir, skip = 43)
#   
#   meta.data <- meta.data[meta.data$`instrument model` == "NA", -7]
#   
#   seurat_brain_1 <- CreateSeuratObject(counts = counts, meta.data = meta.data, min.cells = 5, project = "GSE131928")
#   
#   #save seurat object
#   if (save_object){
#     saveRDS(seurat_brain_1, file = paste0(saveddata.dir, "/seurat_brain_1.RDS"))
#   }
#   
#   return(seurat_brain_1)
# }

### 2.2 GSE135045

download_GSE135045 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE135045"
  location = "Brain"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # meta data
  gse <- getGEO("GSE135045")
  pdata <- pData(gse[[1]])
  
  pdata <- pdata[pdata$source_name_ch1 == "high-grade primary glioma", c(2, 41:44)]
  colnames(pdata) <- vapply(strsplit(colnames(pdata),":"), '[', 1, FUN.VALUE=character(1))
  
  
  link <- pdata$supplementary_file_1
  file.dir = paste0(study.dir, "/", basename(link))
  download.file(url = link, destfile = file.dir)
  
  #expression data
  counts <- fread(paste0(file.dir), data.table = FALSE)
  rownames(counts) <- counts$V1
  counts <- counts[,-1]
  
  #meta data
  df<- do.call("rbind", replicate(pdata,n=ncol(counts),simplify = FALSE))
  rownames(df) <- colnames(counts)
  
  seurat_brain_2 <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5, project = ID)
  
  #save seurat object
  if(save_object){
    saveRDS(seurat_brain_2, file = paste0(saveddata.dir, "/seurat_brain_2.RDS"))
  }
  
  return(seurat_brain_2)
}



# download_GSE135045 <- function(rawdata.dir, saveddata.dir, save_object = T){
#   ID = "GSE135045"
#   location = "Brain"
#   study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
#   
#   # meta data
#   gse <- getGEO("GSE135045")
#   pdata <- pData(gse[[1]])
#   
#   for (accession in rownames(pdata)){
#     #download barcode file
#     link = pdata[accession, "supplementary_file_1"]
#     file.dir = paste0(study.dir, "/", "barcodes.tsv.gz")
#     download.file(url = link, destfile = file.dir)
#     
#     #download feature file
#     link = pdata[accession, "supplementary_file_2"]
#     file.dir = paste0(study.dir, "/", "features.tsv.gz")
#     download.file(url = link, destfile = file.dir)
#     
#     #download matrix file
#     link = pdata[accession, "supplementary_file_3"]
#     file.dir = paste0(study.dir, "/", "matrix.mtx.gz")
#     download.file(url = link, destfile = file.dir)
#     
#     # read expression data
#     data <- Read10X(study.dir)
#     
#     # meta data
#     df <- do.call("rbind", replicate(pdata[accession,c(1,2,43:45)],n=ncol(data),simplify = FALSE))
#     colnames(df) <- vapply(strsplit(colnames(df),":"), '[', 1, FUN.VALUE=character(1))
#     patient_data <- str_split_fixed(str_match(df$title, "\\[([^]]+)")[,2], "-", 3)
#     df$patient <- paste0(patient_data[,1], "-", patient_data[,2])
#     df$fragment <- patient_data[,3]
#     
#     if (!exists("seurat_brain_2")){
#       seurat_brain_2 <- CreateSeuratObject(counts = data, meta.data = df, min.cells = 5)
#     } else {
#       seurat_tmp <- CreateSeuratObject(counts = data, meta.data = df, min.cells = 5)
#       seurat_brain_2 <- merge(seurat_brain_2, seurat_tmp)
#     }
#     
#   }
#   
#   #rename layers
#   seurat_brain_2 <- JoinLayers(seurat_brain_2)
#   seurat_brain_2 <- split(seurat_brain_2, f = seurat_brain_2$geo_accession)
#   
#   #save seurat object
#   if(save_object){
#     saveRDS(seurat_brain_2, file = paste0(saveddata.dir, "/seurat_brain_2.RDS"))
#   }
#   
#   return(seurat_brain_2)
# }



##### Breast ###############

### 3.1 GSE118389

download_GSE118389 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE118389"
  location = "Breast"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # meta data
  gse <- getGEO("GSE118389")
  pdata <- pData(gse[[1]])
  
  #download raw data
  link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118389&format=file&file=GSE118389_counts_rsem.txt.gz"
  file.dir = paste0(study.dir, "/", basename(link))
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
  }
  
  counts <- read.table(file.dir, header = T, sep = "\t")
  
  
  ## potential to include more metadata from phenodata
  
  df<- do.call("rbind", replicate(pdata[,c(1:2,40:42)],n=ncol(counts),simplify = FALSE))
  colnames(df) <- vapply(strsplit(colnames(df),":"), '[', 1, FUN.VALUE=character(1))
  df <- as.data.frame(df %>% group_by(geo_accession) %>% slice_sample(n = 1))
  rownames(df) <- df$title
  
  seurat_breast_1 <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5)
  
  #save seurat object
  if(save_object){
    saveRDS(seurat_breast_1, file = paste0(saveddata.dir, "/seurat_breast_1.RDS"))
  }
  
  return(seurat_breast_1)
}

### 3.2 GSE161529

download_GSE161529 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE161529"
  location = "Breast"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # meta data
  gse <- getGEO("GSE161529")
  pdata <- pData(gse[[1]])
  
  link = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161529/suppl/GSE161529_features.tsv.gz"
  file.dir = paste0(study.dir, "/features.tsv.gz")
  if (!file.exists(file.dir)) {
    download.file(url = link, destfile = file.dir)
  }
  
  for (accession in rownames(pdata)){
    #download barcode file
    link = pdata[accession, "supplementary_file_1"]
    file.dir = paste0(study.dir, "/", "barcodes.tsv.gz")
    download.file(url = link, destfile = file.dir)
    
    #download matrix file
    link = pdata[accession, "supplementary_file_2"]
    file.dir = paste0(study.dir, "/", "matrix.mtx.gz")
    download.file(url = link, destfile = file.dir)
    
    # read expression data
    data <- Read10X(study.dir)
    
    # meta data
    df <- do.call("rbind", replicate(pdata[accession,c(1,2,48:54)],n=ncol(data),simplify = FALSE))
    colnames(df) <- vapply(strsplit(colnames(df),":"), '[', 1, FUN.VALUE=character(1))
    rownames(df) <- Cells(data)
    
    if (!exists("seurat_breast_2")){
      seurat_breast_2 <- CreateSeuratObject(counts = data, meta.data = df, min.cells = 5)
    } else {
      seurat_tmp <- CreateSeuratObject(counts = data, meta.data = df, min.cells = 5)
      seurat_breast_2 <- merge(seurat_breast_2, seurat_tmp)
    }
    
  }
  
  #rename layers
  seurat_breast_2 <- JoinLayers(seurat_breast_2)
  seurat_breast_2 <- split(seurat_breast_2, f = seurat_breast_2$geo_accession)
  
  #save seurat object
  if (save_object){
    saveRDS(seurat_breast_2, file = paste0(saveddata.dir, "/seurat_breast_2.RDS"))
  }
  
  return(seurat_breast_2)
}



##### Ovary ###############

### 4.1 GSE184880

download_GSE184880 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE184880"
  location = "Ovary"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # meta data
  gse <- getGEO(ID)
  list_tumor <- pData(gse[[1]])[pData(gse[[1]])[,45]=="High-grade serous ovarian cancer tissue",]
  
  for (subset in list_tumor[,2]){
    #download barcode data
    link = list_tumor[subset, 'supplementary_file_1']
    file.dir = paste0(study.dir, "/", "barcodes.tsv.gz")
    download.file(url = link, destfile = file.dir)
    
    #download feature data
    link = list_tumor[subset, 'supplementary_file_2']
    file.dir = paste0(study.dir, "/", "features.tsv.gz")
    download.file(url = link, destfile = file.dir)
    
    #download matrix data
    link = list_tumor[subset, 'supplementary_file_3']
    file.dir = paste0(study.dir, "/", "matrix.mtx.gz")
    download.file(url = link, destfile = file.dir)
    
    counts <- Read10X(data.dir = study.dir)
    
    df<- do.call("rbind", replicate(list_tumor[subset,c(1,2, 43:46)],n=ncol(counts),simplify = FALSE))
    colnames(df) <- vapply(strsplit(colnames(df),":"), '[', 1, FUN.VALUE=character(1))
    rownames(df) <- colnames(counts)
    
    if (!exists("seurat_ovary_1")){
      seurat_ovary_1 <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5)
    } else {
      seurat_tmp <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5)
      seurat_ovary_1 <- merge(seurat_ovary_1, seurat_tmp)
    }
    
  }
  
  if(save_object){
    saveRDS(seurat_ovary_1, file = paste0(saveddata.dir, "/seurat_ovary_1.RDS"))
  }
  
  return(seurat_ovary_1)
}

### 4.2 GSE154600

download_GSE154600 <- function(rawdata.dir, saveddata.dir, save_object = T){
  ID = "GSE154600"
  location = "Ovary"
  study.dir <- study.dir.setup(rawdata.dir = rawdata.dir, location = location, accession = ID)
  
  # meta data
  gse <- getGEO(ID)
  list_tumor <- pData(gse[[1]])
  
  for (subset in list_tumor[,2]){
    
    #download barcode data
    link = list_tumor[subset, 'supplementary_file_1']
    file.dir = paste0(study.dir, "/", "barcodes.tsv.gz")
    download.file(url = link, destfile = file.dir)
    
    #download feature data
    link = list_tumor[subset, 'supplementary_file_2']
    file.dir = paste0(study.dir, "/", "features.tsv.gz")
    download.file(url = link, destfile = file.dir)
    
    #download matrix data
    link = list_tumor[subset, 'supplementary_file_3']
    file.dir = paste0(study.dir, "/", "matrix.mtx.gz")
    download.file(url = link, destfile = file.dir)
    
    counts <- Read10X(data.dir = study.dir)
    
    df<- do.call("rbind", replicate(list_tumor[subset,c(1,2, 43:47)],n=ncol(counts),simplify = FALSE))
    colnames(df) <- vapply(strsplit(colnames(df),":"), '[', 1, FUN.VALUE=character(1))
    rownames(df) <- colnames(counts)
    
    if (!exists("seurat_ovary_2")){
      seurat_ovary_2 <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5)
    } else {
      seurat_tmp <- CreateSeuratObject(counts = counts, meta.data = df, min.cells = 5)
      seurat_ovary_2 <- merge(seurat_ovary_2, seurat_tmp)
    }
    
  }
  
  if(save_object){
    saveRDS(seurat_ovary_2, file = paste0(saveddata.dir, "/seurat_ovary_2.RDS"))
  }
  
  return(seurat_ovary_2)
}


####### supplementary functions ################################################

check_seurat_existance <- function(accession = NA, study.code = NA, after_UMAP = F, saveddata.dir = "~/DFS/RawData/403-Master-RawData/BioInf/Changsheng_Li/Saveddata"){
  entry <- accession_code_matching(accession = accession, study.code = study.code)
  accession <- entry$accession
  study.code <- entry$code
  location <- entry$location
  
  sub_code <- round((study.code%%1) * 10)
  if (!after_UMAP) {
    seuratfile.dir <- paste0(saveddata.dir, "/seurat_", tolower(location), "_", sub_code, ".RDS")
  } else {
    seuratfile.dir <- paste0(saveddata.dir, "Seurat_object_integrated_UMAP/", study.code, "_integrated_UMAP.RDS")
  }
  
  if(file.exists(seuratfile.dir)){
    return(seuratfile.dir)
  } else {
    return(NA)
  }
}

accession_code_matching <- function(accession = NA, study.code = NA){
  check_list = data.frame(code = c(1.1, 1.2, 2.1, 2.2, 3.1, 3.2, 4.1, 4.2),
                          accession = c("GSE123902", "GSE131907", "GSM4972211", "GSE135045", "GSE118389", "GSE161529", "GSE184880", "GSE154600"),
                          location = c("Lung", "Lung", "Brain", "Brain", "Breast", "Breast", "Ovary", "Ovary"))
  if (is.na(accession)){
    if(is.na(study.code)){
      stop("Please provide at least one variable")
    } else {
      entry <- check_list[check_list$code == study.code,]
    }
  } else {
    entry <- check_list[check_list$accession == accession,]
    if(!is.na(study.code) && entry$code == study.code){
      stop("The accession and study code do not match")
    }
  }
  if (nrow(entry) == 0){
    stop("Entry not found. Check your accession or study code")
  } else {
    return(entry)
  }
}

study.dir.setup <- function(rawdata.dir, location, accession){
  # temporarily change wd
  on.exit(setwd(rawdata.dir))
  
  # create location folder if not exists
  location.dir <- paste0(rawdata.dir, "/", str_to_title(location))
  if (!file.exists(location.dir)){
    dir.create(location.dir)
  }
  # create study folder if not exists
  study.dir = paste0(location.dir, "/", accession)
  if (!file.exists(study.dir)){
    dir.create(study.dir)
  }
  
  return(study.dir)
}
