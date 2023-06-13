#' Validate input parameters differential expression analysis
#' Internal function called by sc_cell_type_de()
#'
#' Validate that there are no user input errors for the differential expression 
#' analysis - sc.cell.type.de
#' @param SCE SingleCellExperiment object, a specialised S4 class for storing 
#' data from single-cell experiments
#' @param design the design formula of class type `formula`. Equation used to 
#' fit the model- data for the generalised linear model.
#' @param pseudobulk_ID the column name in the SCE object to perform pseudobulk 
#' on, usually the patient identifier. This column is used for grouping in the 
#' pseudobulk approach
#' @param celltype_ID the column name in the SCE object for the cell type 
#' variable 
#' @param y the column name in the SCE object for the return variable e.g. 
#' "diagnosis" - Case or disease
#' @param region The column name in the SCE object for the study region. If 
#' there are multiple regions in the study (for example two brain regions). 
#' Pseudobulk values can be derived separately. Default is "single_region" 
#' which will not split by region.
#' @param coef character specifying which level to investigate for the 
#' differential expression analysis e.g. in a case/control study use "case" if 
#' case is the identifier in the y column to get positive fold changes to 
#' relate to case samples. leave as default value for continuous y. 
#' @param control character specifying which control level for the differential 
#' expression analysis e.g. in a case/control/other study use "control" in the 
#' y column to compare against. NOTE only need to specify if more than two 
#' groups in y, leave as default value for two groups or continuous y.
#' @param pval_adjust_method the adjustment method for the p-value in the 
#' differential expression analysis. Default is benjamini hochberg "BH". See  
#' stats::p.adjust for available options
#' @param adj_pval the adjusted p-value cut-off for the differential expression 
#' analysis, 0-1 range
#' @param folder the folder where the graphs from the differential expression 
#' analysis are saved. Default will create a folder in the current working 
#' directory "sc_cell_type_de_graphs". False will skip plotting.
#' @param rmv_zero_count_genes whether genes with no count values in any cell 
#' should be removed. Default is TRUE
#' @param verbose logical indicating if extra information about the 
#' differential expression analysis should be printed
#' @return NULL
#' 
validate_input_parameters_de<-function(SCE, design, pseudobulk_ID, celltype_ID, 
                                       y, region, coef, control, 
                                       pval_adjust_method, adj_pval, folder, 
                                       rmv_zero_count_genes, verbose){
  if(class(SCE)[1]!="SingleCellExperiment")
    stop("Please input the data as a SingleCellExperiment object")
  if(!is.character(celltype_ID))
    stop(paste0("Please input a character for celltype_ID indicating the ",
                "column holding the cell type information"))
  if(is.null(SCE[[celltype_ID]]))
    stop(paste0("The inputted celltype_ID: ",celltype_ID,
                " is not present in the SCE object, perhaps check for ",
                "spelling is correct"))
  #check design is a formula, formatting is correct and variables exist
  if(!inherits(design,"formula"))
    stop(paste0("Please input a formula for the design variable specifying",
                " the comparison. See examples."))
  design_txt <- paste0(deparse(design,width.cutoff = 500),collapse=',')
  # check if formula contains `~` to ensure user inputting correct format
  if(!grepl( "~", design_txt, fixed = TRUE))
    stop(paste0("Please input a correctly formated formula for the design",
                " variable containing a `~`. See examples."))
  design_txt <- gsub(".*~","",design_txt)
  design_txt <- gsub("^\\s+|\\s+$","",strsplit(design_txt, "[+]")[[1]])
  #check for duplicate entries
  if(length(design_txt)!=length(unique(design_txt)))
    stop("There are duplicate entries in the design formula")
  #check each variable to see if they are in SCE object
  for(i in design_txt){
    if(is.null(SCE[[i]]))
      stop(paste0("The inputted value: ",i,
                  " in the design formula is not present in the SCE ",
                  "object, perhaps check for spelling is correct"))
  }
  #only check y if user is setting it themselves
  if(!is.null(y)){
    if(!is.character(y))
      stop(paste0("Please input a character for y indicating the column ",
                  "holding the response variable information"))
    if(is.null(SCE[[y]]))
      stop(paste0("The inputted y value: ",y,
                  " is not present in the SCE object, perhaps check for ",
                  "spelling is correct"))
  }
  else{
    #if y not specified take last value in design matrix
    y <- design_txt[[length(design_txt)]]
    if(is.null(SCE[[y]]))
      stop(paste0("The inputted y value taken from your formula (the ",
                  "last variable): ",y,
                  " \nis not present in the SCE object, perhaps check ",
                  "the spelling is correct"))
  }
  #only check region if user is setting it themselves
  if(region!="single_region"){
    if(!is.character(region))
      stop(paste0("Please input a character for region indicating the ",
                  "column holding the variable information"))
    if(is.null(SCE[[region]]))
      stop(paste0("The inputted region value: ",region,
                  " is not present in the SCE object, perhaps check for ",
                  "spelling is correct"))
  }
  if(!is.null(coef)){
    if(!coef %in% unique(SCE[[y]]))
      stop(paste0("The inputted coef value: ",coef,
                  " is not present in y ", y,
                  "\nin the SCE object, perhaps check for spelling ",
                  "is correct"))
  }
  if(!is.null(control)){
    if(!control %in% unique(SCE[[y]]))
      stop(paste0("The inputted control value: ",control,
                  " is not present in y ", y,
                  "\nin the SCE object, perhaps check for spelling is",
                  " correct"))
  }
  if(!is.character(pseudobulk_ID))
    stop(paste0("Please input a character for pseudobulk_ID indicating the",
                " column holding the patient identifier"))
  if(is.null(SCE[[pseudobulk_ID]]))
    stop(paste0("The inputted patient_ID value: ",pseudobulk_ID,
                " is not present in the SCE object, ",
                "\nperhaps check for spelling is correct"))
  if(!is.character(pval_adjust_method))
    stop(paste0("Please input a character for the pval_adjust_method ",
                "indicating the method to be used"))
  if(class(adj_pval)[1]!="numeric")
    stop(paste0("Please input a 0-1 range number for the adj_pval ",
                "indicating the cut off of significance for the ",
                "differential expression analysis"))
  if(adj_pval<0|adj_pval>1)
    stop(paste0("Please input a 0-1 range number for the adj_pval ",
                "indicating the cut off of significance for the ",
                "differential expression analysis"))
  #if the user doesn't want to save plots from DE analysis they can pass 
  #in a false for the folder input parameter
  if(!isFALSE(folder))
    if(!is.character(folder))
      stop(paste0("For the folder input variable, please pass a director",
                  "y for the plots to be saved in.",
                  "\nGo with the default or pass in FALSE stop plotting"))
  if(!is.logical(verbose))
    stop("Please input TRUE/FALSE for verbose")
  if(!is.logical(rmv_zero_count_genes))
    stop("Please input TRUE/FALSE for rmv_zero_count_genes")
}

#' Single Cell, Cell Type Differential Expression Analysis
#'
#' Perform differential expression analysis across cell types based on single 
#' cell data with Pseudobulk approach. Returns the results of this differential 
#' expression analysis 
#' @param SCE SingleCellExperiment object, a specialised S4 class for storing 
#' data from single-cell experiments. A location of an R, rds or qs file to be 
#' loaded can also be passed. If using R objects make sure SCE object saved 
#' with the name SCE 
#' @param design the design formula of class type `formula`. Equation used to 
#' fit the model- data for the generalised linear model e.g. 
#' expression ~ sex + pmi + disease.
#' @param pseudobulk_ID the column name in the SCE object to perform pseudobulk 
#' on, usually the patient identifier. This column is used for grouping in the 
#' pseudobulk approach
#' @param celltype_ID the column name in the SCE object for the cell type 
#' variable. This is used to identify celltype specific DEGs. If there is only 
#' one cell type in the analysis, the analysis will be automatically updated 
#' accordingly.  
#' @param y the column name in the SCE object for the return variable e.g. 
#' "diagnosis" - Case or disease. Default is the last variable in the design 
#' formula. y can be discrete (logisitc regression) or continuous (linear 
#' regression)
#' @param region The column name in the SCE object for the study region. If 
#' there are multiple regions in the study (for example two brain regions). 
#' Pseudobulk values can be derived separately. Default is "single_region" 
#' which will not split by region.
#' @param coef character specifying which level to investigate for the 
#' differential expression analysis e.g. in a case/control study use "case" if 
#' case is the identifier in the y column to get positive fold changes to 
#' relate to case samples. leave as default value for continuous y. Default is 
#' NULL.
#' @param control character specifying which control level for the differential 
#' expression analysis e.g. in a case/control/other study use "control" in the 
#' y column to compare against. NOTE only need to specify if more than two 
#' groups in y, leave as default value for two groups or continuous y. Default 
#' is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the 
#' differential expression analysis. Default is benjamini hochberg "BH". See  
#' stats::p.adjust for available options
#' @param adj_pval the adjusted p-value cut-off for the differential expression 
#' analysis, 0-1 range
#' @param folder the folder where the graphs from the differential expression 
#' analysis are saved. Default will create a folder in the current working 
#' directory "sc.cell.type.de.graphs". False will skip plotting.
#' @param verbose logical indicating if extra information about the 
#' differential expression analysis should be printed
#' @return A list containing
#' \itemize{
#'   \item \code{celltype_DEGs}: list with the differentially expressed genes 
#'   (DEGs) for each cell type
#'   \item \code{celltype_all_genes}: list with all genes along with their 
#'   differential expression scores for each cell type
#'   \item \code{celltype_counts}: vector with the counts of cells after QC in 
#'   each cell type
#' }
#' @examples
#'\dontrun{
#'# To use the function input the formula of the comparison you want to make 
#'# along with the names of the pseudobulk ID and celltype ID.
#'# If you want to compare disease and control across cell types, taking into 
#'# account sex use formula = ~ sex + disease 
#'# firstly load  your SCE object (can otherwise specify location)
#' library(qs)
#' library(SingleCellExperiment)
#' SCE <- qread("../data/sce.qs")
#' sc.cell.type.de.return_incl_sex<- sc.cell.type.de(SCE_small,design= ~ sex + pathological_diagnosis_original,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="AD")
#' # If you only want to also account for other variables as well as sex, such as postmortem interval (PMI)  
#' sc.cell.type.de.return_sex_pmi<- sc.cell.type.de(SCE_small,design= ~ sex + PMI + pathological_diagnosis_original,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="AD") 
#' # If you only want to compare disease and control across cell types without any other constraints
#' sc.cell.type.de.return<- sc.cell.type.de(SCE_small,design= ~ pathological_diagnosis_original,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="AD")
#' # A good way to validate the DEG analysis approach is to run a sex comparison
#' # You would expect genes on the sex chromosomes to be DEGs, this can be done:
#' sc.cell.type.de.sex.return<- sc.cell.type.de(SCE,design= ~ sex,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="M")
#' #get gene names of chromosomal DEGs
#' library('biomaRt')
#' mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#' genes <- unique(as.character(our_degs$name))
#' gene_IDs <- getBM(filters= "ensembl_gene_id", 
#'                   attributes= c("ensembl_gene_id","chromosome_name","hgnc_symbol"),
#'                   values = genes, mart= mart,useCache = FALSE)
#' gene_IDs <- as.data.table(gene_IDs)
#' setnames(gene_IDs,"ensembl_gene_id","name")
#' DEGs <- rbindlist(sc.cell.type.de.return$celltype_DEGs,idcol="celltype")
#' setkey(DEGs,name)
#' #append gene names
#' DEGs[, gene_name := gene_IDs[DEGs, on=.(name), x.hgnc_symbol]]
#' DEGs[, chromosome_name := gene_IDs[DEGs, on=.(name), x.chromosome_name]]
#' #Xist should be there, key indicator of correct DEG analysis
#' DEGs[gene_name=="XIST",c("celltype","adj_pval","lfc")]
#' #Present in all cell types
#' # There should also be high number of other DEGs on sex chromosomes, check:
#' our_degs[chromosome_name %in% c("X","Y"),
#'          .N,by=.(celltype,chromosome_name)]
#'}
#' @import ggplot data.table
#' @export
sc_cell_type_de <- function(SCE, design, pseudobulk_ID, celltype_ID, y=NULL,
                            region="single_region", coef=NULL, control=NULL,
                            pval_adjust_method = "BH", adj_pval=0.05,
                            folder="sc_cell_type_de_graphs/", 
                            rmv_zero_count_genes=TRUE, verbose=F){
  #source necessary functions from other scripts
  # source("../scripts/validate_input_parameters_de.R")
  # source("../scripts/make_pseudobulk.R")
  # source("../scripts/de_analysis.R")
  # source("../scripts/plot_de_analysis.R")
  
  #need to load SCE if a directory is passed
  if(class(SCE)[1]=="character"){
    if(!file.exists(SCE))
      stop(paste0("Directory to SCE file doesn't exist, please check ",
                  "your file location"))
    #load the dataset
    if(substr(SCE,nchar(SCE)-2,nchar(SCE))==".qs"){
      #load QS objects  
      SCE <- qs::qread(SCE)
    }else if (substr(SCE,nchar(SCE)-3,nchar(SCE))==".rds"){
      # rds object  
      SCE <- readRDS(SCE)
    }else{
      #normal R object - saved object must also be called SCE for this 
      #to work, will hit error in input check otherwise
      load(SCE)
    }
  }
  
  #sense check inputs
  validate_input_parameters_de(SCE, design, pseudobulk_ID, celltype_ID, y, 
                               region, coef, control, pval_adjust_method, 
                               adj_pval, folder, rmv_zero_count_genes, 
                               verbose)
  
  #get counts of each cell type 
  #counts(SCE) <- as.matrix(counts(SCE))
  counts_celltypes <- SCE[[celltype_ID]]
  counts_celltypes <-as.vector(table(counts_celltypes))
  names(counts_celltypes) <- names(table(SCE[[celltype_ID]]))
  
  #first format formula
  design_txt <- paste0(deparse(design,width.cutoff = 500),collapse=',')
  #make design formula minus anything before ~
  formula <- as.formula(gsub(".*~","~",design_txt))
  #if exists remove everything before ~ this format is necessary for glm_gp()
  design_txt <- gsub(".*~","",design_txt)
  #split design by + to get components to bring forward for pseudobulk
  #add in pseudobulk by column
  pb_columns <- c(gsub("^\\s+|\\s+$","",strsplit(design_txt, "[+]")[[1]]),
                  pseudobulk_ID)
  #if y not specified take last value in design matrix
  if(is.null(y))
    y <- design_txt[[length(design_txt)]]
  #Check if continuous or categorical variable to be modeled
  y_contin <- FALSE
  if(is.numeric(SCE[[y]]))
    y_contin <- TRUE
  
  #Get pseudobulk values
  celltypes <- unique(SCE[[celltype_ID]])
  if(isTRUE(verbose))
    message("Deriving pseudobulk data")
  pb_dat <-
    lapply(celltypes,function(x) 
      make_pseudobulk(SCE[,SCE[[celltype_ID]]==x], 
                      pseudobulk_ID=pseudobulk_ID, 
                      pb_columns=pb_columns))
  names(pb_dat) <- celltypes
  
  #run edgeR LRT DE analysis
  celltype_de <-
    de_analysis(pb_dat,formula,y_name=y,y_contin,coef,control,
                pval_adjust_method, adj_pval, verbose)
  #get sig DEGs for each
  celltype_DEGs <- lapply(celltype_de, function(x) x[x$adj_pval<adj_pval,])
  
  unique_genes <-lapply(celltype_DEGs, function(x) x$name)
  unique_degs <- unique(unlist(unique_genes))
  
  #if no DEGs found break but return the DE analysis scores still
  if(length(unique_degs)==0)
    return(list("celltype_DEGs"=celltype_DEGs,
                "celltype_all_genes"=celltype_de,
                "celltype_counts"=counts_celltypes))
  
  if(isTRUE(verbose))
    message(length(unique_degs)," unique DEGs foundacross all cell types")
  
  celltype_all_genes_dt <- data.table::rbindlist(celltype_de,idcol = T)
  setnames(celltype_all_genes_dt,".id","celltype")
  celltype_DEGs_dt <- data.table::rbindlist(celltype_DEGs,idcol = T)
  setnames(celltype_DEGs_dt,".id","celltype")
  
  
  if(!isFALSE(folder)){
    if(isTRUE(verbose))
      message("Plotting the results of differential expression analysis")
    #make plots for DE analysis
    plot_de_analysis(pb_dat,y,celltype_DEGs_dt,celltype_all_genes_dt,
                     counts_celltypes,folder)
  }
  
  #get gene names for DEGs
  genes <- unique(data.table::rbindlist(celltype_DEGs)$name)
  gene_IDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, 
                                keytype = "GENEID", 
                                columns = c("GENEID","SYMBOL"))
  colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
  gene_IDs <- data.table::as.data.table(gene_IDs)
  data.table::setnames(gene_IDs,"ensembl_gene_id","name")
  #remove any dups in the reference set - two names for one ENSEMBL ID
  gene_IDs <- unique(gene_IDs,by="name")
  celltype_DEGs <-
    lapply(celltype_DEGs,
           function(i){setDT(i);setkey(i,name);
             i[,gene_name:=gene_IDs[i, on=.(name),x.hgnc_symbol]];
             setnames(i,"name","ensembl_id");setnames(i,"gene_name","HGNC")})
  
  #get gene names for all genes
  genes <- unique(data.table::rbindlist(celltype_de)$name)
  gene_IDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, 
                                keytype = "GENEID", 
                                columns = c("GENEID","SYMBOL"))
  colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
  gene_IDs <- data.table::as.data.table(gene_IDs)
  data.table::setnames(gene_IDs,"ensembl_gene_id","name")
  #remove any dups in the reference set - two names for one ENSEMBL ID
  gene_IDs <- unique(gene_IDs,by="name")
  celltype_de <-
    lapply(celltype_de,
           function(i){setDT(i);setkey(i,name);
             i[,gene_name:=gene_IDs[i, on=.(name),x.hgnc_symbol]];
             setnames(i,"name","ensembl_id");setnames(i,"gene_name","HGNC")})
  
  return(list("celltype_DEGs"=celltype_DEGs,
              "celltype_all_genes"=celltype_de,
              "celltype_counts"=counts_celltypes))
}


#' Single Cell, Cell Type Differential Expression Analysis
#'
#' Perform differential expression analysis across cell types based on single 
#' cell data with Pseudobulk approach. Returns the results of this differential 
#' expression analysis 
#' @param SCE SingleCellExperiment object, a specialised S4 class for storing 
#' data from single-cell experiments. A location of an R, rds or qs file to be 
#' loaded can also be passed. If using R objects make sure SCE object saved 
#' with the name SCE 
#' @param design the design formula of class type `formula`. Equation used to 
#' fit the model- data for the generalised linear model e.g. 
#' expression ~ sex + pmi + disease.
#' @param pseudobulk_ID the column name in the SCE object to perform pseudobulk 
#' on, usually the patient identifier. This column is used for grouping in the 
#' pseudobulk approach
#' @param celltype_ID the column name in the SCE object for the cell type 
#' variable. This is used to identify celltype specific DEGs. If there is only 
#' one cell type in the analysis, the analysis will be automatically updated 
#' accordingly.  
#' @param y the column name in the SCE object for the return variable e.g. 
#' "diagnosis" - Case or disease. Default is the last variable in the design 
#' formula. y can be discrete (logisitc regression) or continuous (linear 
#' regression)
#' @param region The column name in the SCE object for the study region. If 
#' there are multiple regions in the study (for example two brain regions). 
#' Pseudobulk values can be derived separately. Default is "single_region" 
#' which will not split by region.
#' @param coef character specifying which level to investigate for the 
#' differential expression analysis e.g. in a case/control study use "case" if 
#' case is the identifier in the y column to get positive fold changes to 
#' relate to case samples. leave as default value for continuous y. Default is 
#' NULL.
#' @param control character specifying which control level for the differential 
#' expression analysis e.g. in a case/control/other study use "control" in the 
#' y column to compare against. NOTE only need to specify if more than two 
#' groups in y, leave as default value for two groups or continuous y. Default 
#' is NULL.
#' @param pval_adjust_method the adjustment method for the p-value in the 
#' differential expression analysis. Default is benjamini hochberg "BH". See  
#' stats::p.adjust for available options
#' @param adj_pval the adjusted p-value cut-off for the differential expression 
#' analysis, 0-1 range
#' @param folder the folder where the graphs from the differential expression 
#' analysis are saved. Default will create a folder in the current working 
#' directory "sc.cell.type.de.graphs". False will skip plotting.
#' @param verbose logical indicating if extra information about the 
#' differential expression analysis should be printed
#' @return A list containing
#' \itemize{
#'   \item \code{celltype_DEGs}: list with the differentially expressed genes 
#'   (DEGs) for each cell type
#'   \item \code{celltype_all_genes}: list with all genes along with their 
#'   differential expression scores for each cell type
#'   \item \code{celltype_counts}: vector with the counts of cells after QC in 
#'   each cell type
#' }
#' @examples
#'\dontrun{
#'# To use the function input the formula of the comparison you want to make 
#'# along with the names of the pseudobulk ID and celltype ID.
#'# If you want to compare disease and control across cell types, taking into 
#'# account sex use formula = ~ sex + disease 
#'# firstly load  your SCE object (can otherwise specify location)
#' library(qs)
#' library(SingleCellExperiment)
#' SCE <- qread("../data/sce.qs")
#' sc.cell.type.de.return_incl_sex<- sc.cell.type.de(SCE_small,design= ~ sex + pathological_diagnosis_original,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="AD")
#' # If you only want to also account for other variables as well as sex, such as postmortem interval (PMI)  
#' sc.cell.type.de.return_sex_pmi<- sc.cell.type.de(SCE_small,design= ~ sex + PMI + pathological_diagnosis_original,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="AD") 
#' # If you only want to compare disease and control across cell types without any other constraints
#' sc.cell.type.de.return<- sc.cell.type.de(SCE_small,design= ~ pathological_diagnosis_original,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="AD")
#' # A good way to validate the DEG analysis approach is to run a sex comparison
#' # You would expect genes on the sex chromosomes to be DEGs, this can be done:
#' sc.cell.type.de.sex.return<- sc.cell.type.de(SCE,design= ~ sex,
#' pseudobulk_ID="sample_id", celltype_ID="allan_celltype",coef="M")
#' #get gene names of chromosomal DEGs
#' library('biomaRt')
#' mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#' genes <- unique(as.character(our_degs$name))
#' gene_IDs <- getBM(filters= "ensembl_gene_id", 
#'                   attributes= c("ensembl_gene_id","chromosome_name","hgnc_symbol"),
#'                   values = genes, mart= mart,useCache = FALSE)
#' gene_IDs <- as.data.table(gene_IDs)
#' setnames(gene_IDs,"ensembl_gene_id","name")
#' DEGs <- rbindlist(sc.cell.type.de.return$celltype_DEGs,idcol="celltype")
#' setkey(DEGs,name)
#' #append gene names
#' DEGs[, gene_name := gene_IDs[DEGs, on=.(name), x.hgnc_symbol]]
#' DEGs[, chromosome_name := gene_IDs[DEGs, on=.(name), x.chromosome_name]]
#' #Xist should be there, key indicator of correct DEG analysis
#' DEGs[gene_name=="XIST",c("celltype","adj_pval","lfc")]
#' #Present in all cell types
#' # There should also be high number of other DEGs on sex chromosomes, check:
#' our_degs[chromosome_name %in% c("X","Y"),
#'          .N,by=.(celltype,chromosome_name)]
#'}
#' @import ggplot data.table
#' @export
sc_cell_type_de <- function(SCE, design, pseudobulk_ID, celltype_ID, y=NULL,
                            region="single_region", coef=NULL, control=NULL,
                            pval_adjust_method = "BH", adj_pval=0.05,
                            folder="sc_cell_type_de_graphs/", 
                            rmv_zero_count_genes=TRUE, verbose=F){
  #source necessary functions from other scripts
  # source("../scripts/validate_input_parameters_de.R")
  # source("../scripts/make_pseudobulk.R")
  # source("../scripts/de_analysis.R")
  # source("../scripts/plot_de_analysis.R")
  
  #need to load SCE if a directory is passed
  if(class(SCE)[1]=="character"){
    if(!file.exists(SCE))
      stop(paste0("Directory to SCE file doesn't exist, please check ",
                  "your file location"))
    #load the dataset
    if(substr(SCE,nchar(SCE)-2,nchar(SCE))==".qs"){
      #load QS objects  
      SCE <- qs::qread(SCE)
    }else if (substr(SCE,nchar(SCE)-3,nchar(SCE))==".rds"){
      # rds object  
      SCE <- readRDS(SCE)
    }else{
      #normal R object - saved object must also be called SCE for this 
      #to work, will hit error in input check otherwise
      load(SCE)
    }
  }
  
  #sense check inputs
  validate_input_parameters_de(SCE, design, pseudobulk_ID, celltype_ID, y, 
                               region, coef, control, pval_adjust_method, 
                               adj_pval, folder, rmv_zero_count_genes, 
                               verbose)
  
  #get counts of each cell type 
  #counts(SCE) <- as.matrix(counts(SCE))
  counts_celltypes <- SCE[[celltype_ID]]
  counts_celltypes <-as.vector(table(counts_celltypes))
  names(counts_celltypes) <- names(table(SCE[[celltype_ID]]))
  
  #first format formula
  design_txt <- paste0(deparse(design,width.cutoff = 500),collapse=',')
  #make design formula minus anything before ~
  formula <- as.formula(gsub(".*~","~",design_txt))
  #if exists remove everything before ~ this format is necessary for glm_gp()
  design_txt <- gsub(".*~","",design_txt)
  #split design by + to get components to bring forward for pseudobulk
  #add in pseudobulk by column
  pb_columns <- c(gsub("^\\s+|\\s+$","",strsplit(design_txt, "[+]")[[1]]),
                  pseudobulk_ID)
  #if y not specified take last value in design matrix
  if(is.null(y))
    y <- design_txt[[length(design_txt)]]
  #Check if continuous or categorical variable to be modeled
  y_contin <- FALSE
  if(is.numeric(SCE[[y]]))
    y_contin <- TRUE
  
  #Get pseudobulk values
  celltypes <- unique(SCE[[celltype_ID]])
  if(isTRUE(verbose))
    message("Deriving pseudobulk data")
  pb_dat <-
    lapply(celltypes,function(x) 
      make_pseudobulk(SCE[,SCE[[celltype_ID]]==x], 
                      pseudobulk_ID=pseudobulk_ID, 
                      pb_columns=pb_columns))
  names(pb_dat) <- celltypes
  
  #run edgeR LRT DE analysis
  celltype_de <-
    de_analysis(pb_dat,formula,y_name=y,y_contin,coef,control,
                pval_adjust_method, adj_pval, verbose)
  #get sig DEGs for each
  celltype_DEGs <- lapply(celltype_de, function(x) x[x$adj_pval<adj_pval,])
  
  unique_genes <-lapply(celltype_DEGs, function(x) x$name)
  unique_degs <- unique(unlist(unique_genes))
  
  #if no DEGs found break but return the DE analysis scores still
  if(length(unique_degs)==0)
    return(list("celltype_DEGs"=celltype_DEGs,
                "celltype_all_genes"=celltype_de,
                "celltype_counts"=counts_celltypes))
  
  if(isTRUE(verbose))
    message(length(unique_degs)," unique DEGs foundacross all cell types")
  
  celltype_all_genes_dt <- data.table::rbindlist(celltype_de,idcol = T)
  setnames(celltype_all_genes_dt,".id","celltype")
  celltype_DEGs_dt <- data.table::rbindlist(celltype_DEGs,idcol = T)
  setnames(celltype_DEGs_dt,".id","celltype")
  
  
  if(!isFALSE(folder)){
    if(isTRUE(verbose))
      message("Plotting the results of differential expression analysis")
    #make plots for DE analysis
    plot_de_analysis(pb_dat,y,celltype_DEGs_dt,celltype_all_genes_dt,
                     counts_celltypes,folder)
  }
  
  #get gene names for DEGs
  genes <- unique(data.table::rbindlist(celltype_DEGs)$name)
  gene_IDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, 
                                keytype = "GENEID", 
                                columns = c("GENEID","SYMBOL"))
  colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
  gene_IDs <- data.table::as.data.table(gene_IDs)
  data.table::setnames(gene_IDs,"ensembl_gene_id","name")
  #remove any dups in the reference set - two names for one ENSEMBL ID
  gene_IDs <- unique(gene_IDs,by="name")
  celltype_DEGs <-
    lapply(celltype_DEGs,
           function(i){setDT(i);setkey(i,name);
             i[,gene_name:=gene_IDs[i, on=.(name),x.hgnc_symbol]];
             setnames(i,"name","ensembl_id");setnames(i,"gene_name","HGNC")})
  
  #get gene names for all genes
  genes <- unique(data.table::rbindlist(celltype_de)$name)
  gene_IDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, 
                                keytype = "GENEID", 
                                columns = c("GENEID","SYMBOL"))
  colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
  gene_IDs <- data.table::as.data.table(gene_IDs)
  data.table::setnames(gene_IDs,"ensembl_gene_id","name")
  #remove any dups in the reference set - two names for one ENSEMBL ID
  gene_IDs <- unique(gene_IDs,by="name")
  celltype_de <-
    lapply(celltype_de,
           function(i){setDT(i);setkey(i,name);
             i[,gene_name:=gene_IDs[i, on=.(name),x.hgnc_symbol]];
             setnames(i,"name","ensembl_id");setnames(i,"gene_name","HGNC")})
  
  return(list("celltype_DEGs"=celltype_DEGs,
              "celltype_all_genes"=celltype_de,
              "celltype_counts"=counts_celltypes))
}



# make_pseudobulk
# Calculate the summed pseudobulk values for an SCE object based on one 
# single cell type only. Ensure to filter SCE to pass one cell type's data.
#' @param SCE SingleCellExperiment object, a specialised S4 class for storing 
#' data from single-cell experiments. 
#' @param design the design formula of class type `formula`. Equation used to 
#' fit the model- data for the generalised linear model.
#' @param pseudobulk_ID the column name in the SCE object to perform pseudobulk 
#' on, usually the patient identifier. This column is used for grouping in the 
#' pseudobulk approach
#' @param pb_columns vector, list of annotation column names in the SCE object 
#' to be returned in annot_pb rolled up to pseudobulk level. Default is NULL 
#' which won't return any information in annot_pb.
#' @param region The column name in the SCE object for the study region. If 
#' there are multiple regions in the study (for example two brain regions). 
#' Pseudobulk values can be derived separately. Default is "single_region" 
#' which will not split by region.
#' @param rmv_zero_count_genes whether genes with no count values in any cell 
#' should be removed. Default is TRUE
#' @return A list containing
#' \itemize{
#'   \item \code{sumDat}: matrix of the summed pseudobulk count values
#'   \item \code{annot_pb}: dataframe of the annotation data from the SCE 
#'   rolled up based on the pseudobulk aggregation.
#' }
#' @export
make_pseudobulk <- function(data,pseudobulk_ID, pb_columns=NULL,
                            region="single_region",rmv_zero_count_genes=TRUE){
    allAnnot <- colData(data)
    if(region=="single_region") #constant value for all samples
        allAnnot[[region]] <- "one_region"
    indvs <- as.character(unique(allAnnot[[pseudobulk_ID]]))
    regions <- as.character(unique(allAnnot[[region]]))
    sumDat <- 
        matrix(0,nrow=dim(data)[1],
               ncol=length(indvs)*length(regions))
    colnames(sumDat) <- rep(indvs,length(regions))
    rownames(sumDat) <- rownames(data)
    count=0
    #for a single cell type, look at a single individual in single brain region
    #create list to hold annot values for each region, individual combination 
    annot_list <- vector(mode="list",length=length(regions)*length(indvs))
    names(annot_list) <- paste0(indvs,"_",regions)
    for(region_i in regions){
        for(indv in indvs){
            whichCells = 
                allAnnot[[pseudobulk_ID]]==indv & allAnnot[[region]]==region_i
            theData <- data[,whichCells]
            count <- count+1
            sumDat[,count] <- Matrix::rowSums(counts(theData))
            colnames(sumDat)[count] = sprintf("%s_%s",indv,region_i)
            #Get annotation data
            print_warning <- FALSE
            if(!is.null(pb_columns)){
                annot_i <- colData(theData)
                #there should only be a single value for each variable since 
                #this is single cell type, brain region and person
                #throw warning if not the case for numeric & aggregate 
                #accordingly. Throw error if categorical
                #first restrict data to just those from the design formula
                annot_i <- annot_i[,names(annot_i) %in% pb_columns]
                #change all factor variables to characters
                facts <- sapply(annot_i, is.factor)
                annot_i[facts] <- lapply(annot_i[facts], as.character)
                #check if values are unique
                cols_sum <- 
                    names(annot_i)[sapply(annot_i, 
                                          function(x) length(unique(x))>1)]
                if(length(cols_sum)>=1){
                    print_warning <- TRUE
                    cols_sum_warn <- cols_sum
                    #get numeric columns and categorical separately to deal with 
                    num_cols <- 
                        unlist(lapply(annot_i[,cols_sum,drop = FALSE], 
                                        is.numeric))  
                    #if any categorical throw error, shouldn't be lower level
                    cat_cols <- 
                        names(annot_i[,cols_sum,drop=FALSE][,!num_cols,drop=FALSE])
                    if(length(cat_cols)>0)
                        stop(paste0(c("Your design for the DE analysis contain",
                                      "s ",length(cat_cols)," categorical",
                                      " variable(s) which don't have a unique",
                                      " value at the celltype-region-individua",
                                      "l level. The non-unique column(s) are: ",
                                      cat_cols)))
                    annot_i_pb <- c(
                        as.list(unique(annot_i[,!(names(annot_i) %in% 
                                                        cols_sum)])),
                        sapply(annot_i[,cols_sum][,num_cols], 
                                function(x) mean(x))
                    )
                    #rearrange order to input order
                    annot_i_pb<-annot_i_pb[names(annot_i)]
                }
                else{
                    #all values unique, just get unique values
                    annot_i_pb <- as.list(unique(annot_i))
                }
                annot_list[[paste0(indv,"_",region_i)]] <- unlist(annot_i_pb)
            }    
        }
    }
    if(print_warning){
        print(paste0("Warning: The following cell level data were aggregated ", 
                     "to brain region, cell type, individual level where ",
                     "the mean for these numeric values will be ",
                     "taken:"))
        print(cols_sum_warn)
    }
    #create annotation dataframe from list
    annot_df <- as.data.frame(do.call(rbind, annot_list))
    annot_df$group_sample <- rownames(annot_df)
    rownames(annot_df) <- NULL
    #remove genes with 0 counts
    if(rmv_zero_count_genes)
        sumDat <- sumDat[rowSums(sumDat)!=0,]
    
    return(list("sumDat"=sumDat,"annot_pb"=annot_df))
}

#' Differential Expression Analysis using edgeR LRT on pseudobulk data
#'
#' @param pb_dat A list containing
#' \itemize{
#'   \item \code{sumDat}: matrix of the summed pseudobulk count values
#'   \item \code{annot_pb}: dataframe of the annotation data from the SCE 
#'   rolled up based on the pseudobulk aggregation.
#' }
#' @param formula the validated design formula of class type `formula`. 
#' Equation used to fit the model- data for the generalised linear model 
#' e.g. ~ sex + pmi + disease.
#' @param y_name the column name in the SCE object for the return variable 
#' e.g. "diagnosis" - Case or disease. y can be discrete (logisitc regression) 
#' or continuous (linear regression)
#' @param y_contin is the variable being modelled continuous e.g. if 
#' case/control then TRUE if level of Tau in AD study then FALSE
#' @param coef character specifying which level to investigate for the 
#' differential expression analysis e.g. in a case/control study use "case" if 
#' case is the identifier in the y column to get positive fold changes to 
#' relate to case samples. Leave as default value for continuous y.
#' @param control character specifying which control level for the differential 
#' expression analysis e.g. in a case/control/other study use "control" in the 
#' y column to compare against. NOTE only need to specify if more than two 
#' groups in y, leave as default value for two groups or continuous y. 
#' @param pval_adjust_method the adjustment method for the p-value in the 
#' differential expression analysis. Default is benjamini hochberg "BH". See  
#' stats::p.adjust for available options
#' @param adj_pval the adjusted p-value cut-off for the differential expression 
#' analysis, 0-1 range
#' @param verbose logical indicating if extra information about the 
#' differential expression analysis should be printed
#' @return A list containing differential expression data (a dataframe) for 
#' each cell type. The dataframe contains log fold change (logFC), log counts 
#' per million (logCPM), log ratio (LR), p-value (PValue), adjusted p-value 
#' (adj_pval)
de_analysis <- function(pb_dat,formula,y_name,y_contin,coef, control,
                        pval_adjust_method, adj_pval,verbose){
  #Use Likelihood ratio test from edgeR, found to have best perf in a recent
  #benchmark: https://www.biorxiv.org/content/10.1101/2021.03.12.435024v1.full
  #Other option for edgeR is quasi-likelihood F-tests
  #Run each cell type through
  DE_ct_rtrn <- vector(mode="list",length=length(names(pb_dat)))
  names(DE_ct_rtrn) <- names(pb_dat)
  for(ct_i in names(pb_dat)){
    if(verbose)
      message("Analysing ",ct_i)
    targets <- pb_dat[[ct_i]]$annot_pb
    if(!y_contin){#adjust levels based on user input if not contin
      #chec if factor inputted, if not convert
      if(!is.factor(targets[[y_name]])){
        targets[[y_name]] <- as.factor(targets[[y_name]])
      }
      #reorder levels so coef is last so won't be intercept and chosen
      old_levels <- levels(targets[[y_name]])
      if(length(old_levels)>2){
        if(!is.null(coef)){#if user inputted value for case samples
          if(is.null(control)){
            message(paste0("WARNING: No control passed so coef may",
                           " not be compared against desired ",
                           "variable"))
            new_levels <- 
              c(old_levels[which(coef!=old_levels)],coef)
          }
          else{
            new_levels <- 
              c(control,
                old_levels[which(coef!=old_levels & 
                                   control!=old_levels)],coef)
          }
        }
        else{
          if(is.null(control)){
            message(paste0("WARNING: No control passed so may not ",
                           "be compared against desired variable"))
            new_levels <- old_levels #do nothing no coef or control
          }
          else{
            new_levels <- 
              c(control,old_levels[which(control!=old_levels)])
          }  
        }
      }
      else{
        if(!is.null(coef)){#if user inputted value for case samples
          new_levels <- c(old_levels[which(coef!=old_levels)],coef)
        }
        else{
          new_levels <- old_levels #do nothing no coef
        }
      }
      #assign the new levels
      targets[[y_name]] <- factor(targets[[y_name]],levels=new_levels)
    }
    # create design
    design <- model.matrix(formula, data = targets)   
    #slight diff approach for contin or cat
    if(y_contin){
      y <- edgeR::DGEList(counts = pb_dat[[ct_i]]$sumDat)
    }
    else{#categorical i.e. case/control
      y <- edgeR::DGEList(counts = pb_dat[[ct_i]]$sumDat,
                          group = targets[[y_name]])
    }
    #Calculate normalization factors to scale the raw library sizes.
    y <- edgeR::calcNormFactors(y,method = 'TMM')
    #Maximizes the negative binomial likelihood to give the estimate of the 
    #common,trended and tagwise dispersions across all tags.
    y <- edgeR::estimateDisp(y, design)
    #LRT
    fit <- edgeR::glmFit(y, design = design)
    #slight diff approach for contin or cat
    if(y_contin){#need to specify coefficient target as not done in design
      test <- edgeR::glmLRT(fit,coef=y_name)
    }
    else{#categorical i.e. case/control, no need
      if(!is.null(coef)){#if user inputted value for case samples
        test <- edgeR::glmLRT(fit,coef=length(new_levels))
      }
      else{#they didn't direction may not be as they expect
        message(paste0("WARNING: No coefficient passed so direction ma",
                       "y not relate to desired variable"))
        test <- edgeR::glmLRT(fit)
      }
    }
    pvals <- test$table
    #add adj p-values
    pvals$adj_pval <- stats::p.adjust(pvals$PValue, 
                                      method = pval_adjust_method)
    pvals$name <- rownames(pvals)
    rownames(pvals) <- NULL
    if(verbose){
      numDEGs <- nrow(pvals[pvals$adj_pval<adj_pval,])
      message(numDEGs," DEGs found")
    }
    DE_ct_rtrn[[ct_i]] <- pvals
  }
  return(DE_ct_rtrn)
}


#' Single Cell Cell Type Differential Expression Analysis Plots
#'
#' Create differential expression analysis plots. Run by sc_cell_type_de() 
#' @param pb_dat A list containing
#' \itemize{
#'   \item \code{sumDat}: matrix of the summed pseudobulk count values
#'   \item \code{annot_pb}: dataframe of the annotation data from the SCE 
#'   rolled up based on the pseudobulk aggregation.
#' }
#' @param y the column name in the SCE object for the return variable e.g. 
#' "diagnosis" - Case or disease. y can be discrete (logisitc regression) or 
#' continuous (linear regression)
#' @param celltype_DEGs_dt data table containing the DEGs for each cell type 
#' with their differential expression data
#' @param celltype_all_genes_dt data table containing the all genes for each 
#' cell type with their differential expression data
#' @param celltype_counts vector with the counts of cells after QC in each cell 
#' type
#' @param folder the folder where the graphs from the differential expression 
#' analysis are saved.
#' @param pal colour pallete to use for plots. Needs a vector of 6 colours 
#' (hex). By default, uses a Wes Anderson colour palette.
#' @return NULL
plot_de_analysis <- function(pb_dat,y,celltype_DEGs_dt,celltype_all_genes_dt,
                             counts_celltypes,folder,
                             pal=c(wesanderson::wes_palette("Royal2"),
                                   wesanderson::wes_palette("Moonrise3")[1])){
  logFC = name = NULL
  
  #if it doesn't exist already make folder for plots
  dir.create(folder,showWarnings = F)
  
  celltype_DEGs_dt[,deg_direction:="Down"]
  celltype_DEGs_dt[logFC>0,deg_direction:="Up"]
  #get top 3 most sig p-values for each cell type
  top_up_down_genes <-celltype_DEGs_dt[
    celltype_DEGs_dt[,.I[adj_pval %in% sort(adj_pval)[1:3]],by=celltype]$V1]
  #plot pseudobulk expression
  #first filter cell type expression to genes of interest
  top_degs_pseudobulk_exp <- pb_dat[unique(top_up_down_genes$celltype)]
  for(ct_i in names(top_degs_pseudobulk_exp)){
    ct_genes <- top_up_down_genes[celltype==ct_i,name]
    #filter and melt to get long DF of genes of interest
    #/sum(pb_dat$Micro$sumDat)
    top_degs_pseudobulk_exp[[ct_i]] <- 
      reshape2::melt(top_degs_pseudobulk_exp[[ct_i]]$sumDat[ct_genes,,drop=FALSE])
  }
  #combine
  top_degs_pseudobulk_exp <- 
    data.table::rbindlist(top_degs_pseudobulk_exp,id="celltype")
  setnames(top_degs_pseudobulk_exp,c("celltype","name","group_sample",
                                     "expression"))
  #add on DEG direction
  top_degs_pseudobulk_exp[top_up_down_genes,deg_direction:=i.deg_direction,
                          on=c("name","celltype")]
  #add in y - first get annotation data in single datatable
  annot_dt <- lapply(pb_dat, function(x) x$annot_pb)
  annot_dt <- data.table::rbindlist(annot_dt,id="celltype")
  top_degs_pseudobulk_exp[annot_dt,phenotype:=get(y),
                          on=c("group_sample","celltype")]
  #add in gene names
  genes <- unique(as.character(top_degs_pseudobulk_exp$name))
  gene_IDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, 
                                keytype = "GENEID", 
                                columns = c("GENEID","SYMBOL"))
  colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
  gene_IDs <- data.table::as.data.table(gene_IDs)
  data.table::setnames(gene_IDs,"ensembl_gene_id","name")
  #remove any dups in the reference set - two names for one ENSEMBL ID
  gene_IDs <- unique(gene_IDs,by="name")
  data.table::setkey(top_degs_pseudobulk_exp,name)
  #append gene names
  top_degs_pseudobulk_exp[, gene_name := gene_IDs
                          [top_degs_pseudobulk_exp, on=.(name), 
                            x.hgnc_symbol]]
  #plot increase size A4
  top_degs_pseudobulk_exp_plot <-
    ggplot(top_degs_pseudobulk_exp[gene_name!=""&!is.na(deg_direction),], 
           aes(x = phenotype, y = expression,colour=deg_direction)) +
    geom_jitter(height=0) +
    stat_summary(fun.data = "mean_cl_normal",
                 #aes(shape="mean"), 
                 colour = "grey",
                 geom = "crossbar",
                 show.legend = T)+
    scale_shape_manual("", values=c("mean"="x"))+
    labs(y= "Sum pseudobulk expression counts (unnormalised)", 
         x = "Phenotype",
         colour="DEG direction") +
    facet_wrap(~ celltype+gene_name, scales = "free_y")+
    theme_cowplot()+
    theme(strip.text.x = element_text(size = 6),
          axis.text = element_text(size=6))+
    #scale_colour_viridis(discrete = T)
    scale_colour_manual(values=pal)
  #save the graph to folder
  suppressMessages(ggsave(path = folder,
                          filename = "Pseudobulk_exp_most_sig_genes.pdf",
                          plot=top_degs_pseudobulk_exp_plot,
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in"))
  suppressMessages(ggsave(path = folder,
                          filename = "Pseudobulk_exp_most_sig_genes.png",
                          plot=top_degs_pseudobulk_exp_plot,
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in"))
  
  #plot number of cells per cell type
  counts_celltypes_dt <- 
    data.table::data.table(celltype=names(counts_celltypes),
                           counts=counts_celltypes)
  cell_counts_plot <-
    ggplot(data=counts_celltypes_dt,
           aes(x=factor(celltype),y = counts, fill=factor(celltype)))+ 
    geom_bar(stat="identity")+
    labs(y= "Number of cells after QC", x = "Cell Type",fill="Cell Type") +
    geom_text(aes(x = celltype, y = counts, label = counts),
              vjust=-0.25,size=3) +
    theme_cowplot()+
    theme(axis.text = element_text(size=6))+
    #scale_fill_viridis(discrete = T)
    scale_fill_manual(values=pal)
  #save the graph to folder
  suppressMessages(ggsave(path = folder,
                          filename = "Cell_counts_after_QC.pdf",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=cell_counts_plot))
  suppressMessages(ggsave(path = folder,
                          filename = "Cell_counts_after_QC.png",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=cell_counts_plot))
  
  #plot DEGs per cell type
  deg_per_cell_type_plot <-
    ggplot(data=celltype_DEGs_dt[,.N,by=celltype],
           aes(x=factor(celltype),y = N, fill=factor(celltype)))+ 
    geom_bar(stat="identity")+
    labs(y= "Number of DEGs identified", x = "Cell Type",fill="Cell Type") +
    geom_text(aes(x = celltype, y = N, label = N),
              vjust=-0.25,size=3) +
    theme_cowplot()+
    theme(axis.text = element_text(size=9))+
    #scale_fill_viridis(discrete = T)
    scale_fill_manual(values=pal)
  #save the graph to folder
  suppressMessages(ggsave(path = folder,
                          filename = "DEGs_per_cell_type.pdf",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=deg_per_cell_type_plot))
  suppressMessages(ggsave(path = folder,
                          filename = "DEGs_per_cell_type.png",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=deg_per_cell_type_plot))
  
  #plot degs as proportion of cell numbers
  degs_prop <- celltype_DEGs_dt[,.N,by=celltype]
  data.table::setkey(degs_prop,celltype)
  degs_prop[names(counts_celltypes),num_cells:=counts_celltypes]
  degs_prop[,prop:=N/num_cells]
  degs_prop[,N_prop:=num_cells/sum(num_cells)]
  
  deg_prop_plot <-
    ggplot(data=degs_prop,
           aes(x=factor(celltype)))+ 
    geom_bar(aes(y = prop,fill=factor(celltype)),stat="identity")+
    labs(y= "Proportion of DEGs identified", x = "Cell Type",
         fill="Cell Type") +
    geom_text(aes(x = factor(celltype), y=prop, 
                  label = scales::percent(prop)),
              vjust=-0.25,size=3) +
    theme_cowplot()+
    theme(axis.text = element_text(size=9))+
    #scale_fill_viridis(discrete = T)
    scale_fill_manual(values=pal)
  #save the graph to folder
  suppressMessages(ggsave(path = folder,
                          filename = "DEGs_proportion_per_cell_type.pdf",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=deg_prop_plot))
  suppressMessages(ggsave(path = folder,
                          filename = "DEGs_proportion_per_cell_type.png",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=deg_prop_plot))
  
  
  #plot DEG direction and effect size
  celltype_DEGs_dt[,deg_direction:="Down"]
  celltype_DEGs_dt[logFC>0,deg_direction:="Up"]
  #count
  count_deg_direction_plot<-
    ggplot(data=celltype_DEGs_dt[,.N,by=.(celltype,deg_direction)],
           aes(x = factor(deg_direction), y = N, fill=factor(celltype)))+ 
    geom_bar(stat="identity") + 
    facet_wrap(~factor(celltype))+
    ggtitle("Cell type counts of DEGs")+
    labs(y= "Log2 Fold Change", x = "DEG direction", fill="Cell Type") +
    theme_cowplot()+
    theme(axis.text = element_text(size=9))+
    #scale_fill_viridis(discrete = T)
    scale_fill_manual(values=pal)
  #save the graph to folder
  suppressMessages(
    ggsave(path = folder,
           filename="DEGs_direction_count_per_cell_type.pdf",
           dpi = 1200,width = 8.65,
           height = 10.0, units ="in",
           plot=count_deg_direction_plot))
  suppressMessages(
    ggsave(path = folder,
           filename="DEGs_direction_count_per_cell_type.png",
           dpi = 1200,width = 8.65,
           height = 10.0, units ="in",
           plot=count_deg_direction_plot))
  
  #plot volcano plots
  #get gene names - for sig genes
  genes <- unique(celltype_all_genes_dt[adj_pval<0.05,]$name)
  gene_IDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= genes, 
                                keytype = "GENEID", 
                                columns = c("GENEID","SYMBOL"))
  colnames(gene_IDs) <- c("ensembl_gene_id","hgnc_symbol")
  gene_IDs <- data.table::as.data.table(gene_IDs)
  data.table::setnames(gene_IDs,"ensembl_gene_id","name")
  #in case any ID's have more than one hgnc name found 
  gene_IDs <- gene_IDs[!duplicated(gene_IDs$name),]
  data.table::setkey(gene_IDs,name)
  data.table::setkey(celltype_all_genes_dt,name)
  #append gene names
  celltype_all_genes_dt[, gene_name := 
                          gene_IDs[celltype_all_genes_dt, x.hgnc_symbol]]
  #add colour identifier
  celltype_all_genes_dt[,colour_ident:="Not Significant"]
  celltype_all_genes_dt[adj_pval<0.05 & logFC<0 ,
                        colour_ident:="Decreased DEG"]
  celltype_all_genes_dt[adj_pval<0.05 & logFC>0 ,
                        colour_ident:="Increased DEG"]
  
  cols <- c("Not Significant" = "grey", "Increased DEG" = "#FDE725FF", 
            "Decreased DEG" = "#440154FF")
  #remove cell types without sig genes
  celltype_all_genes_dt <-
    celltype_all_genes_dt[(celltype %in% 
                             unique(celltype_DEGs_dt$celltype)),]
  
  volcano_plot<-
    ggplot(data=celltype_all_genes_dt,
           aes(x = logFC, y = -log10(adj_pval)))+ #scale FDR by log10
    geom_point(aes(colour=colour_ident),stat="identity",size=.2) + 
    facet_wrap(~factor(celltype),scales = "free")+
    geom_hline(yintercept = -log10(0.05), colour="#990000", 
               linetype="dashed",size=.3) + 
    #geom_vline(xintercept = 0, colour="black", size=.3) + 
    labs(y= "-log10(FDR)", x = "Log2 Fold Change", colour="") +
    theme_cowplot()+
    theme(axis.text = element_text(size=9),axis.text.x=element_blank(),
          legend.text=element_text(size=10))+
    scale_colour_manual(values = cols)+
    geom_text_repel( #give top three genes per cell type by adjusted p value
      data = celltype_all_genes_dt[
        celltype_all_genes_dt[,.I[adj_pval %in% 
                                    sort(adj_pval)[1:3]][1:3],
                              by=celltype]$V1],
      aes(label = gene_name),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  #save the graph to folder - increase size A4
  suppressMessages(ggsave(path = folder,
                          filename = "volcano_plot_degs_cell_types.pdf",
                          plot=volcano_plot,
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in"))
  suppressMessages(ggsave(path = folder,
                          filename = "volcano_plot_degs_cell_types.png",
                          plot=volcano_plot,
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in"))
  
  #plot boxplots of DEGs
  #check median log2 fold change for both up and down regulated genes
  #This will give a sense of directional effect size
  degs_median_lfc<-celltype_DEGs_dt[,median(logFC),by=deg_direction]
  #if only up or down reg DEGs add 0 for other
  if(nrow(degs_median_lfc)<2){
    dir <- c("Up","Down")
    degs_median_lfc<-
      rbind(degs_median_lfc,
            data.table::data.table("deg_direction"=
                                     dir[dir!=degs_median_lfc$deg_direction],
                                   "V1"=0))
    #set order so Up first
    data.table::setorder(degs_median_lfc,-deg_direction)
  }
  
  # Compute boxplot statistics
  calc_boxplot_stat <- function(x) {
    coef <- 1.5
    n <- sum(!is.na(x))
    # calculate quantiles
    stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
    iqr <- diff(stats[c(2, 4)])
    # set whiskers
    outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
    if (any(outliers)) {
      stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
    }
    return(stats)
  }
  
  degs_boxplot_plots<-
    ggplot(data=celltype_DEGs_dt,
           aes(x = factor(deg_direction), y = abs(logFC), 
               fill=factor(celltype)))+ 
    stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
    facet_wrap(~factor(celltype),scales="free")+
    ggtitle("Cell type Log2 Fold Change of DEGs",
            subtitle = paste0("Median LFC of Down and Up regulated DEGs: ",
                              round(degs_median_lfc$V1[[2]]*-1,2),", ",
                              round(degs_median_lfc$V1[[1]],2)))+
    labs(y= "Log2 Fold Change", x = "DEG direction", fill="Cell Type") +
    theme_cowplot()+
    theme( axis.text = element_text(size=9))+
    #scale_fill_viridis(discrete = T)
    scale_fill_manual(values=pal)
  #save the graph to folder
  suppressMessages(ggsave(path = folder,
                          filename = "deg_boxplots_cell_types.pdf",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=degs_boxplot_plots))
  suppressMessages(ggsave(path = folder,
                          filename = "deg_boxplots_cell_types.png",
                          dpi = 1200,width = 8.65,
                          height = 10.0, units ="in",
                          plot=degs_boxplot_plots))
} 

