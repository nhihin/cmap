# The following functions are used for our CMap workflow.

# The following function patches the metadata of an input
# Seurat object to ensure 3.2 month old samples are labelled as
# "young" rather than "aged" in the age_binary column.
# This is only relevant for datasets containing samples
# from the batch tg07.
patchSeuratMetadataAge <- function(seurat){
  library(Seurat)
  library(magrittr)
  library(dplyr)

  seurat@meta.data <- seurat@meta.data %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::mutate(age_binary = case_when(
      seq_batch == "tg07" & age_prep == "3.2" ~ "young",
      TRUE ~ age_binary
    )) %>%
    tibble::column_to_rownames("barcode")

  return(seurat)
}


# The following function adds a column "pool_id" to the
# metadata of a Seurat object, meant for pseudobulk-
# aggregation. The pool_id column is structured in the
# following way:
# [grouping_column]_pool[biological replicate number]_n[number of animals within the pool]
addPoolingColumnToSeuratMetadata <- function(seurat,
                                             grouping_column = "treatment_label",
                                             samplesheet_path ="/mnt/home/alkahest.com/dleone/Documents/Analyses/library_processors/nextera_processors/FACS_annotations/FACS_10x_samplesheet_2022-10-27.tsv"){
  library(Seurat)
  library(readr)
  library(dplyr)
  library(tibble)
  library(magrittr)

  samplesheet <- readr::read_tsv(samplesheet_path) %>%
    dplyr::select(sample_ID, FACS_date, animal_ID)

  seurat@meta.data %<>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(samplesheet, by = c("sample_ID", "FACS_date")) %>%
    tibble::column_to_rownames("barcode")

  animalToPool <- seurat@meta.data %>%
    as.data.frame %>%
    dplyr::select(animal_ID, .data[[grouping_column]], seq_batch) %>%
    dplyr::distinct() %>%
    set_rownames(NULL) %>%
    dplyr::arrange(.data[[grouping_column]])

  nanimals <- animalToPool$animal_ID %>%
    lapply(function(x){
      as.list(strsplit(x[[1]], ","))[[1]] %>% length()
    })

  animalToPool$n_animals <- nanimals

  ngroups <- animalToPool[[grouping_column]] %>% table %>% as.list()

  tmp <- animalToPool %>%
    dplyr::group_by(.data[[grouping_column]]) %>%
    dplyr::summarise(biol_rep=n())

  animalToPool %<>%
    left_join(tmp, by = grouping_column)

  animalToPool$pool_id <- paste0(animalToPool[[grouping_column]], "_",
                                 "pool",
                                 animalToPool$biol_rep,
                                 "_n",
                                 animalToPool$n_animals)

  seurat@meta.data %<>%
    tibble::rownames_to_column("barcode") %>%
    left_join(animalToPool, by = c("animal_ID", grouping_column, "seq_batch")) %>%
    tibble::column_to_rownames("barcode")

  return(seurat)
}


# The following function takes a Seurat object
# and returns a SCE object that has been
# log-normalized and filtered such that only
# genes with a minimum count across all cells
# are retained. Uses a lot of memory.
convertSeuratToSCE <- function(seurat,
                               minCountsToRetainGene = 10){

  counts <- GetAssayData(object = seurat,
                         slot = "counts",
                         assay = "RNA")

  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = seurat@meta.data)

  rm(counts)
  gc()

  # Calculate size factor & normalize by
  # dividing counts for each cell by the
  # size factor.
  # sf <- 2^rnorm(ncol(sce))
  # sf <- sf/mean(sf)
  # normcounts(sce) <- t(t(counts(sce))/sf)
  # logcounts(sce) <- log2(normcounts(sce)+1)
  # Quicker way of doing the above:
  logcounts(sce) <- normalizeCounts(sce,
                                    log=TRUE,
                                     pseudo.count = 1)

  if(minCountsToRetainGene > 0){
    sce <- sce[rowSums(counts(sce) > 1) >= minCountsToRetainGene, ]
    #rm(sf)
    gc()
  } else{
    #rm(sf)
    gc()
  }

  return(sce)
}

# The following function takes an SCE object
# and converts it to an aggregated matrix suitable
# for pseudobulk analyses. Gene counts are summed for
# each sample as defined by `pseudobulk_grouping_column`.
aggregateSCEIntoPseudobulkMatrix <- function(sce,
                                             pseudobulk_grouping_column = "sample_ID"){

  require(Matrix.utils)

  # Named vector of cluster names
  cluster_ids <- purrr::set_names(levels(as.factor(sce$Major_celltype)))

  # Total number of clusters
  nc <- length(cluster_ids)

  # Named vector of sample names
  sample_ids <- purrr::set_names(levels(as.factor(sce[[pseudobulk_grouping_column]])))

  # Total number of samples
  ## 11
  ns <- length(sample_ids)

  # Generate sample level metadata


  ## Turn class "table" into a named vector of cells per sample
  n_cells <- table(sce[[pseudobulk_grouping_column]]) %>%  as.vector()
  names(n_cells) <- names(table(sce[[pseudobulk_grouping_column]]))

  ## Match the named vector with metadata to combine it
  m <- match(names(n_cells), sce[[pseudobulk_grouping_column]])

  ## Create the sample level metadata by selecting specific columns
  ei <- data.frame(colData(sce)[m, ],
                   n_cells, row.names = NULL) %>%
    dplyr::select(.data[[pseudobulk_grouping_column]],  "n_cells")

  groups <- colData(sce)[, c("Major_celltype", pseudobulk_grouping_column)]
  groups[[pseudobulk_grouping_column]] <- factor(groups[[pseudobulk_grouping_column]])

  # Aggregate across cluster-sample groups
  # Each row corresponds to aggregate counts for a cluster-sample combo
  pb <- Matrix.utils::aggregate.Matrix(t(logcounts(sce)),
                                       groupings = groups, fun = "sum")

  return(pb)
}

# The following function takes an aggregated matrix and
# splits it into a list of matrices based on cell type.
# Assumed that running this function after
# `aggregateSCEIntoPseudobulkMatrix()`. Not sure what
# behaviour would be like if running in isolation due to
# hardcoded splitting strings.
splitAggregatedMatrixPerCelltype <- function(mat){
  library(magrittr)

  pb <- mat
  splitf <- sapply(
    stringr::str_split(rownames(pb),
                       pattern = "_",
                       n = 2),
    `[`, 1)

  pb <- split.data.frame(pb,
                         factor(splitf))

  pb <- pb %>%
    lapply(function(x){
      sample_names <-
        rownames(x) %>% gsub(
          x = .,
          pattern = "^.*_(FACS.*)$", replacement = "\\1"
        )
      set_colnames(t(x), sample_names)
    })

  return(pb)
}

#test3 <- splitAggregatedMatrixPerCelltype(test2)

# The following function takes as input a list of
# pseudobulk matrices (assumed for each cell type)
# and further splits these by treatment. It then applies
# z-score to the values (assumed to be log-norm counts),
# followed by doing a weighted average across all replicates,
# if zscore is TRUE.
splitPerCelltypeAggregatedMatrixByGroup <- function(sce,
                                                    mat_list,
                                                    grouping_column = "treatment_label",
                                                    pseudobulk_grouping_column = "sample_ID",
                                                    zscore = TRUE){

  require(cmapR)

  map_df <- colData(sce) %>%
    as.data.frame %>%
    dplyr::group_by(.data[[pseudobulk_grouping_column]], .data[[grouping_column]]) %>%
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(sample_treatment = paste0(.data[[pseudobulk_grouping_column]], "_", .data[[grouping_column]]))

  for(celltype in names(mat_list)){
    oldColnames <- colnames(mat_list[[celltype]])
    newColnames <- data.frame(sample_ID = oldColnames) %>%
      magrittr::set_colnames(pseudobulk_grouping_column) %>%
      dplyr::left_join(map_df, by = pseudobulk_grouping_column) %>%
      dplyr::pull(sample_treatment)

    colnames(mat_list[[celltype]]) <- newColnames
  }

  groupsOfInterest <- unique(map_df[[grouping_column]])
  for(celltype in names(mat_list)){
    tmp_list <- list()

    for(group in groupsOfInterest){
      tmp_list[[group]] <- mat_list[[celltype]] %>%
        as.data.frame %>%
        dplyr::select(contains(group))

      if(ncol(tmp_list[[group]]) == 0){
        tmp_list[[group]] <- NULL
        next
      }

      if(zscore == TRUE){

        capLargeZscores <- function(x){
          x[x > 10] <- 10
          x[x < -10] <- -10
          return(x)
        }

        zscorefunc <- function(x){
          x %>%
            cmapR::robust_zscore() %>%
            capLargeZscores()
        }

        tmp_list[[group]] <-  apply(tmp_list[[group]],
                                    2,
                                    zscorefunc)

        if(ncol(tmp_list[[group]]) > 1){
          tmp_list[[group]] <- cmapR::distil(tmp_list[[group]],
                                             dimension = "col")$values
        }

      }
    }
    if(length(tmp_list) == 0L){
      mat_list[[celltype]] <- NULL
    } else{
      mat_list[[celltype]] <- tmp_list
    }

  }
  return(mat_list)
}

# Absmax, returns the absolute maximum value
# while preserving the sign, important for zscores and
# other directional values.
absmax <- function(x){
  x[which.max( abs(x) )]
}

# The following function takes a data.frame structured as follows:
# rownames are mouse ensembl IDs, and another 1 column of values
# such as counts, z-scores, etc. (can be named anything).
# The mouse ensembl IDs are converted to human entrez IDs,
# and any many to 1 mappings are
# dealt with using the `functionForReplicates`, which can be
# any function, but recommended to be
# "sum", "mean", or "absmax". Any mouse ensembl IDs without
# corresponding human entrez IDs are filtered out.
convertMouseEnsemblToHumanEntrezID <- function(df,
                                               functionForReplicates = "absmax"){

  # Only download files necessary for the mapping df once
  if("mm_ensembl_to_human_entrez.rds" %in% list.files(here("data"))){
    mm_ensembl_to_human_entrez <- readRDS(here("data","mm_ensembl_to_human_entrez.rds"))
  } else{
    require(uatools)
    require(org.Hs.eg.db)
    require(clusterProfiler)
    require(here)

    hs_to_mm_ensembl <- uatools::get_scaffold("hs_to_mm_orthologs")

    human_ensembl_to_entrez  <- clusterProfiler::bitr(
      geneID = hs_to_mm_ensembl$human_gene_id,
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = 'org.Hs.eg.db') %>%
      as.data.frame %>%
      set_colnames(c("human_gene_id", "human_entrezid"))

    mm_ensembl_to_human_entrez <- hs_to_mm_ensembl %>%
      as.data.frame %>%
      left_join(human_ensembl_to_entrez,
                by = "human_gene_id",
                multiple = "all") %>%
      dplyr::filter(!is.na(human_entrezid)) %>%
      dplyr::distinct(human_entrezid, mouse_gene_id)

    mm_ensembl_to_human_entrez %>%
      saveRDS(here("data", "mm_ensembl_to_human_entrez.rds"))
  }

  res <- df %>%
    tibble::rownames_to_column("mouse_gene_id") %>%
    dplyr::left_join(mm_ensembl_to_human_entrez,
                     multiple = "all") %>%
    dplyr::select(-mouse_gene_id)

  res <- plyr::ddply(res,
                     "human_entrezid",
                     plyr::numcolwise(functionForReplicates))

  res <- res %>%
    dplyr::filter(!is.na(human_entrezid)) %>%
    tibble::column_to_rownames("human_entrezid")

  return(res)
}


# Given a named list, rank the list in either
# ascending or descending order and return the
# top _n_ names (but not the values themselves)
top_n_names_from_list <- function(x,
                                  direction = "top",
                                  n = 150){
  ordered_values <- order(-x)
  names_of_ordered_values <- names(x[ordered_values])
  if(direction == "top"){
    res <- head(names_of_ordered_values, n)
  } else if(direction == "bottom"){
    res <- tail(names_of_ordered_values, n)
  } else {
    print("Direction not specified, must be either 'top' or 'bottom'")
  }
  return(res)
}

# # Input: A list of objects from running multiple signatureSearch
# # queries using `gess_lincs()`.
# # `groupingColname` is the name of the new column to be
# # created which will have the values of the names of the
# # list items provided in `ls`.
# # Output: A merged data.frame of all results, with a new
# # column `sig` indicating whether a particular query was
# # significant.
# combineResultsFromMultipleSSQueries <- function(ls,
#                                                 groupingColname = "query_celltype",
#                                                 onlyRetainDrugsWithTargetGenes = TRUE,
#                                                 FDRCutoffForSig = 0.05,
#                                                 NCSCutoffForSig = 0.8){
#   require(signatureSearch)
#   full_table <- ls %>%
#     lapply(result) %>%
#     dplyr::bind_rows(.id = groupingColname)
#
#   if(onlyRetainDrugsWithTargetGenes == TRUE){
#     full_table <- full_table %>%
#       dplyr::filter(!is.na(t_gn_sym))
#   }
#
#   full_table <- full_table %>%
#     dplyr::mutate(sig =   (trend == "up" &
#                              WTCS_FDR < FDRCutoffForSig &
#                              NCSct > NCSCutoffForSig))
#
#
# }

# Over-representation analysis using either
# Fisher's exact test (test = "fishers") or
# hypergeometric test (test = "hypergeometric")
# See https://stats.stackexchange.com/questions/288081/use-fishers-exact-test-or-a-hypergeometric-test
ORAUsingFishersExactTest <- function(full_table,
                                     groupingCol = "indication_class",
                                     test = "hypergeometric"){
  table_split_by_groupingCol <- split(full_table,
                                      f = full_table[[groupingCol]])
  sig_drugs_per_groupingCol <- table_split_by_groupingCol %>%
    lapply(function(x){
      x %>%
        dplyr::filter(sig == TRUE) %>%
        .$pert_iname %>% unique
    })
  sig_drugs_per_groupingCol2 <- data.frame(indication = names(sig_drugs_per_groupingCol))
  sig_drugs_per_groupingCol2$sig_drugs <- table_split_by_groupingCol

  # Prepare contingency matrices
  cont_tables <- table_split_by_groupingCol %>%
    lapply(function(x){
      n_sig_with_indication = sum(x$sig)
      n_notsig_with_indication = nrow(x) - n_sig_with_indication
      n_sig_without_indication = sum(full_table$sig) - n_sig_with_indication
      n_notsig_without_indication = nrow(full_table) - n_notsig_with_indication

      res = matrix(
        c(n_sig_with_indication, n_notsig_with_indication,
          n_sig_without_indication, n_notsig_without_indication),
        ncol = 2, byrow = TRUE
      )

      return(res)
    })

  # Calculate p-values for contingency matrices
  p_vals <- list()
  if(test == "fishers"){
    for(x in names(table_split_by_groupingCol)){
      p_vals[[x]] <- fisher.test(cont_tables[[x]])$p.value
    }
  } else if(test == "hypergeometric"){
    for(x in names(table_split_by_groupingCol)){
      p_vals[[x]] <- fisher.test(cont_tables[[x]],
                                 alternative = "less")$p.value
    }
  }

  res2 <- list(
    cont_matrices = cont_tables,
    p_vals = p_vals
  )

  # Add drug information to the final results
  mappingDf <- readRDS("/mnt/home/alkahest.com/nhin/Projects/cmap_internal/data/mappingTable2.rds")


  p_val_df <- data.frame(
    indication = names(p_vals),
    p = p_vals %>% unlist
  )%>%
    dplyr::mutate(BH_p = p.adjust(p, method="BH")) %>%
    dplyr::mutate(sig = BH_p < 0.05) %>%
    dplyr::arrange(p) %>%
    dplyr::left_join(mappingDf, by = "indication") %>%
    dplyr::rename(drugs_univ = drugs)# %>%
    #dplyr::rowwise() %>%
    #dplyr::filter(drugs = drugs[drugs %in% sig_drugs_per_groupingCol])

  return(p_val_df)
}

# Overlap between groups
overlapGroups <- function (listInput, sort = TRUE) {
  require(UpSetR)
  # listInput could look like this:
  # $one
  # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  # $two
  # [1] "a" "b" "d" "e" "j"
  # $three
  # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
  listInputmat    <- fromList(listInput) == 1
  #     one   two three
  # a  TRUE  TRUE  TRUE
  # b  TRUE  TRUE FALSE
  #...
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
    # attr(,"groups")
    #   one   two three
    # FALSE FALSE  TRUE
    #  f  i
    # 12 13
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

