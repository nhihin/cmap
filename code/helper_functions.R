# The following functions are used for our CMap workflow.

# The following function patches the metadata of an input
# Seurat object to ensure 3.2 month old samples are labelled as
# "young" rather than "aged" in the age_binary column.
# This is only relevant for datasets containing samples
# from the batch tg07.
patchSeuratMetadataAge <- function(seurat) {
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
                                             samplesheet_path = "/mnt/home/alkahest.com/dleone/Documents/Analyses/library_processors/nextera_processors/FACS_annotations/FACS_10x_samplesheet_2022-10-27.tsv") {
  library(Seurat)
  library(readr)
  library(dplyr)
  library(tibble)
  library(magrittr)

  samplesheet <- readr::read_tsv(samplesheet_path) %>%
    dplyr::select(sample_ID,
                  FACS_date,
                  animal_ID)

  seurat@meta.data %<>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(samplesheet,
                     by = c("sample_ID", "FACS_date")) %>%
    tibble::column_to_rownames("barcode")

  animalToPool <- seurat@meta.data %>%
    as.data.frame %>%
    dplyr::select(animal_ID,
                  .data[[grouping_column]],
                  seq_batch) %>%
    dplyr::distinct() %>%
    set_rownames(NULL) %>%
    dplyr::arrange(.data[[grouping_column]])

  nanimals <- animalToPool$animal_ID %>%
    lapply(function(x) {
      as.list(strsplit(x[[1]], ","))[[1]] %>% length()
    })

  animalToPool$n_animals <- nanimals

  ngroups <- animalToPool[[grouping_column]] %>%
    table %>%
    as.list()

  tmp <- animalToPool %>%
    dplyr::group_by(.data[[grouping_column]]) %>%
    dplyr::summarise(biol_rep = n())

  animalToPool %<>%
    left_join(tmp, by = grouping_column)

  animalToPool$pool_id <-
    paste0(animalToPool[[grouping_column]],
           "_",
           "pool",
           animalToPool$biol_rep,
           "_n",
           animalToPool$n_animals)

  seurat@meta.data %<>%
    tibble::rownames_to_column("barcode") %>%
    left_join(animalToPool,
              by = c("animal_ID",
                     grouping_column,
                     "seq_batch")) %>%
    tibble::column_to_rownames("barcode")

  return(seurat)
}


# The following function takes a Seurat object
# and returns a SCE object that has been
# log-normalized and filtered such that only
# genes with a minimum count across all cells
# are retained. Uses a lot of memory.
convertSeuratToSCE <- function(seurat,
                               minCountsToRetainGene = 10) {
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
                                    log = TRUE,
                                    pseudo.count = 1)

  if (minCountsToRetainGene > 0) {
    sce <- sce[rowSums(counts(sce) > 1) >= minCountsToRetainGene,]
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
                                             pseudobulk_grouping_column = "sample_ID") {
  require(Matrix.utils)

  # Named vector of cluster names
  cluster_ids <-
    purrr::set_names(levels(as.factor(sce$Major_celltype)))

  # Total number of clusters
  nc <- length(cluster_ids)

  # Named vector of sample names
  sample_ids <-
    purrr::set_names(levels(as.factor(sce[[pseudobulk_grouping_column]])))

  # Total number of samples
  ## 11
  ns <- length(sample_ids)

  # Generate sample level metadata


  ## Turn class "table" into a named vector of cells per sample
  n_cells <-
    table(sce[[pseudobulk_grouping_column]]) %>%  as.vector()
  names(n_cells) <- names(table(sce[[pseudobulk_grouping_column]]))

  ## Match the named vector with metadata to combine it
  m <- match(names(n_cells), sce[[pseudobulk_grouping_column]])

  ## Create the sample level metadata by selecting specific columns
  ei <- data.frame(colData(sce)[m,],
                   n_cells, row.names = NULL) %>%
    dplyr::select(.data[[pseudobulk_grouping_column]],  "n_cells")

  groups <-
    colData(sce)[, c("Major_celltype", pseudobulk_grouping_column)]
  groups[[pseudobulk_grouping_column]] <-
    factor(groups[[pseudobulk_grouping_column]])

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
splitAggregatedMatrixPerCelltype <- function(mat) {
  library(magrittr)

  pb <- mat
  splitf <- sapply(stringr::str_split(rownames(pb),
                                      pattern = "_",
                                      n = 2),
                   `[`, 1)

  pb <- split.data.frame(pb,
                         factor(splitf))

  pb <- pb %>%
    lapply(function(x) {
      sample_names <-
        rownames(x) %>% gsub(x = .,
                             pattern = "^.*_(FACS.*)$",
                             replacement = "\\1")
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
                                                    zscore = TRUE) {
  require(cmapR)

  map_df <- colData(sce) %>%
    as.data.frame %>%
    dplyr::group_by(.data[[pseudobulk_grouping_column]], .data[[grouping_column]]) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(sample_treatment = paste0(.data[[pseudobulk_grouping_column]], "_", .data[[grouping_column]]))

  for (celltype in names(mat_list)) {
    oldColnames <- colnames(mat_list[[celltype]])
    newColnames <- data.frame(sample_ID = oldColnames) %>%
      magrittr::set_colnames(pseudobulk_grouping_column) %>%
      dplyr::left_join(map_df, by = pseudobulk_grouping_column) %>%
      dplyr::pull(sample_treatment)

    colnames(mat_list[[celltype]]) <- newColnames
  }

  groupsOfInterest <- unique(map_df[[grouping_column]])
  for (celltype in names(mat_list)) {
    tmp_list <- list()

    for (group in groupsOfInterest) {
      tmp_list[[group]] <- mat_list[[celltype]] %>%
        as.data.frame %>%
        dplyr::select(contains(group))

      if (ncol(tmp_list[[group]]) == 0) {
        tmp_list[[group]] <- NULL
        next
      }

      if (zscore == TRUE) {
        capLargeZscores <- function(x) {
          x[x > 10] <- 10
          x[x < -10] <- -10
          return(x)
        }

        zscorefunc <- function(x) {
          x %>%
            cmapR::robust_zscore() %>%
            capLargeZscores()
        }

        tmp_list[[group]] <-  apply(tmp_list[[group]],
                                    2,
                                    zscorefunc)

        if (ncol(tmp_list[[group]]) > 1) {
          tmp_list[[group]] <- cmapR::distil(tmp_list[[group]],
                                             dimension = "col")$values
        }

      }
    }
    if (length(tmp_list) == 0L) {
      mat_list[[celltype]] <- NULL
    } else{
      mat_list[[celltype]] <- tmp_list
    }

  }
  return(mat_list)
}

capLargeZscores <- function(x) {
  x[x > 10] <- 10
  x[x < -10] <- -10
  return(x)
}

zscorefunc <- function(x) {
  x %>%
    cmapR::robust_zscore() %>%
    capLargeZscores()
}


# The above, but reworked for a DGEList object from a
# bulk RNA-seq analysis rather than single-cell.
DGEListToPerGroupZscores <- function(dge,
                                     grouping_column = "group",
                                     groups_to_keep) {
  if (!missing(groups_to_keep)) {
    samplesToKeep <- dge$samples %>%
      dplyr::filter(.data[[grouping_column]] %in% groups_to_keep) %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(.data[[grouping_column]], sample) %>%
      split(x = ., f = .$group)
  } else{
    # Use all groups present within the grouping_column
    samplesToKeep <- dge$samples %>%
      tibble::rownames_to_column("sample") %>%
      dplyr::select(.data[[grouping_column]], sample) %>%
      split(x = ., f = .$group)
  }

  mat_list <- list()
  for (group in names(samplesToKeep)) {
    mat_list[[group]] <- dge
    mat_list[[group]]$samples <-
      mat_list[[group]]$samples[samplesToKeep[[group]]$sample, ]
    mat_list[[group]]$counts <-
      mat_list[[group]]$counts[, samplesToKeep[[group]]$sample]

    tmp_log2counts <- mat_list[[group]] %>%
      edgeR::cpm(log = TRUE)

    # Convert each column (sample) from counts to z-score
    mat_list[[group]] <-  apply(tmp_log2counts,
                                2,
                                zscorefunc)

    if (ncol(mat_list[[group]]) > 1) {
      mat_list[[group]] <- cmapR::distil(mat_list[[group]],
                                         dimension = "col")$values
    }
  }
  return(mat_list)
}


# Absmax, returns the absolute maximum value
# while preserving the sign, important for zscores and
# other directional values.
absmax <- function(x) {
  x[which.max(abs(x))]
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
                                               functionForReplicates = "absmax",
                                               keepMouseGeneID = FALSE) {
  # Only download files necessary for the mapping df once
  if ("mm_ensembl_to_human_entrez.rds" %in% list.files(here("data"))) {
    mm_ensembl_to_human_entrez <-
      readRDS(here("data", "mm_ensembl_to_human_entrez.rds"))
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
      OrgDb = 'org.Hs.eg.db'
    ) %>%
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
                     multiple = "all")# %>%
  #dplyr::select(-mouse_gene_id)

  res2 <- plyr::ddply(res,
                      "human_entrezid",
                      plyr::numcolwise(functionForReplicates))

  if (keepMouseGeneID == TRUE) {
    res2 <- res2 %>%
      left_join(res, by = c("human_entrezid", "total"))
  }

  res2 <- res2 %>%
    dplyr::filter(!is.na(human_entrezid)) %>%
    tibble::column_to_rownames("human_entrezid")

  return(res2)
}


# Given a named list, rank the list in either
# ascending or descending order and return the
# top _n_ names (but not the values themselves)
top_n_names_from_list <- function(x,
                                  direction = "top",
                                  n = 150) {
  ordered_values <- order(-x)
  names_of_ordered_values <- names(x[ordered_values])
  if (direction == "top") {
    res <- head(names_of_ordered_values, n)
  } else if (direction == "bottom") {
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
                                     test = "hypergeometric") {
  table_split_by_groupingCol <- split(full_table,
                                      f = full_table[[groupingCol]])
  sig_drugs_per_groupingCol <- table_split_by_groupingCol %>%
    lapply(function(x) {
      x %>%
        dplyr::filter(sig == TRUE) %>%
        .$pert_iname %>% unique
    })
  sig_drugs_per_groupingCol2 <-
    data.frame(indication = names(sig_drugs_per_groupingCol))
  sig_drugs_per_groupingCol2$sig_drugs <- table_split_by_groupingCol

  # Prepare contingency matrices
  cont_tables <- table_split_by_groupingCol %>%
    lapply(function(x) {
      n_sig_with_indication = sum(x$sig)
      n_notsig_with_indication = nrow(x) - n_sig_with_indication
      n_sig_without_indication = sum(full_table$sig) - n_sig_with_indication
      n_notsig_without_indication = nrow(full_table) - n_notsig_with_indication

      res = matrix(
        c(
          n_sig_with_indication,
          n_notsig_with_indication,
          n_sig_without_indication,
          n_notsig_without_indication
        ),
        ncol = 2,
        byrow = TRUE
      )

      return(res)
    })

  # Calculate p-values for contingency matrices
  p_vals <- list()
  if (test == "fishers") {
    for (x in names(table_split_by_groupingCol)) {
      p_vals[[x]] <- fisher.test(cont_tables[[x]])$p.value
    }
  } else if (test == "hypergeometric") {
    for (x in names(table_split_by_groupingCol)) {
      p_vals[[x]] <- fisher.test(cont_tables[[x]],
                                 alternative = "less")$p.value
    }
  }

  res2 <- list(cont_matrices = cont_tables,
               p_vals = p_vals)

  # Add drug information to the final results - IN PROGRESS - IGNORE
  #mappingDf <- readRDS("/mnt/home/alkahest.com/nhin/Projects/cmap_internal/data/mappingTable2.rds")


  p_val_df <- data.frame(indication = names(p_vals),
                         p = unlist(p_vals)) %>%
    dplyr::mutate(BH_p = p.adjust(p, method = "BH")) %>%
    dplyr::mutate(sig = BH_p < 0.05) %>%
    dplyr::arrange(p) %>%
    dplyr::mutate(indication = tolower(indication)) %>%
    #dplyr::left_join(mappingDf, by = "indication") %>%
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
    currentRow <- listInputunique[i, ]
    myelements <-
      which(apply(listInputmat, 1, function(x)
        all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <-
      myelements
    myelements
    # attr(,"groups")
    #   one   two three
    # FALSE FALSE  TRUE
    #  f  i
    # 12 13
  }
  if (sort) {
    grouplist <-
      grouplist[order(sapply(grouplist, function(x)
        length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

# Convert human ensembl to human entrez
# The following function takes a data.frame structured as follows:
# rownames are mouse ensembl IDs, and another 1 column of values
# such as counts, z-scores, etc. (can be named anything).
# The mouse ensembl IDs are converted to human entrez IDs,
# and any many to 1 mappings are
# dealt with using the `functionForReplicates`, which can be
# any function, but recommended to be
# "sum", "mean", or "absmax". Any mouse ensembl IDs without
# corresponding human entrez IDs are filtered out.
convertHumanEnsemblToHumanEntrezID <- function(df,
                                               functionForReplicates = "absmax",
                                               keepHumanGeneID = FALSE) {
  # Only download files necessary for the mapping df once
  if ("hs_ensembl_to_human_entrez.rds" %in% list.files(here("data"))) {
    hs_ensembl_to_human_entrez <-
      readRDS(here("data", "hs_ensembl_to_human_entrez.rds"))
  } else{
    require(uatools)
    require(org.Hs.eg.db)
    require(clusterProfiler)
    require(here)

    hs_ensembl_to_human_entrez <-
      uatools::get_scaffold("ensembl_hs_gene_data") %>%
      dplyr::select(gene_id, hgnc_entrez_id) %>%
      dplyr::distinct(gene_id, hgnc_entrez_id, .keep_all = TRUE) %>%
      dplyr::rename(human_gene_id = gene_id,
                    human_entrezid = hgnc_entrez_id) %>%
      dplyr::filter(!is.na(human_entrezid))

    hs_ensembl_to_human_entrez %>%
      saveRDS(here("data", "hs_ensembl_to_human_entrez.rds"))
  }

  res <- df %>%
    tibble::rownames_to_column("human_gene_id") %>%
    dplyr::left_join(hs_ensembl_to_human_entrez,
                     multiple = "all")# %>%
  #dplyr::select(-human_gene_id)

  res2 <- plyr::ddply(res,
                      "human_entrezid",
                      plyr::numcolwise(functionForReplicates))

  if (keepHumanGeneID == TRUE) {
    res2 <- res2 %>%
      left_join(res, by = c("human_entrezid", "total"))
  }

  res2 <- res2 %>%
    dplyr::filter(!is.na(human_entrezid)) %>%
    tibble::column_to_rownames("human_entrezid")

  return(res2)
}

# In the below network, nodes are either cell lines,
# genes, and or drugs. The input is the CMap results.
createCMapNetwork <- function(cmap_result_table,
                              output_dir,
                              output_filename,
                              sigOnly = TRUE) {
  if (sigOnly == TRUE) {
    cmap_result_table %<>%
      dplyr::filter(sig == TRUE)
  }

  cmap_result_table %<>%
    dplyr::select(cell,
                  pert_iname,
                  t_gn_sym)   %>%
    tidyr::separate_rows(t_gn_sym, sep = ";") %>%
    dplyr::mutate(t_gn_sym = stringr::str_trim(t_gn_sym)) %>%
    dplyr::filter(!is.na(t_gn_sym),!is.na(pert_iname),!is.na(cell))  %>%
    dplyr::rename(Drug = pert_iname,
                  Gene = t_gn_sym,
                  Cell = cell)

  # Prepare edges tables
  cell2drug <- cmap_result_table %>%
    dplyr::select(Cell, Drug) %>%
    dplyr::filter(!is.na(Cell),!is.na(Drug)) %>%
    set_colnames(c("source", "target"))

  drug2gene <- cmap_result_table %>%
    dplyr::select(Drug, Gene) %>%
    dplyr::filter(!is.na(Drug),!is.na(Gene)) %>%
    set_colnames(c("source", "target"))

  edges <- dplyr::bind_rows(cell2drug,
                            drug2gene)

  nodes1 <- data.frame(id = cell2drug$source %>% unique,
                       type = "cell")
  nodes2 <- data.frame(id = cell2drug$target %>% unique,
                       type = "drug")
  nodes3 <- data.frame(id = drug2gene$source %>% unique,
                       type = "drug")
  nodes4 <- data.frame(id = drug2gene$target %>% unique,
                       type = "gene")

  nodes <- bind_rows(nodes1, nodes2, nodes3, nodes4) %>%
    dplyr::distinct(id, type, .keep_all = TRUE)

  res <- list(nodes = nodes, edges = edges)

  nodes %>%
    readr::write_csv(file.path(
      output_dir,
      paste0("CMap_network_", output_filename, "_nodes.csv")
    ))

  edges %>%
    readr::write_csv(file.path(
      output_dir,
      paste0("CMap_network_", output_filename, "_edges.csv")
    ))

  return(res)
}

# Create a network where nodes are drugs, genes, and/or diseases,
# and edges indicate relationships between them either from the
# CMap or DisGeNet analysis side.
createDrugGeneIndicationNetwork <- function(cmap_result_table,
                                            disgenet_result_table,
                                            cmap_cell_line = "HUVEC",
                                            disease_class_names = c("Nutritional and Metabolic Diseases",
                                                                    "Cardiovascular Diseases"),
                                            custom_class_mapping,
                                            output_dir,
                                            output_filename) {
  # First, get cmap_result_table into tidy form
  # so one gene per row.
  cmap_result_table_tidy <- cmap_result_table %>%
    dplyr::filter(cell == cmap_cell_line) %>%
    dplyr::filter(sig == TRUE) %>%
    dplyr::select(pert,
                  pert_iname,
                  starts_with("WTC"),
                  starts_with("NCS"),
                  MOAss,
                  t_gn_sym,
                  sig)   %>%
    tidyr::separate_rows(t_gn_sym, sep = ";") %>%
    dplyr::mutate(t_gn_sym =  stringr::str_trim(t_gn_sym)) %>%
    dplyr::select(pert_iname, t_gn_sym) %>%
    dplyr::filter(!is.na(t_gn_sym),!is.na(pert_iname))  %>%
    dplyr::rename(Drug = pert_iname,
                  Gene = t_gn_sym)

  # Assumed that cell line already accounted for
  # ie. DisGeNet results are from the gene symbols from
  # cmap_cell_line
  disgenet_result_table_tidy <- disgenet_result_table %>%
    dplyr::filter(sig == TRUE) %>%
    dplyr::filter(disease_class_name %in% disease_class_names) %>%
    dplyr::select(Description, shared_symbol) %>%
    tidyr::separate_rows(shared_symbol, sep = ";") %>%
    dplyr::mutate(shared_symbol =  stringr::str_trim(shared_symbol)) %>%
    dplyr::filter(!is.na(shared_symbol))  %>%
    dplyr::rename(Gene = shared_symbol,
                  Indication = Description)

  if (missing(custom_class_mapping)) {
    indication2class <- disgenet_result_table %>%
      dplyr::filter(sig == TRUE) %>%
      dplyr::filter(disease_class_name %in% disease_class_names) %>%
      dplyr::select(disease_class_name, Description) %>%
      dplyr::rename(Indication = Description,
                    Class = disease_class_name)
  } else {
    indication2class <- custom_class_mapping
  }


  joined_table <- full_join(
    cmap_result_table_tidy,
    disgenet_result_table_tidy,
    by = "Gene",
    multiple = "all"
  )

  drug2gene <- joined_table %>%
    dplyr::select(Drug, Gene) %>%
    dplyr::filter(!is.na(Drug),!is.na(Gene)) %>%
    set_colnames(c("source", "target"))

  gene2indication <- joined_table %>%
    dplyr::select(Gene, Indication) %>%
    dplyr::filter(!is.na(Gene),!is.na(Indication)) %>%
    set_colnames(c("source", "target"))

  edges <- dplyr::bind_rows(drug2gene,
                            gene2indication)

  nodes1 <- data.frame(id = drug2gene$source %>% unique,
                       type = "drug")
  nodes2 <- data.frame(id = drug2gene$target %>% unique,
                       type = "gene")
  nodes3 <- data.frame(id = gene2indication$source %>% unique,
                       type = "gene")
  nodes4 <- data.frame(id = gene2indication$target %>% unique,
                       type = "indication")

  nodes <- bind_rows(nodes1, nodes2, nodes3, nodes4) %>%
    dplyr::distinct(id, type, .keep_all = T) %>%
    left_join(indication2class,
              by = c("id" = "Indication"),
              multiple = "all")

  res <- list(nodes = nodes, edges = edges)

  nodes %>%
    readr::write_csv(file.path(output_dir, paste0(output_filename, "_nodes.csv")))

  edges %>%
    readr::write_csv(file.path(output_dir, paste0(output_filename, "_edges.csv")))

  return(res)

}

# Correlation scatterplots for all pairwise
# combinations of a data.frame
plotPairwiseCorrelationScatterplots <- function(df,
                                                cor_method = "spearman",
                                                gene_id_column = NULL,
                                                coord_x_text = -1,
                                                coord_y_text = 0.8,
                                                output_dir,
                                                output_name) {
  if(colnames(df)[[1]] %in% c("entrez",
                              "entrezid",
                              "entrez_id",
                              "ensembl_gene_id",
                              "ensembl",
                              "gene_id",
                              gene_id_column)){
    combinations_for_scatterplot <-
      colnames(df[-1]) %>%
      combn(m = 2) %>%
      t %>%
      as.data.frame
  } else{
    combinations_for_scatterplot <-
      colnames(df) %>%
      combn(m = 2) %>%
      t %>%
      as.data.frame
  }

  scatterplots_globalzsc <- list()
  for (combination in 1:nrow(combinations_for_scatterplot)) {
    comp1 <- combinations_for_scatterplot[combination, "V1"]
    comp2 <- combinations_for_scatterplot[combination, "V2"]

    cor_coef <- cor(df[[comp1]], df[[comp2]],
                    method = cor_method)

    scatterplots_globalzsc[[combination]] <-
      df %>%
      ggplot(aes(x = .data[[comp1]],
                 y = .data[[comp2]])) +
      geom_point(alpha = 0.2,
                 size = 1) +
      theme(aspect.ratio = 1,
            legend.position = "none") +
      geom_abline(intercept = 0,
                  slope = 1,
                  color = "red") +
      geom_hline(yintercept = 0, color = "black") +
      geom_vline(xintercept = 0, color = "black") +
      geom_smooth(method = "lm") +
      geom_text(
        x = coord_x_text,
        y = coord_y_text,
        label = paste0(cor_method, " ",
                       "coef. = ",
                       round(cor_coef, 2))
      )
  }

  names(scatterplots_globalzsc) <-
    apply(combinations_for_scatterplot,
          1,
          function(x) {
            x %>% paste0(collapse = "___")
          })

  for (combination in names(scatterplots_globalzsc)) {
    scatterplots_globalzsc[[combination]] %>%
      ggsave(
        filename = file.path(output_dir,
                             paste0(output_name, "_",
                                    combination,
                                    "_scatterplot.pdf")),
        plot = .,
        width = 8,
        height = 8,
        units = "in"
      )
  }

  return(scatterplots_globalzsc)
}
