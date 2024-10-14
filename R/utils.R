require('Seurat')

# pretty print time difference
time_diff <- function(start_time, end_time = NULL) {
  if (is.null(end_time)) {
    end_time <- proc.time()['elapsed']
  }
  tmp <- capture.output(lubridate::make_difftime(end_time - start_time))
  gsub(pattern = 'Time difference of ', replacement = '', x = tmp)
}

# given a directory, return the transcriptome count matrices
get_counts <- function(path) {
  h5_files <- list.files(path = path, pattern = '^filtered_feature_bc_matrix\\.h5$', recursive = TRUE, full.names = TRUE)
  ret <- list()
  for (f in h5_files) {
    counts <- Read10X_h5(f)
    # in case of HTO multiplexing, separate samples here
    if (inherits(x = counts, what = 'list') && length(counts) == 2 && names(counts)[2] == 'Custom') {
      ret <- append(ret, get_demuxed_counts(h5_file_path = f, counts))
    }
    if (inherits(x = counts, what = 'dgCMatrix')) {
      ret[[basename(dirname(f))]] <- counts
    }
  }
  return(ret)
}

# For filtering of raw count matrix
filter_count_matrix <- function(counts, 
                                percent.mito.th = 15,
                                z.th.counts = c(-2.5, 2.5),
                                z.th.feature_outlier = c(-3, 3),
                                min.features = 300,
                                mito.pattern = "^MT[^0-9a-zA-Z]+",
                                return.seurat = FALSE,
                                sample_id = NULL,
                                verbose = TRUE) {
  
  s <- CreateSeuratObject(counts = counts)
  
  s[['orig.ident']] <- sample_id
  s[['percent.mito']] <- PercentageFeatureSet(s, pattern = mito.pattern)
  # at this point we can do some initial filtering based on gene expression metrics
  
  # There are four steps
  keep <- 1:ncol(s)
  
  # 1) percent mitochondrial reads
  keep_this <- s$percent.mito <= percent.mito.th
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p1 <- ggplot(s@meta.data, aes(x = percent.mito)) + 
    geom_histogram(binwidth = 1) + 
    geom_vline(xintercept = c(-Inf, percent.mito.th), color = 'red') + 
    xlab('Mitochondrial reads in %') + 
    annotate('label', x = Inf, y = Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2:::.pt, vjust = 1, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]
  
  # 2) number of features
  nFeature_RNA <- s$nFeature_RNA[keep]
  nFeature_RNA_logscaled <- scale(log10(nFeature_RNA))
  nFeature_RNA_logscaled[nFeature_RNA >= min.features] <- scale(log10(nFeature_RNA[nFeature_RNA >= min.features]))
  keep_this <- nFeature_RNA >= min.features & nFeature_RNA_logscaled >= z.th.counts[1] & nFeature_RNA_logscaled <= z.th.counts[2]
  th <- c(min(nFeature_RNA[keep_this]) - .Machine$double.eps*10, 
          max(nFeature_RNA[keep_this]) + .Machine$double.eps*10)
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p2 <- ggplot(s@meta.data[keep, ], aes(x = nFeature_RNA)) +
    geom_histogram(binwidth = 33) + 
    geom_vline(xintercept = th, color = 'red') + 
    xlab('Number of genes detected') + 
    annotate('label', x = Inf, y = Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2:::.pt, vjust = 1, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]
  
  # 3) number of molecules
  nCount_RNA <- s$nCount_RNA[keep]
  nCount_RNA_logscaled <- scale(log10(nCount_RNA))
  keep_this <- nCount_RNA_logscaled >= z.th.counts[1] & nCount_RNA_logscaled <= z.th.counts[2]
  th <- c(min(nCount_RNA[keep_this]) - .Machine$double.eps*10, 
          max(nCount_RNA[keep_this]) + .Machine$double.eps*10)
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p3 <- ggplot(s@meta.data[keep, ], aes(x = nCount_RNA)) +
    geom_histogram(binwidth = 33) + 
    geom_vline(xintercept = th, color = 'red') + 
    xlab('Number of transcripts') + 
    annotate('label', x = Inf, y = Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2:::.pt, vjust = 1, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]
  
  # 4) transcripts vs genes
  md <- s@meta.data[keep, ]
  mod <- loess(log10(nFeature_RNA) ~ log10(nCount_RNA), data = md, span = 1)
  md$nFeature_outlier <- scale(mod$residuals) < z.th.feature_outlier[1] | scale(mod$residuals) > z.th.feature_outlier[2]
  keep_this <- !md$nFeature_outlier
  txt <- sprintf("Remove %d of %d\n(%1.1f%%) cells", sum(!keep_this), length(keep), sum(!keep_this) / length(keep) * 100)
  p4 <- ggplot(md, aes(log10(nCount_RNA), log10(nFeature_RNA), color = nFeature_outlier)) + 
    #coord_trans(x = 'log10', y = 'log10') +
    geom_point(size = 0.5) + scale_color_manual(values = c('grey35', 'red'), guide = 'none') +
    xlab('Transcripts') + ylab('Genes') +
    scale_y_continuous(breaks = log10((1:55)*1000), labels = (1:55)*1000) +
    scale_x_continuous(breaks = log10(10^(0:7)), labels = as.character(10^(0:7))) +
    theme(panel.grid.minor = element_blank()) +
    annotation_logticks() +
    annotate('label', x = Inf, y = -Inf, label = txt, size = theme_get()$text$size * 0.75 / ggplot2:::.pt, vjust = 0, hjust = 1, label.r = unit(0, 'cm'))
  keep <- keep[keep_this]
  
  keep <- colnames(s)[keep]
  txt <- sprintf('Sample %s; Keeping %d of %d cells (%1.2f percent)\n', sample_id, length(keep), ncol(s), length(keep)/ncol(s)*100)
  if (verbose) {
    message(txt)
  }
  
  #p <- plot_grid(p1, p2, p3, p4)
  #p <- add_title(plot_obj = p, title_str = txt, rel_height = 0.08, font_face = 'bold')
  
  s <- s[, keep]
  if (!return.seurat) {
    s <- s$RNA@counts
  }
  return(list(filtered = s, figures = list(p1, p2, p3, p4), fig_title = txt))
}

composition_plots <- function(df, group_var, comp_var, group_name, comp_name,
                              freq_limit = 0, max_guide_lines = 9) {
  df_sum <- data.frame(group = df[, group_var],
                       label = df[, comp_var]) %>%
    group_by(group) %>% 
    mutate(group_size = n()) %>%
    group_by(group, label) %>% 
    summarise(n = n(), freq = n / group_size[1])
  other_labels <- group_by(df_sum, label) %>%
    summarise(min_freq = max(freq)) %>%
    filter(min_freq < freq_limit) %>%
    pull(label)
  df_sum$label_filtered <- as.character(df_sum$label)
  df_sum$label_filtered[df_sum$label %in% other_labels] <- 'other'
  df_sum$label_filtered <- factor(df_sum$label_filtered, levels = c(setdiff(df_sum$label_filtered, 'other'), 'other'))
  
  # manually set colors here to be able to assign a different color to 'other'
  cols <- colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(length(levels(df_sum$label_filtered))-1)
  cols <- c(cols, 'gray90')
  
  p1 <- ggplot(df_sum, aes(group, n, fill = label_filtered)) + 
    geom_bar(stat = 'identity') + 
    ylab('Cell count') + xlab(group_name) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p2 <- ggplot(df_sum, aes(group, freq, fill = label_filtered)) + 
    geom_bar(stat = 'identity') + 
    ylab('Frequency') + xlab(group_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ncol_guide <- ceiling(length(unique(df_sum$label_filtered)) / max_guide_lines)
  ncol_guide <- scales::oob_censor_any(ncol_guide, c(1, 10))
  p <- (p1 + p2) + plot_layout(guides = "collect") & 
    scale_fill_manual(name = comp_name, values = cols) &
    guides(fill=guide_legend(ncol=ncol_guide))
  colnames(df_sum)[1:2] <- c(group_name, comp_name)
  return(list(data = df_sum, figure = p, fig1 = p1, fig2 = p2))
}

# perform azimuth annotation using a reference
# query should be a seurat object 
azimuth_annotation <- function(query, reference_path, reference_column, 
                               query_assay = 'RNA',
                               verbose = FALSE) {
  # Ensure all packages needed are installed and loaded
  
  # Ensure Seurat v4.0 or higher is installed
  if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
    stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
  }
  
  # Ensure glmGamPoi is installed
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      BiocManager::install("glmGamPoi")
    }
  }
  
  # Ensure Azimuth is installed
  if (packageVersion(pkg = "Azimuth") < package_version(x = "0.3.1")) {
    stop("Please install azimuth - remotes::install_github('satijalab/azimuth')", call. = FALSE)
  }
  
  library('Seurat')
  library('Azimuth')
  
  # Load the reference
  reference <- LoadReference(path = reference_path)
  DefaultAssay(query) <- query_assay
  query <- DietSeurat(query, assays = query_assay)
  
  # Preprocess with SCTransform (exactly as in Azimuth web app)
  query <- SCTransform(
    object = query,
    assay = query_assay,
    new.assay.name = "refAssay",
    residual.features = rownames(x = reference$map),
    reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
    method = 'glmGamPoi',
    ncells = 2000,
    n_genes = 2000,
    do.correct.umi = FALSE,
    do.scale = FALSE,
    do.center = TRUE, 
    verbose = verbose
  )
  
  # Find anchors between query and reference
  anchors <- FindTransferAnchors(
    reference = reference$map,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
    dims = 1:50,
    n.trees = 20,
    mapping.score.k = 100, 
    verbose = verbose
  )
  
  # Transferred labels are in metadata columns named "predicted.*"
  # The maximum prediction score is in a metadata column named "predicted.*.score"
  # The prediction scores for each class are in an assay named "prediction.score.*"
  refdata <- lapply(X = reference_column, function(x) {
    reference$map[[x, drop = TRUE]]
  })
  names(x = refdata) <- reference_column
  query <- TransferData(
    reference = reference$map,
    query = query,
    dims = 1:50,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE, 
    verbose = verbose
  )
  
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference$map,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE, 
    verbose = verbose
  )
  
  # Calculate the query neighbors in the reference
  # with respect to the integrated embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference$map[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE,
    verbose = verbose
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = reference$map[[]]
  )
  
  # Project the query to the reference UMAP.
  query[["proj.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = reference$map[["refUMAP"]],
    reduction.key = 'UMAP_',
    verbose = verbose
  )
  
  # Calculate mapping score and add to metadata
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  return_columns <- gsub(pattern = 'xxx', replacement = reference_column, 
                         x = c('predicted.xxx.score', 'predicted.xxx', 'mapping.score'))
  ret <- query@meta.data[, return_columns]
  colnames(ret) <- c('az.annotation.score', 'az.annotation', 'az.mapping.score')
  
  return(list(prediction = ret,
              umap = query[['proj.umap']]))
}

azimuth_plots <- function(s, reduction, annotation, annotation_score, mapping_score) {
  # DimPlot of the query, colored by predicted cell type
  p2 <- DimPlot(object = s, reduction = reduction, group.by = annotation, label = TRUE, repel = TRUE) + NoLegend() + ggtitle('Query projected UMAP')
  
  # Plot the score for the predicted cell type of the query
  p3 <- FeaturePlot(object = s, features = annotation_score, reduction = reduction)
  
  p4 <- VlnPlot(object = s, features = annotation_score, group.by = annotation) + NoLegend()
  
  # Plot the mapping score
  p5 <- FeaturePlot(object = s, features = mapping_score, reduction = reduction)
  
  p6 <- VlnPlot(object = s, features = mapping_score, group.by = annotation) + NoLegend()
  
  return(list(p2, p3, p4, p5, p6))
}

pseudobulk <- function(counts, grouping, transformation = NULL) {
  grouping <- droplevels(as.factor(grouping))
  mat <- sapply(levels(grouping), function(gr) {
    rowSums(counts[, grouping == gr, drop = FALSE])
  })
  colnames(mat) <- levels(grouping)
  if (!is.null(transformation)) {
    if (transformation == 'log_rel') {
      sf_target <- median(colSums(mat))
      message('pseudobulk log_rel transformation target counts ', sf_target)
      mat <- apply(mat, 2, function(x) {
        log(x / sum(x) * sf_target + 1)
      })
    }
    if (transformation == 'logCPM') {
      message('pseudobulk transformation: logCPM')
      mat <- apply(mat, 2, function(x) {
        log(x / sum(x) * 1e6 + 1)
      })
    }
    if (transformation == 'logCP10K') {
      message('pseudobulk transformation: logCP10K')
      mat <- apply(mat, 2, function(x) {
        log(x / sum(x) * 1e4 + 1)
      })
    }
    if (transformation == 'DESeq2counts') {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = data.frame(grp = colnames(mat)), design = ~ grp)
      mat <- DESeq2::counts(dds, normalized=TRUE)
    }
    if (transformation == 'DESeq2logcounts') {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = data.frame(grp = colnames(mat)), design = ~ grp)
      dds <- DESeq2::estimateSizeFactors(dds)
      mat <- log1p(DESeq2::counts(dds, normalized=TRUE))
    }
    if (transformation == 'DESeq2vst') {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = data.frame(grp = colnames(mat)), design = ~ grp)
      mat <- SummarizedExperiment::assay(DESeq2::vst(dds))
    }
  }
  return(mat)
}


run_edgeR <- function(counts, group_labels, sample_labels, pval_clip = 1e-20, verbose = FALSE) {
  # create pseudobulk data
  group_labels <- droplevels(as.factor(group_labels))
  sample_labels <- droplevels(as.factor(sample_labels))
  
  if (verbose) {
    message('edgeR')
    print(table(group_labels, sample_labels))
  }
  
  combinations <- expand.grid(levels(group_labels), levels(sample_labels)) %>%
    mutate(grouping_level = paste(Var1, Var2, sep = '.'))
  grouping <- interaction(group_labels, sample_labels, drop = TRUE)
  pb_counts <- pseudobulk(counts = counts, 
                          grouping = grouping)
  pb_grouping <- combinations$Var1[match(x = colnames(pb_counts), table = combinations$grouping_level)]
  pb_grouping <- droplevels(factor(x = pb_grouping, levels = levels(group_labels)))
  
  # run edgeR
  y <- edgeR::DGEList(counts = pb_counts, group = pb_grouping)
  keep <- edgeR::filterByExpr(y)
  if (verbose) {
    message(sprintf('edgeR: %d genes expressed above threshold', sum(keep)))
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~ pb_grouping)
  y <- edgeR::estimateDisp(y, design)
  et <- edgeR::exactTest(y)
  res <- et$table %>% tibble::rownames_to_column(var = 'gene') %>%
    mutate(FDR = p.adjust(PValue, method = 'fdr')) %>%
    mutate(pval_clipped = scales::oob_squish_any(x = PValue, range = c(pval_clip, Inf))) %>%
    arrange(FDR) 
  return(res)
}

# run edgeR on multiple splits of the data
run_edgeR_split <- function(counts, group_labels, sample_labels, split_factor, top_n = 6, ...) {
  split_factor <- droplevels(as.factor(split_factor))
  res_lst <- list()
  for (split_level in levels(split_factor)) {
    sel <- which(split_factor == split_level)
    res <- run_edgeR(counts[, sel, drop = FALSE], group_labels[sel], sample_labels[sel], ...) %>%
      mutate(split_level = split_level, .before = gene) 
    res_lst[[split_level]] <- res
  }
  res <- do.call(rbind, res_lst)
  res$split_level <- factor(res$split_level, levels = levels(split_factor))
  
  top_markers <- filter(res, FDR < 0.05) %>%
    group_by(split_level, sign(logFC)) %>% 
    filter(rank(FDR, ties.method = "first") <= top_n) 
  
  return(list(res = res, top_markers = top_markers))
}

# normalized entropy from vector of observations; NAs are dropped
norm_entropy <- function(observations) {
  observations <- observations[!is.na(observations)]
  p <- as.numeric(table(observations)) / length(observations)
  n <- length(p)
  if (n == 1) {
    return(0)
  }
  h <- -sum(p*log(p))
  #h_max <- log(n)
  h_max <- log(length(observations))
  return(h / h_max)
}

expression_summary <- function(counts, grouping, transformation, goi) {
  pb_mat <- pseudobulk(counts = counts, 
                       grouping = grouping, 
                       transformation = transformation)
  det_rate <- sctransform:::row_mean_grouped_dgcmatrix(matrix = counts > 0, group = grouping, shuffle = FALSE)
  
  df <- melt(t(pb_mat[goi, ]), 
             varnames = c('grp_id', 'gene'), 
             value.name = 'expression') %>% 
    left_join(melt(t(det_rate[goi, ]),
                   varnames = c('grp_id', 'gene'), 
                   value.name = 'det_rate'))
  return(df)
}
