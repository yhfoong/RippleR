#' Aggregate Gene Counts Across Concentric Radius Groups
#'
#' Calculates gene expression abundance across concentric radius groups
#' generated from tumor-centered spatial regions. For each selected gene,
#' the function:
#' \itemize{
#'   \item Computes total gene counts per radius group.
#'   \item Computes gene counts normalized by the total number of cells
#'         within each radius group.
#'   \item Generates and saves publication-ready plots for both summaries.
#' }
#'
#' Expression values are extracted directly from a Giotto object using either
#' raw transcript counts or normalized expression values.
#'
#' A lightweight Seurat object containing the subsetted expression matrix and
#' metadata is also generated and saved for downstream analyses.
#'
#' @param gobject A Giotto object containing spatial transcriptomics
#'   expression data.
#' @param meta.xenium Data frame containing cell-level metadata. Must contain
#'   at minimum:
#'   \itemize{
#'     \item \code{cell_ID}: unique cell identifiers matching Giotto cell IDs.
#'     \item \code{radius_group_label}: concentric radius group assignments.
#'   }
#' @param select.genes Character vector of gene names to summarize and plot.
#' @param count.type Character string specifying which expression matrix to
#'   use. One of:
#'   \itemize{
#'     \item{"Raw"} Raw transcript counts from the Giotto object.
#'     \item{"Normalized"} Normalized expression values from the Giotto object.
#'   }
#' @param plot1.color Color used for plots showing total gene counts per
#'   radius group.
#' @param plot2.color Color used for plots showing gene counts normalized by
#'   total cell number within each radius group.
#'
#'#' @importFrom forcats fct_reorder
#'
#' @return A list containing:
#' \describe{
#'   \item{seu}{A lightweight Seurat object containing the subsetted expression matrix and metadata.}
#' }
#'
#' @details
#' The function extracts either raw or normalized expression values from the
#' supplied Giotto object and restricts the expression matrix to cells present
#' in \code{meta.xenium}.
#'
#' For each gene in \code{select.genes}, two summary plots are generated:
#' \enumerate{
#'   \item Total gene counts per radius group.
#'   \item Gene counts normalized by the total number of cells within each
#'         radius group.
#' }
#'
#' Plots are automatically exported as PNG files to the current working
#' directory.
#'
#' The lightweight Seurat object is also saved as an RDS file for downstream
#' analyses.
#'
#' @seealso
#' \code{\link{generate_concentric_circles_giotto}},
#' \code{\link{annotate_radius_bins_giotto}},
#' \code{\link{label_peritumoral_distal_giotto}}
#'
#' @export

# Function to average raw/ normalized gene count per radius group
gene_count_radius_giotto <- function (gobject, meta.xenium, select.genes, count.type=c("Raw","Normalized"), plot1.color="darkblue", plot2.color="darkgreen"){
  
  # Safecheck for count.type (defaulted to Raw counts)
  count.type <- match.arg(count.type)
  
  message("Processing Giotto Object")
  
  # Get raw count table from Seurat Object
  if(count.type=="Raw"){
    raw_counts.ext <- Giotto::getExpression(gobject, values = "raw")
    raw_counts_mat <- raw_counts.ext[]
    
  } else {
    raw_counts.ext <- Giotto::getExpression(gobject, values = "normalized")
    raw_counts_mat <- raw_counts.ext[]
  }          
  
  # Check for duplicate
  message(
    "Any duplicated cells: ",
    any(duplicated(meta.xenium$cell_ID))
  )
  
  # Therefore, only get raw_counts cell_ID that are common to metadata
  # raw_counts.sub <- raw_counts[, colnames(raw_counts) %in% meta.xenium$cell_ID]
  raw_counts.sub <- raw_counts_mat[, colnames(raw_counts_mat) %in% meta.xenium$cell_ID]
  
  # Safecheck for select.genes list
  if (is.null(select.genes) || length(select.genes) == 0) {
    
    stop(
      "'select.genes' must contain at least one gene.",
      call. = FALSE
    )
    
  } else if (!all(select.genes %in% rownames(raw_counts_mat))) {
    
    missing.genes <- setdiff(
      select.genes,
      rownames(raw_counts_mat)
    )
    
    stop(
      paste0(
        "The following select.genes were not found in the Giotto dataset: ",
        paste(missing.genes, collapse = ", ")
      ),
      call. = FALSE
    )
    
  }
  
  
  # Can't pivot longer raw_counts.sub because it's a 477 × 916,041 dgCMatrix. That would create ~437 million rows (most zeros) and explode RAM.
  
  # Make sure metadata correspond to Xenium cell_ID
  common.cells <- intersect(
    colnames(raw_counts.sub),
    meta.xenium$cell_ID
  )
  
  raw_counts.sub <- raw_counts.sub[, common.cells]
  
  meta.xenium <- meta.xenium[
    match(common.cells, meta.xenium$cell_ID),
  ]
  
  # Plot raw gene counts by radius groups_____________________________________________
  # Create a lightweight Seurat object containing subsetted raw counts and metadata
  seu <- SeuratObject::CreateSeuratObject(counts = raw_counts.sub, meta.data = meta.xenium)
  
  # Save lightweight Seurat object as RDS
  message("Saving lightweight Seurat object with only gene counts and metadata")
  saveRDS(seu, file=paste0("Seurat_GiottoObj_",count.type,"Counts.RDS"))
  
  
  
  # Get count table from subsetted Seurat Object
  raw_counts.sf <- raw_counts.sub
  
  
  # Automate the plotting and exporting process for gene count by radius group
  message("Plotting gene count by radius group")
  
  plot.gene <- list()
  
  for (i in select.genes){
    # Process GOF one by one
    counts.gof <- tibble::tibble(
      cell_ID = colnames(raw_counts.sf),
      Gene.Count = as.numeric(raw_counts.sf[i, ]),
      Gene = i
    )
    
    # Merge counts.gof with meta.xenium to get radius group info
    counts.merge <- merge(counts.gof, meta.xenium, by="cell_ID")
    
    # Aggregate counts by radius group
    counts.agg <- counts.merge |> dplyr::group_by(radius_group_label) |> dplyr::summarise(Gene.sum = sum(Gene.Count)) |> dplyr::arrange(dplyr::desc(Gene.sum))
    
    
    # Factor relevel radius group by numbering
    counts.agg$radius.num <- stringr::str_extract(counts.agg$radius_group_label, "\\d+$") |> as.integer()
    counts.agg <- counts.agg |> dplyr::mutate(radius_group_label = forcats::fct_reorder(radius_group_label, radius.num, .fun = identity))
    
    # Plot Gene count by radius group and store in list
    plot.gene[[i]] <- ggplot2::ggplot(data = counts.agg , ggplot2::aes(x = radius_group_label, y = Gene.sum, group=1)) + ggplot2::geom_line(color=plot1.color, linewidth=1.5, show.legend=FALSE) + 
      ggplot2::geom_point(size = 3, shape = 21, fill = "white") + ggplot2::ylab(paste0(i, " ", count.type, " Gene Count")) + ggplot2::xlab("Radius Group") + ggplot2::theme_bw() + ggplot2::theme(axis.title=ggplot2::element_text(size=14, color="black"), axis.text.y=ggplot2::element_text(size=14, color="black"), axis.text.x = ggplot2::element_text(size=14, color="black", angle = 90, hjust = 1, vjust=1), axis.title.x= ggplot2::element_blank())
    
    ggplot2::ggsave(plot.gene[[i]], filename = paste0("Gene_count_by_radius_group_", count.type, "_", i, ".png"), width = 10,height = 8,units = "in",dpi = 300)       
    
  }
  
  
  
  # Plot raw gene counts normalized with total cell number per radius group_____________________________________________
  # Read in new subsetted Seurat RDS file
  # Get count table from subsetted Seurat Object
  raw_counts.sf <- raw_counts.sub
  
  # Automate the plotting and exporting process for raw gene count normalized with total cells per radius group
  message("Plotting gene count normalized with total cell number per radius group")
  
  plot.gene <- list()
  
  for (i in select.genes){
    # Process GOF one by one
    counts.gof <- tibble::tibble(
      cell_ID = colnames(raw_counts.sf),
      Gene.Count = as.numeric(raw_counts.sf[i, ]),
      Gene = i
    )
    
    # Merge counts.gof with meta.xenium to get radius group info
    counts.merge <- merge(counts.gof, meta.xenium, by="cell_ID")
    
    # Count total number of cells by radius group
    cell.count.group <- counts.merge |> dplyr::group_by(radius_group_label) |> dplyr::count() |> dplyr::rename(Cell.count.per.group=n)
    cell.count.group$radius.num <- stringr::str_extract(cell.count.group$radius_group_label, "\\d+$") |> as.integer()
    cell.count.group <- cell.count.group |> dplyr::arrange(radius.num)
    
    # Aggregate counts by radius group
    counts.agg <- counts.merge |> dplyr::group_by(radius_group_label) |> dplyr::summarise(Gene.sum = sum(Gene.Count)) |> dplyr::arrange(dplyr::desc(Gene.sum))
    
    # Merge counts.agg with cell count per radius group
    counts.agg.total <- merge(counts.agg, cell.count.group, by="radius_group_label")
    
    # Calculate raw counts normalized with total cell count per radius group
    counts.agg.total$Genes.sum.normalized.total.cells <- counts.agg.total$Gene.sum/ counts.agg.total$Cell.count.per.group
    
    
    # Factor relevel radius group by numbering
    counts.agg.total$radius.num <- stringr::str_extract(counts.agg.total$radius_group_label, "\\d+$") |> as.integer()
    counts.agg.total <- counts.agg.total |> dplyr::mutate(radius_group_label = forcats::fct_reorder(radius_group_label, radius.num, .fun = identity))
    
    
    # Plot Gene count by radius group and store in list
    plot.gene[[i]] <- ggplot2::ggplot(data = counts.agg.total, ggplot2::aes(x = radius_group_label, y = Genes.sum.normalized.total.cells, group=1)) + ggplot2::geom_line(color=plot2.color, linewidth=1.5, show.legend=FALSE) + 
      ggplot2::geom_point(size = 3, shape = 21, fill = "white") + ggplot2::ylab(paste0(i, " ", count.type, " Gene Count/ Total Cells Per Radius Group")) + ggplot2::xlab("Radius Group") + ggplot2::theme_bw() + ggplot2::theme(axis.title=ggplot2::element_text(size=14, color="black"), axis.text.y=ggplot2::element_text(size=14, color="black"), axis.text.x = ggplot2::element_text(size=14, color="black", angle = 90, hjust = 1, vjust=1), axis.title.x= ggplot2::element_blank())
    
    ggplot2::ggsave(plot.gene[[i]], filename = paste0("Gene_count_by_radius_group[normTotalCellsPerGroup]_", count.type, "_", i, ".png"), width = 10,height = 8,units = "in",dpi = 300)       
    
  }
  
  result <- list(seu=seu)
  
  return(invisible(result))
}


