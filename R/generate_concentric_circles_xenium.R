#' Generate Concentric Circles Around a Tumor Region
#'
#' @param xenium.obj A Seurat Xenium object containing cell metadata and
#'   spatial coordinates.
#' @param sample Character string specifying the sample name to analyze.
#'   Must match values stored in `orig.ident`.
#' @param L.annot Character string specifying the metadata column containing
#'   cell annotations.
#' @param Tumor Character string specifying the annotation label identifying
#'   tumor cells.
#' @param conc.radius.start Numeric. Starting radius for the first circle
#'   (default = 200).
#' @param conc.radius.stepSize Numeric. Increment between successive circle
#'   radii (default = 200).
#' @param conc.direction Character string specifying the direction of circle
#'   expansion. One of:
#'   \itemize{
#'     \item{"full"} Full concentric circles.
#'     \item{"right"} Right-facing semicircles.
#'     \item{"left"} Left-facing semicircles.
#'     \item{"top"} Upward-facing semicircles.
#'     \item{"bottom"} Downward-facing semicircles.
#'   }
#' @param save_output Logical. If TRUE, save the results as an RDS file.
#' @param output_file Character string specifying the output RDS filename.
#' 
#' @importFrom rlang .data
#' 
#' @return A list containing:
#' \describe{
#'   \item{xenium.sub}{Metadata subset corresponding to the selected sample.}
#'   \item{tumor_centroid_coords}{Data frame containing tumor centroid coordinates.}
#'   \item{center_point_sfc}{Tumor centroid represented as an sf point object.}
#'   \item{bbox}{Bounding box coordinates of the selected sample.}
#'   \item{max_radius}{Maximum allowable radius constrained by the bounding box.}
#'   \item{radii}{Vector of generated circle radii.}
#'   \item{clip_rect}{Clipping polygon used for directional semicircles. NULL when `conc.direction = "full"`.}
#'   \item{concentric_circles_sf}{sf object containing concentric circles or semicircles.}
#' }
#'
#' @details
#' Spatial coordinates are expected to be stored in metadata columns
#' named `x.coord` and `y.coord`. The tumor centroid is computed as the
#' mean x- and y-coordinate of all cells matching the specified tumor
#' annotation.
#'
#' @export


# Function to generate concentric circles
generate_concentric_circles_xenium <- function (xenium.obj, sample, L.annot="L3.annot", Tumor = "Tumor.1", conc.radius.start = 200, conc.radius.stepSize =200, conc.direction = c("right", "left", "top", "bottom", "full"), save_output = TRUE, output_file = "concentric.circle.outputs.rds")
{
  # SAFE-PROOF: Automatically grabs the first direction string (or defaults to "right")
  conc.direction <- match.arg(conc.direction)
  
  # Extract metadata from gobject
  md <- xenium.obj[[]]
  
  # Ensure the target vector is a clean character string
  md$orig.ident <- as.character(md$orig.ident)
  
  
  # Filter for specific sample
  xenium.sub <- md[md$orig.ident == sample,]
  
  # Get tumor coords
  xenium.tumor <- xenium.sub[xenium.sub[[L.annot]] %in% Tumor,]
  
  
  # Right-side concentric circles around Tumor center:
  # Calculate tumor centroid coordinates
  tumor_centroid_coords <- xenium.tumor |>
    dplyr::summarise(x = mean(x.coord), y = mean(y.coord))
  
  # Convert centroid to sf POINT with empty CRS (planar)
  center_point_sfc <- sf::st_sfc(sf::st_point(c(tumor_centroid_coords$x, tumor_centroid_coords$y)), crs = sf::NA_crs_)
  
  # Get bounding box of the dataset
  bbox <- xenium.sub |>
    dplyr::summarise(min_x = min(x.coord, na.rm = TRUE), max_x = max(x.coord, na.rm = TRUE),
              min_y = min(y.coord, na.rm = TRUE), max_y = max(y.coord, na.rm = TRUE))
  
  # Calculate max allowed radius based on proximity to right, top, and bottom edges
  max_radius <- switch(
    conc.direction,
    
    "right" = min(
      bbox$max_x - tumor_centroid_coords$x,
      tumor_centroid_coords$y - bbox$min_y,
      bbox$max_y - tumor_centroid_coords$y
    ),
    
    "left" = min(
      tumor_centroid_coords$x - bbox$min_x,
      tumor_centroid_coords$y - bbox$min_y,
      bbox$max_y - tumor_centroid_coords$y
    ),
    
    "top" = min(
      bbox$max_y - tumor_centroid_coords$y,
      tumor_centroid_coords$x - bbox$min_x,
      bbox$max_x - tumor_centroid_coords$x
    ),
    
    "bottom" = min(
      tumor_centroid_coords$y - bbox$min_y,
      tumor_centroid_coords$x - bbox$min_x,
      bbox$max_x - tumor_centroid_coords$x
    ),
    
    "full" = min(
      tumor_centroid_coords$x - bbox$min_x,
      bbox$max_x - tumor_centroid_coords$x,
      tumor_centroid_coords$y - bbox$min_y,
      bbox$max_y - tumor_centroid_coords$y
    ),
    
    
    stop("Conc.direction must be one of: right, left, top, bottom, full")
  )
  
  # Define radii within bounding box (every 200 units)
  radii <- seq(conc.radius.start, max_radius, by = conc.radius.stepSize)
  
  # Define clipping rectangle:
  clip_rect <- if (conc.direction == "full") {
    
    NULL
    
  } else 
    
  {switch(
    
    conc.direction,
    
    "right" = {
      sf::st_polygon(list(rbind(
        c(tumor_centroid_coords$x, bbox$min_y),
        c(bbox$max_x, bbox$min_y),
        c(bbox$max_x, bbox$max_y),
        c(tumor_centroid_coords$x, bbox$max_y),
        c(tumor_centroid_coords$x, bbox$min_y)
      )))
    },
    
    "left" = {
      sf::st_polygon(list(rbind(
        c(bbox$min_x, bbox$min_y),
        c(tumor_centroid_coords$x, bbox$min_y),
        c(tumor_centroid_coords$x, bbox$max_y),
        c(bbox$min_x, bbox$max_y),
        c(bbox$min_x, bbox$min_y)
      )))
    },
    
    "top" = {
      sf::st_polygon(list(rbind(
        c(bbox$min_x, tumor_centroid_coords$y),
        c(bbox$max_x, tumor_centroid_coords$y),
        c(bbox$max_x, bbox$max_y),
        c(bbox$min_x, bbox$max_y),
        c(bbox$min_x, tumor_centroid_coords$y)
      )))
    },
    
    "bottom" = {
      sf::st_polygon(list(rbind(
        c(bbox$min_x, bbox$min_y),
        c(bbox$max_x, bbox$min_y),
        c(bbox$max_x, tumor_centroid_coords$y),
        c(bbox$min_x, tumor_centroid_coords$y),
        c(bbox$min_x, bbox$min_y)
      ))) }
    
  ) |>
      sf::st_sfc(crs = sf::NA_crs_)
    
  }
  
  
  # Create semicircle buffers clipped by rectangle
  circle_list <- lapply(radii, function(r) {
    
    circle <- sf::st_buffer(
      center_point_sfc,
      dist = r
    )
    
    if (conc.direction == "full") {
      circle
    } else {
      sf::st_intersection(
        circle,
        clip_rect
      )
    }
    
  })
  
  concentric_circles_sf <- do.call(
    c,
    circle_list
  ) |>
    sf::st_sf(radius = radii)
  
  
  # Save outputs as list
  result <- list(
    xenium.sub=xenium.sub, 
    tumor_centroid_coords=tumor_centroid_coords, 
    center_point_sfc=center_point_sfc,
    bbox=bbox, 
    max_radius=max_radius, 
    radii=radii, 
    clip_rect=clip_rect,
    concentric_circles_sf= concentric_circles_sf
    )
  

  # Save outputs as rds
  if (isTRUE(save_output)) {
    
    saveRDS(
      result,
      file = output_file
    )
    
    message(
      "Saved output as RDS: ",
      normalizePath(output_file)
    )
  }
  
  
  return(invisible(result))
  
}



