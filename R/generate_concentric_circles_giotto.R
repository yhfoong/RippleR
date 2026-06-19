#' Generate Concentric Circles Around a Tumor Region
#'
#' Creates concentric circles or directional semicircles centered on the
#' centroid of a specified tumor population.
#'
#' @param gobject Giotto object.
#' @param sample Sample name.
#' @param x_shift Logical; whether to apply Giotto x-coordinate shift (only needed if joining multiple Giotto objects).
#' @param L.annot Metadata annotation column.
#' @param Tumor Tumor annotation label.
#' @param conc.radius.start Starting radius.
#' @param conc.radius.stepSize Radius increment.
#' @param conc.direction Direction of circle expansion.
#' @param save_output Logical; whether to save output as an RDS file.
#' @param output_file Output RDS filename.
#' 
#' @importFrom rlang .data
#' 
#' @return A list containing:
#' \describe{
#'   \item{xenium.sub}{Filtered metadata}
#'   \item{tumor_centroid_coords}{Tumor centroid coordinates}
#'   \item{center_point_sfc}{Tumor centroid sf point}
#'   \item{bbox}{Bounding box}
#'   \item{max_radius}{Maximum allowable radius}
#'   \item{radii}{Radius vector}
#'   \item{clip_rect}{Clipping polygon}
#'   \item{concentric_circles_sf}{Concentric circle sf object}
#' }
#'
#' @export


# Function to generate concentric circles
generate_concentric_circles_giotto <- function (gobject = gobject, sample ="SLN_stageIII", x_shift=TRUE, L.annot="L3.annot", Tumor = "Tumor.1", conc.radius.start = 200, conc.radius.stepSize =200, conc.direction = c("right", "left", "top", "bottom", "full"), save_output = TRUE, output_file = "concentric.circle.outputs.rds")
{
  # SAFE-PROOF: Automatically grabs the first direction string (or defaults to "right")
  conc.direction <- match.arg(conc.direction)
  
  # Extract metadata from gobject
  md <- as.data.frame(Giotto::pDataDT(gobject))
  
  # Ensure the target vector is a clean character string
  md$orig.ident <- as.character(md$orig.ident)
  
  # Note that when generating multiple Giotto object, we applied x_shift of 27000; y.coord is the same 
  if (isTRUE(x_shift)){
    
    md$x.coord.shift <- md$x.coord + 27000
  } else {md$x.coord.shift <- md$x.coord}
  
  
  # Filter for specific sample
  xenium.sub <- md[md$orig.ident == sample,]

  # Get tumor coords
  xenium.tumor <- xenium.sub[xenium.sub[[L.annot]] %in% Tumor,]

  
  # Right-side concentric circles around Tumor center:
  # Calculate tumor centroid coordinates
  tumor_centroid_coords <- xenium.tumor |>
    dplyr::summarise(x = mean(x.coord.shift), y = mean(y.coord))

  # Convert centroid to sf POINT with empty CRS (planar)
  center_point_sfc <- sf::st_sfc(sf::st_point(c(tumor_centroid_coords$x, tumor_centroid_coords$y)), crs = sf::NA_crs_)

  # Get bounding box of the dataset
  bbox <- xenium.sub |>
    dplyr::summarise(min_x = min(x.coord.shift), max_x = max(x.coord.shift),
              min_y = min(y.coord), max_y = max(y.coord))

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
    xenium.sub = xenium.sub,
    tumor_centroid_coords = tumor_centroid_coords,
    center_point_sfc = center_point_sfc,
    bbox = bbox,
    max_radius = max_radius,
    radii = radii,
    clip_rect = clip_rect,
    concentric_circles_sf = concentric_circles_sf
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



