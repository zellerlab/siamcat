#' @title Plot ordination of a SIAMCAT object
#'
#' @description A thin wrapper around \code{\link[phyloseq]{plot_ordination}}
#' and \code{\link[phyloseq]{ordinate}} with sensible defaults and boilerplate
#' code for SIAMCAT objects. Computes and plots an ordination of the
#' phyloseq object contained within a \link{siamcat-class} object.
#'
#' @usage plot.ordination(siamcat,
#' color.by = NULL, name.color.by = NULL, palette = NULL,
#' font.size = 14, fn.plot = NULL)
#'
#' @param siamcat object of class \link{siamcat-class}
#'
#' @param color.by string, name of a column in the sample metadata to use
#' for coloring points. If \code{NULL} (default), points are not colored.
#' Numeric columns are colored with a continuous color scale; categorical
#' columns are colored with a discrete palette.
#'
#' @param name.color.by string, label for the color legend. If \code{NULL}
#' (default), the value of \code{color.by} is used.
#' 
#' @param shape.by string, name of a column in the sample metadata to use
#' for setting the shape of points. If \code{NULL} (default), points are all of the same shape.
#' Must be a categorical variable.
#'
#' @param name.shape.by string, label for the shape legend. If \code{NULL}
#' (default), the value of \code{shape.by} is used.
#'
#' @param name.color.by string, label for the color legend. If \code{NULL}
#' (default), the value of \code{color.by} is used.
#' 
#' @param rename.values named vector, used to rename the values of the color.by variable.
#' The names of the vector are the original values, and the values of the vector are the new names.
#' If \code{NULL} (default), no renaming is performed.
#' The order of the names in the vector determines the order of the labels in the legend.
#' 
#' @param alpha numeric, transparency of the points
#' 
#' @param palette for continuous \code{color.by}: a ColorBrewer palette name,
#' defaults to \code{"RdBu"}. For categorical \code{color.by}: a vector of
#' valid R colors, defaults to \code{okabe_palette}. If \code{NULL},
#' the appropriate default is used.
#'
#' @param font.size integer, base font size for the plot, defaults to
#' \code{14}
#'
#' @param fn.plot string, filename for the plot (any extension supported by
#' \code{\link[ggplot2]{ggsave}} is allowed). If \code{NULL} (default),
#' the plot is only returned as a ggplot object and not saved to disk.
#'
#' @param width numeric, width of the plot in inches, defaults to 7
#'
#' @param height numeric, height of the plot in inches, defaults to 6
#' 
#' @param title string, title for the plot. If \code{NULL} (default), no title is added.
#'
#' @return Returns the ggplot plot object invisibly
#'
#' @keywords SIAMCAT plot.ordination
#'
#' @export
#'
#' @encoding UTF-8
#'
#' @examples
#' # Example data
#' ordination(siamcat_example)
#'
#' # Simple PCoA with Bray-Curtis dissimilarity
#' plot.ordination(siamcat_example)
#'
#' # Color by a metadata column
#' plot.ordination(siamcat_example, color.by = "disease")
#'
#' # NMDS with Jaccard distance, colored by a continuous variable
#' plot.ordination(siamcat_example,
#'     color.by = "age", name.color.by = "Age (years)")

# called like this to differentiate from phyloseq::plot.ordination
plot.ordination.siamcat <- function(
    siamcat, color.by = NULL, shape.by = NULL, name.color.by = NULL, name.shape.by = NULL, font.size = 14,
    palette = NULL, rename.values = NULL, alpha = 0.6,
    fn.plot = NULL, verbose = 1, width = 7, height = 6, title = NULL
) { 
    if (verbose > 1) message("+++ Starting plot.ordination")

    if(is.null(ordination(siamcat))) {
        stop("This siamcat object does not contain an ordination. Generate one with make.ordination.")
    }

    ord <- ordination(siamcat)$ord
    distance <- ordination(siamcat)$distance
    method <- ordination(siamcat)$method
    xlab <- sprintf("%s1 (%s)", method, distance)
    ylab <- sprintf("%s2 (%s)", method, distance)
    
    if (verbose > 1) message("+++ Extracting metadata")
    meta <- sample_data(siamcat@phyloseq)
    # samples may differ because of dropping zero-abundance samples
    # see make.ordination
    meta <- meta[rownames(meta) %in% rownames(ord$vectors)]
     # randomize order to avoid bias in plotting
    meta <- meta[sample(rownames(meta)), , drop=FALSE]
    temp.phyloseq <- phyloseq(meta=meta)

    if (verbose > 1) message("+++ Generating plot")
    if (!is.null(color.by)) {
        if (! (color.by %in% colnames(meta))) {
            stop("color.by column not found in metadata.")
        }
        if (is.null(name.color.by)) name.color.by <- color.by
    }
    if (!is.null(shape.by) && shape.by %in% colnames(meta)) {
        if (! (shape.by %in% colnames(meta))) {
            stop("shape.by column not found in metadata.")
        }
        if (length(unique(meta[[shape.by]])) > 15) {
            stop("shape.by column has more than 15 unique values.")
        }
        shapes <- 1:length(unique(meta[[shape.by]]))
        set_shapes <- TRUE
        if (is.null(name.shape.by)) name.shape.by <- shape.by
    } else {
        set_shapes <- FALSE
    }

    p <- phyloseq::plot_ordination(temp.phyloseq, ord, color=color.by, shape=shape.by) +
        labs(x=xlab, y=ylab, color=name.color.by, shape=name.shape.by)

    # override the shape
    if (is.null(shape.by)) {
        p$layers[[1]]$aes_params$shape <- 19
    }

    # override alpha
    if (!(is.numeric(alpha) && length(alpha) == 1 && alpha <= 1 && alpha >= 0)) stop ("Value for parameter alpha must be a single float between 0 and 1.")
    p$layers[[1]]$aes_params$alpha <- alpha

    # set color palette
    if (!is.null(color.by)) {
        meta_col <- meta[[color.by]]
        if (is.numeric(meta_col)) {
            if (is.null(palette)) palette <- "RdBu"
            p <- p + scale_color_distiller(palette = palette)
        } else {
            if (is.null(rename.values)) {
                labels <- unique(meta_col)
                names(labels) <- labels
            } else {
                labels <- rename.values
            }
            if (!is.null(palette)) {
                if (length(unique(meta_col)) > length(palette)) {
                    warning("The number of unique values in the color.by column is greater than the length of the provided palette. Discarding the provided palette and using the default palette instead.")
                    palette <- NULL
                }
            }
            if (is.null(palette)) {
                if (length(unique(meta_col)) <= length(okabe_palette)) {
                    palette <- okabe_palette
                } else {
                    palette <- RColorBrewer::brewer.pal(n = length(unique(meta_col)), name = "Set2")
                }
            }
            p <- p + scale_color_manual(values = palette, labels = labels, breaks = names(labels))
        }
    }

    if (set_shapes){
        p <- p + scale_shape_manual(values = shapes)
    }

    p <- p + theme_siamcat(font.size)

    if (!is.null(title)) {
        if (!(is.character(title) && length(title) == 1)) {
            stop("Title must be a single string of length 1.")
        }
        p <- p + ggtitle(title)
    }

    if (verbose > 1) message("+++ Producing plot file")
    if (!is.null(fn.plot)) {
         ggsave(fn.plot, p, bg="white", width=width, height=height)
    }

    return(p)
}