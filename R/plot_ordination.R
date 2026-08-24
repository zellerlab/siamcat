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
#' @param rename.values named vector, used to rename the values of the color.by variable.
#' The names of the vector are the original values, and the values of the vector are the new names.
#' If \code{NULL} (default), no renaming is performed.
#' The order of the names in the vector determines the order of the labels in the legend.
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
    siamcat, color.by = NULL, name.color.by = NULL, palette = NULL, font.size = 14,
    rename.values = NULL, fn.plot = NULL, verbose = 1, width = 7, height = 6,
    title = NULL
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
    temp.phyloseq <- phyloseq(meta=meta)

    if (verbose > 1) message("+++ Generating plot")
    if (!is.null(color.by) && color.by %in% colnames(meta)) {
        if (is.null(name.color.by)) name.color.by <- color.by
        p <- phyloseq::plot_ordination(temp.phyloseq, ord, color=color.by) +
            labs(x=xlab, y=ylab, color=name.color.by)
        # override the shape
        p$layers[[1]]$aes_params$shape <- 1
        meta_col <- meta[[color.by]]
        if (is.numeric(meta_col)) {
            if (is.null(palette)) palette <- "RdBu"
            p <- p + scale_color_distiller(palette = palette)
        } else {
            if (is.null(palette)) palette <- okabe_palette
            if (is.null(rename.values)) {
                labels <- unique(meta_col)
                names(labels) <- labels
            } else {
                labels <- rename.values
            }
            if (length(unique(meta_col)) <= length(okabe_palette)) {
                p <- p + scale_color_manual(values = palette, , labels = labels, breaks = names(labels))
            } else warning("Refusing to set palette with fewer levels than groups in 'group.by'.")
        }
    } else if (!is.null(color.by)) {
        stop("color.by column not found in sample data.")
    } else {
        p <- phyloseq::plot_ordination(siamcat@phyloseq, ord) +
                labs(x=xlab, y=ylab) 
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