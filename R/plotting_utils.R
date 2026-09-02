#!/usr/bin/Rscript
### SIAMCAT - Statistical Inference of Associations between
### Microbial Communities And host phenoTypes R flavor EMBL
### Heidelberg 2012-2018 GNU GPL 3.0


#'@keywords internal ggplot plotting theme used in the package
theme_siamcat <- function(font.size){
    line.size <- 0.5
    half.line <- font.size / 2
    rel.small <- 12 / 14
    small.size <- rel.small * font.size

    theme_grey(base_size = font.size) %+replace%
        theme(
            line = element_line(
                color = "black", linetype = 1,
                lineend = "butt", linewidth = line.size
            ),
            text = element_text(
                size = font.size, hjust = 0.5, vjust = 0.5, angle = 0,
                lineheight = .9, margin = margin(), debug = FALSE
            ),
            axis.line = element_blank(),
            axis.text = element_text(color = "black", size = small.size),
            axis.text.x = element_text(
                margin = margin(t = small.size / 4), vjust = 1
            ),
            axis.text.y = element_text(
                margin = margin(r = small.size / 4), hjust = 1
            ),
            axis.ticks = element_line(color = "black", linewidth = line.size),
            axis.ticks.length = unit(half.line / 2, "pt"),
            axis.title.x = element_text(
                margin = margin(t = half.line / 2), vjust = 1
            ),
            axis.title.y = element_text(
                angle = 90, margin = margin(r = half.line / 2), vjust = 1
            ),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", linewidth = line.size),
            panel.grid = element_blank(),
            legend.background = element_blank(),
            legend.spacing = unit(font.size, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            legend.key = element_blank(),
            legend.key.size = unit(1.1 * font.size, "pt"),
            legend.text = element_text(size = rel(rel.small)),
            legend.title = element_text(hjust = 0),
            legend.box.background = element_blank(),
            legend.box.spacing = unit(font.size, "pt"),
            plot.margin = margin(half.line, half.line, half.line, half.line),
            complete = TRUE
        )
}

#'@keywords internal
theme_siamcat_hist <- function(font.size){
    line.size <- 0.5
    half.line <- font.size / 2
    theme_siamcat(font.size) %+replace%
        theme(
            panel.border = element_blank(),
            axis.line.y = element_line(color = "black", linewidth = line.size, lineend = "square"),
            axis.text.x = element_text(margin = margin(t = half.line)),
            axis.ticks.x = element_blank(),
            plot.title = element_text(
                hjust = 0.5, face = "bold",
                margin = margin(b = half.line)
            )
        )
}

#'@keywords internal
theme_siamcat_forest <- function(font.size){
    line.size <- 0.5
    half.line <- font.size / 2
    theme_siamcat(font.size) %+replace%
        theme(
            axis.title.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.border = element_rect(colour = "black"),
            panel.grid.major.x = element_line(color = "black", linewidth = 0.05),
            legend.position = "none"
        )
}

#'@keywords internal
okabe_palette <- c(
    "black", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

#'@keywords internal
pos_neg_neut_palette <- c("#832424", "#3A3A98", "gray")