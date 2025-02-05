# Interface with snakemake 
library(funr)
library(optparse)
library(knitr)




#Rendering 
library(RColorBrewer)
library(ggplot2)
library(ggnewscale)
library(ggpubr)
library(plotly)
library(cowplot)


#Others 
#library(factoextra)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(data.table)


library(rstatix)
library(ggprism)

library(DT)
library(gridExtra)


Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}



multiplot <- function(..., plotlist = NULL, file, cols = 4, layout = NULL) {
  require(grid)

  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                 ncol = cols, nrow = ceiling(numPlots/cols))
}

if (numPlots == 1) {
print(plots[[1]])

} else {
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

for (i in 1:numPlots) {
  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                  layout.pos.col = matchidx$col))
 }
}
 }


wrap_plots <- function(..., ncol = NULL, nrow = NULL, byrow = NULL,
                       widths = NULL, heights = NULL, guides = NULL,
                       tag_level = NULL, design = NULL) {
  if (is_valid_plot(..1)) {
    plots <- list(...)
  } else if (is.list(..1)) {
    plots <- ..1
  } else {
    stop('Can only wrap ggplot and/or grob objects or a list of them', call. = FALSE)
  }
  if (!all(vapply(plots, is_valid_plot, logical(1)))) stop('Only know how to add ggplots and/or grobs', call. = FALSE)
  if (!is.null(names(plots)) && !is.null(design) && is.character(design)) {
    area_names <- unique(trimws(strsplit(design, '')[[1]]))
    area_names <- sort(setdiff(area_names, c('', '#')))
    if (all(names(plots) %in% area_names)) {
      plot_list <- vector('list', length(area_names))
      names(plot_list) <- area_names
      plot_list[names(plots)] <- plots
      plot_list[vapply(plot_list, is.null, logical(1))] <- list(plot_spacer())
      plots <- plot_list
    }
  }
  Reduce(`+`, plots, init = plot_filler()) + plot_layout(
    ncol = ncol, nrow = nrow, byrow = byrow, widths = widths, heights = heights,
    guides = guides, tag_level = tag_level, design = design
  )
}

#' @importFrom ggplot2 is.ggplot
#' @importFrom grid is.grob
is_valid_plot <- function(x) is.ggplot(x) || is.grob(x)
 
subsetList <- function(myList, elementNames) {
    lapply(elementNames, FUN=function(x) myList[[x]])
}
