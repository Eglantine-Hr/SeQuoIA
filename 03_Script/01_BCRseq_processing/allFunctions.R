## @knitr load_functions

# #####################################################################################
# #####################################################################################
# #####################################################################################
########             Generate a datatable summarizing values                   ######## 
########    For environments (parameters), all values/variables are shown      ######## 
# #####################################################################################
# #####################################################################################
# #####################################################################################


showSimpleDT = function( dataToShow, rownames = TRUE, colnames = "Value")
{
  valuesDF = NULL;
  
  if(is.environment( dataToShow))
  {
    # Extract values from environment as character strings
    envValues = sapply(lapply(dataToShow, function(x) {ifelse(is.null(x), "NULL", x)}), paste, collapse = ", ");
    # Sort them by name and convert to data.frame
    valuesDF = data.frame("Value" = envValues[order(names(envValues))]);
  } else
  {
    valuesDF = dataToShow;
  }
  
  # Create a datatable with collected information
  # Create datatable
  datatable( as.data.frame(valuesDF), 
             class = "compact",
             rownames = rownames,
             colnames = colnames,
             options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                            autoWidth = FALSE,
                            columnDefs = list( # Center all columns
                              list( targets = "_all",
                                    className = 'dt-center')),
                            orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                            ordering = FALSE,
                            paging = FALSE, # Disable pagination (show all)
                            processing = TRUE, 
                            scrollCollapse = TRUE,
                            scroller = TRUE,  # Only load visible data
                            scrollX = TRUE,
                            scrollY = "525px",
                            stateSave = TRUE));
}



############# upper and lower SEM calculation########

lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

################ GEOM_FLAT_VIOLIN_PLOT_DEF##########


## define function of geom_flat_volin
# devtools::install_github(repo = "IndrajeetPatil/ggstatsplot")

"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}

geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }

GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    
    draw_key = draw_key_polygon,
    
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    
    required_aes = c("x", "y")
  )

#################
# Add minor ticks
##################

insert_minor <- function(major_labs, n_minor) {labs <- 
  c( sapply( major_labs, function(x) c(x, rep("", n_minor-1) ) ) )
labs[1:(length(labs)-n_minor+1)]}

####################
# %notin% operator 
####################
`%notin%` <- Negate(`%in%`)

########################################
# Emulate ggplot2 default color palette 
########################################
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}