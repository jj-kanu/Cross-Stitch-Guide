library(imager)
library(tidyverse)
library(tidymodels)
library(sp)
library(scales)
library(cowplot)
#devtools::install_github("sharlagelfand/dmc")
library(dmc)
library(scales)


##FUNCTIONS
process_image <- function(image_file_name, k_list){
  ## process_image converts a picture to cluster_info data that
  ## can be used for other functions.
  ## 
  ## Input:
  ## - image_file_name: The name of an image file in directory
  ##                     or filepath to image file.
  ## - k_list: A vector of numbers determining how many clusters
  ##           will be used for k-means.
  ## 
  ## Output:
  ##   - Returns a list with two indices. The first index contains
  ##     a tibble with a column of kmeans and a column of tidied
  ##     cluster centres with dmc colour info added, each mapped to 
  ##     a row associated with an entry of k_list.
  ##     The second tibble contains the image data for image_file_name
  ##     with only x,y,R,G,B columns.
  ## 
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager')
  ##   k_vec <- c(6,7,8)
  ##   process_image(fpath, k_vec)
  
  im <- imager::load.image(image_file_name)
  tidy_dat <-as.data.frame(im, wide = "c")%>% 
    rename(R = c.1, G = c.2, B = c.3)
  
  dat <-select(tidy_dat,c(-x,-y))
  im_dat <- select(tidy_dat,x,y,R,G,B)
    
    cluster_info <- 
     tibble(k = k_list) %>%
      mutate(
       kclust = map(k, ~kmeans(x=dat, centers = .x, nstart=4)),
      )
    cent_temp<- cluster_info %>% pull(kclust)
    count=1
    for(i in cent_temp){
      cent_temp[[count]] = tidy(cent_temp[[count]])%>% mutate(col =rgb(R,G,B))
      tmp <- cent_temp[[count]] %>% pull(col)
      tmp_list = vector(mode = "list", length = nrow(cent_temp[[count]]))
      tmp_count=1
      for(j in tmp){
        tmp_list[[tmp_count]] = dmc(tmp[[tmp_count]])
        tmp_count = tmp_count+1
      }
      
      tmp_list <- as_tibble_col(tmp_list)
      tmp_list <- tmp_list %>%
        rename(dmc = value)
      
      cent_temp[[count]] <- cent_temp[[count]] %>%
        add_column(tmp_list)
      
      count = count+1
    }
    cent_temp <- as_tibble_col(cent_temp)
    cluster_info <- cluster_info %>%
      add_column(cent_temp)%>%
      rename(centres = value)
    return(list(cluster_info,im_dat))
}

scree_plot <- function(cluster_info){
  ## scree_plot produces a scree plot.
  ##
  ## Input: 
  ## -cluster_info: A list produced from the process_image function
  ##                containing kmeans, tidied clusters, and image data.
  ## 
  ## Output:
  ## -Returns a scree plot based on the kmeans in cluster_info.
  ##
  ## Example: 
  ## cluster_info <- process_image("squid.jpg",c(6:10))
  ## scree_plot_ex <- scree_plot(cluster_info)
  ## scree_plot
  temp_clusts <-cluster_info[[1]]%>%
    mutate(
      glanced =map(kclust, glance)
    )
  
  clusterings <-temp_clusts%>%
    unnest(cols =c(glanced))
  
  ggplot(clusterings,
         aes(k,tot.withinss))+
    labs(x="Number of Clusters",
         y="Within Groups Sum of Squares")+
    geom_line()+
    geom_point()
}

colour_strips <- function(cluster_info){
  ## colour_strips produces colour strips with the DMC colour closest 
  ## to the cluster centre colour.
  ## 
  ## Input:
  ##   -cluster_info: A list produced from the process_image function
  ##                 containing kmeans, tidied clusters, and image data.
  ## Output:
  ##   -Returns a list with a size based on number of potential clusters
  ##   used in process_image(). Each index in list has a colour strip
  ##   made with hex codes for the dmc colours most associated with the
  ##   clustered colours of the original image.
  ## 
  ## Example: 
  ## cluster_info <- process_image("squid.jpg",c(6:10))
  ## strip_ex <- colour_strips(cluster_info)
  ## strip_ex[[2]]
  
  strips = vector(mode = "list", length = nrow(cluster_info[[1]]))
  for (row in 1:nrow(cluster_info[[1]])){
    tmp_strips<-tibble(cluster_info[[1]]$centres[[row]]$dmc)
    hexes = vector(mode = "list", length = nrow(tmp_strips))
    for (tmp_row in 1:nrow(tmp_strips)){
      hexes[[tmp_row]] = tmp_strips[[1]][[tmp_row]]$hex
    }
    t <- tibble(colours = hexes,
                squares = purrr::map(colours, ~ square(.x, 2)))
    strips[[row]] = plot_grid(plotlist = t$squares)
  }
  return(strips)
}

make_pattern <- function(cluster_info, knum, x_size, black_white =FALSE, background_colour = NULL){
  ## make_pattern plots a cross-stitch using a reduced version 
  ## of the original image.
  ## 
  ## Input:
  ##   -cluster_info: The output of process_image().
  ##   -k: The chosen cluster size.
  ##   -x_size: The (approximate) total number of possible stitches in 
  ##            the horizontal direction
  ##   -black_white: (logical) Print the pattern in black and white (TRUE) 
  ##                  or colour (FALSE,default)
  ##   -background_colour: The colour of the background, which should not 
  ##                       be stitched in the pattern. (Default is to not 
  ##                       have a colour)
  ## 
  ## Output:
  ##   -Returns a cross-stitch pattern that can be followed, complete with
  ##    a legend that shows hasthread colour, and a guide grid.
  ## 
  ## Example:
  ## cluster_info <- process_image("squid.jpg",c(6:10))
  ## pattern <- make_pattern(cluster_info, 8, 50)
  ## pattern
  ## #For black and white:
  ## pattern1 <- make_pattern(cluster_info, 8, 50,black_white=TRUE)
  ## pattern1
  ## #To add a background colour:
  ## pattern2 <- make_pattern(cluster_info, 8, 50,background_colour="light grey")
  ## pattern2
  
  row_dat <- filter(cluster_info[[1]], k==knum)
  broomed_dat <-augment(row_dat$kclust[[1]], cluster_info[[2]])%>% rename(cluster = .cluster)
  agg_image <- change_resolution(broomed_dat, x_size)
  
  hexes = vector(mode = "list", length = nrow(row_dat$centres[[1]]))
  dmc_names = vector(mode = "list", length = nrow(row_dat$centres[[1]]))
  for (tmp_row in 1:nrow(row_dat$centres[[1]])){
    hexes[[tmp_row]] = row_dat$centres[[1]]$dmc[[tmp_row]]$hex[[1]]
    dmc_names[[tmp_row]] = paste(row_dat$centres[[1]]$dmc[[tmp_row]]$dmc,
                                 row_dat$centres[[1]]$dmc[[tmp_row]]$name)
  }
  
  gg = ggplot(agg_image,aes(x=x, y = y))
    if(!black_white){
    gg = gg+ geom_point(aes(name= NULL,
                   group = cluster,
                   col = cluster,
                   shape =cluster))+
    scale_colour_manual(name = NULL,
                        breaks = c(1:knum),
                        values = hexes,
                        label = dmc_names
                        )
    }
  else{
    gg = gg+ geom_point(aes(name= NULL,
                            group = cluster,
                            shape =cluster))+
      scale_colour_manual(name = NULL,
                          breaks = c(1:knum),
                          values = hexes,
                          label = dmc_names
      )
  }
    gg = gg + scale_shape_manual(name = NULL,
                       breaks = c(1:knum),
                       values = c(1:knum),
                        label = dmc_names
                       )
    if(is.null(background_colour)){
      gg = gg + theme(
        legend.position="bottom",
        legend.text = element_text(size = 6),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                        colour = "black"), 
        panel.grid.minor = element_line(size = 0.4, linetype = 'solid',
                                        colour = "light grey"), 
        panel.background = element_blank()
      )
    }
    else{
      gg = gg + theme(
        legend.position="bottom",
        legend.text = element_text(size = 6),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_line(size = 0.8, linetype = 'solid',
                                        colour = "black"), 
        panel.grid.minor = element_line(size = 0.4, linetype = 'solid',
                                        colour = "light grey"), 
        panel.background = element_rect(
          fill = background_colour,
          colour = background_colour)
      )
    }
    
    gg= gg+ scale_y_reverse()
    return(gg)
}

##HELPER FUNCTIONS
square <- function(x, label_size) { 
  ggplot()  + 
    coord_fixed(xlim=c(0,1), ylim = c(0,1))+  
    theme_void() + 
    theme(plot.background = element_rect(fill = x)) + 
    geom_text(aes(0.5,0.5),label = x , size = label_size)
}
change_resolution <- function(image_df, x_size)
{
  ## change_resolution(image_df, x_size) subsamples an image to produce
  ## a lower resolution image. Any non-coordinate columns in the data
  ## frame are summarized with their most common value in the larger
  ## grid cell.
  ##
  ## Input:
  ## - image_df: A data frame in wide format. The x-coordinate column MUST
  ##             be named 'x' and the y-coordinate column MUST be named 'y'.
  ##             Further columns have no naming restrictions.
  ## - x_size:   The number of cells in the x-direction. The number of cells
  ##             in the vertical direction will be computed to maintain the 
  ##             perspective. There is no guarantee that the exact number
  ##             of cells in the x-direction is x_size
  ##
  ## Output:
  ## - A data frame with the same column names as image_df, but with fewer 
  ##   entries that corresponds to the reduced resolution image.
  ##
  ## Example:
  ##   library(imager)
  ##   library(dplyr)
  ##   fpath <- system.file('extdata/Leonardo_Birds.jpg',package='imager') 
  ##   im <- load.image(fpath)
  ##   im_dat<- as.data.frame(im,wide = "c") %>% rename(R = c.1, G = c.2, B = c.3) %>%
  ##            select(x,y,R,G,B)
  ##   agg_image <- change_resolution(im_dat, 50)
  
  if(!require(sp)) {
    stop("The sp packages must be installed. Run install.packages(\"sp\") and then try again.")
  }
  if(!require(dplyr)) {
    stop("The dplyr packages must be installed. Run install.packages(\"dplyr\") and then try again.")
  }
  
  sp_dat <- image_df 
  gridded(sp_dat) = ~x+y
  
  persp = (gridparameters(sp_dat)$cells.dim[2]/gridparameters(sp_dat)$cells.dim[1])
  y_size = floor(x_size*persp)
  orig_x_size = gridparameters(sp_dat)$cells.dim[1]
  orig_y_size = gridparameters(sp_dat)$cells.dim[2]
  
  x_res = ceiling(orig_x_size/x_size)
  y_res = ceiling(orig_y_size/y_size)
  
  gt = GridTopology(c(0.5,0.5), c(x_res, y_res),
                    c(floor(orig_x_size/x_res), floor(orig_y_size/y_res)))
  SG = SpatialGrid(gt)
  agg = aggregate(sp_dat, SG, function(x) names(which.max(table(x)))[1] )
  agg@grid@cellsize <- c(1,1)
  df <- agg %>% as.data.frame %>% rename(x = s1, y = s2)  %>% select(colnames(image_df))
  
  return(df)
  
}

## Testing
# process_image("sponge.png",c(6,7,8))
# 
# make_pattern(cluster_info, 6, 50)
# 
# kclust <-kmeans(dat, centers =   7, nstart = 8)
# centres <-tidy(kclust)
# 
# centres <- centres%>% mutate(col =rgb(R,G,B))
# centres


