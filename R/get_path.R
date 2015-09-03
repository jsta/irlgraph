#'@name get_path
#'@title Generate a shortest path raster from a specific graph node
#'@export
#'@importFrom igraph shortest.paths
#'@param egraph2
#'@param snode
#'@param csurf

get_path <- function(egraph2, snode, csurf){
  
  #ADD A CHECK TO MAKE SURE SNODE EXISTS IN GRID!!!
  
  spaths <- igraph::shortest.paths(egraph2, v = snode)#171
  
  result <- as(csurf, "RasterLayer")
  result[] <- NA
  result[] <- spaths
  result
}
