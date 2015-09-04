#'@name acc_path
#'@title Generate a shortest path raster from a specific graph node
#'@export
#'@importFrom igraph shortest.paths
#'@param egraph2
#'@param snode
#'@param scoord
#'@param csurf

acc_path <- function(graph, snode = NULL, scoord = NULL, costsurf){
  
  #ADD A CHECK TO MAKE SURE SNODE EXISTS IN GRID!!!
  
  if(is.null(snode) & all(is.null(scoord))){
    stop("Must supply either a starting node or coordinates")
  }
  
  if(is.null(snode)){
    snode <- raster::cellFromXY(costsurf, scoord)
  }
  
  #snode <- 2
  spaths <- igraph::shortest.paths(graph, v = snode)
  
  result <- as(costsurf, "RasterLayer")
  result[] <- NA
  result[] <- spaths
  result
}
