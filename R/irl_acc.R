#'@name irl_acc
#'@title Generate an accumulated cost surface from an irregular landscape graph
#'@export
#'@param dm matrix cost raster
#'@param poicoords matrix
#'@param grainprop numeric
#'@description This high level function allows for generation of an accumulated cost surface from a single function call. 
#'@examples \dontrun{
#'set.seed(123) #make reproducible
#'
#'dm <- as.matrix(read.delim(system.file(
#'"extdata/etherington20120code/cost-surface20x20.txt", package = "irlgraph"),
#'skip = 6, na.strings = "-9999", header = FALSE, sep = " "))
#'
#'costsurf <- raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0,
#' xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2]) #neccessary to set resolution
#'costsurf[] <- dm
#'
#'poicoords <- matrix(c(10.5, 10.5), ncol = 2)
#'
#'result <- irl_acc(dm, poicoords, grainprop = 0.25, costsurf = costsurf, scoord = c(10.5, 10.5))
#'}

irl_acc <- function(dm, poicoords = NA, cutoff = 0, grainprop = 0.25, costsurf, scoord, snode = NULL, irregular = TRUE){
  
  graph <- irl_graph(dm, poicoords = poicoords, grainprop = grainprop, irregular = irregular, cutoff = cutoff)
  result<-acc_path(graph = graph$graph, scoord = scoord, snode = snode, costsurf = costsurf)
  
  impute_na(result, costsurf, graph)
}