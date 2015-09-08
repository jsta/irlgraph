#'@name impute_na
#'@title Impute NA cells according to nearest neighbor
#'@export
#'@importFrom spatstat as.ppp nncross
#'@importFrom raster xyFromCell cellFromXY
#'@param csurf RasterLayer
#'@param cells numeric cell numbers
#'@param nullcells numeric cell numbers
#'@param result RasterLayer
#'@param cellcoords matrix xy coordinates

impute_na <- function(result, csurf, graph){
  
  cells<-graph$cells
  infcells <- cells[which(result[cells]==Inf)]
  
  nonnullcells <- (1:length(csurf))[!(1:length(csurf) %in% graph$cells)]
  nonnullcells <- nonnullcells[!(nonnullcells %in% graph$nullcells)]
  nonnullcells <- c(nonnullcells, infcells)
  
  cells <- cells[!cells %in% nonnullcells]
  
  cellcoords <- raster::xyFromCell(csurf, cells)
  nonnullcellcoords <- raster::xyFromCell(csurf, nonnullcells)
  
  fcells <- spatstat::as.ppp(cellcoords,c(0,dim(csurf)[1],0,dim(csurf)[2]))
  icells <- spatstat::as.ppp(nonnullcellcoords,c(0,dim(csurf)[1],0,dim(csurf)[2]))
  ncross <- spatstat::nncross(icells,fcells, what = "which")
  
  result[nonnullcells] <- result[raster::cellFromXY(csurf, cellcoords[ncross,])]
  result
}
