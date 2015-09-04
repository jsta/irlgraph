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

impute_na <- function(csurf, cells, nullcells, result, cellcoords){
  nonnullcells <- (1:length(csurf))[!(1:length(csurf) %in% cells)]
  nonnullcells <- nonnullcells[!(nonnullcells %in% nullcells)]
  
  fcells <- spatstat::as.ppp(raster::xyFromCell(csurf,cells),c(0,20,0,20))
  icells <- spatstat::as.ppp(raster::xyFromCell(csurf,nonnullcells),c(0,20,0,20))
  ncross <- spatstat::nncross(icells,fcells, what = "which")
  result[nonnullcells] <- result[raster::cellFromXY(csurf, cellcoords[ncross,])]
  result
}
