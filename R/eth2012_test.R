# library(sp)
# library(tripack)
# library(igraph)

dm<-as.matrix(read.delim("inst/extdata/etherington20120code/cost-surface20x20.txt",skip = 6, na.strings = "-9999", header = FALSE, sep = " "))

r<-raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
r[]<-dm
csurf<-r

get_cells<-function(dm,grainprop){
  
  csurf<-raster::raster(dm)
  
  nullcells  <- NA
  limitcells <- NA
  vipcells   <- NA
  
  xmax<-dim(csurf)[1]
  ymax<-dim(csurf)[2]
  
  allcells<-raster::rasterToPoints(is.na(csurf))
  nullcells<-which(is.na(raster::extract(csurf,allcells[,1:2])))
  
  limitcells<-c(raster::cellFromRow(csurf,c(1,ymax)),raster::cellFromCol(csurf,c(1,xmax)))
  limitcells<-limitcells[!(limitcells %in% nullcells)]
  limitcells<-c(limitcells,which(raster::extract(raster::boundaries(csurf),allcells[,1:2])==1))
  
  vipcells<-which(raster::extract(raster::boundaries(csurf,classes = TRUE, type="outer"),allcells[,1:2])==1)
  vipcells<-vipcells[!(vipcells %in% limitcells)]
  vipcells<-vipcells[!(vipcells %in% nullcells)]
  
  uncells<-raster::cellFromXY(csurf, allcells)[!raster::cellFromXY(csurf, allcells) %in% c(nullcells, vipcells, limitcells)]
  set.seed(123)
  graincells<-sample(uncells, length(uncells)*grainprop)

  #   sp::plot(csurf)
  #   points(raster::xyFromCell(csurf,nullcells), col = "green")
  #   points(raster::xyFromCell(csurf,limitcells), col = "red")
  #   points(raster::xyFromCell(csurf,vipcells), col = "black")
  #   points(raster::xyFromCell(csurf,graincells), col = "blue")
  
  cells<-c(limitcells, vipcells, graincells)
  cells<-cells[!duplicated(cells)]
  cells<-cells[order(cells)]
  
  return(list(cells =cells, nullcells = nullcells))
}

cells<-get_cells(dm,grainprop = 0.25)
nullcells<-cells$nullcells
nullcoords<-raster::xyFromCell(csurf,nullcells)[,1:2]
cells<-cells$cells
cellcoords<-raster::xyFromCell(csurf,cells)[,1:2]
allcoords<-rbind(cellcoords,nullcoords)

cellcoords[2,2]<-cellcoords[2,2] + 0.01 #avoids jitter error
allcoords[2,2]<-allcoords[2,2] + 0.01 #avoids jitter error
dtri<-tripack::tri.mesh(cellcoords[,1],cellcoords[,2])

#code block adapted from stackoverflow (search term: how can i extract distance between the points after a delaunay triangulation with deldir in R?)
dneigh<-tripack::neighbours(dtri)
orig<-rep(c(1:length(dneigh)), sapply(dneigh, length))
neigh<-unlist(dneigh)
keep<-!(orig %in% nullcells) & !(neigh %in% nullcells)
orig<-orig[keep]
neigh<-neigh[keep]
elist<-cbind(orig, neigh)

egraph<-igraph::graph_from_edgelist(elist, directed = TRUE)
spaths<-igraph::shortest.paths(egraph, v = 118)

#sp::plot(csurf)

result<-as(csurf, "RasterLayer")
result[cells]<-spaths
result[!(1:length(csurf) %in% cells)]<-NA
#sp::plot(result)
#points(cellcoords[118,1],cellcoords[118,2])

# result[361:370]

nonnullcells<-(1:length(csurf))[!(1:length(csurf) %in% cells)]
nonnullcells<-nonnullcells[!(nonnullcells %in% nullcells)]
infcells<-which(result[]==Inf)
nonnullcells<-c(nonnullcells,infcells)
cellcoords<-cellcoords[!(cells %in% infcells),]
cells<-cells[!(cells %in% infcells)]

#   sp::plot(csurf)
#   points(raster::xyFromCell(csurf,nullcells), col = "green")
#   points(raster::xyFromCell(csurf,cells), col = "blue")
#   points(raster::xyFromCell(csurf,nonnullcells), col = "orange")

#imputation==============================================================================#
fcells<-spatstat::as.ppp(raster::xyFromCell(csurf,cells),c(0,20,0,20))
icells<-spatstat::as.ppp(raster::xyFromCell(csurf,nonnullcells),c(0,20,0,20))
ncross<-spatstat::nncross(icells,fcells, what = "which")

#result[nonnullcells]<-result[ncross]
result[nonnullcells]<-result[cellFromXY(csurf, cellcoords[ncross,])]

plot(result)
points(cellcoords[118,1],cellcoords[118,2])
result[161:180]
