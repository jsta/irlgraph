# library(sp)
# library(tripack)
# library(igraph)

dm<-as.matrix(read.delim("inst/extdata/etherington20120code/cost-surface20x20.txt",skip = 6, na.strings = "-9999", header = FALSE, sep = " "))

csurf<-raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
csurf[]<-dm

#dm<-dm[5:10,5:10]

irgraph<-function(dm){
  csurf<-raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
  csurf[]<-dm
  
  #get cells=========================================================#
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
  allcells<-c(cells,nullcells)
  allcoords<-rbind(cellcoords,nullcoords)
  
  #create graph======================================================#
  create_tri<-function(cellcoords){
    cellcoords[2,2]<-cellcoords[2,2] + 0.01 #avoids jitter error
    tripack::tri.mesh(cellcoords[,1],cellcoords[,2])
  }
  
  dtri<-create_tri(allcoords)
  
  #create adjacency matrix===========================================#
  
  create_adjmat<-function(dtri,cells){
    #code block adapted from stackoverflow (search term: how can i extract distance between the points after a delaunay triangulation with deldir in R?)
    dneigh<-tripack::neighbours(dtri)
    orig<-allcells[rep(c(1:length(dneigh)), sapply(dneigh, length))]
    neigh<-allcells[unlist(dneigh)]
    keep<-!(orig %in% nullcells) & !(neigh %in% nullcells)
    
    orig<-orig[keep]
    neigh<-neigh[keep]
    elist<-cbind(orig, neigh)
    
    if(length(unique(elist[,1]))!=length(cells)){
      warning("Mismatch between edgelist and nonnull cells")
    }
    
    datavals<-cbind(raster::getValues(csurf)[elist[,1]], raster::getValues(csurf)[elist[,2]])
    datavals<-apply(datavals, 1, function(x) mean(x, na.rm = TRUE))
    
    adjmat<-Matrix::Matrix(0, nrow = raster::ncell(csurf), ncol =  raster::ncell(csurf))
    adjmat[elist]<-as.vector(1/datavals)
    adjmat<-Matrix::forceSymmetric(adjmat)
    adjmat
  }
  
  adjmat<-create_adjmat(dtri,cells)
  
  #egraph2<-igraph::graph.adjacency(adjmat, mode = "directed", weighted = TRUE)
  egraph2<-igraph::graph.adjacency(adjmat, mode = "undirected", weighted = TRUE)
  igraph::E(egraph2)$weight <- 1/igraph::E(egraph2)$weight 
  egraph2
}

egraph2<-irgraph(dm=dm)

#calc from specific origin=========================================#
#spaths<-igraph::shortest.paths(egraph2, v = 118)
get_path<-function(egraph2, snode, csurf){
  spaths<-igraph::shortest.paths(egraph2, v = snode)#171
  #igraph::shortest.paths(egraph2, v = 2, to = 9)
  
  result<-as(csurf, "RasterLayer")
  result[]<-NA
  result[]<-spaths
  result
}
result<-get_path(egraph2, 2, csurf)

#imputation========================================================#
impute_na<-function(csurf, cells, result, cellcoords){
  nonnullcells<-(1:length(csurf))[!(1:length(csurf) %in% cells)]
  nonnullcells<-nonnullcells[!(nonnullcells %in% nullcells)]
  
  fcells<-spatstat::as.ppp(raster::xyFromCell(csurf,cells),c(0,20,0,20))
  icells<-spatstat::as.ppp(raster::xyFromCell(csurf,nonnullcells),c(0,20,0,20))
  ncross<-spatstat::nncross(icells,fcells, what = "which")
  result[nonnullcells]<-result[raster::cellFromXY(csurf, cellcoords[ncross,])]
  result
}

sp::plot(impute_na(csurf, cells, result, cellcoords))
#==================================================================#

# sp::plot(csurf)
# points(raster::xyFromCell(csurf,nullcells), col = "green")
# points(raster::xyFromCell(csurf,cells), col = "blue")
# points(raster::xyFromCell(csurf,nonnullcells), col = "orange")

# plot(dtri)

# sp::plot(result)
# points(cellcoords[118,1],cellcoords[118,2])
# result[161:180]

#tripack - igraph regular graph====================================#

allcoords<-xyFromCell(csurf,1:ncell(csurf))




#gdistance regular graph===========================================#
x<-gdistance::transition(csurf, function(x) 1/mean(x), directions=8)
#fromCoords<-xyFromCell(csurf,cells[118])
fromCells<-cells[2]

#source("/home/jose/R/scripts/gdistance/R/internal-functions.R")
tr <- transitionMatrix(x)
tr <- rBind(tr,rep(0,nrow(tr)))
tr <- cBind(tr,rep(0,nrow(tr)))

startNode <- nrow(tr) #extra node to serve as origin
adjP <- cbind(rep(startNode, times=length(fromCells)), fromCells)
tr[adjP] <- Inf
adjacencyGraph <- graph.adjacency(tr, mode="directed", weighted=TRUE)

E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight		
shortestPaths <- shortest.paths(adjacencyGraph, v=startNode)[-startNode]

result <- as(x, "RasterLayer")
result <- setValues(result, shortestPaths)	
plot(result)




