#'@name irgraph
#'@title Generate an irregular landscape graph
#'@param dm matrix cost raster
#'@import sp
#'@importFrom Matrix Matrix
#'@importFrom igraph graph.adjacency E
#'@importFrom raster raster extract boundaries cellFromCol cellFromRow cellFromXY xyFromCell
#'@importFrom deldir deldir
#'@export
#'@examples
#'
#'dm<-matrix(c( 1,  1,  1,  1,  1,
#'  NA, NA, NA, NA, 1, 
#'  1,  1,  1,  1,  1,  
#'  1,  1,  1,  1,  1,  
#'  1,  1,  1,  1,  1
#'  ), ncol = 5, byrow = TRUE)
#'  
#'irgraph(dm)

irgraph<-function(dm){
  csurf<-raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
  csurf[]<-dm
  
  #get cells=========================================================#
  get_cells<-function(dm,grainprop){
    
    csurf<-raster::raster(dm)
    
    xmax<-dim(csurf)[1]
    ymax<-dim(csurf)[2]
    
    allcells<-raster::rasterToPoints(is.na(csurf))
    nullcells<-which(is.na(raster::extract(csurf,allcells[,1:2])))
    
    limitcells<-c(raster::cellFromRow(csurf,c(1,ymax)),raster::cellFromCol(csurf,c(1,xmax)))
    limitcells<-limitcells[!(limitcells %in% nullcells)]
    limitcells<-c(limitcells,which(raster::extract(raster::boundaries(csurf),allcells[,1:2])==1))
    limitcells<-limitcells[!duplicated(limitcells)]
    
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
    
    return(list(cells =cells, nullcells = nullcells, limitcells = limitcells))
  }
  
  cells<-get_cells(dm,grainprop = 0.25)
  nullcells<-cells$nullcells
  nullcoords<-raster::xyFromCell(csurf,nullcells)[,1:2]
  cells<-cells$cells
  cellcoords<-raster::xyFromCell(csurf,cells)[,1:2]
  allcells<-c(cells,nullcells)
  allcells<-allcells[order(allcells)]
  allcoords<-rbind(cellcoords,nullcoords)
  allcoords<-allcoords[order(-allcoords[,2], allcoords[,1]),]
  
  #create graph====================================================#
  create_tri<-function(cellcoords){
    #cellcoords<-allcoords
    deldir::deldir(cellcoords[,1], cellcoords[,2])
  }
  
  dtri<-create_tri(allcoords)
  
  
  #create adjacency matrix===========================================#
  
  create_adjmat<-function(dtri,cells){
    orig<-allcells[dtri$delsgs[,5]]
    neigh<-allcells[dtri$delsgs[,6]]
    keep<-!(orig %in% nullcells) & !(neigh %in% nullcells)
    
    orig<-orig[keep]
    neigh<-neigh[keep]
    elist<-cbind(orig, neigh)
    
    #if(length(unique(elist[,1]))!=length(cells)){
    #  warning("Mismatch between edgelist and nonnull cells")
    #}
    
    datavals<-cbind(raster::getValues(csurf)[elist[,1]], raster::getValues(csurf)[elist[,2]])
    datavals<-apply(datavals, 1, function(x) mean(x, na.rm = TRUE))
    
    adjmat<-Matrix::Matrix(0, nrow = raster::ncell(csurf), ncol =  raster::ncell(csurf))
    adjmat[elist]<-as.vector(1/datavals)
    #adjmat<-Matrix::forceSymmetric(adjmat)
    #as(adjmat, "symmetricMatrix")
    adjmat
  }
  
  adjmat<-create_adjmat(dtri,cells)
  
  #egraph2<-igraph::graph.adjacency(adjmat, mode = "directed", weighted = TRUE)
  egraph2<-igraph::graph.adjacency(adjmat, mode = "directed", weighted = TRUE)
  igraph::E(egraph2)$weight <- 1/igraph::E(egraph2)$weight 
  list(cells = cells, cellcoords = cellcoords, nullcells = nullcells, egraph2 = egraph2)
}

