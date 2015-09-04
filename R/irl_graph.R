#'@name irl_graph
#'@title Generate an irregular landscape graph
#'@param dm matrix cost raster
#'@param poicoords matrix
#'@param grainprop numeric
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
#'plot(raster(dm))
#'text(raster(dm), which(raster(dm)[] == 1))    
#'
#'graph <- irl_graph(dm)
#'plot(graph$graph)

irl_graph <- function(dm, poicoords = NA, grainprop = 0.25){
  csurf <- raster::raster(nrows=dim(dm)[1], ncols=dim(dm)[2], resolution=1, xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
  csurf[] <- dm
  
  get_cells <- function(dm, grainprop, poicoords){
    
    csurf <- raster::raster(nrows=dim(dm)[1], ncols=dim(dm)[2], resolution=1, xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
    csurf[] <- dm
    
    xmax <- dim(csurf)[1]
    ymax <- dim(csurf)[2]
    
    allcells <- raster::rasterToPoints(is.na(csurf))
    nullcells <- which(is.na(raster::extract(csurf,allcells[,1:2])))
    limitcells <- c(raster::cellFromRow(csurf,c(1,ymax)),raster::cellFromCol(csurf,c(1,xmax)))
    limitcells <- limitcells[!(limitcells %in% nullcells)]
    limitcells <- c(limitcells,which(raster::extract(raster::boundaries(csurf),allcells[,1:2]) == 1))
    limitcells <- limitcells[!duplicated(limitcells)]
    
    vipcells <- which(raster::extract(raster::boundaries(csurf,classes = TRUE, type="outer"),allcells[,1:2]) == 1)
    vipcells <- vipcells[!(vipcells %in% limitcells)]
    vipcells <- vipcells[!(vipcells %in% nullcells)]
    
    uncells <- raster::cellFromXY(csurf, allcells)[!raster::cellFromXY(csurf, allcells) %in% c(nullcells, vipcells, limitcells)]
    #set.seed(123)
    graincells <- sample(uncells, length(uncells) * grainprop)
    
    poicells<-NULL
    if(length(poicoords)>1){
    poicells <- apply(poicoords, 1, function(x) raster::cellFromXY(csurf, x))
    }
    
    cells <- c(limitcells, vipcells, graincells, poicells)
    cells <- cells[!duplicated(cells)]
    cells <- cells[order(cells)]
    
    return(list(cells =cells, nullcells = nullcells, limitcells = limitcells, vipcells = vipcells, graincells = graincells, poicells = poicells))
  }
  
  cells <- get_cells(dm = dm, poicoords = poicoords, grainprop = grainprop)
  poicells <- cells$poicells
  limitcells <- cells$limitcells
  vipcells <- cells$vipcells
  graincells <- cells$graincells
  nullcells <- cells$nullcells
  nullcoords <- raster::xyFromCell(csurf, nullcells)[,1:2]
  cells <- cells$cells
  cellcoords <- raster::xyFromCell(csurf, cells)[,1:2]
  
  allcells <- c(cells, nullcells)
  allcells <- allcells[order(allcells)]
  allcoords <- rbind(cellcoords, nullcoords)
  allcoords <- allcoords[order(-allcoords[,2], allcoords[,1]),]
  
  #create graph====================================================#
  create_tri <- function(cellcoords){
    #cellcoords<-allcoords
    deldir::deldir(cellcoords[,1], cellcoords[,2], suppressMsge = TRUE)
  }
  
  dtri <- create_tri(allcoords)
  
  #create adjacency matrix===========================================#
  create_adjmat <- function(dtri, cells){
    
    orig <- allcells[dtri$delsgs[,5]]
    neigh <- allcells[dtri$delsgs[,6]]
    keep <- !(orig %in% nullcells) & !(neigh %in% nullcells)
    
    orig <- orig[keep]
    neigh <- neigh[keep]
    elist <- cbind(orig, neigh)
    
    edgelengths<-rep(1, nrow(elist))
    celldiff <- elist[,1] - elist [,2]
    longedges<-which(!celldiff %in% c(1, dim(csurf)[1]))#dim(csurf)[1] - 1, dim(csurf)[1]  + 1
    longedgecoords_from<-raster::xyFromCell(csurf, elist[longedges,1])
    longedgecoords_to<-raster::xyFromCell(csurf, elist[longedges,2])
    
    elength_from_coords <- function(from, to){
      a <- from[1] - to[1]
      b <- from[2] - to[2]
      sqrt(a^2 + b^2)
    }
    
    edgelengths[longedges] <- apply(cbind(longedgecoords_from, longedgecoords_to), 1, function(x) elength_from_coords(x[1:2], x[3:4]))
    
    #if(length(unique(elist[,1]))!=length(cells)){
    #  warning("Mismatch between edgelist and nonnull cells")
    #}
    
    datavals <- cbind(raster::getValues(csurf)[elist[,1]], raster::getValues(csurf)[elist[,2]])
    datavals <- apply(datavals, 1, function(x) mean(x, na.rm = TRUE))
    
    adjmat <- Matrix::Matrix(0, nrow = raster::ncell(csurf), ncol =  raster::ncell(csurf))
    adjmat[elist] <- as.vector(1 / (datavals * edgelengths))
    #adjmat<-Matrix::forceSymmetric(adjmat)
    #as(adjmat, "symmetricMatrix")
    adjmat
  }
  
  adjmat <- create_adjmat(dtri, cells)
  
  egraph2 <- igraph::graph.adjacency(adjmat, mode = "directed", weighted = TRUE)
  igraph::E(egraph2)$weight <- 1 / igraph::E(egraph2)$weight 
  
  list(cells = cells, cellcoords = cellcoords, nullcells = nullcells, graph = egraph2, tri = dtri, limitcells = limitcells, vipcells = vipcells, graincells = graincells, poicells = poicells)
}
