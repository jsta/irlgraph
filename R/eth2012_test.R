#==================================================================#

# dm<-as.matrix(read.delim("inst/extdata/etherington20120code/cost-surface20x20.txt",skip = 6, na.strings = "-9999", header = FALSE, sep = " "))
# 
# egraph2<-irgraph(dm=dm)
# 
# result<-get_path(egraph2$egraph2, 171, csurf)
# 
# sp::plot(impute_na(csurf, egraph2$cells, egraph2$nullcells, result, egraph2$cellcoords))
# points(raster::xyFromCell(csurf,191))

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
# 
# allcoords<-xyFromCell(csurf,1:ncell(csurf))




#gdistance regular graph===========================================#
# x<-gdistance::transition(csurf, function(x) 1/mean(x), directions=8)
# #fromCoords<-xyFromCell(csurf,cells[118])
# fromCells<-cells[2]
# 
# #source("/home/jose/R/scripts/gdistance/R/internal-functions.R")
# tr <- transitionMatrix(x)
# tr <- rBind(tr,rep(0,nrow(tr)))
# tr <- cBind(tr,rep(0,nrow(tr)))
# 
# startNode <- nrow(tr) #extra node to serve as origin
# adjP <- cbind(rep(startNode, times=length(fromCells)), fromCells)
# tr[adjP] <- Inf
# adjacencyGraph <- graph.adjacency(tr, mode="directed", weighted=TRUE)
# 
# E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight		
# shortestPaths <- shortest.paths(adjacencyGraph, v=startNode)[-startNode]
# 
# result <- as(x, "RasterLayer")
# result <- setValues(result, shortestPaths)	
# plot(result)
