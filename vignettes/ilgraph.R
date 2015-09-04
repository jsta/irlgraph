## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github("jsta/ilgraph")

## ----eval=TRUE, fig.align='center', fig.height=4, fig.width=6------------
library(ilgraph)

dm<-as.matrix(read.delim(system.file("extdata/etherington20120code/cost-surface20x20.txt",
package = "ilgraph"), skip = 6, na.strings = "-9999", header = FALSE, sep = " "))

egraph2<-irgraph(dm=dm)

costsurf <- raster::raster(nrows=dim(dm)[1],ncols=dim(dm)[2],resolution=1,xmn=0, xmx = dim(dm)[1], ymn = 0, ymx = dim(dm)[2])
costsurf[] <- dm


sp::plot(costsurf)
legend("top",inset  = c(0, -0.15), legend = c("Points of interest", "Very important", "Landscape limit", "Null data", "Ecological grain"), col = viridis::viridis(5), pch = 1, horiz = TRUE, xpd = TRUE, cex = 0.7)

points(raster::xyFromCell(costsurf,egraph2$nullcells), col = viridis::viridis(5)[4])
points(raster::xyFromCell(costsurf,egraph2$cells), col = viridis::viridis(5)[1])

#plot(dtri)


## ----eval=FALSE----------------------------------------------------------
#  result<-get_path(egraph2$egraph2, snode = 171, csurf)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

