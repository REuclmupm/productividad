## First attempt using Parallel.

library(raster)
library(ncdf4)
library(rasterVis)
library(rgdal)
library(maps)
library(maptools)
library(parallel)

## Calculate the annual productivity using parallel function

##########################################################################
## 1. DATA
#########################################################################

load('data/boundaries.Rdata')
load('data/linea.Rdata')
load('data/mascaraClustersSat.Rdata')
load('data/cmsafRasterdailyRadiation.Rdata')

## Creo un raster con valores medios mensuales para que sea más pequeño

idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-30'), 'day')
SISS <- setZ(SISS, idx)

month <- function(x) as.numeric(format(x, '%m'))

SISmm <- zApply(SISS, by=month, fun='mean') ## medias mensuales de irradiancia

## paso a irradiacion Wh/m2 para poder emplearlo en el calculo de la GEf

SISmm <- SISmm*24

########################################################################
## 2. Gef
#######################################################################

###############################################################
## mi prueba
#####################

## Esta es la función que quiero hacer el paralelo

fooGef <- function(g0, idx){
    n <- length(g0)
    lat <- g0[1]
    g0d <- list(file = zoo(data.frame(G0 = g0[2:n]),  idx),
                lat = lat)
    Prod <- calc(g0d, fun=function(x){
                         prod <- prodGCPV(lat=lat, dataRad= g0d, keep.night=FALSE)
                         result <- as.data.frameY(prod)[c('Eac', 'Edc', 'Yf')] ##the results are yearly values
                         as.numeric(result)
                     }
                 )
                 ## para sacar los valores anuales. Yf es productividad.
                 Prod_Yf <- as.numeric(as.data.frame(subset(Prod, 'layer.3')))
                 return(Prod_Yf)
}



##################

## Este es el código que hace en paralelo la función que se le expecifica en FUN. Para las funciones mean, sum etc funciona, pero da problemas con la que yo he difinido. Creo que el problema está en que lo que devuelve la función fooGef es un raster y lo que hay al final del codigo de paralelo es para poner en formato raster algo que estaba como data.frame.

fooParallel <- function(data, filename="", nodes=detectCores(), blocks=6,...){
    ## latitude values as a new raster
    y <- init(data, v='y')
    idx <- getZ(data)

    bs <- blockSize(data, minblocks=blocks*nodes)


fooGef <- function(g0){
    n <- length(g0)
    lat <- g0[1]
    g0d <- list(file = zoo(data.frame(G0 = g0[2:n]),  idx),
                lat = lat)
    Prod <- lapply(g0d, FUN=function(y) calc(y, fun=function(x){
                         prod <- prodGCPV(lat=lat, dataRad= g0d, keep.night=FALSE)
                         result <- as.data.frameY(prod)[c('Eac', 'Edc', 'Yf')] ##the results are yearly values
                         as.numeric(result)
                     })
                 )
                 ## para sacar los valores anuales. Yf es productividad.
                 Prod_Yf <- as.numeric(as.data.frame(subset(Prod, 'layer.3')))
                 return(Prod_Yf)
}

    
  ## List with the indices of blocks for each node
    iCluster <- splitIndices(bs$n, nodes)
    resCl <- mclapply(iCluster,
                      ## Each node receives an element of iCluster, a set of indices
                      function(icl){
                          resList <- lapply(icl, function(i){
                              ## An element of icl is the number of block to
                              ## be read by each node
                              vals <- getValues(data, bs$row[i], bs$nrows[i])
                              lat <- getValues(y, bs$row[i], bs$nrows[i])
                              vals <- cbind(lat, vals)
                              cat(i, ':', range(lat), '\n')
                              res0 <- try(apply(vals, MARGIN=1L, FUN=fooGef))
                              cat(i, ':', range(res0), '\n')
                              if (inherits(res0, 'try-error')) res0 <- NA
                              else  res0
                          })
                          do.call(c, resList)
                    }, mc.cores = nodes)                                                                                              
  ## The result of mclapply is a list with as many elements as nodes
  ## Each element of the list is a matrix with 3 columns (resCl0)
  ## corresponding to a block as defined by bs.
  resCl <- do.call(c, resCl)  
    
  out <- stack(data) 
  out <- setValues(out, resCl)
  if (filename!='') out <- writeRaster(out, filename=filename)
  out
}

prueba <- fooParallel(SISmm)
