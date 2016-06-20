## First attempt using Parallel.

library(raster)
library(ncdf4)
library(rasterVis)
library(rgdal)
library(maps)
library(maptools)
library(parallel)
library(solaR)

## Calculate the annual productivity using parallel function

##########################################################################
## 1. DATA
#########################################################################

load('data/boundaries.Rdata')
load('data/linea.Rdata')
load('data/mascaraClustersSat.Rdata')
load('data/cmsafRasterdailyRadiation.Rdata')

## Creo un raster con valores medios mensuales para que sea más pequeño

SISS <- SISS*24
idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-30'), 'day')
SISS <- setZ(SISS, idx)

month <- function(x) as.numeric(format(x, '%m'))

SISmm <- zApply(SISS, by=month, fun='mean') ## medias mensuales de irradiancia

## paso a irradiacion Wh/m2 para poder emplearlo en el calculo de la GEf

########################################################################
## 2. ProdCGPV
#######################################################################

###############################################################
## mi prueba
#####################

## Productividad anual. Dependiendo del tipo de seguidor en la función fooProd dentro de fooParallel tendré que cambiar el modeTrk de prodCGPV

## Este es el código que hace en paralelo la función que se le expecifica en FUN.

fooParallel <- function(data, filename="", modeTrk = 'fixed', nodes=detectCores(), blocks=6,...){
    ## latitude values as a new raster
    y <- init(data, v='y')
    idx <- getZ(data)

    bs <- blockSize(data, minblocks=blocks*nodes)
    
    fooProd <- function(g0){
        n <- length(g0)
        lat <- g0[1]
        Prod <- prodGCPV(lat= lat,
                         dataRad= list(G0dm=g0[2:n]),
                         keep.night=FALSE, modeTrk = modeTrk)
        result <- as.data.frameY(Prod)[c('Yf')] ##the results are yearly values
        result <- as.numeric(result) ## para sacar los valores anuales. Yf es productividad.
        return(result)
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
                              res0 <- try(apply(vals, MARGIN=1L, FUN=fooProd))
                              cat(i, ':', range(res0), '\n')
                              if (inherits(res0, 'try-error')) res0 <- NA
                              else  res0
                          })
                          do.call(c, resList)
                      }, mc.cores = nodes)                                                                                              
    ## The result of mclapply is a list with as many elements as nodes
    ## Each element of the list is a matrix with 1 columns (resCl0)
    ## corresponding to a block as defined by bs.
    resCl <- do.call(c, resCl)  
    
    out <- raster(data) 
    out <- setValues(out, resCl)
    if (filename!='') out <- writeRaster(out, filename=filename)
    out
}

prueba <- fooParallel(SISmm)

#############################################################
## 20 years productivity
############################################################

## Quiero calcular la productividad anual para cada uno de los años del periodo.

library(zoo)
## Medias mensuales por cada año
SISmm240 <- zApply(SISS, by=as.yearmon, fun='mean') ## SISS ya esta multiplicado por 24.
## Índice temporal (año-mes)
idxSISmm <- getZ(SISmm240)
## Trocea el índice temporal en grupos definidos por el año. El
## resultado es una lista con tantos elementos como años.
idLayers <- split(1:nlayers(SISmm240),
                  year(idxSISmm)) 
## Ahora recorremos esta lista con lapply, aplicando a cada grupo
## (año) la función fooParallel, eligiendo previamente las capas que
## corresponden a ese año
yProdFixed <- lapply(idLayers, FUN = function(idx)
{
    SISmm <- subset(SISmm240, idx)
    fooParallel(SISmm, modeTrk = 'fixed')
})
## El resultado es una lista con 20 rasters, con el resultado del
## cálculo de productividad para cada año del periodo. En fooParallell
## cambiamos modeTrk para obtener los resultados para cada tipo
## de seguidor.
yProdFixed <- stack(yProdFixed)
writeRaster(yProdFixed, filename='YearlyProductivity20_fixed')

## Seguimiento eje horizontal
yProdHoriz <- lapply(idLayers, FUN = function(idx)
{
    SISmm <- subset(SISmm240, idx)
    fooParallel(SISmm, modeTrk = 'horiz')
})
yProdHoriz <- stack(yProdHoriz)
writeRaster(yProdHoriz, filename='YearlyProductivity20_horiz')

## Seguimiento Doble Eje
yProdTwo <- lapply(idLayers, FUN = function(idx)
{
    SISmm <- subset(SISmm240, idx)
    fooParallel(SISmm, modeTrk = 'two')
})
yProdTwo <- stack(yProdTwo)
writeRaster(yProdTwo, filename='YearlyProductivity20_two')
