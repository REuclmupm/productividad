## Cálculo de la productividad anual media.

library(raster)
library(rasterVis)
library(maps)
library(ncdf4)
library(maptools)
library(mapdata)
library(rgdal)
library(parallel)
library(solaR)

###############################################################################
## 1. DATOS
##############################################################################
load('data/boundaries.Rdata')
load('data/linea.Rdata')
load('data/mascarasClusters.Rdata')
load('data/cmsafRasterdailyRadiation.Rdata') 

## Creo un raster con valores medios mensuales para que sea más pequeño

idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-30'), 'day')
SISS <- setZ(SISS, idx)

month <- function(x) as.numeric(format(x, '%m'))

SISmm <- zApply(SISS, by=month, fun='mean') ## medias mensuales de irradiancia

## paso a irradiacion Wh/m2 para poder emplearlo en el calculo de la GEf

SISmm <- SISmm*24

#############################################################################
## 2. GEF
############################################################################

## A continuación calculo la irradiancia efectiva para el plano del generador sin utilizar la paralelización. 

## añado al raster una capa con la latitud

latLayer <- init(SISmm, v='y') ## extraigo la latitud

foo <- function(x, ...)
{
    gef <- calcGef(lat = x[1], dataRad = list(G0dm = x[2:13]))
    result <- as.data.frameY(gef)[c('Gefd', 'Befd', 'Defd')] ##the results are yearly values
    as.numeric(result)
}


gefS <- calc(stack(latLayer, SISmm), foo, overwrite=TRUE)
names(gefS) <- c('Gefd', 'Befd', 'Defd')##Three layers

levelplot(subset(gefS, 'Gefd')) +
    layer(sp.lines(linea))

#############################################################################
## 3. PRODUCTIVIDAD ANUAL PARA CADA CELDA. (por tipo de seguidor)
#############################################################################

## Calculo de la productividad si en cada celda hubiera un generador

## Simplifico. Produtividad anual:

PR <- 0.74

## GEF que sale de la funcion calcGef tiene unidades de kWh/m2 y es la irradiancia efectiva anual.

gefS <- subset(gefS, 'Gefd')
productividad <- PR*(gefS)

levelplot(mask(productividad, boundaries_sp),
          margin=FALSE,
          main='PRODUCTIVIDAD ANUAL kWh/kWp') +
    layer(sp.lines(linea))

## Aplico la mascara de clusters para la productividad y calculo la "media de productividad" por cluster para comparar". 

prod_clusters <- function(x, y){
    ## x es raster con datos de radiación efectiva anual, y es el cluster que quiero analizar.
    prod <- mask(x,mascarasClusters[[y]])
    prod_media <- cellStats(prod, stat='mean')
    sol <- list(prod, prod_media)
    return(sol)
}

prodBycluster <- list(1:6)
for (i in 1:6) prodBycluster[[i]] <-prod_clusters(productividad, i)

mediaBycluster <- c(1:6)
for (i in 1:6) mediaBycluster[i] <- prodBycluster[[i]][[2]]

#######################################################################
## Calculo de productividad con seguidores a DOBLE EJE:
#######################################################################


foo_two<- function(x, ...){
gef <- calcGef(lat=x[1], dataRad=list(G0dm=x[2:13]), modeTrk='two')
result <- as.data.frameY(gef)[c('Gefd', 'Befd', 'Defd')] ##the results are yearly values
as.numeric(result)
}

gefS_two<- calc(stack(latLayer, SISmm), foo_two, overwrite=TRUE)
names(gefS_two)=c('Gefd', 'Befd', 'Defd')##Three layers

gef_two<-subset(gefS_two, 'Gefd')

productividad_two <- PR*(gef_two)

prodBycluster_two<- list(1:6)
for (i in 1:6) prodBycluster_two[[i]] <-prod_clusters(productividad_two, i)

mediaBycluster_two<- c(1:6)
for (i in 1:6) mediaBycluster_two[i] <- prodBycluster_two[[i]][[2]]


#############################################################################
## Calculo de productividad con seguidores N-S:
#############################################################################


foo_ns<- function(x, ...){
gef <- calcGef(lat=x[1], dataRad=list(G0dm=x[2:13]), modeTrk='horiz')
result <- as.data.frameY(gef)[c('Gefd', 'Befd', 'Defd')] ##the results are yearly values
as.numeric(result)
}

gefS_ns<- calc(stack(latLayer, SISmm), foo_ns, overwrite=TRUE) 
names(gefS_ns)=c('Gefd', 'Befd', 'Defd')##Three layers

gef_ns<-subset(gefS_ns, 'Gefd')
 
productividad_ns <- PR*(gef_ns)

prodBycluster_ns<- list(1:6)
for (i in 1:6) prodBycluster_ns[[i]] <-prod_clusters(productividad_ns, i)

mediaBycluster_ns<- c(1:6)
for (i in 1:6) mediaBycluster_ns[i] <- prodBycluster_ns[[i]][[2]]

##########################################################################
## USO prodGCPV
#########################################################################

## El calculo de la radiación efectiva se hace con las funciones foo que hemos definido antes.

## 1. seguidor estático

foo <- function(x, ...){
prod <- prodGCPV(lat=x[1], dataRad=list(G0dm=x[2:13]), keep.night=FALSE)
result <- as.data.frameY(prod)[c('Eac', 'Edc', 'Yf')] ##the results are yearly values
as.numeric(result)
} ## para sacar los valores anuales. Yf es productividad.


Prod  <- calc(stack(latLayer, SISmm), foo, overwrite=TRUE)
Prod_Yf <- subset(Prod, 'layer.3')

prodBycluster <- list(1:6)
for (i in 1:6) prodBycluster[[i]] <-prod_clusters(Prod_Yf,i)

mediaBycluster<- c(1:6)
for (i in 1:6) mediaBycluster[i] <- prodBycluster[[i]][[2]]


## 2. Seguidor N-S

foo_ns<- function(x, ...){
prod <- prodGCPV(lat=x[1], dataRad=list(G0dm=x[2:13]), modeTrk='horiz')
result <- as.data.frameY(prod)[c('Eac', 'Edc', 'Yf')] ##the results are yearly values
as.numeric(result)
}

Prod_ns <- calc(stack(latLayer, SISmm), foo_ns, overwrite=TRUE)
Prod_ns_Yf <- subset(Prod_ns, 'layer.3')

prodBycluster_ns<- list(1:6)
for (i in 1:6) prodBycluster_ns[[i]] <-prod_clusters(Prod_ns_Yf,i)

mediaBycluster_ns<- c(1:6)
for (i in 1:6) mediaBycluster_ns[i] <- prodBycluster_ns[[i]][[2]]


## 3. Seguidor Doble eje

foo_two<- function(x, ...){
prod <- prodGCPV(lat=x[1], dataRad=list(G0dm=x[2:13]), modeTrk='two')
result <- as.data.frameY(prod)[c('Eac', 'Edc', 'Yf')] ##the results are yearly values
as.numeric(result)
}

Prod_two <- calc(stack(latLayer, SISmm), foo_two, overwrite=TRUE)
Prod_two_Yf <- subset(Prod_two, 'layer.3')

prodBycluster_two<- list(1:6)
for (i in 1:6) prodBycluster_two[[i]] <-prod_clusters(Prod_two_Yf,i)

mediaBycluster_two<- c(1:6)
for (i in 1:6) mediaBycluster_two[i] <- prodBycluster_two[[i]][[2]]

#####################################################################################
## GRAFICAS
#####################################################################################

## Productividad anual:

anual <- stack(Prod_Yf, Prod_ns_Yf, Prod_two_Yf)
names(anual) <- c("estático", "norte-sur", "doble-eje")

levelplot(mask(anual, boundaries_sp), margin=FALSE, main='PRODUCTIVIDAD ANUAL kWh/kWp')+layer(sp.lines(linea))

dev.copy(png, file='productividad_anual_sat.png')
dev.off()

## matriz de medias de productividad por cluster:

byClusterMatrix <- cbind(mediaBycluster, mediaBycluster_ns, mediaBycluster_two)
rownames(byClusterMatrix) <- c("C1", "C2", "C3", "C4", "C5", "C6")

library(lattice)
barchart(byClusterMatrix, stack=F, xlab='kWh/Kwp', main='Productividad anual media por cluster PROMES')
dev.copy(png, file='productividad_anual_clusters_sat.png')
dev.off()
