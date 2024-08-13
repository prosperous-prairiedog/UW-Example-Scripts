library(sp)
library(rgdal)
library(sf)
library(tmap)
library(spdep)
library(raster)
library(fasterize)
library(tidyverse)
library(snow)
library(RColorBrewer)
library(rasterVis)    
library(ggplot2)      
library(colorspace)
library(dplyr)
library(tidyr)
library(readr)
library(maptools)
library(MMWRweek)
library(plotrix)
library(tidygeocoder)

### load Tarrant county test data
setwd("C:/Users/josmcmil/OneDrive - Texas Tech University/JRM_Tarrant County Aedes/Tarrant County")
Dat <- as_tibble(read.csv("Tarrant test data.csv"))

#convert to date and epi week
Dat2 <- Dat %>%
  mutate(Date.2 = as.Date(Collection_Date,"%m/%d/%Y")) %>%
  mutate(Week = MMWRweek(Date.2)[[2]],
         Year = MMWRweek(Date.2)[[1]])

syn.fun <- function(YR, TRP) { 
  Dat.3 <- Dat2 %>% filter(Year==YR, Trap_Method==TRP)

###collection effort per address
Effort <- tapply(Dat.3$Quantity, list(Dat.3$Date.2, Dat.3$Address, Dat.3$Week), sum)
Effort[is.na(Effort)] <- 0
dim <- dim(Effort)

Trp.ngt.wk <- matrix(nrow=dim[2], ncol=dim[3])
for (i in 1:dim[2]) {
  for (j in 1:dim[3]) {
    Trp.ngt.wk[i,j] = sum(Effort[,i,j]>0)
  }	}

######      correcting for trapping effort
Sm <- tapply(Dat.3$Quantity, list(Dat.3$Address, Dat.3$Week), sum, na.rm=T)
Sm[is.na(Sm)] <- 0

Sm.c <- matrix(nrow=dim[2], ncol=dim[3])
for (i in 1:dim[2]) {
  for (j in 1:dim[3]) {
    Sm.c[i,j] <- Sm[i,j]/Trp.ngt.wk[i,j]
    }
  }
Sm.c[is.na(Sm.c)] <- 0

### as a dataframe
library(reshape)
Sm.c.dat <- as.data.frame(melt(Sm.c))
colnames(Sm.c.dat) <- c("Site","Week","Cr.Collection")

CRD <- Dat.3 %>% group_by(Address) %>%
  summarise(lat = mean(Lat), lng = mean(Long))

### synchrony 
library(synchrony)
Dat <- as_tibble(Sm.c.dat) %>%
  mutate(Lat = rep(c(CRD$lat), dim[3]), Long = rep(c(CRD$lng), dim[3])) %>%
  na.omit() 
Dat.2 <- as.data.frame(subset(Dat,select=c("Lat","Long","Week","Cr.Collection")))
Dat.2w <- reshape(data=Dat.2,idvar=c("Lat","Long"),timevar=c("Week"),
                  direction= "wide")
Crrlgrm <-vario(data=Dat.2w,
                type= "kendall",extent=0.9,nrands=1000,
                is.centered=TRUE,quiet=TRUE)
par(mar=c(4.5, 4.5, 2, 1))
p <- plot(Crrlgrm, rug=TRUE, ci=TRUE, ylim=c(-.1, .1), 
          xlab="Distance (km)",main=paste(YR,TRP,"Aedes aegypti synchrony - Tarrant County", sep=" "))
return(p)
}

par(mfrow=c(6,1))
syn.fun(YR=2017, TRP="GRAVID")
syn.fun(YR=2018, TRP="GRAVID")
syn.fun(YR=2019, TRP="GRAVID")
syn.fun(YR=2020, TRP="GRAVID")
syn.fun(YR=2021, TRP="GRAVID")
syn.fun(YR=2022, TRP="GRAVID")

#checking which are core sites
tab <- tapply(Dat2$Quantity, list(Dat2$Address, Dat2$Year), sum)
tab2 <- cbind(tab, rowSums(tab>0, na.rm=T))
sum(tab2[,5]==1)

### 130 address sampled only 1 year
### 334 address sampled 2 years
### 16 address sampled 3 years
### 11 sampled 4 years
#### 2016
###collection effort per address
names(Dat.16)
Effort.16 <- tapply(Dat.16$Quantity, list(Dat.16$Date.2, Dat.16$Address, Dat.16$Week), sum)
Effort.16[is.na(Effort.16)] <- 0

Trp.ngt.wk.16 <- matrix(nrow=181, ncol=33)
for (i in 1:181) {
  for (j in 1:33) {
    Trp.ngt.wk.16[i,j] = sum(Effort.16[,i,j]>0)
  }	}
rownames(Trp.ngt.wk.16) <- levels(as.factor(Dat.16$Address))
colnames(Trp.ngt.wk.16) <- levels(as.factor(Dat.16$Week))

######      correcting for trapping effort
Sm.16 <- tapply(Dat.16$Quantity, list(Dat.16$Address, Dat.16$Week), sum, na.rm=T)
Sm.16[is.na(Sm.16)] <- 0

Sm.16.c <- matrix(nrow=181, ncol=33)
for (i in 1:181) {
  for (j in 1:33) {
    Sm.16.c[i,j] <- Sm.16[i,j]/Trp.ngt.wk.16[i,j]
  }
}

dimnames(Sm.16.c)[[1]] <- levels(as.factor(Dat.16$Address))
dimnames(Sm.16.c)[[2]] <- levels(as.factor(Dat.16$Week))

### as a dataframe
library(reshape)
Sm.16.c.dat <- as.data.frame(melt(Sm.16.c))
colnames(Sm.16.c.dat) <- c("Address","Week","Cr.Collection")

Dat.16.v2 <- Dat.16 %>%
  group_by(Address) %>%
  summarise(Lat.2 = mean(Lat, na.rm=T), Long.2 = mean(Long, na.rm=T))


### synchrony 2016
library(synchrony)
Dat <- as_tibble(Sm.16.c.dat) %>%
  mutate(Lat = rep(c(Dat.16.v2$Lat.2), 33), Long = rep(c(Dat.16.v2$Long.2), 33)) %>%
  na.omit() 
Dat.2 <- as.data.frame(subset(Dat,select=c("Lat","Long","Week","Cr.Collection")))
Dat.2w <- reshape(data=Dat.2,idvar=c("Lat","Long"),timevar=c("Week"),
                  direction= "wide")

### tough to figure out if there is spatial autocorrelation from these plots... too may NAs
Crrlgrm <-vario(data=Dat.2w,
                type= "kendall",extent=0.9,nrands=1000,
                is.centered=TRUE,quiet=TRUE)
par(mar=c(4.5, 4.5, 2, 1))
plot(Crrlgrm, rug=TRUE, ci=TRUE, ylim=c(-.3, .3), 
     xlab="Distance (km)",main="Aedes aegypti (females only) collection synchrony - 2016 Tarrant County",
     cex.lab=2, cex.axis=2)

#### 2017
###collection effort per ZIPCODE
Effort.17 <- tapply(Dat.17$Quantity, list(Dat.17$Date.2, Dat.17$Address, Dat.17$Week), sum)
Effort.17[is.na(Effort.17)] <- 0
dim(Effort.17)

Trp.ngt.wk.17 <- matrix(nrow=235, ncol=34)
for (i in 1:235) {
  for (j in 1:34) {
    Trp.ngt.wk.17[i,j] = sum(Effort.17[,i,j]>0)
  }	}
rownames(Trp.ngt.wk.17) <- levels(as.factor(Dat.17$Address))
colnames(Trp.ngt.wk.17) <- levels(as.factor(Dat.17$Week))

######      correcting for trapping effort
Sm.17 <- tapply(Dat.17$Quantity, list(Dat.17$Address, Dat.17$Week), sum, na.rm=T)
Sm.17[is.na(Sm.17)] <- 0

Sm.17.c <- matrix(nrow=235, ncol=34)
for (i in 1:235) {
  for (j in 1:34) {
    Sm.17.c[i,j] <- Sm.17[i,j]/Trp.ngt.wk.17[i,j]
  }
}
Sm.17.c[is.na(Sm.17.c)] <- 0
dimnames(Sm.17.c)[[1]] <- levels(as.factor(Dat.17$Address))
dimnames(Sm.17.c)[[2]] <- levels(as.factor(Dat.17$Week))

### as a dataframe
library(reshape)
Sm.17.c.dat <- as.data.frame(melt(Sm.17.c))
colnames(Sm.17.c.dat) <- c("Address","Week","Cr.Collection")

Dat.17.v2 <- Dat.17 %>%
  group_by(Address) %>%
  summarise(Lat.2 = mean(Lat, na.rm=T), Long.2 = mean(Long, na.rm=T))
dim(Dat.17.v2)

### synchrony 
library(synchrony)
Dat <- as_tibble(Sm.17.c.dat) %>%
  mutate(Lat = rep(c(Dat.17.v2$Lat.2), 34), Long = rep(c(Dat.17.v2$Long.2), 34)) %>%
  na.omit() 
Dat.2 <- as.data.frame(subset(Dat,select=c("Lat","Long","Week","Cr.Collection")))
Dat.2w <- reshape(data=Dat.2,idvar=c("Lat","Long"),timevar=c("Week"),
                  direction= "wide")
Dat.2w[is.na(Dat.2w)] <- 0
Crrlgrm <-vario(data=Dat.2w,
                type= "kendall",extent=0.9,nrands=1000,
                is.centered=TRUE,quiet=TRUE)
par(mar=c(4.5, 4.5, 2, 1))
plot(Crrlgrm, rug=TRUE, ci=TRUE, ylim=c(-.075, .075), xlab="Distance (km)",main="2017 Tarrant County - based on ZIP")
#### PROB NEED TO USE ANNUAL VALUES???
