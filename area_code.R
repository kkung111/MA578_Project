#Kelly Kung
#12/10/17

##########################################################################
# 
# setting up the model and variables for precipitation data
#
##########################################################################

setwd("/Users/kkung/Documents/MA 578/Project")
library(ncdf4)
weather<-nc_open("pnwrain.50km.daily.4994.nc")
#details about the dataset
print(weather)
weatherDat<-ncvar_get(weather, attributes(weather$var)$names[1])
weatherDat
dim(weatherDat) 

#create the location groups
areaGroup<-array(NA, dim = dim(weatherDat))
latlong<-t(weatherDat[,,1])
areaGroup[c(1:9),c(1:6),]<-1
areaGroup[c(1:8),c(7:9),]<-1
areaGroup[c(1:7),c(10:14),]<-1
areaGroup[c(1:6),c(15:17),]<-1
areaGroup[is.na(areaGroup)]<-2
t(areaGroup[,,1])
