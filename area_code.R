#Kelly Kung
#12/10/17

##########################################################################
# 
# setting up the model and variables for precipitation data
# area, season, years, and month arrays
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
areaGroup[c(1:9),c(1:6),]<-1
areaGroup[c(1:8),c(7:9),]<-1
areaGroup[c(1:7),c(10:14),]<-1
areaGroup[c(1:6),c(15:17),]<-1
areaGroup[is.na(areaGroup)]<-2
t(areaGroup[,,1])

#code the years, months, and cities
#code up the number of days for each month and season
daysinMonth<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
spr<-(31 - 19) + 30 + 31 + 20
sum<-(30 - 20) + 31 + 31 + 21
aut<-(30 - 21) + 31 + 30 + 20
win1<-31 + 28 + 19
win2<-(31 - 20)

#set up the arrays
years<-array(NA, dim = dim(weatherDat))
seasons<-array(NA, dim = dim(weatherDat)) 
months<-array(NA, dim = dim(weatherDat))

#start the indices
year = 1
count = 1
diff = dim(weatherDat)[3] #this keeps track of our while loop

while(diff>0){
  #a sequence of number to keep track of the days in a year
  elts<-count:(count + 364)
  #note: 272 = 16*17
  monthDays<-c(rep(1, 31*272), rep(2, 28*272), rep(3, 31*272), rep(4, 30*272), 
               rep(5, 31*272), rep(6, 30*272), rep(7, 31*272), rep(8, 31*272),
               rep(9, 30*272), rep(10, 31*272), rep(11, 30*272), rep(12, 31*272))
  seasonDays<-c(rep("win", win1*272), rep("spr", spr*272), rep("sum", sum*272), 
                rep("aut", aut*272), rep("win", win2*272))
  #16801 is the length of all the days so we just cut off the fill ins if it's too long
  if(TRUE %in% (elts>16801)){
    index<-which(elts == 16801)
    elts<-elts[1:index]
    monthDays<-monthDays[1:index]
    seasonDays<-seasonDays[1:index]
  }
  #fill in the arrays 
  years[,,elts]<-year
  months[,,elts]<-monthDays
  seasons[,,elts]<-seasonDays
  
  #increment for the while loop
  count = count + 365
  year = year + 1
  diff = dim(weatherDat)[3] - count
}

#years--> correlation, time series component A(1), banded--> after 2 years no correlation because correlation
#decreases as years go on; enough data --> choose some priors that are more influential
#model selection--> g prior, elner's 