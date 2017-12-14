#Kelly Kung
#12/15/17
#define rain/no rain

#setwd("/Users/kkung/Documents/GitHub/MA578_Project")
source("area_code.R")

#look to see how the aggregate precipitation is for one year
sumDay1<-c()
sumDay2<-c()
#for(i in 1:365){
#  sumDay1<-c(sumDay1, sum(weatherDat[,,i][areaGroup==1], na.rm=TRUE))
#  sumDay2<-c(sumDay2, sum(weatherDat[,,i][areaGroup==2], na.rm=TRUE))
#}

#plot(sumDay1, type= "l")
#table(sumDay1)
#plot(sumDay2, type = "l")


#median(sumDay1)
#sum(sumDay1<56.1)/length(sumDay1)

#define the dry seasons as this for now
wcutoff<-9.9#quantile(sumDay1, .3)
ecutoff<-23.2#quantile(sumDay2, .5)

#code up the definitions of rainvnorain

#set up a list to put each 
weatherDatE<-weatherDatW<-weatherDat
weatherDatW[areaGroup==2]<-NA
weatherDatE[areaGroup==1]<-NA
eastYear<-list()
westYear<-list()
for( i in 1:46) {
  westYear[[i]] <- weatherDatW[,,seq(365*(i-1)+1,365*i)]
  eastYear[[i]] <- weatherDatE[,,seq(365*(i-1)+1,365*i)]
}
dim(westYear[[1]])

sumGroup<-function(rain){
  rainFall<-apply(rain, 3, sum, na.rm=TRUE)
  return(rainFall)
}
sumWest<-lapply(westYear, sumGroup)
sumEast<-lapply(eastYear, sumGroup)

rainLevelW<-list()
rainLevelE<-list()
for(i in 1:46){
  rainLevelW[[i]]<-sumWest[[i]]>wcutoff
  rainLevelE[[i]]<-sumEast[[i]]>ecutoff
}
