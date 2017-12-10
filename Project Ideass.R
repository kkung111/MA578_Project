#Kelly Kung
#11/26/17
#just initial looking at MS 578 data

setwd("/Users/kkung/Documents/MA 578/Project")
library(ncdf4)
weather<-nc_open("pnwrain.50km.daily.4994.nc")
#details about the dataset
print(weather)
weatherDat<-ncvar_get(weather, attributes(weather$var)$names[1])
weatherDat
dim(weatherDat) #latitude (17), longitude (16), time (16801 in days since 1949)
                #each cross section shows the mm/day amount of rainfall 
                #32767 missing data points

#look at the distribution of the rainfall -- pretty left skewed
hist(weatherDat)
hist(weatherDat[,,1])
hist(weatherDat[,,2])
hist(weatherDat[3,,])

#idea for weather: 
#1) take a subset of data, like a couple of days before and then fit a model to predict the amount of rainfall
#2) use surrounding areas to predict the amount of rainfall in certain area
#3) impute missing data for certain days/area since there are many missing values
#4) compare rainfall in different areas, like a t-test

#--------------------------------------------------------------------------------------------------

players<-read.csv("players.txt")
colnames(players)
head(players)

regSeason<-read.csv("player_regular_season_career.txt")
colnames(regSeason)

coaches<-read.csv("coaches_season.txt")
head(coaches)
#ideas for bball:
#1) are coaches or top player more important for wins during that year? -- maybe separate models for coach and player
#2) predict winners using some stats
#3) can we identify who will be in allstar games base on players stats in the season--> lin model maybe?
#4) 