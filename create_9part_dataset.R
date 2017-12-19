# Split data into 9 parts

#grid2

#num of sections
area_count <- 9
years <- 46

tf_data2 <- list(rep(0,9))
for( i in seq(1,area_count)) {
  tf_data2[[i]] <- list(rep(0,years))
  for( j in seq(1, years)){
    #locs <- array(rep(grid2 == i,365),dim = c(16,17,365))
    locs <- (grid2 == i)
   this_year <- tf_weather_data[,, seq(j*365-364,j*365)]
  tf_data2[[i]][[j]] <- apply(this_year,3,function(x){sum(x&locs,na.rm=T)>(.5*sum(locs)) })
  }
}

MrainLevel <- list(rep(0,9))
for( i in seq(1,9)) {
MrainLevel[[i]] <- matrix(unlist(lapply(tf_data[[i]],CalculateBounds)),ncol=2,byrow=T)
write.matrix(MrainLevel[[i]],file=paste("Area",i,"MinMaxDataV2.tsv"),sep="\t")
}

rain_levels_area <- MrainLevel
