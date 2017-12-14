

input <- sample(c(0,1),365,replace = T, prob = c(.4,.6))

sample <- weatherDat[7,5,1:365]
sample <- list(rep(0,45))
for(i in seq(1,45)) {
  sample[[i]] <-  (weatherDat[5,7,seq(i*365-364,i*365)]>0)
}

CalculateBounds <- function(input) {
best_score <- -1
best_loc <- c(-1,-1)
for(i in seq(101,263)) {
  for(j in seq(i+1,304)) {
    score <- sum(input[seq(1,i)]) + sum(1-input[seq(i+1,j)]) + sum(input[seq(j+1,365)])
    if(score > best_score) {
      best_score <- score
      best_loc <- c(i,j)
    }
    if(score == best_score & (best_loc[2]- best_loc[1]< j-i)) {
      best_score<-score
      best_loc <- c(i,j)
    }
  }
}
return(best_loc)
}
 bounds <- lapply(sample,CalculateBounds)

