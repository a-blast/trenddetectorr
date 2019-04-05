# Austin Armstrong 2018
# Atomata grow data statistics

library(tidyverse)
library(jsonlite)
library(data.table)
library(magrittr)
library(broom)
library(wmtsa)
library(tuneR)

getModuleNames <- function(filename) {
  raw_json <- read_json(filename)
  moduleNames <- raw_json$reads %>% names
  return (moduleNames)
}

cleanRawJSONResponse <-function(filename,moduleName) {
  raw_json <- read_json(filename)
  moduleRecords <- raw_json$reads[[moduleName]]
  sensorNames <- raw_json$reads[[moduleName]] %>% names
  raw_df <- map(sensorNames,function(name) moduleRecords %>% pluck(name) %>% as.double %>% tibble) %>% reduce(bind_cols)
  colnames(raw_df) <- sensorNames
  return (raw_df)
}

# define the stream to evaluate change points

IncrementalChangePointDetector <- function(time.series, max.exponential.basis,min.num.points,delta,smooth.data=TRUE){
  # make sure the time.series is long enough to play ball
  if(dim(time.series)[1] < (2*min.num.points)) {
    return()
  }
  
  if(smooth.data){
    smoothed.fit <- with(time.series, smooth.spline(time, dat))
    time.series <- tibble(dat=smoothed.fit$y,time=smoothed.fit$x)
  }
  
  FindLikelihoodCriteria <- function(time.series.subset,delta){
    time.series.subset$time <- time.series.subset$time - min(time.series.subset$time)
    se <- lm(dat~poly(time,max.exponential.basis), data=time.series.subset) %>% 
      residuals() %>% 
      map(function(x) x^2) %>% 
      reduce(`+`)
    return(se)
  }
  
  FindSplits <- function(){
    candidate.optimal.likelihood.criteria <- Inf
    points.to.check.seq <- min.num.points:(dim(time.series)[1]-min.num.points-1)
    print(points.to.check.seq)
    for(time.point in points.to.check.seq) {
      left.hand.split <- FindLikelihoodCriteria(time.series %>% slice(1:time.point))
      right.hand.split <-  FindLikelihoodCriteria(time.series %>% slice((time.point+1):(dim(time.series)[1])))
      likelihood.criteria <- left.hand.split + right.hand.split
      if (likelihood.criteria < candidate.optimal.likelihood.criteria){
        candidate.optimal.likelihood.criteria <- likelihood.criteria
        candidate.split <- time.point
      }
    }
    return(list(optimal.likelihood.criteria=candidate.optimal.likelihood.criteria,
                split=candidate.split))
  }
  
  split.likelihood.criteria <- FindSplits()
  no.split.likelihood.criteria <- FindLikelihoodCriteria(time.series)
  
  print(c("DIFF:",(no.split.likelihood.criteria - split.likelihood.criteria$optimal.likelihood.criteria)/no.split.likelihood.criteria, delta))
  if((no.split.likelihood.criteria - split.likelihood.criteria$optimal.likelihood.criteria)/no.split.likelihood.criteria > delta){
    print(split.likelihood.criteria)
    return(split.likelihood.criteria$split)
  } else {
    return(0)
  }
}

incremental <- aisdank_df %>% 
  filter(name=="V1C0001") %>% 
  unnest() %>% 
  filter(time >= lastWeek*1.00027) %>%
  filter(time<=max(time-1000000)) %>%
  select(temp, time) %>%
  set_colnames(c("dat","time"))

change.points <- c()
split.point<-0

for(point in 15:length(incremental$dat)){
  if(is.null(change.points)){
    split.point <- IncrementalChangePointDetector(incremental %>% slice(1:point),1,6,0)
  }else{
    split.point <- IncrementalChangePointDetector(incremental %>% slice((change.points[length(change.points)]):point),1,6,0.3)
    
    if(split.point != 0){
      split.point <- change.points[length(change.points)] + split.point
    }
  }
  print(c("CHACKsplit.point",split.point))
  if (split.point != 0){
    print("NOT NULL")
    change.points <- c(change.points,split.point)
  }
}

frame.to.plot %>%
  #pull(hum) %>%
  #plot(type="l")
  mutate(splitTime = ifelse(rowname %in% change.points,time,NA)) %>%
  #mutate(fit.point = ifelse(rowname %in% fit.point.bounds$point,fit.point.bounds$fit.hum[]))
  #mutate(rownameNum = as.numeric(rowname)) %>%
  ggplot +
  geom_line(aes(x=time, y=temp)) +
  #geom_line(aes(x=time, y=c(diff(time),NA)))
  #geom_point(aes(x=time, y=fit.hum)) +
  geom_vline(aes(xintercept=splitTime))


plot(dat ~ time, data=incremental, type="l")
smoothed <- with(incremental, smooth.spline(time, dat))
smoothed$y %>% plot(type="l")
smoothed$x
incremental %>% smooth.spline() %>% plot()

detectChangePointsBatch_Stream <- function(timeSeries,maxExponentialBasis,maxNumPoints){
  
  newChangePoint <- c()
  changePoints <- c()
  candidates <- c()
  newSegmentSplits <- c()
  
  ### initialize the stream ###
  streamInit <- function() {
    # initialize a df that contains an entry for the split point and the new likelihood criteria
    # define the entries of the new change point
    newChangePoint <- findAndUpdateCandidates(c(1,length(timeSeries$time),"left"), set.change.point = TRUE)
    print(newChangePoint)
    #candidates <<- data.frame(split=c(1,length(timeSeries)),optimalLikelihoodCriteria=c(NA,NA))
    newSegmentSplits <<- list(c(1,newChangePoint$split,"left"),c(newChangePoint$split+1,length(timeSeries$time),"right"))
  }
  
  findAndUpdateCandidates <- function(ts.bounds, set.change.point=FALSE){
    candidateOptimalLikelihoodCriteria <- Inf
    
    print(c("RANGE CHECK!",ts.bounds[1:2])) 
    
    lowBound <- as.numeric(ts.bounds[1])
    highBound <- as.numeric(ts.bounds[2])
    
    if((highBound-lowBound)/maxNumPoints <= 2){
      print("TS BOUNDS TOOSMALL")
      return()
    }
   
    #print(range(maxNumPoints:(length(tS)-maxNumPoints-1)))
    points.to.check.seq <- (lowBound+maxNumPoints):(highBound-maxNumPoints-1)
    historical <- c()
    for (time.point in points.to.check.seq) {
      leftHandSplit <- findLikelihoodCriteria(timeSeries %>% slice(lowBound:time.point),highBound)
      rightHandSplit <- findLikelihoodCriteria(timeSeries %>% slice((time.point+1):highBound),highBound)
      likelyhood_criteria <-  leftHandSplit + rightHandSplit
      historical <- c(historical,likelyhood_criteria)
      
      if (likelyhood_criteria < candidateOptimalLikelihoodCriteria) {
        print(c("NEW_LIKELIHOOD: ", likelyhood_criteria, "Point:", time.point))
        candidateSplit <- time.point
        candidateOptimalLikelihoodCriteria <- likelyhood_criteria
        optimalLeftHandSplit <- leftHandSplit
        optimalRightHandSplit <- rightHandSplit
      }
      
    
    }
    
    
    
    print(c("TESTTHIS!!!!!!",time.point+1,highBound))
    
    if(candidateSplit < maxNumPoints){
      print("TS TOOSMALL")
      return()
    }
    
    print(c("CANDIDATE SPLIT",candidateSplit))
    
    if (ts.bounds[3] == "left"){
      parent.point <- highBound
      if (is.null(candidates)){
        parent.loss <- findLikelihoodCriteria(timeSeries %>% slice(lowBound:highBound),0)
      } else {
        parent.loss <- candidates %>% 
          filter(split == parent.point) %>%
          pull(left.loss)
      }
    } else {
      parent.point <- lowBound-1
      if (is.null(candidates)){
        parent.loss <- findLikelihoodCriteria(timeSeries %>% slice(lowBound:highBound),0)
      } else {
        parent.loss <- candidates %>% 
          filter(split == parent.point) %>%
          pull(right.loss)
        }
    }
    
    candidates <<- candidates %>% bind_rows(data.frame(split=candidateSplit,
                                                       split.time=timeSeries$time[candidateSplit],
                                                       left.loss=optimalLeftHandSplit,
                                                       right.loss=optimalRightHandSplit,
                                                       parent=parent.point,
                                                       parent.side=ts.bounds[3],
                                                       parent.side.loss=parent.loss,
                                                       accepted.change.point=set.change.point))
    
    return(candidates)
    
  }
  
  findLikelihoodCriteria <- function(timeSeriesSubset,max){
    timeSeriesSubset$time <- timeSeriesSubset$time - min(timeSeriesSubset$time)
    se <- lm(dat~poly(time,1), data=timeSeriesSubset) %>% residuals() %>% map(function(x) x^2) %>% reduce(`+`)
    #se <- lm(dat~time, data=timeSeriesSubset) %>% residuals() %>% map(function(x) x^2) %>% reduce(`+`)
    return(se)
  }
 
  GetNeighbors <- function(point) {
    print(c("IN GETNEIGHBORS:, CHANGE POINTS",candidates))
    
    
    
    neighbor.low <- candidates %>%
      filter(accepted.change.point) %>%
      filter(split < point) %>%
      pull(split) %>%
      max()
    
    neighbor.high <- candidates %>%
      filter(accepted.change.point) %>%
      filter(split > point) %>%
      pull(split) %>%
      min()
    
    if(!is.finite(neighbor.low)) {
      neighbor.low <- 1
    }
    if(!is.finite(neighbor.high)) {
      neighbor.high <- length(timeSeries$time)
    }
    return(list(low=neighbor.low,high=neighbor.high))
  }
   
  setNewSegmentSplits <- function(optimalPointSplit) {
    #get the neighbor points around the new optimal point
    
    neighbor.points <- GetNeighbors(optimalPointSplit)
    
    print(c("NEIGHBORS:",neighbor.points))
    lowSegment <- c(neighbor.points$low,optimalPointSplit,"left")
    highSegment <- c(optimalPointSplit+1,neighbor.points$high,"right")

    print(c("HighSegmentBounds:",optimalPointSplit+1,length(timeSeries)))    
    
    print("CHANGEPOINTS")
    print(changePoints)
    print("NEIGHBORPOINTS")
    print(c(neighbor.points))
    print(c("NEWSEGMENTS",lowSegment,highSegment))
    # Define the new relevant time series sections to explore
    newSegmentSplits <<- list(lowSegment, highSegment)
  }
  
  nextChangePoint <- function() {
    print(map(newSegmentSplits,print))
    map(newSegmentSplits,findAndUpdateCandidates)
    # pick the change point with the lowest loss
    print(c("C_PRE:",candidates))
    
    optimalPoint <- candidates %>% 
      filter(!accepted.change.point) %>%
      mutate(relative.scores = (left.loss+right.loss)-parent.side.loss ) %>%
      filter(relative.scores == min(relative.scores))
    
    optimalPoint$accepted.change.point[1] <- TRUE
    
    candidates <<- candidates %>%
      filter(split != optimalPoint$split[1]) %>%
      bind_rows(optimalPoint)
    # remove it from the candidates
    print(c("C_POST:",candidates))
    # add it to the chage points
    # update the new segments of the time series
    setNewSegmentSplits(optimalPoint$split)
    
    return(list(points=candidates,funct=nextChangePoint))
    
  }
  
  streamInit()
  #print(map(newSegmentSplits,length))
  
  return(nextChangePoint())
  
}



stream <- aisdank_df %>% 
  filter(name=="V1C0001") %>% 
  unnest() %>% 
  filter(time >= lastWeek*1.00027) %>%
  filter(time<=max(time-1000000)) %>%
  select(temp, time) %>%
  set_colnames(c("dat","time"))
  #Rev() %>%
  detectChangePointsBatch_Stream(1,6) %>%
  eventLoop()

eventLoop <- function(stream){
  #print(stream$points %>% filter(!is.na(relative.scores)) %>% pull(relative.scores) %>% max())
  while(
    (stream$points %>% filter(!is.na(relative.scores)) %>% pull(relative.scores) %>% max()) <= -1) {
    #print(as.tibble(stream$points))
    stream <- stream$funct()
  }
  return(stream)
}

subseq <- aisdank_df %>% 
  filter(name=="V1C0001") %>% 
  unnest() %>% 
  filter(time >= lastWeek*1.00027) %>%
  filter(time<=max(time-1000000)) %>%
  select(temp, time) %>%
  set_colnames(c("dat","time")) %>%
  slice(1:169)

findLikelihoodCriteria_debug <- function(timeSeriesSubset){
  #ts_df <- data.frame(series=timeSeriesSubset, time=1:length(timeSeriesSubset))
  se <- lm(dat~poly(time,3), data=timeSeriesSubset) %>% residuals() %>% map(function(x) x^2) %>% reduce(`+`)
  timeSeriesSubset$time <- timeSeriesSubset$time - min(timeSeriesSubset$time)
  #se <- lm(dat~time, data=timeSeriesSubset) %>% residuals() %>% map(function(x) x^2) %>% reduce(`+`)
  return(se)
}

ans <- c()

for(point in 6:162) {
  leftHandSplit_debug <- findLikelihoodCriteria_debug(subseq %>% slice(1:point))
  rightHandSplit_debug <- findLikelihoodCriteria_debug(subseq %>% slice(point+1:169))
  likelyhood_criteria <-  leftHandSplit + rightHandSplit
  ans <- c(ans,likelyhood_criteria)
}

without.split <- findLikelihoodCriteria(subseq)

subseq %>% pull(dat) %>% plot(type="l")
plot(ans,type="l")
abline(h=without.split)


ts <- stream <- aisdank_df %>% 
  filter(name=="V1C0001") %>% 
  unnest() %>% 
  filter(time >= lastWeek*1.00027) %>%
  filter(time<=max(time-1000000)) %>%
  pull(hum)

stream <- stream$funct()

accepted <- stream$points %>% filter(accepted.change.point)
points.to.map <- c(0,sort(accepted$split),length(ts))
points.to.map <- mapply(c, points.to.map[-length(points.to.map)], points.to.map[-1], SIMPLIFY = FALSE)

GetLMCords <- function(range.pair) {
  ts.to.fit <- ts[range.pair[1]:range.pair[2]] %>%
    as.data.frame() %>%
    rownames_to_column()
  model <- lm(. ~ as.numeric(rowname), data=ts.to.fit)
  cords <- map(range.pair, function(point) {
    as.numeric(model$coefficients[1] + (model$coefficients[2] * point))
  })
  return(cords)
}

fit.point.bounds <- unlist(map(points.to.map,GetLMCords)) %>%
  cbind(unlist(points.to.map))

colnames(fit.point.bounds) <- c("fit.hum","rowname")
fit.point.bounds <- as.data.frame(fit.point.bounds)
fit.point.bounds$rowname <- fit.point.bounds$rowname+1


frame.to.plot <- 
  aisdank_df %>%
  filter(name=="V1C0001") %>%
  unnest() %>%
  filter(time>=lastWeek*1.00027) %>%
  filter(time<=max(time-1000000)) %>%
  rownames_to_column()

frame.to.plot$rowname <- as.numeric(frame.to.plot$rowname)

frame.to.plot <- frame.to.plot %>% left_join(fit.point.bounds,"rowname")


frame.to.plot %>%
  #pull(hum) %>%
  #plot(type="l")
  mutate(splitTime = ifelse(rowname %in% accepted$split,time,NA)) %>%
  #mutate(fit.point = ifelse(rowname %in% fit.point.bounds$point,fit.point.bounds$fit.hum[]))
  #mutate(rownameNum = as.numeric(rowname)) %>%
  ggplot +
  geom_line(aes(x=time, y=temp)) +
  #geom_line(aes(x=time, y=c(diff(time),NA)))
  #geom_point(aes(x=time, y=fit.hum)) +
  geom_vline(aes(xintercept=splitTime))

diff(frame.to.plot$time) %>%
  as.tibble() %>%
  filter(value > 1) %>%
  filter(value < 310000) %$%
  hist(.$value, breaks=1000)

stream
moduleNames <- getModuleNames('~/Atomata/dataLab/AIsDank_raw.json')

aisdank_df <- map(moduleNames,
                  function(name) cleanRawJSONResponse('~/Atomata/dataLab/AIsDank_raw.json',name)
                  %>%
                  add_column(name))%>% reduce(bind_rows)

aisdank_df <- aisdank_df %>% group_by(name) %>% nest()

map()

stream$points$split + 3000

timeRange <- aisdank_df %>%
  filter(name=="V1C0001") %>%
  unnest()%>%
  select('time') %>%
  range()

halfway <- ((timeRange %>% diff()) / 2) + timeRange[1]
lastWeek <- timeRange[2] - 6.048e+8

aisdank_df %>%
  filter(name=="V1C0001") %>%
  unnest()%>%
  


ariam_df <- cleanRawJSONResponse('AIsDank_raw.json','AriaM')

ariam_df %>%
  ggplot +
  aes(x=time,y=co2)+
  geom_point()
 


test <- aisdank_df %>%
  filter(name=="V1C0001") %>%
  unnest() %>%
  filter(time>=lastWeek*1.00025) %>%
  filter(time<=max(time-100000000)) %>%
  rownames_to_column() %>%
  mutate(splitTime = ifelse(rowname %in% accepted$split,time,NA)) %>%
  filter(as.numeric(rowname) >= 251)

test$rowname <- as.numeric(test$rowname)

lm(temp ~ rowname, data=test) %>% residuals() %>% map(function(x) x^2) %>% reduce(`+`)

  do(data.frame(tidy(lm(time ~ temp, data=.)))) %>%
  residuals() %>% map(function(x) x^2) %>% reduce(`+`)
  unnest()
  mutate(Intercept=coef(model)[1], Slope=coef(model)[2]) %>%
  select(-model)

test <- test %>% as.tibble() %>%
  rownames_to_column() %>%
  set_colnames(c("time","dat"))
  
test$time <- as.numeric(test$time)

test.stream <- test %>%
  detectChangePointsBatch_Stream(1,5)

accepted <- test.stream$points %>% filter(accepted.change.point)

test.stream <- test.stream$funct()

test %>%
  mutate(splitTime = ifelse(time %in% accepted$split,time,NA)) %>%
  ggplot +
  geom_line(aes(x=time, y=dat)) +
  geom_vline(aes(xintercept=splitTime))
  ggplot
 
