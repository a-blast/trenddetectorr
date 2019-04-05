library(tidyverse)
library(memoise)
library(rredis)
library(digest)
library(purrr)

BatchChangePointsStream <- function(time.series,
                                    max.exponential.basis=1,
                                    max.num.points){

    candidates <- c()
    new.segment.splits <- c()

### initialize the stream ###
    StreamInit <- function() {
        initial.range <- c(1,length(time.series$time),"left")

        first.change.point <- FindAndUpdateCandidates(initial.range, 
                                                      set.change.point = TRUE)

        left.split.range <- c(1,first.change.point$split,"left")
        right.split.range <- c(first.change.point$split+1,
                               length(time.series$time),"right")

        new.segment.splits <<- list(left.split.range,right.split.range)
    }

    FindAndUpdateCandidates <- function(ts.bounds, set.change.point=FALSE){
        candidate.optimal.likelihood.criteria <- Inf

        low.bound <- as.numeric(ts.bounds[1])
        high.bound <- as.numeric(ts.bounds[2])

        if((high.bound-low.bound)/max.num.points <= 2){
            #print("TS BOUNDS TOOSMALL")
            return()
        }

        points.to.check.seq <- (low.bound+max.num.points):
            (high.bound-max.num.points-1)

        for (time.point in points.to.check.seq) {
            left.hand.split <- FindLikelihoodCriteria(
                time.series %>% 
                slice(low.bound:time.point)
            )

            right.hand.split <- FindLikelihoodCriteria(
                time.series %>% 
                slice((time.point+1):high.bound)
            )

            likelyhood.criteria <-  left.hand.split + right.hand.split

            if (likelyhood.criteria < candidate.optimal.likelihood.criteria) {
                candidate.split <- time.point
                candidate.optimal.likelihood.criteria <- likelyhood.criteria
                optimalleft.hand.split <- left.hand.split
                optimalright.hand.split <- right.hand.split
            }


        }

        if (ts.bounds[3] == "left"){
            parent.point <- high.bound
            if (is.null(candidates)){
                parent.loss <- FindLikelihoodCriteria(
                    time.series %>% 
                    slice(low.bound:high.bound)
                )
            } else {
                parent.loss <- candidates %>% 
                    filter(split == parent.point) %>%
                    pull(left.loss)
            }
        } else {
            parent.point <- low.bound-1
            if (is.null(candidates)){
                parent.loss <- FindLikelihoodCriteria(
                    time.series %>% 
                    slice(low.bound:high.bound)
                )
            } else {
                parent.loss <- candidates %>% 
                    filter(split == parent.point) %>%
                    pull(right.loss)
            }
        }

        candidates <<- candidates %>% 
            bind_rows(data.frame(split=candidate.split,
                                 split.time=time.series$time[candidate.split],
                                 left.loss=optimalleft.hand.split,
                                 right.loss=optimalright.hand.split,
                                 parent=parent.point,
                                 parent.side=ts.bounds[3],
                                 parent.side.loss=parent.loss,
                                 accepted.change.point=set.change.point))

        return(candidates)

    }

    FindLikelihoodCriteria <- function(time.series.subset){
        time.series.subset$time <- time.series.subset$time - min(time.series.subset$time)
        se <- lm(dat~poly(time,max.exponential.basis), data=time.series.subset) %>% 
            residuals() %>% 
            map(function(x) x^2) %>% 
            reduce(`+`)
        return(se)
    }

    GetNeighbors <- function(point) {

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
            neighbor.high <- length(time.series$time)
        }
        return(list(low=neighbor.low,high=neighbor.high))
    }

    SetNewSegmentSplits <- function(optimal.point.split) {
                                        #get the neighbor points around the new optimal point

        neighbor.points <- GetNeighbors(optimal.point.split)

        lowSegment <- c(neighbor.points$low,optimal.point.split,"left")
        highSegment <- c(optimal.point.split+1,neighbor.points$high,"right")
        new.segment.splits <<- list(lowSegment, highSegment)
    }

    NextChangePoint <- function() {

        map(new.segment.splits,FindAndUpdateCandidates)

        ## pick the change point with the lowest loss
        optimal.point <- candidates %>% 
            filter(!accepted.change.point) %>%
            mutate(relative.scores = (left.loss+right.loss)-parent.side.loss ) %>%
            filter(relative.scores == min(relative.scores))

        optimal.point$accepted.change.point[1] <- TRUE

        candidates <<- candidates %>%
            filter(split != optimal.point$split[1]) %>%
            bind_rows(optimal.point)

        ## add it to the chage points
        ## update the new segments of the time series
        SetNewSegmentSplits(optimal.point$split)

        return(list(points=candidates,funct=NextChangePoint))

    }

    StreamInit()
    return(NextChangePoint())

}


ParBatchChangePointsStream <- function(time.series,
                                    max.exponential.basis=1,
                                    max.num.points,
                                    start.point=1){

    redisConnect()
    candidates <- c()
    new.segment.splits <- c()

### initialize the stream ###
    StreamInit <- function() {
        initial.range <- c(1,length(time.series$time),"left")

        first.change.point <- FindAndUpdateCandidates(initial.range, 
                                                      set.change.point = TRUE)

        left.split.range <- c(1,first.change.point$split,"left")
        right.split.range <- c(first.change.point$split+1,
                               length(time.series$time),"right")

        new.segment.splits <<- list(left.split.range,right.split.range)
    }

    FindAndUpdateCandidates <- function(ts.bounds, set.change.point=FALSE){
        candidate.optimal.likelihood.criteria <- Inf

        low.bound <- as.numeric(ts.bounds[1])
        high.bound <- as.numeric(ts.bounds[2])

        if((high.bound-low.bound)/max.num.points <= 2){
            return()
        }

        points.to.check.seq <- (low.bound+max.num.points):
            (high.bound-max.num.points-1)

        for (time.point in points.to.check.seq) {

            left.hand.digest <- digest((low.bound+start.point):(time.point+start.point))
            right.hand.digest <- digest((time.point+start.point+1):(high.bound+start.point))

            left.hand.split <- redisGet(paste(left.hand.digest,"test1",sep=":"))

            if(is.null(left.hand.split)){
                left.hand.split <- MemFindLikelihoodCriteria(
                    time.series %>% 
                    slice(low.bound:time.point)
                )
                redisSet(paste(left.hand.digest,"test2",sep=":"),left.hand.split)
            }#else{print("leftGet")}

            right.hand.split <- redisGet(right.hand.digest)

            if(is.null(right.hand.split)){
                right.hand.split <- MemFindLikelihoodCriteria(
                    time.series %>% 
                    slice((time.point+1):high.bound)
                )
                redisSet(paste(right.hand.digest,"test2",sep=":"),right.hand.split)
            }#else{print("rightGet")}

            likelyhood.criteria <-  left.hand.split + right.hand.split

            if (likelyhood.criteria < candidate.optimal.likelihood.criteria) {
                candidate.split <- time.point
                candidate.optimal.likelihood.criteria <- likelyhood.criteria
                optimalleft.hand.split <- left.hand.split
                optimalright.hand.split <- right.hand.split
            }


        }

        if (ts.bounds[3] == "left"){
            parent.point <- high.bound
            if (is.null(candidates)){
                parent.loss <- FindLikelihoodCriteria(
                    time.series %>% 
                    slice(low.bound:high.bound)
                )
            } else {
                parent.loss <- candidates %>% 
                    filter(split == parent.point) %>%
                    pull(left.loss)
            }
        } else {
            parent.point <- low.bound-1
            if (is.null(candidates)){
                parent.loss <- FindLikelihoodCriteria(
                    time.series %>% 
                    slice(low.bound:high.bound)
                )
            } else {
                parent.loss <- candidates %>% 
                    filter(split == parent.point) %>%
                    pull(right.loss)
            }
        }

        newCandidate <-data.frame(split=candidate.split,
                                  split.time=time.series$time[candidate.split],
                                  left.loss=optimalleft.hand.split,
                                  right.loss=optimalright.hand.split,
                                  parent=parent.point,
                                  parent.side=ts.bounds[3],
                                  parent.side.loss=parent.loss,
                                  accepted.change.point=set.change.point)

        candidates <<- candidates %>% 
            bind_rows(newCandidate)

        return(newCandidate)

    }

    FindLikelihoodCriteria <- function(time.series.subset){
        time.series.subset$time <- time.series.subset$time - min(time.series.subset$time)
        se <- lm(dat~poly(time,max.exponential.basis), data=time.series.subset) %>% 
            residuals() %>% 
            map(function(x) x^2) %>% 
            reduce(`+`)
        return(se)
    }

    MemFindLikelihoodCriteria <- memoise(FindLikelihoodCriteria)

    GetNeighbors <- function(point) {

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
            neighbor.high <- length(time.series$time)
        }
        return(list(low=neighbor.low,high=neighbor.high))
    }

    SetNewSegmentSplits <- function(optimal.point.split) {
                                        #get the neighbor points around the new optimal point

        neighbor.points <- GetNeighbors(optimal.point.split)

        lowSegment <- c(neighbor.points$low,optimal.point.split,"left")
        highSegment <- c(optimal.point.split+1,neighbor.points$high,"right")
        new.segment.splits <<- list(lowSegment, highSegment)
    }

    NextChangePoint <- function() {
        
        segment.digests <- map(new.segment.splits,
                               function(segment) digest(c(FindAndUpdateCandidates,
                                                          as.numeric(segment[1:2])+start.point,
                                                          "temp")))

        map(segment.digests,redisGet) %>%
            list(segment.digests,new.segment.splits) %>%
            pmap(function(redisOut,digest,segment) {
                if(is.null(redisOut)){
                    newCandidate <- FindAndUpdateCandidates(segment)
                    redisSet(digest,newCandidate)
                } else {
                    candidates <<- candidates %>%
                        bind_rows(redisOut)
                }
            })

        
        ## pick the change point with the lowest loss
        optimal.point <- candidates %>% 
            filter(!accepted.change.point) %>%
            mutate(relative.scores = (left.loss+right.loss)-parent.side.loss ) %>%
            filter(relative.scores == min(relative.scores))

        optimal.point$accepted.change.point[1] <- TRUE

        candidates <<- candidates %>%
            filter(split != optimal.point$split[1]) %>%
            bind_rows(optimal.point)

        ## add it to the chage points
        ## update the new segments of the time series
        SetNewSegmentSplits(optimal.point$split)

        return(list(points=candidates,funct=NextChangePoint))

    }

    StreamInit()
    return(NextChangePoint())

}

