library(foreach)
library(doParallel)
library(parallel)

stream <- V1C0001.df %>%
    filter(time >= lastWeek*1.00027) %>%
    filter(time<=max(time-1000000)) %>%
    select(temp, time) %>%
    set_colnames(c("dat","time")) %>%
    BatchChangePointsStream(1,6) %>%
    loop()

stream$points %>%
    filter(!is.na(relative.scores)) %>%
    rownames_to_column() %>%
    select(rowname,relative.scores) %>%
    plot(type="l")



PlotChangePoints <- function(stream.points){
    accepted <- stream.points %>% filter(accepted.change.point)

    plt <- V1C0001.df %>%
        filter(time >= lastWeek*1.00027) %>%
        filter(time<=max(time-1000000)) %>%
        select(temp, time) %>%
        set_colnames(c("dat","time")) %>%
        rownames_to_column() %>%
        mutate(rowname=as.numeric(rowname)) %>%
        mutate(splitTime = ifelse(rowname %in% accepted$split,time,NA)) %>%
        ggplot +
        geom_line(aes(y=dat,x=time)) +
        geom_vline(aes(xintercept=splitTime))

    print(plt)
}

loop <- function(stream){
    while(
        stream$points %>% filter(!is.na(relative.scores)) %>% pull(relative.scores) %>% max() <= -1){
            stream <- stream$funct()
            #PlotChangePoints(stream$points)
            #print(stream$points %>% filter(!is.na(relative.scores)) %>% pull(relative.scores) %>% max())
            #invisible(readline(prompt="Press [enter] to continue"))
        }
    return(stream)
}



no_cores <- detectCores() -1

registerDoParallel(no_cores)
stopImplicitCluster()
V1C0001.df <- aisdank_df %>% filter(name=="V1C0001") %>% unnest() %>% select(-name)

V1C0001.df

out <- foreach(start = 1:(288*2),#(dim(V1C0001.df)[1]- 288),
               .combine = bind_rows)  %dopar%
    GetSubStream(start)

accepted <- out %>%
    filter(accepted.change.point) %>%
    mutate(split2=split.time)

frame.to.plot <- V1C0001.df %>% head(n=1000) %>% rownames_to_column

frame.to.plot %>%
    full_join(accepted,by=c("rowname"="split")) %>%
    filter(time < 1.53941e+12) %>%
    filter(time > 1.53940e+12) %>%
    ggplot +
    geom_line(aes(x=time,y=scale(temp)))+
    geom_vline(aes(xintercept=split2),alpha=1/3)

frame.to.plot %>%
    full_join(accepted,by=c("rowname"="split")) %>%
    filter(time < 1.53941e+12) %>%
    filter(!is.na(split2)) %>%
    filter(time > 1.53940e+12) %>%
    select(-c(temp,hum,co2,Lux,vpd,split.time,left.loss,right.loss,
              parent,parent.side,parent.side.loss,accepted.change.point)) %>%
    group_by(split2) %>%
    count()

GetSubStream <- function(start){
    redisConnect()

    print(start)
    sub.stream.df <- V1C0001.df %>%
        slice(start:(start+288)) %>%
        select(temp,time) %>%
        set_colnames(c("dat","time"))

    sub.stream.digest <- digest(c(sub.stream.df,ParBatchChangePointsStream,1,6))
    sub.stream.points <- redisGet(sub.stream.digest)

    if(is.null(sub.stream.points)){
        sub.stream <- sub.stream.df %>%
            ParBatchChangePointsStream(1,6,start) %>%
            loop()
        sub.stream.points <- sub.stream$points %>% mutate(start.point=start)
        redisSet(sub.stream.digest,sub.stream.points)
    }

    sub.stream.points$split <- sub.stream.points$split + start

    return(sub.stream.points)
}

stopCluster(cl)


parLapply(cl, 2:4,
          function(exponent)
              2^exponent)

aisdank_df %>%
    filter(name=="V1C0001") %>%
    unnest() %>%
    select(temp,time)

    filter(time >= lastWeek*1.00027) %>%
    filter(time<=max(time-1000000)) %>%
    select(temp, time) %>%
    set_colnames(c("dat","time"))
