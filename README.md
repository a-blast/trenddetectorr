# trenddetectorr

Break a time series into linear model trends in a hierarchichal & realtime fashion to a selected minnimum loss criterion.

Currently incomplete as a package, but runnable as funcitons. See workspace.R, where I've been playing with calling the trenddetector functions I've written. Got a little wild with Closures, so the code needs to be refactored a bit for better unit testability.

I ultimately wrote this to play with the implementation of realtime sytems & get a feeling for implementing larger projects in R. 
###BUT:
If I were to do it over again for utilization in a production system I would write it in python (or Clojure!!). I love R for prototyping ML models & data wrangling, but it can be rather unwieldy & slow for a production system.
