setwd("~/Downloads/temp_BJ_revision/Final_push/Push2/Figure4")
library(StochKit2R) 
warning("StochKit2R is buggy in multicore mode. Doesn't tally all trajectories when saving histograms.")
library(tidyverse)
library(reshape2)

files <- list.files() %>%
  subset(., endsWith(., ".xml"))

map(files, function(x) ssa(modelFile = x,
                           outputDir = x %>% paste0("Kinetic_", .),
                           time = 36000,
                           realizations = 10^5,
                           intervals = 2,
                           keepHistograms = TRUE,
                           bins = 1000,
                           p = 1,
                           force = TRUE))


