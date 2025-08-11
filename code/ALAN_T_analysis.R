library(ATNr)
library(tidyverse)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

res <- read.csv("../output/results_combined.csv")


ggplot(filter(res, temperature == "15"), aes(x = light, y = pers_tot, color = as.factor(light.effect))) +
  geom_smooth(method = "lm")
