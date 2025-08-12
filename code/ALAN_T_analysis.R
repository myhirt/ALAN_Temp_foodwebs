library(ATNr)
library(tidyverse)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

res <- read.csv("../output/ALAN_T_results.csv")


ggplot(filter(res, temperature == "15"), aes(x = light, y = tot_ext, color = as.factor(light.effect))) +
  geom_smooth(method = "lm")
