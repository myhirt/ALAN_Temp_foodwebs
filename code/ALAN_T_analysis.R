library(ATNr)
library(tidyverse)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

res <- read.csv("../output/ALAN_T_results.csv")


ggplot(filter(res, temperature == "22"), aes(x = light, y = tot_ext, color = as.factor(light.effect))) +
  geom_smooth(method = "lm")


res.array <- read.csv("../output/results_combined.csv")

ggplot(res.array, aes(x = light, y = pers_tot, color = as.factor(temperature))) +
  geom_smooth(method = "lm") +
  facet_wrap(~light.effect) +
  theme_classic()

ggplot(filter(res.array, temperature == "15"), aes(x = light, y = pers_tot, color = as.factor(light.effect))) +
  geom_smooth(method = "lm") +
  theme_classic()
