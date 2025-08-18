library(ATNr)
library(tidyverse)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

res <- read.csv("../output/ALAN_T_results.csv")


ggplot(filter(res, temperature == "15"), aes(x = light, y = tot_ext, color = as.factor(light.effect))) +
  geom_smooth(method = "lm")


res.array <- read.csv("../output/results_combined.csv")

ggplot(filter(res.array, temperature == "20"), aes(x = light, y = pers_tot, color = as.factor(light.effect))) +
  geom_smooth(method = "lm") +
  theme_classic()

ggplot(res.array, aes(x = connectance, y = pers_tot, color = as.factor(temperature))) +
  geom_smooth(method = "lm") +
  theme_classic()
