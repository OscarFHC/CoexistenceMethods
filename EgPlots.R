### This is the R scripts for making 

library(deSolve)
library(ggplot)
library(tidyverse)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - alpha11*x - alpha21*y)
    dy = y*ry*(1 - alpha21*x - alpha22*y)
    return(list(c(dx, dy)))
  })
}

Pars <- c(rx = 0.1, ry = 0.08, alpha11 = 0.8, alpha21 = 0.6, alpha21 = 0.6, alpha22 = 0.9)
State <- c(x = 0.1, y = 0.1)
Time <- seq(0, 150, by = 1)

out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time)) %>%
  mutate(x_rand = rnorm(151, 0, 0.05)) %>%
  mutate(y_rand = rnorm(151, 0, 0.05)) %>%
  mutate(x_dat = x + x_rand) %>%
  mutate(y_dat = y + y_rand)

out %>%
  ggplot() + 
    geom_point(aes(x = time, y = x_dat), size = 1) + 
    geom_line(aes(x = time, y = x), linetype = 2, size = 2)
