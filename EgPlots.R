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
  mutate(x_dat = x + rnorm(1, 0, 1)) %>%
  mutate(y_dat = y + rnorm(1, 0, 1))

out