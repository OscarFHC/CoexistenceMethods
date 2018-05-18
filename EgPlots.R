##### This is the R scripts for making example plots in the project of multiple species coexistence methods/models
###################################################################################################################
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
#################################################################################
##### The LV model ##############################################################
#################################################################################
##### The pre-determined model ######################################
LVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx1 = x1*rx*(1 - a11*x1)
    dy1 = y1*ry*(1 - a22*y1)
    dx = x*rx*(1 - a11*x - a21*y)
    dy = y*ry*(1 - a21*x - a22*y)
    return(list(c(dx1, dy1, dx, dy)))
  })
}
##### The pre-determined model ######################################
##### Parameter values for the model ################################
Pars <- c(rx = 0.1, ry = 0.05, a11 = 0.8, a12 = 0.6, a21 = 0.6, a22 = 1.5)
State <- c(x1 = 0.1, y1 = 0.1, x = 0.1, y = 0.1)
Time <- seq(0, 600, by = 1)
##### Parameter values for the model ################################
##### creating the fake data set ####################################
out <- as.data.frame(ode(func = LVmod, y = State, parms = Pars, times = Time, hmax=0.1, maxsteps=10000)) %>%
  mutate(x1_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(y1_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(x_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(y_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(x1_dat = x + x_rand) %>%
  mutate(y1_dat = y + y_rand) %>%
  mutate(x_dat = x + x_rand) %>%
  mutate(y_dat = y + y_rand)
##### creating the fake data set ####################################

##### Fitting parameters
####### for Mono culture_x ####################################  
xstart <- c(x = 0.1)
times <- out$time
  
LotVmod_x <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - alpha11*x)
    return(list(c(dx)))
  })
}

nll <- function (par) {
  y <- ode(func = LotVmod_x, y = xstart, parms = par, times = times, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out$x_dat, mean = y[,2], sd = 0.1, log = TRUE))
}

fit_x <- optim(par = c(rx = 0.1, alpha11 = 0.8), fn = nll)
####### for Mono culture_x ####################################
####### for Mono culture_y ####################################
ystart <- c(y = 0.1)
times <- out$time

LotVmod_y <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dy = y*ry*(1 - alpha22*y)
    return(list(c(dy)))
  })
}

nll <- function (par) {
  y <- ode(func = LotVmod_y, y = ystart, parms = par, times = times, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out$y_dat, mean = y[,2], sd = 0.1, log = TRUE))
}

fit_y <- optim(par = c(ry = 0.1, alpha22 = 0.8), fn = nll)
####### for Mono culture_y ####################################
####### for bi-culture ########################################
xstart <- c(x = 0.1, y = 0.1)
times <- out$time

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - alpha11*x - alpha12*y)
    dy = y*ry*(1 - alpha21*x - alpha22*y)
    return(list(c(dx, dy)))
  })
}

nll <- function (par) {
  y <- ode(func = LotVmod, y = xstart, parms = par, times = times, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out$x_dat, mean = y[,2], sd = 0.1, log = TRUE) + 
       dnorm(x = out$y_dat, mean = y[,3], sd = 0.1, log = TRUE))
}

LV_fit_two <- optim(par = c(rx = 0.1, ry = 0.08, alpha11 = 0.8, alpha12 = 0.6, alpha21 = 0.6, alpha22 = 0.9),
                 fn = nll)
####### for bi-culture ########################################

####### use fitted parameters to plot a fitted line #################
fitted_Param <- c(rx = LV_fit_two$par[[1]], ry = LV_fit_two$par[[2]], a11 = LV_fit_two$par[[3]], 
                  a12 = LV_fit_two$par[[4]], a21 = LV_fit_two$par[[5]], a22 = LV_fit_two$par[[6]])
State <- c(x1 = 0.01, y1 = 0.01, x = 0.1, y = 0.1)
Time <- seq(0, 600, by = 1)

fitted_Curve <- as.data.frame(ode(func = LVmod, y = State, parms = fitted_Param, times = Time,
                                  hmax = 0.1, maxsteps = 10000))

out <- out %>%
  subset(select = c("time", "x_dat", "y_dat")) %>%
  gather(key = "sp", value = "den", -time)

LV_plot = ggplot() + 
  geom_point(data = out, aes(x = time, y = den, shape = factor(sp)), size = 1) + 
  scale_shape_manual(values = c(1, 19), labels = c("Species 1", "Species 2"),
                     guide_legend("")) +
  geom_line(data = fitted_Curve, aes(x = time, y = x), linetype = "dashed", size = 1.5) +
  geom_line(data = fitted_Curve, aes(x = time, y = y), linetype = "dashed", size = 1.5) + 
  labs(x = "Time", y = "Population density") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 12, b = 12, l = 12, "pt"),
        legend.margin = margin(0, 0, 0, 0, "pt"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Fig2_LV.pdf", 
       plot = LV_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
####### use fitted parameters to plot a fitted line #################
#################################################################################
##### The LV model ##############################################################
#################################################################################

#################################################################################
##### The sensitivity method ####################################################
#################################################################################
##### Create fake data set 
####### calculate the carrying capacity of competiting sp_i ##########
LVmod_x1 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx1 = x1*rx1*(1 - a11*x1)
    return(list(c(dx1)))
  })
}

Pars_x1 <- c(rx1 = 0.1, a11 = 0.8)
State <- c(x1 = 0.01)
Time <- seq(0, 600, by = 1)

out_x1 <- as.data.frame(ode(func = LVmod_x1, y = State, parms = Pars_x1, times = Time, hmax=0.1, maxsteps=10000)) %>%
  mutate(x1_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(x1_dat = x1 + x1_rand) 
Kx <- mean(out_x1[c((length(Time)-100):length(Time)), "x1_dat"])

nll_x <- function (par) {
  y <- ode(func = LVmod_x1, y = State, parms = par, times = Time, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out_x1$x1_dat, mean = y[,2], sd = 0.1, log = TRUE))
}

fit_x <- optim(par = c(rx1 = 0.1, a11 = 0.8), fn = nll_x)
fitted_Curve_x1 <- as.data.frame(ode(func = LVmod_x1, y = State, parms = fit_x$par, times = Time,
                                     hmax = 0.1, maxsteps = 10000))
####### calculate the carrying capacity of competiting sp_i ##########
####### calculate the carrying capacity of competiting sp_j ##########
LVmod_y1 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dy1 = y1*ry1*(1 - a22*y1)
    return(list(c(dy1)))
  })
}

Pars_y <- c(ry1 = 0.05, a22 = 1.5)
State <- c(y1 = 0.01)
Time <- seq(0, 600, by = 1)

out_y1 <- as.data.frame(ode(func = LVmod_y1, y = State, parms = Pars_y, times = Time, hmax=0.1, maxsteps=10000)) %>%
  mutate(y1_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(y1_dat = y1 + y1_rand) 
Ky <- mean(out_y1[c((length(Time)-100):length(Time)), "y1_dat"])

nll_y <- function (par) {
  y <- ode(func = LVmod_y1, y = State, parms = par, times = Time, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out_y1$y1_dat, mean = y[,2], sd = 0.1, log = TRUE))
}

fit_y <- optim(par = c(ry1 = 0.1, a22 = 0.8), fn = nll_y)
fitted_Curve_y1 <- as.data.frame(ode(func = LVmod_y1, y = State, parms = fit_y$par, times = Time,
                                     hmax = 0.1, maxsteps = 10000))
####### calculate the carrying capacity of competiting sp_j ##########

##### Estimate invasion per capita growth rate for i invading j
####### create data for sp_i invading sp_j ##########################
LVmod_x <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - a11*x - a12*y)
    return(list(c(dx)))
  })
}

Pars <- c(rx = 0.1, a11 = 0.8, a12 = 0.6, y = Ky)
State <- c(x = 0.01)
Time <- seq(0, 600, by = 1)

out <- as.data.frame(ode(func = LVmod_x, y = State, parms = Pars, times = Time, hmax = 0.1, maxsteps = 10000)) %>%
  mutate(x_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(x_dat = x + x_rand)
####### create data for sp_i alone and sp_i invading sp_j ############
####### estimate invading growth rate for sp_i invading sp_j #########
out_x <- out %>%
  subset(time %in% c(0:50)) %>%
  subset(select = c("time", "x_dat"))

Expmod_b <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx
    return(list(c(dx)))
  })
}

Time_x <- out_x$time

nll_x <- function (x) {
  par <- c(rx = x)
  y <- ode(func = Expmod_b, y = State, parms = par, times = Time_x, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out_x$x_dat, mean = y[,2], sd = 0.1, log = TRUE))
}

fit_sen_x <- optimize(nll_x, c(0, 1), tol = 0.0001)

fitted_Curve_x <- as.data.frame(ode(func = Expmod_b, y = State, parms = c(rx = fit_sen_x$minimum), times = Time_x,
                                    hmax = 0.1, maxsteps = 10000))

a12 <- (fit_x$par["rx1"] - fit_sen_x$minimum)/fit_x$par["rx1"]
a12_correct <- (fit_x$par["rx1"] - fit_sen_x$minimum)/fit_x$par["rx1"]/Ky
####### estimate invading growth rate for sp_i invading sp_j #########
####### use fitted parameters to plot a fitted line ##################
x_dat <- out_x1 %>%
  full_join(out, by = "time") %>%
  subset(select = c("time", "x1_dat", "x_dat")) %>%
  gather(key = "sp", value = "den", -time) %>%
  subset(den != "NA" & time < 125)

x_fit <- fitted_Curve_x1 %>%
  full_join(fitted_Curve_x, by = "time") %>%
  gather(key = "sp", value = "den_fit", -time) %>%
  subset(den_fit != "NA" & time < 125)

x_Sen_plot = ggplot() + 
  geom_point(data = x_dat, aes(x = time, y = den, shape = factor(sp)), size = 1) + 
  scale_shape_manual(values = c(1, 19), labels = c("invading", "alone"),
                     guide_legend("")) +
  geom_line(data = x_fit, aes(x = time, y = den_fit, linetype = factor(sp)), size = 1.5, show.legend = F) +
  scale_linetype_manual(values = c("dashed", "dashed"), labels = c("invading", "alone"),
                        guide_legend("")) +
  labs(title = expression("Species " * italic(i) * " invading species " * italic(j)), x =  "", y =  "" ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 6, r = 6, b = 0, l = 36, "pt"),
        plot.title = element_text(size = 18),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))
####### use fitted parameters to plot a fitted line ##################

##### Estimate invasion per capita growth rate for j invading i
####### create data for sp_j invading sp_i ##########################
LVmod_y <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dy = y*ry*(1 - a22*y - a21*x)
    return(list(c(dy)))
  })
}

Pars <- c(ry = 0.05, a22 = 1.5, a21 = 0.6, x = Kx)
State <- c(y = 0.01)
Time <- seq(0, 1200, by = 1)

out <- as.data.frame(ode(func = LVmod_y, y = State, parms = Pars, times = Time, hmax = 0.1, maxsteps = 10000)) %>%
  mutate(y_rand = rnorm(length(Time), 0, 0.05)) %>%
  mutate(y_dat = y + y_rand)
####### create data for sp_j alone and sp_i invading sp_j ############
####### estimate invading growth rate for sp_j invading sp_i #########
out_y <- out %>%
  subset(time %in% c(0:250)) %>%
  subset(select = c("time", "y_dat"))

Expmod_b <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dy = y*ry
    return(list(c(dy)))
  })
}

Times_y <- out_y$time

nll_y <- function (y) {
  par <- c(ry = y)
  y <- ode(func = Expmod_b, y = State, parms = par, times = Times_y, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = out_y$y_dat, mean = y[,2], sd = 0.1, log = TRUE))
}

fit_sen_y <- optimize(nll_y, c(0, 1), tol = 0.0001)

State <- c(y = 0.01)
Time <- seq(0, nrow(out_y), by = 1)

fitted_Curve_y <- as.data.frame(ode(func = Expmod_b, y = State, parms = c(ry = fit_sen_y$minimum), times = Times_y,
                                    hmax = 0.1, maxsteps = 10000))

a21 <- (fit_y$par["ry1"] - fit_sen_y$minimum)/fit_y$par["ry1"]
a21_correct <- (fit_y$par["ry1"] - fit_sen_y$minimum)/fit_y$par["ry1"]/Kx
####### estimate invading growth rate for sp_j invading sp_i #########
####### use fitted parameters to plot a fitted line ##################
y_dat <- out_y1 %>%
  full_join(out, by = "time") %>%
  subset(select = c("time", "y1_dat", "y_dat")) %>%
  gather(key = "sp", value = "den", -time) %>%
  subset(den != "NA" & time < 300)

y_fit <- fitted_Curve_y1 %>%
  full_join(fitted_Curve_y, by = "time") %>%
  gather(key = "sp", value = "den_fit", -time) %>%
  subset(den_fit != "NA" & time < 300)

y_Sen_plot = ggplot() + 
  geom_point(data = y_dat, aes(x = time, y = den, shape = factor(sp)), size = 1) + 
  scale_shape_manual(values = c(1, 19), labels = c("invading", "alone"),
                     guide_legend("")) +
  geom_line(data = y_fit, aes(x = time, y = den_fit, linetype = factor(sp)), size = 1.5, show.legend = F) +
  scale_linetype_manual(values = c("dashed", "dashed"), labels = c("invading", "alone"),
                        guide_legend("")) +
  labs(title = expression("Species " * italic(j) * " invading species " * italic(i)), x =  "", y =  "" ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 0, r = 6, b = 6, l = 36, "pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))
####### use fitted parameters to plot a fitted line ##################

####### Combine two invasion growth plots ############################
Sen_plot <- plot_grid(x_Sen_plot + theme(legend.position="none"), y_Sen_plot + theme(legend.position="none"), 
                    labels = c("a.", "b."), ncol = 1, align = 'v')
shared_legend <- get_legend(x_Sen_plot + theme(legend.position = c(0.5, 0.58)))
Sen_plot <- plot_grid(Sen_plot, shared_legend, nrow = 1, rel_widths = c(3, 0.3))

Sen_plot <- ggdraw(Sen_plot) + 
  draw_label(bquote("Time"), x = 0.52, y = 0.02, size = 18) + 
  draw_label(bquote("Population density"), angle = 90, x = 0.04, y = 0.52, size = 18)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Fig3_Sen.eps", 
       plot = Sen_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
####### use fitted parameters to plot a fitted line ##################
#################################################################################
##### The sensitivity method ####################################################
#################################################################################

#################################################################################
##### The NFD method ############################################################
#################################################################################
##### The pre-determined model ######################################
LVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx1 = x1*rx*(1 - a11*x1)
    dy1 = y1*ry*(1 - a22*y1)
    dx = x*rx*(1 - a11*x - a21*y)
    dy = y*ry*(1 - a21*x - a22*y)
    return(list(c(dx1, dy1, dx, dy)))
  })
}
##### The pre-determined model ######################################
##### Parameter values for the model ################################
Pars <- c(rx = 0.1, ry = 0.05, a11 = 0.8, a12 = 0.6, a21 = 0.6, a22 = 1.5)
State <- c(x1 = 0.1, y1 = 0.1, x = 0.1, y = 0.1)
Time <- seq(0, 600, by = 1)
##### Parameter values for the model ################################
##### creating the fake data set ####################################
out <- as.data.frame(ode(func = LVmod, y = State, parms = Pars, times = Time, hmax=0.1, maxsteps=10000))
out %>% ggplot() + 
  geom_point(aes(x = time, y = x1), shape = 1, col = "blue") + 
  geom_point(aes(x = time, y = y1), shape = 1, col = "red") + 
  geom_point(aes(x = time, y = x), shape = 19, col = "blue") + 
  geom_point(aes(x = time, y = y), shape = 19, col = "red")
##### creating the fake data set ####################################
##### creating NFD for the species_i ################################
LV_NFD_x <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - a11*x - a12*y)
    return(list(c(dx)))
  })
}
freq <- seq(from = 0.001, to = 1 - 0.001, by = 0.001)
NFD_x_dat <- data.frame("freq" = freq,
                        "Gr" = rep(0, length(freq)),
                        "InGr" = rep(0, length(freq)))
for (i in c(1:length(NFD_x_dat[,"freq"]))){
  fre_y <- 1 - NFD_x_dat[i,"freq"]
  fre_x <- NFD_x_dat[i,"freq"]
  Pars <- c(rx = 0.1, a11 = 0.8, a12 = 0.6, y = fre_y)
  State <- c(x = fre_x)
  Time <- seq(0, 1, by = 1)
  out <- as.data.frame(ode(func = LV_NFD_x, y = State, parms = Pars, times = Time, hmax=0.1, maxsteps=10000))
  NFD_x_dat[i,"Gr"] <- out[2,"x"] - out[1,"x"]
  NFD_x_dat[i,"InGr"] <- NFD_x_dat[i,"Gr"]/State
}

NFD_x_plot <- NFD_x_dat %>%
  ggplot() + 
  geom_line(aes(x = freq, y = InGr)) +
  labs(x = expression("Frequency (%) of species " * italic(i)), y =  "" ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 0, r = 0, b = 6, l = 24, "pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))

##### creating NFD for the species_i ################################
##### creating NFD for the species_j ################################
LV_NFD_y <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dy = y*ry*(1 - a22*y - a21*x)
    return(list(c(dy)))
  })
}
freq <- seq(from = 0.001, to = 1 - 0.001, by = 0.001)
NFD_y_dat <- data.frame("freq" = freq,
                        "Gr" = rep(0, length(freq)),
                        "InGr" = rep(0, length(freq)))
for (i in c(1:length(NFD_y_dat[,"freq"]))){
  fre_x <- 1 - NFD_y_dat[i,"freq"]
  fre_y <- NFD_y_dat[i,"freq"]
  Pars <- c(ry = 0.05, a22 = 1.5, a21 = 0.6, x = fre_x)
  State <- c(y = fre_y)
  Time <- seq(0, 1, by = 1)
  out <- as.data.frame(ode(func = LV_NFD_y, y = State, parms = Pars, times = Time, hmax=0.1, maxsteps=10000))
  NFD_y_dat[i,"Gr"] <- out[2,"y"] - out[1,"y"]
  NFD_y_dat[i,"InGr"] <- NFD_y_dat[i,"Gr"]/State
}

NFD_y_plot <- NFD_y_dat %>%
  ggplot() + 
  geom_line(aes(x = freq, y = InGr)) +
  labs(x = expression("Frequency (%) of species " * italic(j)), y =  "" ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 0, r = 6, b = 6, l = 6, "pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))

##### creating NFD for the species_j ################################
##### combine two NFD plots #########################################
NFD_plot <- plot_grid(NFD_x_plot, NFD_y_plot, labels = c("a.", "b."), nrow = 1, align = 'h')

NFD_plot <- ggdraw(NFD_plot) + 
  #draw_label(bquote("Frequency"), x = 0.52, y = 0.02, size = 18) + 
  draw_label(expression(italic("per capita") * "growth rate"), angle = 90, x = 0.02, y = 0.52, size = 18)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Fig4_NFD.pdf", 
       plot = NFD_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
##### combine two NFD plots #########################################

#################################################################################
##### The NFD method ############################################################
#################################################################################

