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
####### Combine two invasion growth plots ############################
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
##### creating NFD data #############################################
LV_NFD <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - a11*x - a12*y)
    dy = y*ry*(1 - a21*x - a22*y)
    return(list(c(dx, dy)))
  })
}

NFD_dat <- expand.grid(freq = seq(from = 0.001, to = 1 - 0.001, by = 0.001),
                         B = seq(from = 0.1, to = 2, by = 0.1)) %>%
  mutate(InGr_x = 0) %>%
  mutate(InGr_y = 0)
  
for (i in c(1:nrow(NFD_x_dat))){
  B <- NFD_x_dat$"B"[i]
  den_x <- B * (NFD_x_dat$freq[i])
  den_y <- B * (1 - NFD_x_dat$freq[i])
  Pars <- c(rx = 0.1, a11 = 0.8, a12 = 0.6, ry = 0.05, a21 = 0.6, a22 = 1.5)
  NFD_dat$InGr_x[i] <- Pars[[1]] * (1 - Pars[[2]] * den_x - Pars[[3]] * den_y)
  NFD_dat$InGr_y[i] <- Pars[[4]] * (1 - Pars[[5]] * den_x - Pars[[6]] * den_y) 
}
##### creating NFD data #############################################
##### making NFD plot ###############################################
NFD_x_plot <- NFD_dat %>%
  subset(B %in% c(0.1, 1, 2)) %>% 
  ggplot() + 
  geom_line(aes(x = freq, y = InGr_x, linetype = factor(B)), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "longdash"), labels = c("0.1", "1", "2"),
                        guide_legend("Community \nbiomass (B)")) + 
  labs(x = expression("Frequency (%) of species " * italic(i)), y = "" ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 6, r = 0, b = 6, l = 36, "pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))

NFD_y_plot <- NFD_dat %>%
  subset(B %in% c(0.1, 1, 2)) %>% 
  ggplot() + 
  geom_line(aes(x = 1 - freq, y = InGr_y, linetype = factor(B)), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "longdash"), labels = c("0.1", "1", "2"),
                        guide_legend("Community \nbiomass (B)")) + 
  labs(x = expression("Frequency (%) of species " * italic(j)), y = "" ) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 6, r = 6, b = 6, l = 0, "pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))
##### making NFD plot ###############################################
####### Make 3D plot, but it is not satisfying enough thouth ########
# install.packages("plotly")
# library(plotly)
# NFD_x_dat <- list(freq = seq(from = 0.001, to = 1 - 0.001, by = 0.001),
#                   B = seq(from = 0.1, to = 100, by = 0.1),
#                   InGr = matrix(0, length(seq(from = 0.001, to = 1 - 0.001, by = 0.001)), 
#                                 length(seq(from = 0.1, to = 100, by = 0.1))))
# 
# for (i in c(1:length(NFD_x_dat$freq))){
#   for (j in c(1:length(NFD_x_dat$B))){
#     B <- NFD_x_dat$B[j]
#     den_y <- B * (1 - NFD_x_dat$freq[i])
#     den_x <- B * (NFD_x_dat$freq[i])
#     
#     Pars <- c(rx = 0.1, a11 = 0.8, a12 = 0.6)
#     NFD_x_dat$InGr[i,j] <- Pars[[1]] * (1 - Pars[[2]] * den_x - Pars[[3]] * den_y)
#   }
# }
# plot_ly(x = NFD_x_dat$freq, y = NFD_x_dat$B, z = NFD_x_dat$InGr) %>% add_surface()
####### Make 3D plot, but it is not satisfying enough thouth ########

##### combine two NFD plots #########################################
NFD_plot <- plot_grid(NFD_x_plot + theme(legend.position="none"), 
                      NFD_y_plot + theme(legend.position="none"), 
                      labels = c("a.", "b."), nrow = 1, align = 'h')
shared_legend <- get_legend(NFD_x_plot + theme(legend.position = c(0.5, 0.58)))
NFD_plot <- plot_grid(NFD_plot, shared_legend, nrow = 1, rel_widths = c(3, 0.3))
NFD_plot <- ggdraw(NFD_plot) + 
  draw_label(expression(italic("per capita") * "growth rate"), angle = 90, x = 0.04, y = 0.52, size = 18)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Fig4_NFD.pdf", 
       plot = NFD_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
##### combine two NFD plots #########################################

#################################################################################
##### The NFD method ############################################################
#################################################################################

#################################################################################
##### MacArthur's consumer resource model #######################################
#################################################################################

set.seed(20180521)
Consumption_dat <- data.frame("res" = seq(from = 1, to = 10, by = 1),
                              "cx" = rnorm(10, 10, 1), 
                              "cy" = rnorm(10, 12, 1))

MCR_plot <- ggplot() + 
  geom_point(data = Consumption_dat, aes(x = cy, y = cx)) +
  geom_abline(intercept = 0, slope = 1) +
  #geom_line(data = Consumption_dat, aes(x = cy, y = predict(lm(cx ~ cy, data = Consumption_dat)))) +
  scale_x_continuous(limits = c(8, 14), expand = c(0, 0)) +
  scale_y_continuous(limits = c(8, 14), expand = c(0, 0)) + 
  labs(x = expression("Consumption of species " * italic(i) * " on resource " * italic(l)),
       y = expression("Consumption of species " * italic(j) * " on resource " * italic(l))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 18),
        plot.margin = margin(t = 6, r = 0, b = 6, l = 0, "pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 6, l = 0, "pt")),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 6, b = 0, l = 6, "pt")))
MCR_plot <- ggdraw(MCR_plot) + 
  draw_label(expression("1:1 line"), x = 0.95, y = 0.9, size = 18)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Fig5_MCR.jpeg", 
       plot = MCR_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
#################################################################################
##### MacArthur's consumer resource model #######################################
#################################################################################

#################################################################################
##### Tilman's consumer resource model ##########################################
#################################################################################

TCR_dat <- data.frame("R1_v" = seq(from = 0.3, to = 2, length = 100),
                      "R1_h" = seq(from = 0.8, to = 2, length = 100),
                      "R2_v" = seq(from = 0.7, to = 2, length = 100),
                      "R2_h" = seq(from = 0.5, to = 2, length = 100))

TCR_plot <- TCR_dat %>%
  ggplot() + 
    geom_line(aes(x = 0.8, y = R1_v)) + 
    geom_line(aes(x = R1_h, y = 0.3)) + 
    geom_line(aes(x = 0.5, y = R2_v)) + 
    geom_line(aes(x = R2_h, y = 0.7)) +
    #scale_x_continuous(limits = c(0, 2), expand = c(0.5, 0.5)) +
    #scale_y_continuous(limits = c(0, 2), expand = c(0.5, 0.5)) + 
    labs(x = expression("Density of resource " * italic(i)),
         y = expression("Density of resource " * italic(j))) + 
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size = 18),
          plot.margin = margin(t = 6, r = 0, b = 6, l = 0, "pt"),
          legend.margin = margin(0, 0, 0, 0, "cm"),
          legend.key = element_rect(color = "white", fill = "white"), 
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.line = element_line(colour = "black"), 
          axis.text = element_text(size = 14), 
          axis.title.x = element_text(size = 18, margin = margin(t = 12, r = 0, b = 6, l = 0, "pt")),
          axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 6, b = 0, l = 6, "pt")))

TCR_plot <- ggdraw(TCR_plot) + 
  draw_label(expression("Zero growth isocline of Species " * italic(i)), x = 0.85, y = 0.15, size = 12) + 
  draw_label(expression("Zero growth isocline of Species " * italic(j)), x = 0.85, y = 0.34, size = 12)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Fig6_TCR.eps", 
       plot = TCR_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
#################################################################################
##### Tilman's consumer resource model ##########################################
#################################################################################
