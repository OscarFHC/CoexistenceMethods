##### This is the R scripts for making example plots in the project of multiple species coexistence methods/models
###################################################################################################################
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
library(ggpubr)

###################################################################################################################
##### This section is to estimate parameter values for LV model from empirical data ###############################
###################################################################################################################
### Read in data
dat_Mono_ori <- read.table(file = "D:/Manuscript/CoexistenceMethods_Figs/Chl_Clo_Mono.csv", sep = ",", header = TRUE, fill = TRUE)
dat_Bi_ori <- read.table(file = "D:/Manuscript/CoexistenceMethods_Figs/Chl_Clo_Bi.csv", sep = ",", header = TRUE, fill = TRUE)

### Preping a function to use boostrapping for standard error calculation
NonParam <- function(den){ # This function does not log-transform data before doing permutation
  set.seed(1032)
  means = replicate(5000, mean(sample(den, size=length(den), replace=TRUE), na.rm = TRUE))
  sd(means, na.rm = TRUE)
}

### Organize data; use the mean of the three replicates as the empirical data for parameter fitting
### I scaled down the actual numbers (/1000000) to avoid possible numerical errors when fitting parameters
dat_Mono <- dat_Mono_ori %>%
  group_by(tp) %>%
  summarize(Chl_m = mean(chl/1000000, na.rm = TRUE),
            Clo_m = mean(clo/1000000, na.rm = TRUE),
            Chl_se = NonParam(chl/1000000),
            Clo_se = NonParam(clo/1000000))
dat_Bi <- dat_Bi_ori %>%
  group_by(tp) %>%
  summarize(Chl_m = mean(chl/1000000, na.rm = TRUE),
            Clo_m = mean(clo/1000000, na.rm = TRUE),
            Chl_se = NonParam(chl/1000000),
            Clo_se = NonParam(clo/1000000))

### Creating the model
LV_Mono <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - a11*x)
    return(list(c(dx)))
  })
}
LV_Bi <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(Kx - a11*x - a12*y)
    dy = y*ry*(Ky - a21*x - a22*y)
    return(list(c(dx, dy)))
  })
}

### set up initial values and time step for growth curve fitting
Init_Chl <- c(x = dat_Mono$Chl_m[1])
Init_Clo <- c(x = dat_Mono$Clo_m[1])
Time <- seq(1, nrow(dat_Mono), by = 1)

### Creating the negative log likelihood function to be minimized
nll_x_Mono <- function (param) {
  y <- ode(func = LV_Mono, y = Init_Chl, parms = param, times = Time, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = dat_Mono$Chl_m, mean = y[,2], sd = 0.1, log = TRUE))
}
nll_y_Mono <- function (param) {
  y <- ode(func = LV_Mono, y = Init_Clo, parms = param, times = Time, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = dat_Mono$Clo_m, mean = y[,2], sd = 0.1, log = TRUE))
}

### Use "optim" function to minimize negative log likelihood (maximize likelihood)
Par_x_init <- c(rx = 0.45, Kx = 2, a11 = 1)
fit_x <- optim(par = Par_x_init, fn = nll_x_Mono)
fit_y <- optim(par = Par_x_init, fn = nll_y_Mono)

### Set up initial values (the same as the empirical measurement) and time steps
Init <- c(x = dat_Bi$Chl_m[1], y = dat_Bi$Clo_m[1])
Time <- seq(1, nrow(dat_Bi), by = 1)
names(fit_y$par) <- c("ry", "Ky", "a22") # Rename the parameters to match the parameter name in the "LV_bi" function 

nll_Bi <- function (par) {
  y <- ode(func = LV_Bi, y = Init, parms = c(fit_x$par, fit_y$par, par), times = Time, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = dat_Bi$Chl_m, mean = y[,2], sd = 0.1, log = TRUE) + 
         dnorm(x = dat_Bi$Clo_m, mean = y[,3], sd = 0.1, log = TRUE))
}

### Use "optim" function to minimize negative log likelihood (maximize likelihood)
fit_Bi <- optim(par = c(a12 = 0.5, a21 = 0.5), fn = nll_Bi)
###################################################################################################################
##### This section is to estimate parameter values for LV model from empirical data ###############################
###################################################################################################################



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
Pars <- c(rx = 0.1, ry = 0.08, a11 = 0.8, a12 = 0.6, a21 = 0.6, a22 = 1.2)
State <- c(x1 = 0.01, y1 = 0.01, x = 0.01, y = 0.01)
Time <- seq(0, 100, by = 1)
##### Parameter values for the model ################################
##### creating the fake data set ####################################
Growth_Curve <- as.data.frame(ode(func = LVmod, y = State, parms = Pars, times = Time, hmax = 0.1, maxsteps = 10000))
##### creating the fake data set ####################################
##### use fitted parameters to plot a fitted line ###################
X1 <- ggplot() + 
  geom_line(data = Growth_Curve, aes(x = time, y = x1), color = "#000099", size = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  labs(title = expression(paste("Species " * italic(i) * " monoculture")), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 52), "pt"),
        plot.title = element_text(size = 20, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank())

Y1 <- ggplot() +
  geom_line(data = Growth_Curve, aes(x = time, y = y), color = "#FF6600", size = 2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  labs(title = expression(paste("Species " * italic(j) * " monoculture")), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 24), "pt"),
        plot.title = element_text(size = 20, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank())

Bi <- ggplot() + 
  geom_line(data = Growth_Curve, aes(x = time, y = x), color = "#000099", size = 2) +
  geom_line(data = Growth_Curve, aes(x = time, y = y), color = "#FF6600", size = 2) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) + 
  labs(title = expression(paste("Species " * italic(i) * " and " * italic(j) * " biculture")),
       x = expression(paste("Time")),
       y = "") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 52), "pt"),
        plot.title = element_text(size = 20, vjust = 0.3),
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12, r = 0, b = 12, l = 0)))

plot_up <- plot_grid(X1, Y1, labels = c("A.", "B."), ncol = 2)
Plot_Bi1 <- plot_grid(plot_up, Bi, labels = c("", "C."), ncol = 1) %>%
  ggdraw() + 
  draw_label(expression(paste("Population density")), angle = 90, x = 0.04, y = 0.5, size = 24)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Ver3/Fig1_LV.pdf", 
       plot = Plot_Bi1, width = 35, height = 24, units = c("cm"), dpi = 600)

##### use fitted parameters to plot a fitted line ###################
#################################################################################
##### The LV model ##############################################################
#################################################################################

#################################################################################
##### The sensitivity method ####################################################
#################################################################################
####### Creating the model to predict growth curves  ###############
LVmod_sen <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx1 = x1*rx*(1 - a11*x1)
    dy1 = y1*ry*(1 - a22*y1)
    dx = x*rx*(1 - a11*x - a21*(1/a22))
    dy = y*ry*(1 - a21*(1/a11) - a22*y)
    return(list(c(dx1, dy1, dx, dy)))
  })
}
####### Creating the model to predict growth curves  ###############
##### Parameter values for the model ################################
Pars <- c(rx = 0.1, ry = 0.08, a11 = 0.8, a12 = 0.6, a21 = 0.6, a22 = 1.2)
State <- c(x1 = 0.01, y1 = 0.01, x = 0.01, y = 0.01)
Time <- seq(0, 300, by = 1)
##### Parameter values for the model ################################
##### creating the fake data set ####################################
Growth_Sen <- as.data.frame(ode(func = LVmod_sen, y = State, parms = Pars, times = Time, hmax = 0.1, maxsteps = 10000))
##### creating the fake data set ####################################
##### use fitted parameters to plot a fitted line ###################
x_to_y <- ggplot() + 
  geom_line(data = Growth_Sen[1:250,], aes(x = time, y = y1), color = "#FF6600", size = 2) +
  geom_line(data = Growth_Sen[1:100,], aes(x = time + 150, y = x), color = "#000099", size = 2) +
  geom_segment(aes(x = 10, y = 0.05, xend = 25, yend = 0.35), color = "black", size = 2, arrow = arrow()) + 
  geom_segment(aes(x = 160, y = 0.05, xend = 175, yend = 0.35), color = "black", size = 2, arrow = arrow()) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 250)) + 
  labs(title = expression(paste("Species ", italic(i), " invading species ", italic(j))), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 52), "pt"),
        plot.title = element_text(size = 20, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank())

y_to_x <- ggplot() + 
  geom_line(data = Growth_Sen[1:250,], aes(x = time, y = x1), color = "#000099", size = 2) +
  geom_line(data = Growth_Sen[1:100,], aes(x = time + 150, y = y), color = "#FF6600", size = 2) +
  geom_segment(aes(x = 10, y = 0.05, xend = 25, yend = 0.75), color = "black", size = 2, arrow = arrow()) + 
  geom_segment(aes(x = 160, y = 0.05, xend = 175, yend = 0.35), color = "black", size = 2, arrow = arrow()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 250)) + 
  labs(title = expression(paste("Species ", italic(j), " invading species ", italic(i))), 
       x = expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 52), "pt"),
        plot.title = element_text(size = 20, vjust = 0.3),
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12, r = 0, b = 12, l = 0)))
##### use fitted parameters to plot a fitted line ####################
####### Combine two invasion growth plots ############################
uy <- expression(atop("Growth rate of species", paste(italic(j), " growing alone (", mu[j], ")")))
uxy <- expression(atop("Growth rate of species", paste(italic(i), " invading species (", mu[j], ")")))
ux <- expression(atop("Growth rate of species", paste(italic(i), " growing alone (", mu[i], ")")))
uyx <- expression(atop("Growth rate of species", paste(italic(j), " invading species (", mu[i], ")")))
Sen_plot <- plot_grid(x_to_y , y_to_x, labels = c("A.", "B."), ncol = 1) %>% 
  ggdraw() + 
    draw_label(bquote("Population density"), angle = 90, x = 0.04, y = 0.5, size = 24) + 
    draw_label(uy, x = 0.175, y = 0.69, size = 20) +
    draw_label(uxy, x = 0.7, y = 0.69, size = 20) +
    draw_label(ux, x = 0.175, y = 0.35, size = 20) +
    draw_label(uyx, x = 0.7, y = 0.225, size = 20)
ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Ver3/Fig2_Sen.pdf", 
       plot = Sen_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
####### Combine two invasion growth plots ############################
#################################################################################
##### The sensitivity method ####################################################
#################################################################################

#################################################################################
##### The NFD method ############################################################
#################################################################################
##### creating NFD data #############################################
FD_dat <- expand.grid(freq = seq(from = 0.001, to = 1 - 0.001, by = 0.001),
                      B = seq(from = 0.1, to = 2, by = 0.1)) %>%
          mutate(InGr_x = 0) %>%
          mutate(InGr_y = 0)
Pars <- c(rx = 0.1, ry = 0.08, a11 = 0.8, a12 = 0.6, a21 = 0.6, a22 = 1.2)  

for (i in c(1:nrow(FD_dat))){
  B <- FD_dat$"B"[i]
  den_x <- B * (FD_dat$freq[i])
  den_y <- B * (1 - FD_dat$freq[i])
  FD_dat$InGr_x[i] <- Pars[[1]] * (1 - Pars[[3]] * den_x - Pars[[4]] * den_y)
  FD_dat$InGr_y[i] <- Pars[[2]] * (1 - Pars[[5]] * den_x - Pars[[6]] * den_y) 
}
##### creating NFD data #############################################
##### making NFD plot ###############################################
FD_x_plot1 <- FD_dat %>% # This plot is just for the legend
  subset(B %in% c(0.1, 0.5, 1)) %>% 
  ggplot() + 
  geom_line(aes(x = freq, y = InGr_x, linetype = factor(B)), size = 1) + 
  geom_segment(aes(x = 0.4, y = -0.036, xend = 0.6, yend = -0.036), color = "red", size = 2) + 
  geom_segment(aes(x = 0.6, y = -0.036, xend = 0.6, yend = -0.044), color = "red", size = 2) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"), labels = c("0.1", "0.5", "1"),
                        guide_legend("Community \nbiomass (B)")) + 
  labs(x = expression("Frequency (%) of species " * italic(i)), 
       y = expression(paste(italic("per capita"), "growth rate of species ", italic(i)))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 52), "pt"),
        plot.title = element_text(size = 20, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12, r = 0, b = 12, l = 0)),
        axis.title.y = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 8, b = 0, l = 0)))

FD_x_plot <- FD_dat %>%
  subset(B %in% c(0.1, 0.5, 1)) %>% 
  ggplot() + 
  geom_line(aes(x = freq, y = InGr_x, linetype = factor(B)), size = 1, color = "#000099") + 
  geom_segment(aes(x = 0.4, y = 0.032, xend = 0.6, yend = 0.032), color = "red", size = 2) + 
  geom_segment(aes(x = 0.6, y = 0.032, xend = 0.6, yend = 0.028), color = "red", size = 2) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"), labels = c("0.1", "0.5", "1"),
                        guide_legend("Community \nbiomass (B)")) + 
  labs(x = expression("Frequency (%) of species " * italic(i)), 
       y = expression(paste(italic("per capita"), "growth rate of species ", italic(i)))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 48), "pt"),
        plot.title = element_text(size = 20, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12, r = 0, b = 12, l = 0)),
        axis.title.y = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 8, b = 0, l = 0)))

FD_y_plot <- FD_dat %>%
  subset(B %in% c(0.1, 0.5, 1)) %>% 
  ggplot() + 
  geom_line(aes(x = 1 - freq, y = InGr_y, linetype = factor(B)), size = 1, color = "#FF6600") + 
  geom_segment(aes(x = 0.4, y = 0.0128, xend = 0.6, yend = 0.0128), color = "red", size = 2) + 
  geom_segment(aes(x = 0.6, y = 0.0128, xend = 0.6, yend = 0.0032), color = "red", size = 2) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"), labels = c("0.1", "0.5", "1"),
                        guide_legend("Community \nbiomass (B)")) + 
  labs(x = expression("Frequency (%) of species " * italic(j)), 
       y = expression(paste(italic("per capita"), "growth rate of species ", italic(j))), size = 18) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 12, 0, 36), "pt"),
        plot.title = element_text(size = 20, vjust = 0.4),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12, r = 0, b = 12, l = 0)),
        axis.title.y = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 8, b = 0, l = 0)))

##### making NFD plot ###############################################
##### combine two NFD plots #########################################
FD_plot <- plot_grid(FD_x_plot + theme(legend.position="none"), 
                     FD_y_plot + theme(legend.position="none"), 
                     labels = c("A.", "B."), nrow = 1, align = 'h')
shared_legend <- get_legend(FD_x_plot1 + theme(legend.position = c(0.5, 0.58), 
                                               legend.title = element_text(size = 20),
                                               legend.text = element_text(size = 20)))
FD_plot <- plot_grid(FD_plot, shared_legend, nrow = 1, rel_widths = c(3, 0.4)) %>%
  ggdraw() + 
  draw_label(expression("Slope = NFD " != alpha), x = 0.345, y = 0.3, size = 20, colour = "red") + 
  draw_label(expression("Slope = NFD " != alpha), x = 0.78, y = 0.425, size = 20, colour = "red")

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Ver2/Fig3_FD.tiff", 
       plot = FD_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
##### combine two NFD plots #########################################

#################################################################################
##### The NFD method ############################################################
#################################################################################

#################################################################################
##### MacArthur's consumer resource model #######################################
#################################################################################
set.seed(20180521)
Cons <- data.frame("cx" = rnorm(10, 10, 1),
                   "cy" = rnorm(10, 12, 1))
Cons_dat <- data.frame("res" = seq(from = 0, to = 20, by = 1),
                       "con" = seq(from = 0, to = 20, by = 1)/4)
box = data.frame(xmin = c(Cons$cy[10] - 0.15, 8.22), xmax = c(Cons$cy[10] + 0.15, 10.45),
                 ymin = c(Cons$cx[10] - 0.15, 11.51), ymax = c(Cons$cx[10] + 0.15, 13.9))

MCR_out <- ggplot() + 
  geom_point(data = Cons, aes(x = cy, y = cx), size = 3) +
  geom_abline(intercept = 0, slope = 1, size = 1.5) +
  geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),  
            colour = "black", fill = "white", alpha = 0, size = 2) +
  geom_segment(aes(x = box$xmin[1], y = box$ymin[1], xend = box$xmin[2], yend = box$ymin[2]), color = "black", size = 2) +
  geom_segment(aes(x = box$xmax[1], y = box$ymax[1], xend = box$xmax[2], yend = box$ymax[2]), color = "black", size = 2) + 
  scale_x_continuous(limits = c(8, 14), expand = c(0, 0)) +
  scale_y_continuous(limits = c(8, 14), expand = c(0, 0)) + 
  labs(x = expression("Consumption of species " * italic(i) * " on resource " * italic(l) * "(" * italic(C[il]) * ")"),
       y = expression("Consumption of species " * italic(j) * " on resource " * italic(l) * "(" * italic(C[jl]) * ")")) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 20),
        plot.margin = margin(t = 12, r = 24, b = 12, l = 12, "pt"),
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 12, r = 0, b = 12, l = 0, "pt")),
        axis.title.y = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 12, b = 0, l = 12, "pt")))

MCR_in <- ggplot() + 
  geom_line(data = Cons_dat, aes(x = res, y = con), size = 2) + 
  scale_x_continuous(limits = c(0, 21), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0)) + 
  labs(x = expression("Density of resource " * italic(l)),
       y = expression(atop("Resource " * italic(l) * " eaten ", "per species " * italic(i) * "per time"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 12, r = 12, b = 12, l = 12, "pt"),
        aspect.ratio = 0.6,
        axis.line = element_line(colour = "black"), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 12, r = 0, b = 0, l = 0, "pt")),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(t = 0, r = 12, b = 0, l = 0, "pt")))

MCR_plot <- ggdraw(MCR_out) +
  draw_label(expression("1:1 line"), x = 0.95, y = 0.9, size = 20) + 
  draw_plot(MCR_in, x = 0.02, y = 0.629, width = 0.5, height = 0.3)
  
ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Ver2/Fig4_MCR.eps", 
       plot = MCR_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
#################################################################################
##### MacArthur's consumer resource model #######################################
#################################################################################

#################################################################################
##### Tilman's consumer resource model ##########################################
#################################################################################
rxR1_f <- function(R1){(rx * R1)/(R1 + Kx1) - D}
rxR2_f <- function(R2){(rx * R2)/(R2 + Kx2) - D}
ryR1_f <- function(R1){(ry * R1)/(R1 + Ky1) - D}
ryR2_f <- function(R2){(ry * R2)/(R2 + Ky2) - D}

rx <- 1.2
ry <- 1
Kx1 <- 12
Kx2 <- 25
Ky1 <- 20
Ky2 <- 12
D <- 0.2

r_dat <- data.frame("R1" = seq(from  = 0, to = 100, length = 1000),
                    "R2" = seq(from  = 0, to = 100, length = 1000)) %>%
  mutate(rxR1 = rxR1_f(R1),
         rxR2 = rxR2_f(R2),
         ryR1 = ryR1_f(R1),
         ryR2 = ryR2_f(R2))

rxR1_p <- r_dat %>% 
  subset(rxR1>0) %>%
  ggplot() +
    geom_line(aes(x = R1, y = rxR1), size = 2, color = "#000099") + 
    geom_abline(intercept = 1.0, slope = 0, size = 1.5, linetype = "dashed") +
    geom_segment(x = -30, y = 0.4, xend = 12, yend = 0.4, size = 1.5, linetype = "dashed") + 
    geom_segment(x = 12, y = 0.4, xend = 12, yend = -0.3, size = 1.5, linetype = "dashed") + 
    draw_label(expression(paste(italic(r[i]), " - D" )), x = -19, y = 1.08, size = 20) +
    draw_label(expression(paste(frac(1, 2), " ",  italic(r[i]), " - D" )), x = -16, y = 0.54, size = 20) +
    scale_x_continuous(limits = c(-30, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0)) + 
    labs(x = expression("Density of resource " * italic(i)),
         y = " ") + 
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size = 18),
          plot.margin = margin(t = 6, r = 12, b = 0, l = 52, "pt"),
          aspect.ratio = 0.6,
          axis.line = element_line(colour = "black"), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 36, r = 0, b = 0, l = 0, "pt")))

rxR2_p <- r_dat %>% 
  subset(rxR2>0) %>%
  ggplot() +
    geom_line(aes(x = R2, y = rxR2), size = 2, color = "#000099") + 
    geom_abline(intercept = 1.0, slope = 0, size = 1.5, linetype = "dashed") +
    geom_segment(x = -30, y = 0.4, xend = 25, yend = 0.4, size = 1.5, linetype = "dashed") + 
    geom_segment(x = 25, y = 0.4, xend = 25, yend = -0.3, size = 1.5, linetype = "dashed") + 
    draw_label(expression(paste(italic(r[i]), " - D" )), x = -19, y = 1.08, size = 20) +
    draw_label(expression(paste(frac(1, 2), " ",  italic(r[i]), " - D" )), x = -16, y = 0.54, size = 20) +
    
    scale_x_continuous(limits = c(-30, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0)) + 
    labs(x = expression("Density of resource " * italic(j)),
         y = "") + 
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size = 18),
          plot.margin = margin(t = 6, r = 12, b = 0, l = 52, "pt"),
          aspect.ratio = 0.6,
          axis.line = element_line(colour = "black"), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 36, r = 0, b = 0, l = 0, "pt")))

ryR1_p <- r_dat %>% 
  subset(ryR1>0) %>%
  ggplot() +
    geom_line(aes(x = R1, y = ryR1), size = 2, color = "#FF6600") + 
    geom_abline(intercept = 0.8, slope = 0, size = 1.5, linetype = "dashed") +
    geom_segment(x = -30, y = 0.3, xend = 20, yend = 0.3, size = 1.5, linetype = "dashed") + 
    geom_segment(x = 20, y = 0.3, xend = 20, yend = -0.3, size = 1.5, linetype = "dashed") + 
    draw_label(expression(paste(italic(r[j]), " - D" )), x = -19, y = .88, size = 20) +
    draw_label(expression(paste(frac(1, 2), " ",  italic(r[j]), " - D" )), x = -16, y = 0.44, size = 20) +
    scale_x_continuous(limits = c(-30, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0)) + 
    labs(x = expression("Density of resource " * italic(i)),
         y = "") + 
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size = 18),
          plot.margin = margin(t = 6, r = 12, b = 0, l = 52, "pt"),
          aspect.ratio = 0.6,
          axis.line = element_line(colour = "black"), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 36, r = 0, b = 0, l = 0, "pt")))

ryR2_p <- r_dat %>% 
  subset(ryR2>0) %>%
  ggplot() +
    geom_line(aes(x = R2, y = ryR2), size = 2, color = "#FF6600") + 
    geom_abline(intercept = 0.8, slope = 0, size = 1.5, linetype = "dashed") +
    geom_segment(x = -30, y = 0.3, xend = 12, yend = 0.3, size = 1.5, linetype = "dashed") + 
    geom_segment(x = 12, y = 0.3, xend = 12, yend = -0.3, size = 1.5, linetype = "dashed") + 
    draw_label(expression(paste(italic(r[j]), " - D" )), x = -19, y = .88, size = 20) +
    draw_label(expression(paste(frac(1, 2), " ",  italic(r[j]), " - D" )), x = -16, y = 0.44, size = 20) +
    scale_x_continuous(limits = c(-30, 100), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0)) + 
    labs(x = expression("Density of resource " * italic(j)),
         y = "") + 
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.title = element_text(size = 18),
          plot.margin = margin(t = 6, r = 12, b = 0, l = 52, "pt"),
          aspect.ratio = 0.6,
          axis.line = element_line(colour = "black"), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(size = 24, face = "bold", margin = margin(t = 36, r = 0, b = 0, l = 0, "pt")))

TCR_plot <- plot_grid(rxR1_p, ryR1_p, rxR2_p, ryR2_p, labels=c("A.", "B.", "C.", "D."), ncol = 2) %>%
  ggdraw() + 
    draw_label(expression(italic("per capita") * " growth rate of species " * italic(i)), angle = 90,  x = 0.04, y = 0.5, size = 24) + 
    draw_label(expression(italic("per capita") * " growth rate of species " * italic(j)), angle = 90,  x = 0.54, y = 0.5, size = 24) + 
    draw_label(expression(italic(R[ii]^"*")), x = 0.171, y = 0.58, size = 20) +  
    draw_label(expression(italic(R[ji]^"*")), x = 0.679, y = 0.58, size = 20) + 
    draw_label(expression(italic(R[ij]^"*")), x = 0.179, y = 0.08, size = 20) + 
    draw_label(expression(italic(R[jj]^"*")), x = 0.672, y = 0.08, size = 20) + 
    draw_label(expression(italic(K[ii])), x = 0.204, y = 0.579, size = 20) +  
    draw_label(expression(italic(K[ji])), x = 0.731, y = 0.579, size = 20) + 
    draw_label(expression(italic(K[ij])), x = 0.247, y = 0.079, size = 20) + 
    draw_label(expression(italic(K[jj])), x = 0.705, y = 0.079, size = 20)

ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/Ver3/Fig5_TCR.pdf", 
       plot = TCR_plot, width = 35, height = 24, units = c("cm"), dpi = 600)
#################################################################################
##### Tilman's consumer resource model ##########################################
#################################################################################
