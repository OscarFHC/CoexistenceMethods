##### This is the R scripts for making example plots in the project of multiple species coexistence methods/models
###################################################################################################################
library(deSolve)
library(ggplot2)
library(cowplot)
library(tidyverse)
#################################################################################

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

### Creating the one sp LV model
LV_Mono <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(1 - a11*x)
    return(list(c(dx)))
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

### Creating predicted value
LV_x <- as.data.frame(ode(func = LV_Mono, y = Init_Chl, parms = c(fit_x$par), times = Time, hmax = 0.1, maxsteps = 10000))
LV_y <- as.data.frame(ode(func = LV_Mono, y = Init_Clo, parms = c(fit_y$par), times = Time, hmax = 0.1, maxsteps = 10000))

### Now fix the parameter values estimated from mono-cultures to estimate parameters in bi-culture
### Organize data; use the mean of the three replicates as the empirical data for parameter fitting
### I scaled down the actual numbers (/1000000) to avoid possible numerical errors when fitting parameters
dat_Bi <- dat_Bi_ori %>%
  group_by(tp) %>%
  summarize(Chl_m = mean(chl/1000000, na.rm = TRUE),
            Clo_m = mean(clo/1000000, na.rm = TRUE),
            Chl_se = NonParam(chl/1000000),
            Clo_se = NonParam(clo/1000000))

### prepare two species LV model
LV_Bi <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*rx*(Kx - a11*x - a12*y)
    dy = y*ry*(Ky - a21*x - a22*y)
    return(list(c(dx, dy)))
  })
}

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

### Creating predicted value
LV_Bi <- as.data.frame(ode(func = LV_Bi, y = Init, 
                           parms = c(fit_x$par, fit_y$par, fit_Bi$par), times = Time, hmax = 0.1, maxsteps = 10000))

### Calculating R2 for each species in mono-culture and bi culture
SST_x <- sum((dat_Mono$Chl_m - mean(dat_Mono$Chl_m))^2)
SSR_1_x <- sum((dat_Mono$Chl_m - LV_x$x)^2)
R2_1_x <- 1 - SSR_1_x/SST_x

SST_y <- sum((dat_Mono$Clo_m - mean(dat_Mono$Clo_m))^2)
SSR_1_y <- sum((dat_Mono$Clo_m - LV_y$x)^2)
R2_1_y <- 1 - SSR_1_y/SST_y

SST_Bi <- sum((dat_Bi$Chl_m - mean(dat_Bi$Chl_m))^2 + (dat_Bi$Clo_m - mean(dat_Bi$Clo_m))^2)
SSR_1_Bi <- sum((dat_Bi$Chl_m - LV_Bi$x)^2 + (dat_Bi$Clo_m - LV_Bi$y)^2)
R2_1_Bi <- 1 - SSR_1_Bi/SST_Bi
### Overlay empirical data with fitted LV model predictions
X1 <- ggplot() + 
  geom_point(data = dat_Mono, aes(x = tp, y = Chl_m), color = "black") + 
  geom_errorbar(data = dat_Mono, aes(x = tp, ymin = Chl_m - Chl_se, ymax = Chl_m + Chl_se), width = 0.1) +
  geom_line(data = LV_x, aes(x = time, y = x), color = "black") +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 2.5)) + 
  labs(title = expression(paste("Mono-culture of ", italic(Chlorella))), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 42), "pt"),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14))

Y1 <- ggplot() + 
  geom_point(data = dat_Mono, aes(x = tp, y = Clo_m), color = "blue") + 
  geom_errorbar(data = dat_Mono, aes(x = tp, ymin = Clo_m - Clo_se, ymax = Clo_m + Clo_se), width = 0.1, color = "blue") +
  geom_line(data = LV_y, aes(x = time, y = x), color = "blue") +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 2.5)) + 
  labs(title = expression(paste("Mono-culture of ", italic(Closteriopsis))), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14))

Bi_1 <- ggplot() + 
  geom_point(data = dat_Bi, aes(x = tp, y = Chl_m), color = "black") + 
  geom_errorbar(data = dat_Bi, aes(x = tp, ymin = Chl_m - Chl_se, ymax = Chl_m + Clo_se), width = 0.1) +
  geom_line(data = LV_Bi, aes(x = time, y = x), color = "black") +
  geom_point(data = dat_Bi, aes(x = tp, y = Clo_m), color = "blue") + 
  geom_errorbar(data = dat_Bi, aes(x = tp, ymin = Clo_m - Clo_se, ymax = Clo_m + Clo_se), width = 0.1, color = "blue") +
  geom_line(data = LV_Bi, aes(x = time, y = y), color = "blue") + 
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 2.5)) + 
  labs(title = expression(paste("Bi-culture")),
       x = expression(paste("Time")),
       y = "") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 42), "pt"),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 4, r = 0, b = 12, l = 24)))

plot_up <- plot_grid(X1, Y1, labels = c("A.", "B."), ncol = 2, align = 'h')
Plot_Bi1 <- plot_grid(plot_up, Bi_1, labels = c("", "C."), ncol = 1) %>%
  ggdraw() + 
    draw_label(expression(paste("Density ( ", 10^6, "/mL)")), 
               angle = 90, x = 0.04, y = 0.5, size = 18) + 
    draw_label("error bars represent standard error of the mean", 
               x = 0.85, y = 0.02, size = 14)
ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/FigX_Chl_Clo_1.jpeg", 
       plot = Plot_Bi1, width = 35, height = 24, units = c("cm"), dpi = 600)

#################################################################################
#################################################################################
### estimate parameters direcely from two species culture

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

### prepare LV model
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

### Set up initial values (the same as the empirical measurement) and time steps
Init <- c(x = dat_Bi$Chl_m[1], y = dat_Bi$Clo_m[1])
Time <- seq(0, nrow(dat_Bi), by = 1)

nll_Bi <- function (par) {
  y <- ode(func = LV_Bi, y = Init, parms = par, times = Time, hmax = 0.1, maxsteps = 10000)
  -sum(dnorm(x = dat_Bi$Chl_m, mean = y[,2], sd = 0.1, log = TRUE) + 
       dnorm(x = dat_Bi$Clo_m, mean = y[,3], sd = 0.1, log = TRUE))
}

### Use "optim" function to minimize negative log likelihood (maximize likelihood)
fit_Bi_2 <- optim(par = c(rx = 0.45, ry = 0.45,
                          Kx = 2, Ky = 2, 
                          a11 = 1, a12 = 1, a21 = 1, a22 = 1), fn = nll_Bi)

### Creating predicted value
Init_Chl <- c(x = dat_Mono$Chl_m[1])
Init_Clo <- c(x = dat_Mono$Clo_m[1])
Time <- seq(1, nrow(dat_Mono), by = 1)
LV_x_2 <- as.data.frame(ode(func = LV_Mono, y = Init_Chl, parms = c(fit_Bi_2$par), times = Time, hmax = 0.1, maxsteps = 10000))
LV_y_2 <- as.data.frame(ode(func = LV_Mono, y = Init_Clo, parms = c(fit_Bi_2$par), times = Time, hmax = 0.1, maxsteps = 10000))

Init <- c(x = dat_Bi$Chl_m[1], y = dat_Bi$Clo_m[1])
Time <- seq(1, nrow(dat_Bi), by = 1)
LV_Bi_2 <- as.data.frame(ode(func = LV_Bi, y = Init, parms = fit_Bi_2$par, times = Time, hmax = 0.1, maxsteps = 10000))

### Calculating R2 for each species in mono-culture and bi culture
SST_x <- sum((dat_Mono$Chl_m - mean(dat_Mono$Chl_m))^2)
SSR_2_x <- sum((dat_Mono$Chl_m - LV_x_2$x)^2)
R2_2_x <- 1 - SSR_2_x/SST_x

SST_y <- sum((dat_Mono$Clo_m - mean(dat_Mono$Clo_m))^2)
SSR_2_y <- sum((dat_Mono$Clo_m - LV_y_2$x)^2)
R2_2_y <- 1 - SSR_2_y/SST_y

SST_Bi <- sum((dat_Bi$Chl_m - mean(dat_Bi$Chl_m))^2 + (dat_Bi$Clo_m - mean(dat_Bi$Clo_m))^2)
SSR_2_Bi <- sum((dat_Bi$Chl_m - LV_Bi_2$x)^2 + (dat_Bi$Clo_m - LV_Bi_2$y)^2)
R2_2_Bi <- 1 - SSR_2_Bi/SST_Bi
### Overlay empirical data with fitted LV model predictions
X2 <- ggplot() + 
  geom_point(data = dat_Mono, aes(x = tp, y = Chl_m), color = "black") + 
  geom_errorbar(data = dat_Mono, aes(x = tp, ymin = Chl_m - Chl_se, ymax = Chl_m + Chl_se), width = 0.1) +
  geom_line(data = LV_x_2, aes(x = time, y = x), color = "black") +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 2.5)) + 
  labs(title = expression(paste("Mono-culture of ", italic(Chlorella))), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 42), "pt"),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14))

Y2 <- ggplot() + 
  geom_point(data = dat_Mono, aes(x = tp, y = Clo_m), color = "blue") + 
  geom_errorbar(data = dat_Mono, aes(x = tp, ymin = Clo_m - Clo_se, ymax = Clo_m + Clo_se), width = 0.1, color = "blue") +
  geom_line(data = LV_y_2, aes(x = time, y = x), color = "blue") +
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 2.5)) + 
  labs(title = expression(paste("Mono-culture of ", italic(Closteriopsis))), 
       x = "", #expression(paste("Time")),
       y = "") + #expression(paste("Density ( ", 10^6, "/mL)"))) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        plot.title = element_text(size = 14, hjust = 0.5, vjust = 0.2),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14))

Bi_2 <- ggplot() + 
  geom_point(data = dat_Bi, aes(x = tp, y = Chl_m), color = "black") + 
  geom_errorbar(data = dat_Bi, aes(x = tp, ymin = Chl_m - Chl_se, ymax = Chl_m + Clo_se), width = 0.1) +
  geom_line(data = LV_Bi_2, aes(x = time, y = x), color = "black") +
  geom_point(data = dat_Bi, aes(x = tp, y = Clo_m), color = "blue") + 
  geom_errorbar(data = dat_Bi, aes(x = tp, ymin = Clo_m - Clo_se, ymax = Clo_m + Clo_se), width = 0.1, color = "blue") +
  geom_line(data = LV_Bi_2, aes(x = time, y = y), color = "blue") + 
  scale_y_continuous(expand = c(0, 0.1), limits = c(0, 2.5)) + 
  labs(title = expression(paste("Bi-culture")),
       x = expression(paste("Time")),
       y = "") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 42), "pt"),
        plot.title = element_text(size = 14, hjust = 0.5),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14), 
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 4, r = 0, b = 12, l = 24)))

plot_up <- plot_grid(X2, Y2, labels = c("A.", "B."), ncol = 2, align = 'h')
Plot_Bi1 <- plot_grid(plot_up, Bi_2, labels = c("", "C."), ncol = 1) %>%
  ggdraw() + 
  draw_label(expression(paste("Density ( ", 10^6, "/mL)")), 
             angle = 90, x = 0.04, y = 0.5, size = 18) + 
  draw_label("error bars represent standard error of the mean", 
             x = 0.85, y = 0.02, size = 14)
ggsave(filename = "D:/Manuscript/CoexistenceMethods_Figs/FigX_Chl_Clo_2.jpeg", 
       plot = Plot_Bi1, width = 35, height = 24, units = c("cm"), dpi = 600)
