rm(list = ls(all = TRUE))                                                                  fit$value
source("models_desolveC.R")  ## load models
require("stats4")
require(twoSppLv)
source("nll_poisson.R") ## load neg log lik

DATA <- read.csv ("Competition Model Fit_DT_plus one.csv") 
modelday <- DATA [,1]

#the counter determines the bottle number +3 (columns before start) for monoculture data to be fit
bottle = 51
biculture = bottle - 24
countera = bottle + 3 + biculture - 1
counterb = countera  + 1
y= 13          #overall time series length, e.g. for plotting
ya= 13         #number of points included for time series of species A
yb= 13         #number of days included for time series of species B

######### Generate deterministic time series given parameterss for two species LV competition

PredictC <- function(ICa, ICb, ra, rb, ka, kb, aab, aba, d){
	#	Return 2-dim 22 step ts with IC as first row
  
  TS_length = 22 # length of time series (TS)
	q = c(ra, rb, ka, kb, aab, aba, d)
	out = numeric(2 * (1 + TS_length))  # length of TS +1 for 0 condition
	dim(out) = c(1 + TS_length, 2)
	out[1,1] = ICa
	out[1,2] = ICb
	for(i in 1:TS_length){out [i+1, ] = M_2(out[i,] + c(1,1), q)}
  return(out[modelday, ])    # returns densities for points in the TS for which data is available
 	}

#### locally defined logl function, uses DATA
NLL = function (ICa, ICb, ra, rb, ka, kb, aab, aba, d) {
  
		## likelihood function for a run of data with two species
		
    predict1 = PredictC (ICa, ICb, ra, rb, ka, kb, aab, aba, d)
    negLogLik = nll(DATA[1:ya, countera], predict1[1:ya, 1]) + nll(DATA[1:yb, counterb], predict1[1:yb, 2])
		return(negLogLik)
	}	
	

	####################################################
	##### now we MINIMIZE negaive log likelihood
	####################################################

#Fixed parameters need to be entered in manually from the estimated monoculture parameters


ra= 0.55
rb= 0.31
ka= 1184418
kb =   113211
ICa = 1
ICb = 300
aab = -10
aba =  0

preds<-   PredictC(ICa, ICb, ra, rb, ka, kb, aab, aba, d = 0.1)
matplot(preds)
matpoints(DATA[1:y, countera:counterb])



      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"= aba)
      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)
      fit<-mle(NLL, start = parmsstart, method = "SANN", fixed = parmsfixed)
      

      ICa = as.numeric(coef(fit)[1])
      ICb = as.numeric(coef(fit)[2])
      aab = as.numeric(coef(fit)[7])
      aba =  as.numeric(coef(fit)[8])

      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" =aab, "aba"= aba)
      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)
      fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed)
      
      
      ICa = as.numeric(coef(fit)[1])
      ICb = as.numeric(coef(fit)[2])
      aab = as.numeric(coef(fit)[7])
      aba =  as.numeric(coef(fit)[8])

      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" =aab, "aba"= aba)
      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)
      fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed)
      
            ICa = as.numeric(coef(fit)[1])
      ICb = as.numeric(coef(fit)[2])
      aab = as.numeric(coef(fit)[7])
      aba =  as.numeric(coef(fit)[8])

      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" =aab, "aba"= aba)
      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)
      fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed)
      
#      
#      ICa = as.numeric(coef(fit)[1])
#      ICb = as.numeric(coef(fit)[2])
#      aab = as.numeric(coef(fit)[7])
#      aba =  as.numeric(coef(fit)[8])
#
#      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" =aab, "aba"= aba)
#      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)
#      fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed)
      
#finalparms<-c(100000, 200, 0.45, 0.35, ka, kb, 0, 0.09, 0.01)

finalparms = c(ICa= as.numeric(coef(fit)[1]), ICb = as.numeric(coef(fit)[2]), ra = as.numeric(coef(fit)[3]), rb = as.numeric(coef(fit)[4]),  ka = as.numeric(coef(fit)[5]), kb = as.numeric(coef(fit)[6]), aab = as.numeric(coef(fit)[7]), aba = as.numeric(coef(fit)[8]), d = 0.1)
finalparms = as.numeric(finalparms)

 fish<- function(x) {
	t = x+1
  parms = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[7], finalparms[8], finalparms[9]))
	out = numeric(2 * t)
	dim(out) = c(t, 2)
	out[1, ]=NA
  out[2, 1] = as.numeric(finalparms[1])
  out[2, 2] = as.numeric(finalparms[2])
	for(i in 3:t){out[i, ] = M_2(c(out[i-1, 1], out[i-1, 2]), parms)}
	return(out)
  }

par(mar = c(5, 5, 4, 5))
plot(1:24, fish(23)[ ,1], pch = NA, xlab = "time (days)",
  ylab = "Density (cells/mL)", xaxt = 'n', ylim = c(min(DATA[1:y, countera:counterb], na.rm=T)- 1000,
  max(DATA[1:y, countera: counterb], na.rm=T) + 1000))

lines(1:24,fish(23)[,1],col="red")
points(DATA[1:y,1]+1,DATA[1:y,countera],col="red")

par(new = T)
plot(1:24, fish(23)[,2], pch = NA, xlab = "time (days)", ylab = "Density (cells/mL)",
  xaxt = 'n', yaxt = 'n', ylim = c(min(DATA[1:y, counterb], na.rm = T), max(DATA[1:y, counterb], na.rm = T) + 1000))

points(DATA[1:y, 1] + 1, DATA[1:y, counterb], col="green")
aty <- seq(par("yaxp")[1], par("yaxp")[2], (par("yaxp")[2] - par("yaxp")[1])/ par("yaxp")[3])
lines(1:24,fish(23)[,2], col="green")
axis(4, at = aty, labels = format(aty, scientific = FALSE))
text(24, max(DATA[1:y, counterb] ,na.rm = T) +500, toString(bottle))

atx <- seq(par("xaxp")[1], par("xaxp")[2], (par("xaxp")[2] - par("xaxp")[1])/ par("xaxp")[3])
axis(1,at=atx,labels = c(10,20,30,40))


mean<-mean(c(DATA[1:ya,countera], DATA[1:yb,counterb]))
RSQ = function(ICa, ICb, ra, rb, ka, kb, aab, aba, d) {
		predictf = PredictC(ICa, ICb, ra, rb, ka, kb, aab, aba, d)
    RSS<-sum(c((DATA[1:ya, countera] - predictf[1:ya, 1])^2, (DATA[1:yb, counterb] - predictf[1:yb, 2])^2))
    rsq<-1-(RSS/sum(c((DATA[1:ya, countera]-mean)^2,(DATA[1:yb,counterb]- mean)^2)))
    print(rsq)
	}


RSQ(as.numeric(finalparms[1]),as.numeric(finalparms[2]),as.numeric(finalparms[3]),
  as.numeric(finalparms[4]),as.numeric(finalparms[5]),as.numeric(finalparms[6]),
  as.numeric(finalparms[7]),as.numeric(finalparms[8]),as.numeric(finalparms[9]))
  print(finalparms)

