rm(list = ls(all = TRUE))                                                                 
source("models_desolveC.R")  ## load models
require("stats4")
require(twoSppLv)
source("nll_poisson.R") ## load neg log lik

DATA <- read.csv ("Competition Model Fit_DT_plus one.csv") 
modelday <- DATA [,1]

#the counter determines the bottle number +3 (columns before start) for monoculture data to be fit
bottle = 65
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
	for(i in 1:TS_length){out [i+1, ] = M_2(out[i,]+c(1,1), q)}
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

ra= 0.40
rb= 0.31
ka=  3840
kb = 113211 

#ICa = 100
#ICb = 300
#aab = 0.02
#aba = 20 

AAB<- seq (0.01, 0.05, by = 0.01)
ABA<- seq (22, 27, by = 1)
ICA<- c (50, 100, 200, 300)
ICB<- c (100, 200, 300, 400) 

tic<-Sys.time()
count=0
goodfit <- matrix (nrow = length(AAB)*length(ABA)*length(ICA)*length(ICB), ncol = 9)
    for (j in 1:length(AAB))    {
        aab = AAB[j]

    for (v in 1:length (ABA))   {
         aba = ABA [v]

    for (b in 1: length (ICA))  {
        ICa = ICA[b]

    for (u in 1: length(ICB))   {
        ICb = ICB[u]

count = count + 1
                 
      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"= aba)
      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)  
      tryCatch( {fit<-mle(NLL, start = parmsstart, method = "SANN", fixed = parmsfixed)}, error = function (e) { print("test") } ) 
               
goodfit[count, 1] = aab
goodfit[count, 2] = aba
goodfit[count, 3] = ICa
goodfit[count, 4] = ICb
goodfit[count, 5] = coef(fit)[7]
goodfit[count, 6] = coef(fit)[8]
goodfit[count, 7] = coef(fit)[1]
goodfit[count, 8] = coef(fit)[2]
goodfit[count, 9] = AIC(fit)

print(count)

        }
    }
  }
}

toc1<-Sys.time()
goodfit = data.frame(goodfit)
colnames(goodfit)=c("aab","aba","ICa","ICb","aabest","abaest","ICaest","ICbest","AIC")
goodfitord <- goodfit[order(goodfit$AIC),]


########## Nelder-Mead run # 1

sample = 10
count=0
goodfit1 <- matrix (nrow = sample, ncol = 9)
    for (j in 1:sample)    {
        aab = goodfitord$aabest[j]
        aba = goodfitord$abaest[j]
        ICa = goodfitord$ICaest[j]
        ICb = goodfitord$ICbest[j]

parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"=aba)
parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)  
try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)

goodfit1[j, 1] = aab
goodfit1[j, 2] = aba
goodfit1[j, 3] = ICa
goodfit1[j, 4] = ICb
goodfit1[j, 5] = coef(fit)[7]
goodfit1[j, 6] = coef(fit)[8]
goodfit1[j, 7] = coef(fit)[1]
goodfit1[j, 8] = coef(fit)[2]
goodfit1[j, 9] = AIC(fit)

print(j)
}



goodfit1 = data.frame(goodfit1)
colnames(goodfit1)= c("aab","aba","ICa","ICb","aabest","abaest","ICaest","ICbest","AIC")

########## Nelder-Mead run # 2

goodfit2 <- matrix (nrow = sample, ncol = 9)
    for (j in 1:sample)    {
        aab = goodfit1$aabest[j]
        aba = goodfit1$abaest[j]
        ICa = goodfit1$ICaest[j]
        ICb = goodfit1$ICbest[j]

parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"=aba)
parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)  
try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)

goodfit2[j, 1] = aab
goodfit2[j, 2] = aba
goodfit2[j, 3] = ICa
goodfit2[j, 4] = ICb
goodfit2[j, 5] = coef(fit)[7]
goodfit2[j, 6] = coef(fit)[8]
goodfit2[j, 7] = coef(fit)[1]
goodfit2[j, 8] = coef(fit)[2]
goodfit2[j, 9] = AIC(fit)

print(j)
}

##########   Nelder-Mead run # 3

goodfit2 = data.frame(goodfit2)
colnames(goodfit2)=c("aab","aba","ICa","ICb","aabest","abaest","ICaest","ICbest","AIC")

goodfit3 <- matrix (nrow = sample, ncol = 9)
    for (j in 1:sample)    {
        aab = goodfit2$aabest[j]
        aba = goodfit2$abaest[j]
        ICa = goodfit2$ICaest[j]
        ICb = goodfit2$ICbest[j]

parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"=aba)
parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)  
try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)

goodfit3[j, 1] = aab
goodfit3[j, 2] = aba
goodfit3[j, 3] = ICa
goodfit3[j, 4] = ICb
goodfit3[j, 5] = coef(fit)[7]
goodfit3[j, 6] = coef(fit)[8]
goodfit3[j, 7] = coef(fit)[1]
goodfit3[j, 8] = coef(fit)[2]
goodfit3[j, 9] = AIC(fit)

print(j)
}

##########   Nelder-Mead run # 4

goodfit3 = data.frame(goodfit3)
colnames(goodfit3)=c("aab","aba","ICa","ICb","aabest","abaest","ICaest","ICbest","AIC")

goodfit4 <- matrix (nrow = sample, ncol = 9)
    for (j in 1:sample)    {
        aab = goodfit3$aabest[j]
        aba = goodfit3$abaest[j]
        ICa = goodfit3$ICaest[j]
        ICb = goodfit3$ICbest[j]

parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"=aba)
parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)  
try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)

goodfit4[j, 1] = aab
goodfit4[j, 2] = aba
goodfit4[j, 3] = ICa
goodfit4[j, 4] = ICb
goodfit4[j, 5] = coef(fit)[7]
goodfit4[j, 6] = coef(fit)[8]
goodfit4[j, 7] = coef(fit)[1]
goodfit4[j, 8] = coef(fit)[2]
goodfit4[j, 9] = AIC(fit)

print(j)
}

##########   Nelder-Mead run # 5

goodfit4 = data.frame(goodfit4)
colnames(goodfit4)=c("aab","aba","ICa","ICb","aabest","abaest","ICaest","ICbest","AIC")

goodfit5 <- matrix (nrow = sample, ncol = 9)
    for (j in 1:sample)    {
        aab = goodfit4$aabest[j]
        aba = goodfit4$abaest[j]
        ICa = goodfit4$ICaest[j]
        ICb = goodfit4$ICbest[j]

parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"=aba)
parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)  
try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)

goodfit5[j, 1] = aab
goodfit5[j, 2] = aba
goodfit5[j, 3] = ICa
goodfit5[j, 4] = ICb
goodfit5[j, 5] = coef(fit)[7]
goodfit5[j, 6] = coef(fit)[8]
goodfit5[j, 7] = coef(fit)[1]
goodfit5[j, 8] = coef(fit)[2]
goodfit5[j, 9] = AIC(fit)

print(j)
}

goodfit5 = data.frame(goodfit5)
colnames(goodfit5)=c("aab","aba","ICa","ICb","aabest","abaest","ICaest","ICbest","AIC")
goodfit5ord <- goodfit5[order(goodfit5$AIC),]

toc2<-Sys.time()
### Timer
tictoc1 = tic-toc1
tictoc1
tictoc2=tic-toc2
tictoc2
############ GRAPHING

finalparms = c(ICa = goodfit5ord$ICaest[1], ICb = goodfit5ord$ICbest[1],
 ra = ra, rb = rb, ka = ka, kb = kb, aab = goodfit5ord$aabest[1], aba = goodfit5ord$abaest[1], d = 0.1)
        
finalparms = as.numeric(finalparms)

 fish<- function(x) {
	t = x+1
  parms = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[7], finalparms[8], finalparms[9]))
	out = numeric(2 * t)
	dim(out) = c(t, 2)
	out[1, ]=NA
  out[2, 1] = as.numeric(finalparms[1])
  out[2, 2] = as.numeric(finalparms[2])
	for(i in 3:t){out[i, ] = M_2( out[i-1,]+c(1,1), parms)}        
	return(out)
  }


par(mar = c(5, 5, 4, 5))
plot(1:24, fish(23)[ ,1], pch = NA, xlab = "time (days)",
  ylab = "Density (cells/mL)", xaxt = 'n', ylim = c(0, max(DATA[1:y, countera], na.rm=T) + 1000))

lines(1:24,fish(23)[,1],col="red")
points(DATA[1:y,1]+1,DATA[1:y,countera],col="red")

par(new = T)
plot(1:24, fish(23)[,2], pch = NA, xlab = "time (days)", ylab = "Density (cells/mL)",
  xaxt = 'n', yaxt = 'n', ylim = c(0, max(DATA[1:y, counterb], na.rm = T) + 1000))

aty <- seq(par("yaxp")[1], par("yaxp")[2], (par("yaxp")[2] - par("yaxp")[1])/ par("yaxp")[3])

points(DATA[1:y, 1] + 1, DATA[1:y, counterb], col="green")
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

AIC = goodfit5ord$AIC[1]
AIC