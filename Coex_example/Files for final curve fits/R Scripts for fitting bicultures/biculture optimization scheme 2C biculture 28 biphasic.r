rm(list = ls(all = TRUE))                                                               
source("models_desolveC.R")  ## load models
require("stats4")
require(twoSppLv)
source("nll_poisson.R") ## load neg log lik

DATA <- read.csv ("Competition Model Fit_DT_plus one.csv")
modelday <- DATA [,1]

#the counter determines the bottle number +3 (columns before start) for monoculture data to be fit
bottle = 108
biculture = bottle - 24
countera = bottle + 3 + biculture - 1
counterb = countera  + 1
y = 13          #overall time series length, e.g. for plotting
ya = 13         #number of points included for time series of species A
yb = 13         #number of days included for time series of species B

######### Generate deterministic time series given parameterss for two species LV competition                              

PredictC <- function(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b){
	#	Return 2-dim 22 step ts with IC as first row

  TS_length = 22 # length of time series (TS)
	q1 = c(ra, rb, ka, kb, aab1, aba1, d)
	q2 = c(ra, rb, ka, kb, aab2, aba2, d)
	out = numeric(2 * (1 + TS_length))  # length of TS +1 for 0 condition
	dim(out) = c(1 + TS_length, 2)
	out[1,1] = ICa
	out[1,2] = ICb
	for(i in 1:TS_length)
  { if (i < b) {out [i+1, ] = M_2(out[i,] + c(1,1), q1)}
  else { out [i+1, ] = M_2(out[i,] + c(1,1), q2) }
  }
  return(out[modelday, ])    # returns densities for points in the TS for which data is available
 	}

#### locally defined logl function, uses DATA
NLL = function (ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b) {

		## likelihood function for a run of data with two species

    predict1 = PredictC (ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b)
    negLogLik = nll(DATA[1:ya, countera], predict1[1:ya, 1]) + nll(DATA[1:yb, counterb], predict1[1:yb, 2])
		return(negLogLik)
	}

	#################


ra= 0.78
rb= 0.37
ka= 1797
kb= 36697

ICa =  200
ICb =  500
aab1 =  -1.6
aba1 = 2.5
aab2 = -0.6
aba2 = 0
b = 16

finalparms = c(ICa = ICa, ICb = ICb, ra = ra, rb = rb, ka = ka, kb = kb, aab1 = aab1, aba1 = aba1, aab2 = aab2, aba2 = aba2, d = 0.1, b = b)

finalparms = as.numeric(finalparms)

fish<- function(x) {
	t = x+1
  parms1 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[7], finalparms[8], finalparms[11]))
  parms2 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[9], finalparms[10], finalparms[11]))
	out = numeric(2 * t)
	dim(out) = c(t, 2)
	b = as.numeric(finalparms[12])
	out[1, ]=NA
  out[2, 1] = as.numeric(finalparms[1])
  out[2, 2] = as.numeric(finalparms[2])
  for(i in 3:t){
  if (i < b) {out [i, ] = M_2(out[i-1,] + c(1,1), parms1)}
  else { out [i, ] = M_2(out[i-1,] + c(1,1), parms2) }
  }
	return(out)
  }

par(mar = c(5, 5, 4, 5))
plot(1:24, fish(23)[ ,1], pch = NA, xlab = "time (days)",
  ylab = "Density (cells/mL)", xaxt = 'n', ylim = c(0,
  max(DATA[1:y, countera: counterb], na.rm=T) + 1000))

lines(1:24,fish(23)[,1],col="red")
points(DATA[1:y,1]+1,DATA[1:y,countera],col="red")

par(new = T)
plot(1:24, fish(23)[,2], pch = NA, xlab = "time (days)", ylab = "Density (cells/mL)",
  xaxt = 'n', yaxt = 'n', ylim = c(0, max(DATA[1:y, counterb], na.rm = T) + 1000))

points(DATA[1:y, 1] + 1, DATA[1:y, counterb], col="green")
aty <- seq(par("yaxp")[1], par("yaxp")[2], (par("yaxp")[2] - par("yaxp")[1])/ par("yaxp")[3])
lines(1:24,fish(23)[,2], col="green")
axis(4, at = aty, labels = format(aty, scientific = FALSE))
text(24, max(DATA[1:y, counterb] ,na.rm = T) +500, toString(bottle))


parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
fit<-mle(NLL, start = parmsstart, method = "SANN", fixed = parmsfixed)


ICa = as.numeric(coef(fit)[1])
ICb = as.numeric(coef(fit)[2])
aab1 = as.numeric(coef(fit)[7])
aba1 =  as.numeric(coef(fit)[8])
aab2 = as.numeric(coef(fit)[9])
aba2 =  as.numeric(coef(fit)[10])
#
#parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
#parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
#fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed)



ICA <- c(200, 400)
ICB <- c(500, 1000)
AAB1 <-  c(-1, -2, -4, -8)
ABA1 <-  c(2, 3, 4)
AAB2 <- c(-0.5, -0.6, -0.7)
ABA2 <- c(-0.4, 0, 0.4)
B <- c(14, 15, 16, 17, 18)


tic<-Sys.time()
count=0
goodfit <- matrix (nrow = length(AAB1)*length(ABA1)*length(AAB2)*length(ABA2)*length(ICA)*length(ICB)*length(B), ncol = 14)
    for (j in 1:length(AAB1))    {
        aab1 = AAB1[j]

    for (v in 1:length (ABA1))   {
         aba1 = ABA1 [v]
         
    for (g in 1:length(AAB2))    {
        aab2 = AAB2[g]

    for (f in 1:length(ABA2))    {
        aba2 = ABA2[f]

    for (b in 1: length (ICA))  {
        ICa = ICA[b]

    for (u in 1: length(ICB))   {
        ICb = ICB[u]
        
    for (e in 1: length(B))   {
        b = B[e]

count = count + 1

      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
      tryCatch( {fit<-mle(NLL, start = parmsstart, method = "SANN", fixed = parmsfixed)}, error = function (e) { print("test") } )

goodfit[count, 1] = aab1
goodfit[count, 2] = aba1
goodfit[count, 3] = aab2
goodfit[count, 4] = aba2
goodfit[count, 5] = ICa
goodfit[count, 6] = ICb
goodfit[count, 7] = coef(fit)[7]
goodfit[count, 8] = coef(fit)[8]
goodfit[count, 9] = coef(fit)[9]
goodfit[count, 10] = coef(fit)[10]
goodfit[count, 11] = coef(fit)[1]
goodfit[count, 12] = coef(fit)[2]
goodfit[count, 13] = coef(fit)[12]
goodfit[count, 14] = AIC(fit)

print(count)

        }}}
        }
    }
  }
}

toc1<-Sys.time()
goodfit = data.frame(goodfit)
colnames(goodfit)=c("aab1","aba1","aab2","aba2","ICa","ICb","aabest1","abaest1","aabest2","abaest2","ICaest","ICbest", "b", "AIC")
goodfitord <- goodfit[order(goodfit$AIC),]


########### Nelder-Mead run # 1
#
#sample = 10
#count=0
#goodfit1 <- matrix (nrow = sample, ncol = 14)
#    for (j in 1:sample)    {
#        aab1 = goodfitord$aabest1[j]
#        aba1 = goodfitord$abaest1[j]
#        aab2 = goodfitord$aabest2[j]
#        aba2 = goodfitord$abaest2[j]
#        ICa = goodfitord$ICaest[j]
#        ICb = goodfitord$ICbest[j]
#        b = goodfitord$b[j]
#
#      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
#      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
#      try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)
#
#goodfit1[count, 1] = aab1
#goodfit1[count, 2] = aba1
#goodfit1[count, 3] = aab2
#goodfit1[count, 4] = aba2
#goodfit1[count, 5] = ICa
#goodfit1[count, 6] = ICb
#goodfit1[count, 7] = coef(fit)[7]
#goodfit1[count, 8] = coef(fit)[8]
#goodfit1[count, 9] = coef(fit)[9]
#goodfit1[count, 10] = coef(fit)[10]
#goodfit1[count, 11] = coef(fit)[1]
#goodfit1[count, 12] = coef(fit)[2]
#goodfit1[count, 13] = coef(fit)[12]
#goodfit1[count, 14] = AIC(fit)
#
#print(j)
#}
#
#
#
#goodfit1 = data.frame(goodfit1)
#colnames(goodfit1)=c("aab1","aba1","aab2","aba2","ICa","ICb","aabest1","abaest1","aabest2","abaest2","ICaest","ICbest", "b", "AIC")
#
########### Nelder-Mead run # 2
#
#goodfit2 <- matrix (nrow = sample, ncol = 14)
#    for (j in 1:sample)    {
#        aab1 = goodfit1$aabest1[j]
#        aba1 = goodfit1$abaest1[j]
#        aab2 = goodfit1$aabest2[j]
#        aba2 = goodfit1$abaest2[j]
#        ICa = goodfit1$ICaest[j]
#        ICb = goodfit1$ICbest[j]
#        b = goodfit1$b[j]
#
#      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
#      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
#      try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)
#
#goodfit2[count, 1] = aab1
#goodfit2[count, 2] = aba1
#goodfit2[count, 3] = aab2
#goodfit2[count, 4] = aba2
#goodfit2[count, 5] = ICa
#goodfit2[count, 6] = ICb
#goodfit2[count, 7] = coef(fit)[7]
#goodfit2[count, 8] = coef(fit)[8]
#goodfit2[count, 9] = coef(fit)[9]
#goodfit2[count, 10] = coef(fit)[10]
#goodfit2[count, 11] = coef(fit)[1]
#goodfit2[count, 12] = coef(fit)[2]
#goodfit2[count, 13] = coef(fit)[12]
#goodfit2[count, 14] = AIC(fit)
#
#print(j)
#}
#
###########   Nelder-Mead run # 3
#
#goodfit2 = data.frame(goodfit2)
#colnames(goodfit2)=c("aab1","aba1","aab2","aba2","ICa","ICb","aabest1","abaest1","aabest2","abaest2","ICaest","ICbest", "b", "AIC")
#
#goodfit3 <- matrix (nrow = sample, ncol = 14)
#    for (j in 1:sample)    {
#        aab1 = goodfit2$aabest1[j]
#        aba1 = goodfit2$abaest1[j]
#        aab2 = goodfit2$aabest2[j]
#        aba2 = goodfit2$abaest2[j]
#        ICa = goodfit2$ICaest[j]
#        ICb = goodfit2$ICbest[j]
#        b = goodfit2$b[j]
#
#      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
#      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
#      try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)
#
#goodfit3[count, 1] = aab1
#goodfit3[count, 2] = aba1
#goodfit3[count, 3] = aab2
#goodfit3[count, 4] = aba2
#goodfit3[count, 5] = ICa
#goodfit3[count, 6] = ICb
#goodfit3[count, 7] = coef(fit)[7]
#goodfit3[count, 8] = coef(fit)[8]
#goodfit3[count, 9] = coef(fit)[9]
#goodfit3[count, 10] = coef(fit)[10]
#goodfit3[count, 11] = coef(fit)[1]
#goodfit3[count, 12] = coef(fit)[2]
#goodfit3[count, 13] = coef(fit)[12]
#goodfit3[count, 14] = AIC(fit)
#
#
#print(j)
#}
#
###########   Nelder-Mead run # 4
#
#goodfit3 = data.frame(goodfit3)
#colnames(goodfit3)=c("aab1","aba1","aab2","aba2","ICa","ICb","aabest1","abaest1","aabest2","abaest2","ICaest","ICbest", "b", "AIC")
#
#goodfit4 <- matrix (nrow = sample, ncol = 14)
#    for (j in 1:sample)    {
#    
#        aab1 = goodfit3$aabest1[j]
#        aba1 = goodfit3$abaest1[j]
#        aab2 = goodfit3$aabest2[j]
#        aba2 = goodfit3$abaest2[j]
#        ICa = goodfit3$ICaest[j]
#        ICb = goodfit3$ICbest[j]
#        b = goodfit3$b[j]
#        
#      parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab1" = aab1, "aba1"= aba1, "aab2" = aab2, "aba2"= aba2)
#      parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1, "b" = b)
#      try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)
#
#goodfit4[count, 1] = aab1
#goodfit4[count, 2] = aba1
#goodfit4[count, 3] = aab2
#goodfit4[count, 4] = aba2
#goodfit4[count, 5] = ICa
#goodfit4[count, 6] = ICb
#goodfit4[count, 7] = coef(fit)[7]
#goodfit4[count, 8] = coef(fit)[8]
#goodfit4[count, 9] = coef(fit)[9]
#goodfit4[count, 10] = coef(fit)[10]
#goodfit4[count, 11] = coef(fit)[1]
#goodfit4[count, 12] = coef(fit)[2]
#goodfit4[count, 13] = coef(fit)[12]
#goodfit4[count, 14] = AIC(fit)
#
#
#print(j)
#}
#
###########   Nelder-Mead run # 5
#
#goodfit4 = data.frame(goodfit4)
#colnames(goodfit4)=c("aab1","aba1","aab2","aba2","ICa","ICb","aabest1","abaest1","aabest2","abaest2","ICaest","ICbest", "b", "AIC")
#
#goodfit5 <- matrix (nrow = sample, ncol = 14)
#    for (j in 1:sample)    {
#    
#        aab1 =  goodfit4$aabest1[j]
#        aba1 =  goodfit4$abaest1[j]
#        aab2 =  goodfit4$aabest2[j]
#        aba2 =  goodfit4$abaest2[j]
#        ICa =   goodfit4$ICaest[j]
#        ICb =   goodfit4$ICbest[j]
#        b = goodfit3$b[j]
#
#parmsstart = list ("ICa" = ICa, "ICb" = ICb, "aab" = aab, "aba"=aba)
#parmsfixed = list ("ra" = ra, "rb" = rb, "ka" = ka, "kb" = kb, "d" = 0.1)
#try(fit<-mle(NLL, start = parmsstart, method = "Nelder-Mead", fixed = parmsfixed), silent = T)
#
#goodfit5[count, 1] = aab1
#goodfit5[count, 2] = aba1
#goodfit5[count, 3] = aab2
#goodfit5[count, 4] = aba2
#goodfit5[count, 5] = ICa
#goodfit5[count, 6] = ICb
#goodfit5[count, 7] = coef(fit)[7]
#goodfit5[count, 8] = coef(fit)[8]
#goodfit5[count, 9] = coef(fit)[9]
#goodfit5[count, 10] = coef(fit)[10]
#goodfit5[count, 11] = coef(fit)[1]
#goodfit5[count, 12] = coef(fit)[2]
#goodfit5[count, 13] = coef(fit)[12]
#goodfit5[count, 14] = AIC(fit)
#
#print(j)
#}
#
#
#goodfit5 = data.frame(goodfit5)
#colnames(goodfit5)=c("aab1","aba1","aab2","aba2","ICa","ICb","aabest1","abaest1","aabest2","abaest2","ICaest","ICbest", "b", "AIC")

goodfitord <- goodfit[order(goodfit5$AIC),]

toc2<-Sys.time()
### Timer
tictoc1 = tic-toc1
tictoc1
tictoc2=tic-toc2
tictoc2
############ GRAPHING

finalparms = c(ICa = goodfitord$ICaest[1], ICb = goodfitord$ICbest[1],
 ra = ra, rb = rb, ka = ka, kb = kb, aab1 = goodfitord$aabest1[1], aba1 = goodfitord$abaest1[1],
 aab2 = goodfitord$aabest2[1], aba2 = goodfitord$abaest2[1], d = 0.1, b = goodfitord$b[1])

finalparms = as.numeric(finalparms)

fish<- function(x) {
	t = x+1
	b = as.numeric(finalparms[12])
  parms1 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[7], finalparms[8], finalparms[11]))
  parms2 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[9], finalparms[10], finalparms[11]))
	out = numeric(2 * t)
	dim(out) = c(t, 2)
	out[1, ]=NA
  out[2, 1] = as.numeric(finalparms[1])
  out[2, 2] = as.numeric(finalparms[2])
  for(i in 3:t){
  if (i < b) {out [i, ] = M_2(out[i-1,] + c(1,1), parms1)}
  else { out [i, ] = M_2(out[i-1,] + c(1,1), parms2) }
  }
	return(out)
  }

par(mar = c(5, 5, 4, 5))
plot(1:24, fish(23)[ ,1], pch = NA, xlab = "time (days)",
  ylab = "Density (cells/mL)", xaxt = 'n', ylim = c(0,
  max(DATA[1:y, countera: counterb], na.rm=T) + 1000))

lines(1:24,fish(23)[,1],col="red")
points(DATA[1:y,1]+1,DATA[1:y,countera],col="red")

par(new = T)
plot(1:24, fish(23)[,2], pch = NA, xlab = "time (days)", ylab = "Density (cells/mL)",
  xaxt = 'n', yaxt = 'n', ylim = c(0, max(DATA[1:y, counterb], na.rm = T) + 1000))

points(DATA[1:y, 1] + 1, DATA[1:y, counterb], col="green")
aty <- seq(par("yaxp")[1], par("yaxp")[2], (par("yaxp")[2] - par("yaxp")[1])/ par("yaxp")[3])
lines(1:24,fish(23)[,2], col="green")
axis(4, at = aty, labels = format(aty, scientific = FALSE))
text(24, max(DATA[1:y, counterb] ,na.rm = T) +500, toString(bottle))




mean<-mean(c(DATA[1:ya,countera], DATA[1:yb,counterb]))
meana<- mean(DATA[1:ya,countera])
meanb<-mean(DATA[1:ya,counterb]) 

RSQt = function(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b) {
		predictf = PredictC(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b)
    RSS<-sum(c((DATA[1:ya, countera] - predictf[1:ya, 1])^2, (DATA[1:yb, counterb] - predictf[1:yb, 2])^2))
    rsq<-1-(RSS/sum(c((DATA[1:ya, countera]-mean)^2,(DATA[1:yb,counterb]- mean)^2)))
    print(rsq)
	}
	
RSQs = function(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b) {
		predictf = PredictC(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b)
    RSS<-sum(c((DATA[1:ya, countera] - predictf[1:ya, 1])^2, (DATA[1:yb, counterb] - predictf[1:yb, 2])^2))
    rsq<-1-(RSS/sum(c((DATA[1:ya, countera]-meana)^2,(DATA[1:yb,counterb]- meanb)^2)))
    print(rsq)
	}



RSQt(as.numeric(finalparms[1]),as.numeric(finalparms[2]),as.numeric(finalparms[3]),
  as.numeric(finalparms[4]),as.numeric(finalparms[5]),as.numeric(finalparms[6]),
  as.numeric(finalparms[7]),as.numeric(finalparms[8]),as.numeric(finalparms[9]),
  as.numeric(finalparms[10]),as.numeric(finalparms[11]),as.numeric(finalparms[12]))
  
RSQs(as.numeric(finalparms[1]),as.numeric(finalparms[2]),as.numeric(finalparms[3]),
  as.numeric(finalparms[4]),as.numeric(finalparms[5]),as.numeric(finalparms[6]),
  as.numeric(finalparms[7]),as.numeric(finalparms[8]),as.numeric(finalparms[9]),
  as.numeric(finalparms[10]),as.numeric(finalparms[11]),as.numeric(finalparms[12]))
  
  print(finalparms)

AIC = goodfitord$AIC[1]
AIC

#write.csv(goodfitord, "goodfitord 28 biphasic.csv")
b