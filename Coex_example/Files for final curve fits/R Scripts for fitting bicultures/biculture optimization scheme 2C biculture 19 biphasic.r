rm(list = ls(all = TRUE))                                                                 
source("models_desolve.R")  ## load models
require("stats4")
#require(twoSppLv)
source("nll_poisson.R") ## load neg log lik

DATA <- read.csv ("Competition Model Fit_DT_plus one.csv")
modelday <- DATA [,1]

#the counter determines the bottle number +3 (columns before start) for monoculture data to be fit
bottle = 79
biculture = bottle - 24
countera = bottle + 3 + biculture - 1
counterb = countera  + 1
y= 13          #overall time series length, e.g. for plotting
ya= 13         #number of points included for time series of species A
yb= 13         #number of days included for time series of species B


Predict <- function(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b){
#Return 2-dim 22 step ts with IC as first row

  TS_length = 22 # length of time series (TS)
	q1 = c(ra, rb, ka, kb, aab1, aba1, d)
	q2 = c(ra, rb, ka, kb, aab2, aba2, d)
	out = numeric(2 * (1 + TS_length))  # length of TS +1 for 0 condition
	dim(out) = c(1 + TS_length, 2)
	out[1,1] = ICa
	out[1,2] = ICb
	for(i in 1:TS_length){
  if (i < b) {out [i+1,] =  M_2(out[i,] + c(1,1), q1)  }
  else {  out [i+1, ] = M_2(out[i, ] + c(1,1), q2)}
  }
  return(out[modelday, ])    # returns densities for points in the TS for which data is available
 	}

#### locally defined logl function, uses DATA
NLL = function (ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b) {

		## likelihood function for a run of data with two species

    predict1 = Predict (ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b)
    negLogLik = nll(DATA[1:ya, countera], predict1[1:ya, 1]) + nll(DATA[1:yb, counterb], predict1[1:yb, 2])
		return(negLogLik)
	}

	####################################################
	##### now we MINIMIZE negaive log likelihood
	####################################################

#Fixed parameters need to be entered in manually from the estimated monoculture parameters

ra=   0.31
rb=   0.79
ka=   113211
kb =  78798
#
#ICa = 2000
#ICb = 200
#aab1 = 1
#aba1 =  2
#aab2 = 0.4
#aba2 =  -6
#b=7

ICA<- c(2000)
ICB<- c(200)
AAB1<- c(1, 4)
ABA1<- c(2, 8)
AAB2<- c(0.4, 0.8)
ABA2<- c(-3, -9)
B <- c(7, 8, 9, 10, 11)
#
#
#finalparms = c(ICa = ICa, ICb = ICb, ra = ra, rb = rb, ka = ka, kb = kb, aab1 = aab1, aba1 = aba1, aab2 = aab2, aba2 = aba2, d = 0.1, b = b)
#
#finalparms = as.numeric(finalparms)
#
#fish<- function(x) {
#	t = x+1
#  parms1 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[7], finalparms[8], finalparms[11]))
#  parms2 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[9], finalparms[10], finalparms[11]))
#	out = numeric(2 * t)
#	dim(out) = c(t, 2)
#	out[1, ]=NA
#  out[2, 1] = as.numeric(finalparms[1])
#  out[2, 2] = as.numeric(finalparms[2])
#  for(i in 3:t){
#  if (i < b) { out [i, ] =  M_2(out[i-1,] + c(1,1), parms1) }
#  else { out [i, ] = M_2(out[i-1, ] + c(1,1), parms2)}
#  }
#	return(out)
#  }
#
#par(mar = c(5, 5, 4, 5))
#plot(1:24, fish(23)[ ,1], pch = NA, xlab = "time (days)",
#  ylab = "Density (cells/mL)", xaxt = 'n', ylim = c(0,
#  max(DATA[1:y, countera: counterb], na.rm=T) + 1000))
#
#lines(1:24,fish(23)[,1],col="red")
#points(DATA[1:y,1]+1,DATA[1:y,countera],col="red")
#
#par(new = T)
#plot(1:24, fish(23)[,2], pch = NA, xlab = "time (days)", ylab = "Density (cells/mL)",
#  xaxt = 'n', yaxt = 'n', ylim = c(0, max(DATA[1:y, counterb], na.rm = T) + 1000))
#
#points(DATA[1:y, 1] + 1, DATA[1:y, counterb], col="green")
#aty <- seq(par("yaxp")[1], par("yaxp")[2], (par("yaxp")[2] - par("yaxp")[1])/ par("yaxp")[3])
#lines(1:24,fish(23)[,2], col="green")
#axis(4, at = aty, labels = format(aty, scientific = FALSE))
#text(24, max(DATA[1:y, counterb] ,na.rm = T) +500, toString(bottle))
#atx <- seq(par("xaxp")[1], par("xaxp")[2], (par("xaxp")[2] - par("xaxp")[1])/ par("xaxp")[3])
#axis(1,at=atx,labels = c(10,20,30,40))
#
#########################

tic<-Sys.time()
count=0
goodfit <- matrix (nrow = length(AAB1)*length(ABA1)*length(AAB2)*length(ABA2)*length(ICA)*length(ICB)*length(B), ncol = 14)
    for (a in 1:length(AAB1))    {
        aab1 = AAB1[a]

    for (b in 1:length (ABA1))   {
         aba1 = ABA1 [b]
         
    for (j in 1:length (AAB2))   {
         aab2 = AAB2 [j]
    
    for (d in 1:length (ABA2))   {
         aba2 = ABA2 [d]

    for (e in 1: length (ICA))  {
        ICa = ICA[e]

    for (f in 1: length(ICB))   {
        ICb = ICB[f]
        
    for (g in 1: length(B))    {
        b = B[g]

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
goodfit[count, 13] = b
goodfit[count, 14] = AIC(fit)

print(count) }
          }
        }
      }
    }
  }
}

toc1<-Sys.time()
goodfit = data.frame(goodfit)
colnames(goodfit)=c("aab1","aba1","aab2","aba2","ICa","ICb","aab1est","aba1est","aab2est","aba2est","ICaest","ICbest", "b", "AIC")
goodfitord <- goodfit[order(goodfit$AIC),]

################################################

### Timer
tictoc1 = tic-toc1
tictoc1

        aab1 = goodfitord$aab1est[j]
        aba1 = goodfitord$aba1est[j]
        aab2 = goodfitord$aab2est[j]
        aba2 = goodfitord$aba2est[j]
        ICa = goodfitord$ICaest[j]
        ICb = goodfitord$ICbest[j] 
        b = goodfitord$b[j]


finalparms = c(ICa = ICa, ICb = ICb, ra = ra, rb = rb, ka = ka, kb = kb, aab1 = aab1, aba1 = aba1, aab2 = aab2, aba2 = aba2, d = 0.1, b = b)

finalparms = as.numeric(finalparms)

fish<- function(x) {
	t = x+1
  parms1 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[7], finalparms[8], finalparms[11]))
  parms2 = as.numeric(c(finalparms[3], finalparms[4], finalparms[5], finalparms[6], finalparms[9], finalparms[10], finalparms[11]))
	out = numeric(2 * t)
	dim(out) = c(t, 2)
	out[1, ]=NA
  out[2, 1] = as.numeric(finalparms[1])
  out[2, 2] = as.numeric(finalparms[2])
  for(i in 3:t){
  if (i < b) { out [i, ] =  M_2(out[i-1,] + c(1,1), parms1) }
  else { out [i, ] = M_2(out[i-1, ] + c(1,1), parms2)}
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
atx <- seq(par("xaxp")[1], par("xaxp")[2], (par("xaxp")[2] - par("xaxp")[1])/ par("xaxp")[3])
axis(1,at=atx,labels = c(10,20,30,40))


mean<-mean(c(DATA[1:ya,countera], DATA[1:yb,counterb]))
RSQ = function(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b) {
		predictf = Predict(ICa, ICb, ra, rb, ka, kb, aab1, aba1, aab2, aba2, d, b)
    RSS<-sum(c((DATA[1:ya, countera] - predictf[1:ya, 1])^2, (DATA[1:yb, counterb] - predictf[1:yb, 2])^2))
    rsq<-1-(RSS/sum(c((DATA[1:ya, countera]-mean)^2,(DATA[1:yb,counterb]- mean)^2)))
    print(rsq)
	}


RSQ(as.numeric(finalparms[1]),as.numeric(finalparms[2]),as.numeric(finalparms[3]),
  as.numeric(finalparms[4]),as.numeric(finalparms[5]),as.numeric(finalparms[6]),
  as.numeric(finalparms[7]),as.numeric(finalparms[8]),as.numeric(finalparms[9]),
  as.numeric(finalparms[10]),as.numeric(finalparms[11]),as.numeric(finalparms[12]))
  print(finalparms)

AIC = goodfit5ord$AIC[1]
AIC

write.csv(goodfitord, "Bottle 79 goodfitord.csv")
