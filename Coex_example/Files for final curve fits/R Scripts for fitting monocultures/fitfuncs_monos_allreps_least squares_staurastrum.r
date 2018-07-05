rm(list = ls(all = TRUE))         
source("models.R")  ## load models
require("stats4")

##the counter determines the bottle number +3 (columns before start) for monoculture data to be fit
species = 7
counter = species+3
y=17 #number of data points to use
jump=3

  DATA<-read.csv("Competition monocultures.csv")

modelday<-DATA[1:y,1]

#########generate deterministic ts given parms  for single species logistic growth
Predict<-function(IC, r, k, d){

	#	return 1-dim 22 step ts with IC as first row
	TS_length = 22 # length of ts -1 for 0 condition
 	q = c(r, k, d)
  out = numeric((1+TS_length))
	dim(out) = c(1,1+TS_length)
	out[1] = IC
	for(i in 1:TS_length){out[i+1]= M_1(out[i],q)}
	return(out[modelday])
	}


##### locally defined LS function, uses DATA
 LS = function(par,d) {
		##likelihood function for a run of data (with either one)  
    IC = par[1]
    r = par[2]
    k = par[3]
    d = d 
      predict = Predict(IC,r,k,d)
      least = sum(c((DATA[1:y,jump*species+1]-predict)^2, (DATA[1:y,jump*species+3]-predict)^2))                                  
    return(least)
	}


	####################################################
	##### now we MINIMIZE the Least Squares
	####################################################

maximum =c(max(DATA[1:y,jump*species+1]),    max(DATA[1:y,jump*species+3]))               
        
IC= mean(c(DATA[1,jump*species+1], DATA[1,jump*species+3]))                   
meank<-mean(maximum)

s=seq(0.8, 1.2, by=0.1)
rchange=0.1  
count=0
goodfit<-matrix(nrow=5,ncol=7)
for (b in 1:5){
      for (c in 1:7) {
      p=s[c]                                                                                          
      r = rchange * b
      k = p * meank
      d = 0.1       
      tryCatch ( {fit<- optim (fn = LS, par = c( IC, r, k), "d"=d, method = "SANN")      }, error = function (e) { print("test") } )   
      goodfit[b,c]= fit$value 
      count = count + 1
      print(count)
              }
  }

rstart<-which(goodfit==min(goodfit,na.rm=TRUE),arr.ind=TRUE)[1]*rchange
kstart<-s[which(goodfit==min(goodfit,na.rm=TRUE),arr.ind=TRUE)[2]]*meank
       
fit<- optim (fn = LS, par = c(IC,  rstart, kstart),  d = 0.1, method = "SANN")

#######################
####################### Graphing the model and data
#######################

finalparms=c(IC= fit$par[1], r=fit$par[2], k=fit$par[3], d=0.1)

fish<-function(x){
	t=x+1
  parms = as.numeric( c(finalparms[2], finalparms[3], finalparms[4]) )
	out = numeric(t)
	out[1]=NA
  out[2]=finalparms[1]
	for(i in 3:t){out[i]= M_1(out[i-1],parms)}
	return(out)
  }

plotmax<-max(c(max(DATA[1:y,jump*species+1]),  max(DATA[1:y,jump*species+3])))

plot(1:24,fish(23),pch=NA,xlab="time (days)",ylab="Density (cells/mL)",xaxt='n',ylim=c(0,max(plotmax)+0.1*plotmax))


lines(1:24,fish(23))
axis(1,1:24,seq(2,48,2))
points(DATA[1:y,1]+1,DATA[1:y,jump*species+1],col=2)

points(DATA[1:y,1]+1,DATA[1:y,jump*species+3],col=4)
points(DATA[1:y,1]+1,rep(50,y),pch=3)
text(19,plotmax+0.1*plotmax,"Staurastrum punctulatum")


mean<-mean(c(DATA[1:y,jump*species+1], DATA[1:y,jump*species+3]) )                 

RSQ = function(IC,r,k,d) {
		predict = Predict(IC,r,k,d)
    RSS<-sum(c((DATA[1:y,jump*species+1]-predict)^2, (DATA[1:y,jump*species+3]-predict)^2))                                    
    rsq<-1-(RSS/sum(c((DATA[1:y,jump*species+1]-mean)^2, (DATA[1:y,jump*species+3]-mean)^2 )))                               
    print(rsq)
	}


RSQ(as.numeric(finalparms[1]),as.numeric(finalparms[2]),as.numeric(finalparms[3]),as.numeric(finalparms[4]))
print(finalparms)
############################################################