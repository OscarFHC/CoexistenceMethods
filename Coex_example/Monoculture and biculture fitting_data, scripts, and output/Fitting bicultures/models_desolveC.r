require("deSolve")

#Load the ODE solver


### A discrete time map 
######



DT_all <- function(N, d){
	
## d is the every other day 'instantaneous' dilution
	
#  N = A or c(A, B)  
  
out = N * (1 - d)
	
return (out)	
	
}
                      



### M2: Lotka-Volterra Competition between species
###

M_2 = function(IC, q){
	
## We want the state of the system after 2 days.
	
# IC = c(starting A, starting B)	
	
# q = c(ra, rb, ka, kb, aab, aba, d) -- the LAST value in the parameter vector is dilution rate
	
  

d = q[length(q)]
	
log.ic = log(DT_all(IC, d))
	
start = c(A = log.ic [1], B = log.ic [2])
  
params = c (ra = q[1], rb = q[2], ka = q[3], kb = q[4], aab = q[5], aba = q[6])
  
times=seq(0,2,by=1)
  
sol <- ode(
y=start,
 times = times, parms = params,
 func="M_2", dllname="twoSppLv", initfunc="initmod", nout=1, rtol=10^-8, atol=10^-8)
output<-exp(sol[3,2:3])

	
return(round(output))
}	