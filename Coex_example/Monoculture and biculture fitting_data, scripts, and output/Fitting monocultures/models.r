require("odesolve")		# Load the ODE solver
                          


#### A discrete time map 
###
DT_all<-function(N, d){
	## d is the every other day 'instantaneous' dilution
	
#  N = V or c(A, B)
	
out = N*(1-d)
	
return (out)	
	
}




### A hybridsystem of ODEs corresponding to each model

### Each model is a function that accepts three arguments:

# time points for solutions (t),

# initial conditions of variables (ic),
# and parameter values (q)



### M1: Exponential growth 
###
M_1 = function(IC,q){
	
# the ASSUMPTION is that we want the state of the system
	
# at 2 timesteps forward (i.e. 1 timestep is a day and we solve for every two days)
	
# IC = A   (average initial density for this species)
	
# q = c(r, k, d), d is the dilution rate

	
d = q[length(q)]
	
start=DT_all(IC, d)
	
Tfinal = 2 # 2 days 
	
times = c(0,Tfinal)
	
output = lsoda(start, times, CT_1, q)  
#entire output of lsoda	
	
end = output[2,2]
	
return(end)
	
}
	
	



### The ODEs for Model 1 


CT_1 = function(t, IC, q) {
 			
	
# Function requires timepoints for solns, initial conditions, and parameter values
   	
# ic = A	
  
  
A = IC[1]; # Initial conditions for 
	
	
# q = c(r, k)
	r = q[1]; k = q[2]
 
dA = A*r*((k - A)/k)
                  	
	
list(c(dA))		
	
}


### End of Model 1
###

