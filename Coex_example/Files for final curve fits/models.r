
#rm(list = ls(all = TRUE)) 	# Clear the memory
require("odesolve")		# Load the ODE solver
                          
#########################################################################
#### A discrete time map 
#########################################################################

DT_all<-function(N, d){
	## d is the every other day 'instantaneous' dilution
	#  N = V or c(A, B)
	out = N*(1-d)
	return (out)	
	}


#########################################################################
#### A hybridsystem of ODEs corresponding to each model
#########################################################################

#### Each model is a function that accepts three arguments:
# time points for solutions (t),
# initial conditions of variables (ic),
# and parameter values (q)

###########################################################################
### M1: Exponential growth 
###########################################################################

M_1 = function(IC,q){
	# the ASSUMPTION is that we want the state of the system
	# at 2 timesteps forward (i.e. 1 timestep is a day and we solve for every two days)
	# IC = A   (average initial density for this species)
	# q = c(r, k, d), d is the dilution rate

	d = q[length(q)]
	start=DT_all(IC, d)
	Tfinal = 2 # 2 days 
	times = c(0,Tfinal)
	output = lsoda(start, times, CT_1, q)  #entire output of lsoda	
	end = output[2,2]
	return(end)
	}
	
	
######## The ODEs for Model 1 
CT_1 = function(t, IC, q) { 			
	# Function requires timepoints for solns, initial conditions, and parameter values
   	# ic = A	
  
  A = IC[1]; # Initial conditions for 
	
	# q = c(r, k)
	r = q[1]; k = q[2]
  		
	dA = A*r*((k - A)/k)
                  	
	list(c(dA))		
	}

###########################################################################
### End of Model 1
###########################################################################
                      
###########################################################################
### M2: Lotka-Volterra Competition between species
###########################################################################

M_2 = function(ic,q){
	## the ASSUMPTION is that we want the state of the system
	# ic = c(A, B)	
	# q = c(ra, rb, ka, kb ... , d) -- the LAST value in the parameter vector is dilution rate
	d = q[length(q)]
	
  start=DT_all(ic, d)
	Tfinal = 2 # days between dilutions
	times = c(0,Tfinal)
	output = lsoda(start, times, CT_2, q)  #entire output of lsoda	
  end = output[2,2:3]
	return(end)
	}
            
######## The ODEs for Model 1 
CT_2 = function(t, ic, q) { 			
	# Function requires timepoints for solns, initial conditions, and parameter values
  # ic = c(A, B)	
	# q = c(ra,rb,ka,kb,aab,aba,d)

  A = ic[1]; B = ic[2]                 	
	ra = q[1]; rb = q[2]; ka = q[3]; kb = q[4]; aab=q[5]; aba=q[6]
	
	dA = A*ra*(1 -((A+aab*B)/ka))	
	dB = B*rb*(1-((B+aba*A)/kb))
	list(c(dA,dB))		
	}

###########################################################################
### End of Model 2
###########################################################################