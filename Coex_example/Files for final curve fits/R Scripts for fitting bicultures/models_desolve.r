rm(list = ls(all = TRUE)) 	# Clear the memory
require("deSolve")		# Load the ODE solver
                          
#########################################################################
#### A discrete time map 
#########################################################################

DT_all <- function(N, d){
	## d is the every other day 'instantaneous' dilution
	#  N = A or c(A, B)
	out = N * (1 - d)
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
	# IC = Starting A	
	# q = c(r, k, d) -- the LAST value in the parameter vector is dilution rate
	
  d = q[length(q)]
	Tfinal = 2 # 2 days 
	times = seq(0,100, by=1)
	#times = c(0, Tfinal)
	start=c(A = DT_all(IC, d))
	parms=c(r = q[1], k = q[2], d = q[3])
  output = ode(y = start, times = times, func = CT_1, parms = parms)
	end = as.numeric(output [2,2])
	return(end)
	}
	
	
######## The ODEs for Model 1 
CT_1 = function(t, start, parms) { 			
  with(as.list(c(start, parms)),{
	# Function requires timepoints for solns, state variables (A), and parameter values
			
	dA = A * r * ((k - A) / k)
                  	
	list(dA)		
		}    )
	}


###########################################################################
### End of Model 1
###########################################################################
                      
###########################################################################
### M2: Lotka-Volterra Competition between species
###########################################################################

M_2 = function(IC, q){
	## We want the state of the system after 2 days.
	# IC = c(starting A, starting B)	
	# q = c(ra, rb, ka, kb, aab, aba, d) -- the LAST value in the parameter vector is dilution rate
	
  d = q[length(q)]
	start = c(A = DT_all(IC, d)[1], B = DT_all (IC, d)[2])

	Tfinal = 2 # days between dilutions
	times = c(0, Tfinal) 
	#times = seq(0,100, by=1)
	#start = IC
  parms = c (ra = q[1], rb = q[2], ka = q[3], kb = q[4], aab = q[5], aba = q[6], d =q [7])
  output = ode(y = start, times = times, func = CT_2, parms = parms, rtol=10^-8, atol=10^-8)  #entire output of ode
  end = as.numeric(output[2, 2:3])   
	return(end)
	#return(output)
	}
            
######## The ODEs for Model 2 
CT_2 = function(t, start, q) { 
  with(as.list(c(start, q)), {			
  
	# Function requires timepoints for solutions, state variables (A and B), and parameters
  # state = c(A, B)	
	# parameters = c(ra, rb, ka, kb, aab, aba, d)
	
	dA = A * ra * ((ka - A - aab * B) / ka)	
	dB = B * rb * ((kb - B - aba * A) / kb)
	list(c(dA, dB))		
	}    )
	}