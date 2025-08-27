NEURON {				
	POINT_PROCESS GradGABAa		
	POINTER PreActiv
	RANGE C, g, gmax, Erev, timestep, tau, g_inf, vref, thres, slop
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Erev	 = -70	 (mV)			
	gmax	 = 0	 (uS)			
	tau 	 = 3 	 (ms)

	g_inf 	 = 0 	 (umho)
	vref     = 1    (mV)                    
	thres    = -45  (mV) 
	slop	 = 0.2			 	
}

ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C				
	PreActiv			
}

INITIAL {
	C = 0
	g = 0
}

BREAKPOINT { LOCAL delta_g
	  
	C = 1/(1+exp(4*slop*(thres-PreActiv)/vref))  
					   	     

	g_inf = gmax* C	
	delta_g = (dt/2) * (g_inf-g)/tau      
					      

	g = g + delta_g	

	

	i = (g*(v-Erev))
	
	
	
	
	
	
	
}