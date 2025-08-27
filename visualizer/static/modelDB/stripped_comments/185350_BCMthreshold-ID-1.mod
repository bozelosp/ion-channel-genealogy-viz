NEURON {
	POINT_PROCESS BCMthreshold
	GLOBAL d0, p0, scount0, scounttau, alpha, alpha_scount
	GLOBAL d, p, tspike
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	
	d0		
	p0		

	scounttau  		
	alpha			
	
	scount0 		
	boltzman = 0	

}

ASSIGNED {
	v (mV)
	flagOLD
	
	alpha_scount	
	d			
	p			
	boltzfactor 	
	pf				
	output			
	tspike (ms)
}
STATE {
	scount			
}
INITIAL {
	
	tspike = -100000
	



	net_send(0, 1)
	flagOLD = 1
	scount = scount0
	d = d0		
	p = p0
	alpha_scount = alpha*scount
	boltzfactor = exp( - 1.0 / scounttau)
	if(boltzman == 0) {pf = 1.0 / (1 - boltzfactor)} else {pf = 1.0}
	output = 0
}

BREAKPOINT { 
	
	scount = scount  + output/pf
	SOLVE state METHOD cnexp								
  	pf = pf * boltzfactor + 1.0
  	alpha_scount = alpha * scount							
	if (alpha_scount > 0) {p = p0/alpha_scount} else {p = p}
	d = d0*alpha_scount
    output = 0
}

DERIVATIVE state{
	scount' = -scount / pf
}







NET_RECEIVE(w) {		
	INITIAL {w=0}




	if (flag == 0) {
		output = 0 			
	} else if (flag == 2) { 	
		tspike = t				
		output = 1				
		
	} else { 

		








		WATCH (v > 0) 2 
						
	}
}