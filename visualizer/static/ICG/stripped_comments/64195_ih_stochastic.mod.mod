UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 			(mV)
	eh  		(mV) 		
	ghdbar=0.00015 		(S/cm2)		
	gamma=680e-15		(S)		
	seed
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT i
	RANGE ghdbar,N,N_open,eh
}

STATE {
	l
}

ASSIGNED {
	i (mA/cm2)
	dt			(ms)
	area			(um2)
	N			
	N_open			
}

INITIAL {								
	N=abs(((1e-8*area*ghdbar)/gamma)+0.5)				
	N_open=abs(N*alpha(v)/(alpha(v)+beta(v)))
	l=0								

	set_seed(seed)
}


BREAKPOINT {
	SOLVE states METHOD cnexp					
	i = ((N_open*gamma)/(1e-8*area))*(v-eh)			
}


FUNCTION alpha(v(mV)) {
	
	alpha = 6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)			
}

FUNCTION beta(v(mV)) {
	
	beta = 193*exp(v/33.1)			
}


DERIVATIVE states {     						
	l' =  l			    					
	noise()
}

PROCEDURE noise() {
	LOCAL N_close,N_open_merk,N_close_merk,a,b,prob_open,prob_close
	a=alpha(v)
	b=beta(v)
	N_open_merk=N_open
	N_close_merk=N-N_open

	
	
	prob_open=1-exp(-dt*a/1000)					
	FROM ii=1 TO N_close_merk {
		if (scop_random()<= prob_open)	{			
			N_open=N_open+1
		}
	}

	
	prob_close=1-exp(-dt*b/1000)	
	FROM ii=1 TO N_open_merk {
		if (scop_random()<= prob_close)	{		
			N_open=N_open-1
		}
	}
}