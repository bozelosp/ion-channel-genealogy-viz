UNITS {
        (nA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
  POINT_PROCESS KCa
  NONSPECIFIC_CURRENT ik
  RANGE dgkbar, egk, ctau, thresh, gk
}


PARAMETER {
  	dgkbar = 0.1 (mS/cm2) <0,1e9>
	egk   = -90 (mV)   
  	ctau =  50   (ms)
  	thresh = -20 (mV)
}


STATE {
        c
}

ASSIGNED {
  	v (mV)
	area (um2) 
 	dgk (uS) 
  	gk (uS) 
	ik (nA)
}

BREAKPOINT {
	SOLVE state METHOD cnexp 
	gk = dgk*c
	ik = gk*(v - egk) 

}

INITIAL {
  	dgk = dgkbar*area*1e-5 
   	gk = 0 
	c = 0
	net_send(0, 1)
}

DERIVATIVE state {     
	c' = -c/ctau
}

NET_RECEIVE (null) {
	if (flag==1) { 
		WATCH (v > thresh) 2 
	}

	if (flag==2) { 
		WATCH (v < thresh) 3 
        
		
   	}
	
	if (flag==3) { 
		c = c+1 
        at_time(t)
		net_send(0, 1)
   	}
}