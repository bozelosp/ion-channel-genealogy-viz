NEURON {
    POINT_PROCESS chirp
    RANGE amp, t1, dur, Finit, beta, ctype
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp) 
    (mV) = (millivolt)
}

PARAMETER {
	amp	=	1.0	(nA)	
	t1 = 	1000 (ms)	
	dur =	20000 (ms)	
	Finit =	0.05 (Hz)	
	beta =	0.24 (Hz/s)	
	ctype =	2	(1)		
}

ASSIGNED {
	i	(nA)			
}


BREAKPOINT {
	LOCAL tc,  pi, uf, mstos
	

	mstos = 1000 (ms/s)	
	uf = 1 (s)			

	pi=3.14159265358979323846
	
	
    if ((t < t1) || (t > t1+dur)) {  
    	i = 0
    } else {
    	tc=(t-t1)/mstos
    	if (ctype == 0) {
    		i = -amp*sin(2*pi*Finit*tc)
    	} else if (ctype == 1) {
    		i = -amp*sin(pi*beta*(tc+Finit/beta)^2)
    	} else if (ctype == 2) {
    		i = -amp*sin(2*pi*Finit*tc*exp(tc*beta*uf))
    	}
    }
}