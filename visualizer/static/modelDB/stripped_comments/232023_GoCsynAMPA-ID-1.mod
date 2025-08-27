NEURON {
	POINT_PROCESS GoCsynAMPA
	NONSPECIFIC_CURRENT i
	RANGE g,gmax, Cdur,Erev
	RANGE tau_1, tau_2
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(nS) = (nanosiemens)
	(um) = (micrometer)
	PI   = (pi)(1)
}

PARAMETER {
	
	tau1		= .4		(/ms/mM)
        tau2		= 3		(/ms)
        gmax		= 2000 		(pS)
	Cdur		= 0.3		(ms)
	Erev		= 0		(mV)
	Tmax		= 1  (mM)
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(pS)		
	r1		(/ms)
	T		(mM)
	ton		(ms)
	w		(1)		

}

STATE {
	C
	o
}


INITIAL {
	C		=	1
	o		=	0
	T		=       0	(mM)
	ton		=  	-1	(ms)
	w		=	1	(1)
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
	g = w*gmax*o
	i = (1e-6) * g * (v - Erev) 
}


KINETIC kstates {
	
	r1 = 1/tau1 * T
	~ C  <-> o	(r1,1/tau2)
	CONSERVE C+o = 1
}


NET_RECEIVE(weight, on, nspike, t0 (ms), tsyn (ms)) {
	INITIAL {
		nspike = 1
		tsyn=t
	}
   	if (flag == 0) {
		w =weight
		nspike = nspike + 1
		if (!on) {
			ton = t
			t0 = t
			on = 1
		        T = Tmax
			tsyn = t
		}
		net_send(Cdur,nspike)
    	}
	if(flag==nspike){
		T = 0
		on = 0
	}
}