NEURON {
	POINT_PROCESS Synapse
	NONSPECIFIC_CURRENT i
	RANGE g,gmax, Cdur,Erev
	RANGE tau1,tau2
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
        gmax		= 1150 		(pS)
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

}

STATE {
	C
	O
}


INITIAL {
	C		=	1
	O		=	0
	T		=	0 	(mM)
	ton		=  	-1   (ms)
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
	g =gmax * O
	i = (1e-6) * g * (v - Erev) 
}


KINETIC kstates {
	
	r1 = 1/tau1 * T
	~ C  <-> O	(r1,1/tau2)
	CONSERVE C+O = 1
}


NET_RECEIVE(weight, on, nspike, t0 (ms), tsyn (ms)) {
	INITIAL {
		nspike = 1
		tsyn=t
	}
   	if (flag == 0) {
		nspike = nspike + 1
		if (!on) {
			ton=t
			t0=t
			on=1
		        T=Tmax
			tsyn=t
		}
		net_send(Cdur,nspike)
    	}
	if(flag==nspike){
		T = 0
		on = 0
	}
}