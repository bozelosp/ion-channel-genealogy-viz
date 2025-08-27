NEURON {
	POINT_PROCESS gaba
	RANGE tau1, tau2
	RANGE erev, g, i, q
    RANGE damod, maxMod, level, max2, lev2
	
	NONSPECIFIC_CURRENT i
}


UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	erev    = -60.0     (mV)
	tau1    =   0.5     (ms)    
    tau2    =   7.5     (ms)    
    q       =   2       
    
    damod       = 0
    maxMod      = 1
    max2        = 1
    level       = 0
    lev2        = 0
}


ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
}


STATE {
	A (uS)
	B (uS)
}


INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A       = 0
	B       = 0
	tp      = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor  = -exp(-tp/tau1) + exp(-tp/tau2)
	factor  = 1/factor
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	
	g = (B - A) * modulation(maxMod,max2,level,lev2)
	i = g * (v - erev)
}


DERIVATIVE state {
	A' = -A/tau1*q
	B' = -B/tau2*q
}


NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
}


FUNCTION modulation(m1,m2,l1,l2) {
    
    
    modulation = 1 + damod * ( (m1-1)*l1 + (m2-1)*l2 )
    if (modulation < 0) {
        modulation = 0
    } 
}