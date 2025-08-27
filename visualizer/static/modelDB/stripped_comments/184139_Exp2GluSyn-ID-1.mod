NEURON {
	POINT_PROCESS GluSyn
	RANGE ntar, e, i, mg, mgshift, tau1, tau2, tau3, tauD, tauF, f, Pb, Gnmda, Gampa
	NONSPECIFIC_CURRENT i, inmda, iampa
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	tau1 = 0.5 (ms)			  
	tau2 = 4 (ms) <1e-9,1e9>  
	tau3 = 42 (ms) <1e-9,1e9> 
	e = 0.0	(mV)
	mg = 1 (mM) 			  
	sh = 0 (mV)
	ntar = .3	(1) < 0, 1 >  
	
	f = 2 (1) < 0, 1e9 >    
	tauF = 100 (ms) < 1e-9, 1e9 > 
	tauD = 500 (ms) < 1e-9, 1e9 > 
	Pb = 0.3 (1) < 0, 1 >	
}

ASSIGNED {
	v (mV)
	i (nA)
	inmda (nA)
	iampa (nA)
	Gnmda (uS)
	factor
}

STATE {
	Anmda (uS)
	Bnmda (uS)
	Aampa (uS)
	Gampa (uS)
}

INITIAL {
	LOCAL tp
	if (tau2/tau3 > .9999) {
		tau2 = .9999*tau3
	}
	Aampa = 0
	Gampa = 0
	Anmda = 0
	Bnmda = 0
	tp = (tau2*tau3)/(tau3 - tau2) * log(tau3/tau2)
	factor = -exp(-tp/tau2) + exp(-tp/tau3)
	factor = 1/factor
	mgblock(v)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	Gnmda = Bnmda - Anmda
	inmda = Gnmda*(v - e)*mgblock(v)
	iampa = Gampa *(v - e)
	i = inmda + iampa
}

DERIVATIVE state {
	Aampa' = - Aampa/tau1
	Gampa' = Aampa/tau1 - Gampa/tau1	
	
	Anmda' = -Anmda/tau2
	Bnmda' = -Bnmda/tau3
}

NET_RECEIVE(weight (uS), F, D, tsyn (ms)) {
	INITIAL {
		
		F = 1
		D = 1
		tsyn = t
    }
	F = 1 + (F-1)*exp(-(t - tsyn)/tauF)
	D = 1 - (1-D)*exp(-(t - tsyn)/tauD)
	tsyn = t 	
	
	Aampa = Aampa + weight*exp(1)*F*D*Pb
	Anmda = Anmda + ntar*weight*factor*F*D*Pb
	Bnmda = Bnmda + ntar*weight*factor*F*D*Pb

	F = F + f  
	D = D*(1 - Pb*F)	
}


FUNCTION mgblock(v) {
	mgblock = 1 / (1 + exp(0.062 * -(v-sh)) * (mg / 3.57))	
}