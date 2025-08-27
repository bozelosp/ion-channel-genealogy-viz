NEURON {
	POINT_PROCESS Exp3NMDA
	NONSPECIFIC_CURRENT i
	RANGE tau1, tau2, v0_tau2, st_tau2, tau3, tauV, e, i, gVI, st_gVD, v0_gVD, Mg, K0, delta, wf
	GLOBAL inf
	THREADSAFE
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
	(S)  = (siemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(J)  = (joules)
}

PARAMETER {

	tau1 = 8.8		(ms)	<1e-9,1e9>	
	tau2 = 500		(ms)
	v0_tau2 = 161.11	(mV)	
	st_tau2 =0.30342 (ms/mV)	

	tauV = 7		(ms)	<1e-9,1e9>	
							
							
							
							
							
	st_gVD = 0.007	(1/mV)	
	v0_gVD = -100	(mV)	
	gVI = 1			(uS)	
	Q10 = 1.52				
	T0 = 26			(degC)	
	celsius 		(degC)	

	Mg = 1			(mM)	
	K0 = 4.1		(mM)	
	delta = 0.8 	(1)		

	e = -0.7		(mV)	
}

CONSTANT {
	T = 273.16	(degC)
	F = 9.648e4	(coul)	
	R = 8.315	(J/degC)
	z = 2		(1)		
}

ASSIGNED {
	v		(mV)
	dt		(ms)
	i		(nA)
	g		(uS)
	factor
	wf
	inf		(uS)
	tau		(ms)





}

STATE {
	A
	B
	C
	gVD (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	
	
	tau = tauV * Q10^((T0 - celsius)/10(degC))

	gVD = 0
	wf = 1
	Mgblock(v)
	rates(v)




}

BREAKPOINT {
	SOLVE state METHOD runge  

	i = (B - A)*(gVI + gVD)*Mgblock(v)*(v - e)
}

DERIVATIVE state {
	rates(v)
	A' = -A/tau1
	B' = -B/tau2
	
	gVD' = (B/wf)*(inf-gVD)/tau
	
}

NET_RECEIVE(weight) {
	wf = weight*factor
	A = A + wf
	B = B + wf
}

FUNCTION Mgblock(v(mV)) {
	
	Mgblock = 1 / (1 + (Mg/K0)*exp((0.001)*(-z)*delta*F*v/R/(T+celsius)))
}

PROCEDURE rates(v (mV)) { 
	inf = (v - v0_gVD) * st_gVD * gVI
}