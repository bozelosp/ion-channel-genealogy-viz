NEURON {
	POINT_PROCESS Exp5NMDA2
	NONSPECIFIC_CURRENT i
	RANGE tau1, tau2_0, a2, b2, wtau2, tau3_0, a3, b3, tauV, e, i, gVI, gVDst, gVDv0, Mg, K0, delta, tp, wf, tau_D, d
	GLOBAL inf, tau2, tau3
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

	tau1 = 1.69		(ms)	<1e-9,1e9>	
										

	tau2_0 = 3.97	(ms)
	a2 = 0.70		(ms)
	b2 = 0.0243		(1/mV)
	wtau2= 0.65		<1e-9,1> 
	
	

	tau3_0 = 41.62	(ms)
	a3 = 34.69		(ms)
	b3 = 0.01		(1/mV)
	
	
	
	Q10_tau1 = 2.2			
	Q10_tau2 = 3.68			
	Q10_tau3 = 2.65			
	T0_tau	 = 35	(degC)	
	
	
	tp = 30			(ms)	
							

	
	
	d = 0.2 (1) < 0, 1 >     
	tau_D = 2500 (ms) < 1e-9, 1e9 >

	tauV = 7		(ms)	<1e-9,1e9>	
							
							
							
							
							
	gVDst = 0.007	(1/mV)	
	gVDv0 = -100	(mV)	
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
	q10_tau2
	q10_tau3
	inf		(uS)
	tau		(ms)
	tau2	(ms)
	tau3	(ms)
	wtau3
}

STATE {
	A		
	B		
	C		
	D		
	gVD (uS)
}

INITIAL { 
	Mgblock(v)
	
	tau1 = tau1 * Q10_tau1^((T0_tau - celsius)/10(degC))
	q10_tau2 = Q10_tau2^((T0_tau - celsius)/10(degC))
	q10_tau3 = Q10_tau3^((T0_tau - celsius)/10(degC))
	
	tau  = tauV * Q10^((T0 - celsius)/10(degC))
	
	rates(v)
	wtau3 = 1 - wtau2
	
	
	
	factor = -exp(-tp/tau1) + wtau2*exp(-tp/tau2) + wtau3*exp(-tp/tau3)
	factor = 1/factor

	A = 0
	B = 0
	C = 0
	D = 1
	gVD = 0
	wf = 1
}

BREAKPOINT {
	SOLVE state METHOD derivimplicit 

	i = (wtau3*C + wtau2*B - A)*(gVI + gVD)*Mgblock(v)*(v - e)
}

DERIVATIVE state {
	rates(v)
	A' = -(A-wf*(1-d)*D')/tau1
	B' = -(B-wf*(1-d)*D')/tau2
	C' = -(C-wf*(1-d)*D')/tau3
	D' = (1-D)/tau_D
	
	gVD' = ((wtau3*C + wtau2*B)/wf)*(inf-gVD)/tau
}

NET_RECEIVE(weight) {
	wf = weight*factor*D
	A = A + wf
	B = B + wf
	C = C + wf

	D = D * d
	wf = weight*factor
}

FUNCTION Mgblock(v(mV)) {
	
	Mgblock = 1 / (1 + (Mg/K0)*exp((0.001)*(-z)*delta*F*v/R/(T+celsius)))
}

PROCEDURE rates(v (mV)) { 
	inf = (v - gVDv0) * gVDst * gVI
	
	tau2 = (tau2_0 + a2*(1-exp(-b2*v)))*q10_tau2
	tau3 = (tau3_0 + a3*(1-exp(-b3*v)))*q10_tau3
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	if (tau2/tau3 > .9999) {
		tau2 = .9999*tau3
	}
}