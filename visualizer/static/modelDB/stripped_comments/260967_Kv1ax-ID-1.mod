NEURON {
	
	SUFFIX kv1ax
	USEION k READ ek WRITE ik
	RANGE gkv1, ikv1, gbar
	GLOBAL ninf, ntau
	GLOBAL h1inf, h1tau
	GLOBAL h2inf, h2tau
	GLOBAL h3inf, h3tau
	THREADSAFE
}

UNITS {
	
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

CONSTANT {


	a1 = 47.7	
	a2 = 41.59	
	b1 = 58.78	
	c1 = 13.28	
	c2 = 48.57	

	a1H1 = 0.347e-03	
	a2H1 = 12.4	
	b2H1 = 61.1	
	c1H1 = 7.28	
	c2H1 = 2.76	

	a12H2fact = 1.60e-02	

	propH2 = 0.82	

	sToMs = 0.001
}

PARAMETER {
	
	v (mV)


	vShift = 10	(mV)
	vShift_inact = -15
	gbar = 11 (pS/um2)   	
	temp = 35.3	(degC)		
	q10  = 3				
	q10h  = 3				
	celsius					(degC)
}

ASSIGNED {
 	ik (mA/cm2) 
 	ikv1 (mA/cm2) 
	gkv1  (mho/cm2)
	ek (mV)

	ninf
	ntau (ms)
	nalpha (1/ms)
	nbeta (1/ms)

	h1inf
	h1tau (ms)
	h1alpha (1/ms)
	h1beta (1/ms)

	h2inf
	h2tau (ms)
	h2alpha (1/ms)
	h2beta (1/ms)

	h3inf
	h3tau (ms)
	h3alpha (1/ms)
	h3beta (1/ms)


	tadj
	tadjh
}

STATE { 
	n
	h1
	h2
}

INITIAL {


	rates(v)
	n = ninf
	h1 = h1inf
	h2 = h2inf
}

BREAKPOINT {
	
	SOLVE states METHOD cnexp
	gkv1 = gbar * ( (1-propH2)*n^8*h1 + propH2*n^8*h2 )
	ikv1 = (1e-4)*gkv1 * (v - ek)
	ik = ikv1
}

DERIVATIVE states {
	
	rates(v)
	n' = (ninf-n)/ntau
	h1' = (h1inf-h1)/h1tau 
	h2' = (h2inf-h2)/h2tau 
}

PROCEDURE rates(v (mV)) {
	
	tadj = q10^((celsius - temp)/10)    
	tadjh = q10h^((celsius - temp)/10)	
	nalpha = tadj*nalphafkt(v-vShift)
	nbeta = tadj*nbetafkt(v-vShift)
	ninf = nalpha/(nalpha+nbeta) 
	ntau = 1/(nalpha + nbeta)

	h1alpha = tadjh*h1alphafkt(v-vShift-vShift_inact)
	h1beta = tadjh*h1betafkt(v-vShift-vShift_inact)
	h1inf = h1alpha/(h1alpha+h1beta) 
	h1tau = 1/(h1alpha + h1beta)

	h2alpha = tadjh*h2alphafkt(v-vShift-vShift_inact)
	h2beta = tadjh*h2betafkt(v-vShift-vShift_inact)
	h2inf = h2alpha/(h2alpha+h2beta) 
	h2tau = 1/(h2alpha + h2beta)
}

FUNCTION nalphafkt(v (mV)) (1/ms) {
	nalphafkt = sToMs * a1*(-(v+b1))/( exp(-(v+b1)/c1) -1)
}

FUNCTION nbetafkt(v (mV)) (1/ms) {
	nbetafkt = sToMs * a2*exp(-(v)/c2)
}

FUNCTION h1alphafkt(v (mV)) (1/ms) {
	h1alphafkt = sToMs * a1H1*exp(-(v)/c1H1)
}

FUNCTION h1betafkt(v (mV)) (1/ms) {
	h1betafkt = sToMs * a2H1/(exp(-(v+b2H1)/c2H1)+1)
}

FUNCTION h2alphafkt(v (mV)) (1/ms) {
	h2alphafkt = sToMs * a12H2fact*a1H1*exp(-(v)/c1H1)
}

FUNCTION h2betafkt(v (mV)) (1/ms) {
	h2betafkt = sToMs * a12H2fact*a2H1/(exp(-(v+b2H1)/c2H1)+1)
}