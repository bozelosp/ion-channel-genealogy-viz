NEURON {
	POINT_PROCESS gaba

	USEION cl READ ecl WRITE icl VALENCE -1
        USEION hco3 READ ehco3 WRITE ihco3 VALENCE -1

	RANGE tau1, tau2, g
	RANGE P, i

	RANGE icl, ihco3, ehco3, e
	GLOBAL total
}

UNITS {
	(mA)    = (milliamp)
	(nA)    = (nanoamp)
	(mV)    = (millivolt)
	(uS)  = (micromho)
	(mM)    = (milli/liter)
	F = (faraday) (coulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {
	tau1	=.1	(ms)	<1e-9,1e9>
	tau2	= 80	(ms)	<1e-9,1e9>

	P       = 0.18		

	celsius = 31    (degC)
}

ASSIGNED {
	v	(mV)		

	icl	(nA)		
	ihco3	(nA)		
	i	(nA)		
				
	g 	(uS)		
				
	factor
	total	(uS)

	ecl	(mV)		
	ehco3	(mV)		

	e	(mV)		
}

STATE {
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
	e = P/(1+P)*ehco3 + 1/(1+P)*ecl
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	g = B - A

	icl = 1/(1+P)*g*(v-ecl)

	ihco3 = P/(1+P)*g*(v-ehco3)
	i = icl + ihco3
	e = P/(1+P)*ehco3 + P/(1+P)*ecl

}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (uS)) {
	A = A + weight*factor
	B = B + weight*factor
	total = total+weight
}