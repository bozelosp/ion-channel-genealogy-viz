NEURON {
	SUFFIX kir2
	USEION k READ ek WRITE ik
	RANGE gbar, g, i, mu
	GLOBAL ninf, ntau
	
}

UNITS {
	(mA) = (milliamp)
	(uA) = (microamp)
	(mV) = (millivolt)
	(mS) = (millimho)
}

PARAMETER {
	gbar = 1.2	(mS/cm2)	<0,1e9>
	ek = -90	(mV)
	vh = -111	(mV)	
	ve = -11	(mV)	
	ntauconst = 0.1	(ms)	
	mu = 1
}

ASSIGNED {
	v	(mV)
	g	(mho/cm2)
	i	(uA/cm2)	
	ik	(mA/cm2)
	ninf	(1)
	ntau	(ms)
	
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = mu*(0.001)*gbar*n
	ik = g*(v - ek)     
	i = (1000)*ik
}

INITIAL {
	rates(v)
	n = ninf
}

DERIVATIVE states { 
	rates(v)
	n' = (ninf-n)/ntau
}




PROCEDURE rates(v(mV)) {
UNITSOFF
	
	ntau = ntauconst
	ninf = 1/(1 + exp(-(v - vh)/ve))
}
UNITSON