NEURON {
	SUFFIX lva
	
	USEION ca WRITE  ica
	RANGE Erev,g, gbar, i
	RANGE k, tee, alpha_1, alpha_2
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 166e-6	(S/cm2) < 0, 1e9 > 
					   
	Erev = 80 (mV)	
			
}

ASSIGNED {
	ica (mA/cm2)
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	k
	tee	
	alpha_1	
	alpha_2	
}

STATE {	m h d }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3 * h
	ica = g * (v - Erev)
	i = ica	
}

INITIAL {
	LOCAL C, E
	
	
	rates(v)
	m = alpham(v)/(alpham(v) + betam(v))
	
	C =  beta1(v) /  alpha1(v)
	E =  alpha2(v) /  beta2(v)
	h = E / (E * C + E + C)
	d = 1 - (1 + C) * h
}

DERIVATIVE states{ 
	rates(v)
	m' = alpham(v) * (1 - m) - betam(v) * m
	h' = alpha1(v) * (1 - h - d) - beta1(v) * h
	d' =  beta2(v) * (1 - h - d) - alpha2(v) * d
}

FUNCTION alpham(Vm (mV)) (/ms) {
	UNITSOFF
	alpham = 3.3 /(1.7 + exp(-(Vm + 28.8)/13.5))
	UNITSON
}

FUNCTION betam(Vm (mV)) (/ms) {
	UNITSOFF
	betam =  3.3 * exp(-(Vm + 63)/7.8)/(1.7 + exp(-(Vm + 28.8)/13.5))
	UNITSON
}

FUNCTION alpha1(Vm (mV)) (/ms) {
	UNITSOFF
	alpha1 = alpha_1
	UNITSON
}

FUNCTION beta1(Vm (mV)) (/ms) {
	UNITSOFF
	beta1 =  k * alpha_1
	UNITSON
}

FUNCTION alpha2(Vm (mV)) (/ms) {
	UNITSOFF
	alpha2 = alpha_2
	UNITSON
}

FUNCTION beta2(Vm (mV)) (/ms) {
	UNITSOFF
	beta2 =  k * alpha_2
	UNITSON
}

PROCEDURE rates(Vm(mV)) {
	k = (0.25 + exp((Vm + 83.5)/6.3))^0.5 - 0.5
	tee = 240.0 / (1 + exp((Vm + 37.4)/30))
	alpha_1 = 2.5 / (tee*(1 + k))	
	alpha_2 = 2.5 * exp(-(Vm + 160.3)/17.8)
}