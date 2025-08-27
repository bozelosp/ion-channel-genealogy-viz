NEURON {
	SUFFIX lva
	
	USEION ca WRITE  ica
	RANGE Erev,g, gbar, i
	RANGE k, alpha_1, alpha_2, beta_1, beta_2, V_s
	GLOBAL mytaum, myminf
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 0.4e-3	(S/cm2) < 0, 1e9 > 
	Erev = 120 (mV)	
	
	V_s = 0 (mV)	
			
}

ASSIGNED {
	ica (mA/cm2)
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	k
	alpha_1 (1)
	alpha_2	(1)
	beta_1 (1)
	beta_2 (1)
	mytaum (ms)
	myminf (1)
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
	m = minf(v)
	
	C =  beta_1 / alpha_1
	E =  alpha_2 / beta_2
	h = E / (E * C + E + C)
	d = 1 - (1 + C) * h
}

DERIVATIVE states{ 
	rates(v)
	m' = (minf(v) - m)/taum(v)		
	h' = alpha_1 * (1 - h - d) - beta_1 * h
	d' =  beta_2 * (1 - h - d) - alpha_2 * d
}

FUNCTION minf(Vm (mV)) (1) {
	minf = 1.0 / (1.0 + exp(-(Vm + V_s + 63)/7.8))
	myminf = minf
}

FUNCTION taum(Vm (mV)) (ms) {
	taum = (1.7 + exp( -(Vm + V_s + 28.8)/13.5 )) / (1.0 + exp( -(Vm + V_s + 63)/7.8) )
	mytaum = taum
}

PROCEDURE rates(Vm(mV)) { LOCAL tau_2
	k = (0.25 + exp((Vm + V_s + 83.5)/6.3))^0.5 - 0.5
	tau_2 = 240.0 / (1 + exp((Vm + V_s + 37.4)/30)) 
	alpha_1 = exp( -(Vm + V_s +160.3)/17.8 )	
	beta_1 = k * alpha_1
	alpha_2 = 1.0 / ( tau_2 * (1.0 + k) )
	beta_2 = k * alpha_2
}