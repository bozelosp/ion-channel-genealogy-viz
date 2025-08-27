NEURON {
	SUFFIX nav1p8
	NONSPECIFIC_CURRENT i
	RANGE gbar, ena, localtemp
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.013 
	ena=60 (mV)
	localtemp = 37

	A_am8 = 3.83 (/ms) 
	B_am8 = 2.58 (mV)
	C_am8 = -11.47 (mV)

	A_ah8 = 0.010152 (/ms) 
	B_ah8 = 105 (mV)
	C_ah8 = 46.33 (mV)

	A_bm8 = 6.894 (/ms) 
	B_bm8 = 61.2 (mV)
	C_bm8 = 19.8 (mV)

	A_bh8 = 0.822853 (/ms)   
	B_bh8 = -21.8 (mV)
	C_bh8 = -11.998 (mV)
}

ASSIGNED {
	v	(mV) 
	i	(mA/cm2)
	g	(S/cm2)
	tau_h	(ms)
	tau_m	(ms)
	minf
	hinf
	Q10gate  
	Q10cond  
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = Q10cond * gbar * m^3 * h
	i = g * (v-ena)
}

INITIAL {
	
	m = alpham(v)/(alpham(v)+betam(v))
	h = alphah(v)/(alphah(v)+betah(v))
	Q10gate = 3^((localtemp-25(degC))/10 (degC))
	Q10cond = 1^((localtemp-25(degC))/10 (degC))
}

DERIVATIVE states {
	rates(v)
	m' = Q10gate * (minf - m)/tau_m
	h' = Q10gate * (hinf - h)/tau_h
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_am8/(1+exp((Vm+B_am8)/C_am8))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	alphah=A_ah8*exp(-(Vm+B_ah8)/C_ah8)
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bm8/(1+exp((Vm+B_bm8)/C_bm8))
}

FUNCTION betah(Vm (mV)) (/ms) {
	betah=A_bh8/(1+exp((Vm+B_bh8)/C_bh8))
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))
	minf = alpham(Vm) * tau_m

	tau_h = 0.8 / (alphah(Vm) + betah(Vm))   
	hinf = alphah(Vm) * tau_h
}