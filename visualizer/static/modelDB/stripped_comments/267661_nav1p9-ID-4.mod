NEURON {
	SUFFIX nav1p9
	NONSPECIFIC_CURRENT i
	RANGE gbar, ena, ActShift, localtemp
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.0003 
	ena=60 (mV)
	localtemp = 37
	ActShift

	A_am9 = 1.548 (/ms) 
	B_am9 = -11.01 (mV)
	C_am9 = -14.871 (mV)

	A_ah9 = 0.2574 (/ms) 
	B_ah9 = 63.264 (mV)
	C_ah9 = 3.7193 (mV)

	A_bm9 = 8.685 (/ms) 
	B_bm9 = 112.4 (mV) 
	C_bm9 = 22.9 (mV)

	A_bh9 = 0.53984 (/ms)   
	B_bh9 = 0.27853 (mV)
	C_bh9 = -9.0933 (mV)
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
	g = Q10cond * gbar * m * h
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
	alpham=A_am9/(1+exp((Vm+B_am9 + ActShift)/C_am9))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	alphah=A_ah9/(1+exp((Vm+B_ah9)/C_ah9))
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bm9/(1+exp((Vm+B_bm9 + ActShift)/C_bm9))
}

FUNCTION betah(Vm (mV)) (/ms) {
	betah=A_bh9/(1+exp((Vm+B_bh9)/C_bh9))
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))
	minf = alpham(Vm) * tau_m

	tau_h = 1.0 / (alphah(Vm) + betah(Vm))
	hinf = alphah(Vm) * tau_h
}