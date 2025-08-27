NEURON {
	SUFFIX nap
	NONSPECIFIC_CURRENT i
	RANGE gbar, ena
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.00005
	ena=60(mV)

	A_amp = 17.235 (/ms) 
	B_amp = 27.58 (mV)
	C_amp = -11.47 (mV)

	A_bmp = 17.235 (/ms) 
	B_bmp = 86.2 (mV)
	C_bmp = 19.8 (mV)
}

ASSIGNED {
	v	(mV) 
	i	(mA/cm2)
	g	(S/cm2)
	tau_m	(ms)
	minf
	hinf
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3
	i = g * (v-ena)
}

INITIAL {
	
	m = alpham(v)/(alpham(v)+betam(v))
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
}

FUNCTION alpham(Vm (mV)) (/ms) {
	alpham=A_amp/(1+exp((Vm+B_amp)/C_amp))
}

FUNCTION betam(Vm (mV)) (/ms) {
	betam=A_bmp/(1+exp((Vm+B_bmp)/C_bmp))
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_m = 1.0 / (alpham(Vm) + betam(Vm))
	minf = alpham(Vm) * tau_m
}