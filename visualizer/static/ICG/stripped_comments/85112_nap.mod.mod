NEURON {
	SUFFIX nap
	
	USEION na READ ena WRITE ina
        RANGE gbar
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 3.6e-8 
	

	A_amp = 17.235 (/ms) 
	B_amp = 27.58 (mV)
	C_amp = -11.47 (mV)

	A_bmp = 17.235 (/ms) 
	B_bmp = 86.2 (mV)
	C_bmp = 19.8 (mV)
}

ASSIGNED {
        ena (mV)
	v	(mV) 
	ina	(mA/cm2)
	g	(S/cm2)
	tau_m	(ms)
	minf
	hinf
}

STATE { m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3
	ina = g * (v-ena)
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