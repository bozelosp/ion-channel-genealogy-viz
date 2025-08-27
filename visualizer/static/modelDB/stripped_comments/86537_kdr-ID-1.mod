NEURON {
	SUFFIX kdr
	NONSPECIFIC_CURRENT i
	RANGE gbar, ek
	RANGE tau_n, n
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.0021 (S/cm2)
	
      ek = -92.34 (mV) 

	A_anF = 0.001265 (/ms) 
	B_anF = 14.273 (mV) 
	C_anF = 10 (mV) 


	A_bnF = 0.125 (/ms)
	B_bnF = 55 (mV)
	C_bnF = -2.5 (mV)









}

ASSIGNED {
	v	(mV) 
	i	(mA/cm2)
	g	(S/cm2)
	tau_n	(ms)
	ninf
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n
	i = g * (v-ek)
}

INITIAL {
	
	n = alphan(v)/(alphan(v)+betan(v))
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/tau_n
}

FUNCTION alphan(Vm (mV)) (/ms) {
	if (-Vm-B_anF != 0) {
		alphan=A_anF*(Vm+B_anF)/(1-exp((-Vm-B_anF)/C_anF))
	} else {
		alphan=A_anF*C_anF
	}
}

FUNCTION betan(Vm (mV)) (/ms) {
	betan=A_bnF*exp((Vm+B_bnF)/C_bnF)
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_n = 1.0 / (alphan(Vm) + betan(Vm))

      ninf = 1/(1+exp((Vm+14.62)/-18.38)) 
}