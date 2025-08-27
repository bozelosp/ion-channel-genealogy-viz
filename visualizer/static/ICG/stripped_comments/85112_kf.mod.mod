NEURON {
	SUFFIX kf
	
	USEION k READ ek WRITE ik
        RANGE gbar
	RANGE tau_n, ninf
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 3e-6 
	


	A_anF = 0.00798 (/ms) 
	B_anF = 72.2 (mV)
	C_anF = 1.1 (mV)

	A_bnF = 0.0142 (/ms) 
	B_bnF = 55 (mV)
	C_bnF = 10.5 (mV)









}

ASSIGNED {
        ek (mV)
	v	(mV) 
	ik	(mA/cm2)
	g	(S/cm2)
	tau_n	(ms)
	ninf
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * n^4
	ik = g * (v-ek)
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
	if (Vm+B_bnF != 0) {
		betan=A_bnF*(-B_bnF-Vm)/(1-exp((Vm+B_bnF)/C_bnF))
	} else {
		betan=A_bnF*C_bnF
	}
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_n = 1.0 / (alphan(Vm) + betan(Vm))
	ninf = alphan(Vm) * tau_n
}