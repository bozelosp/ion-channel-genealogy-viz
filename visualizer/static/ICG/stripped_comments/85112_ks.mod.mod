NEURON {
	SUFFIX ks
	
	USEION k READ ek WRITE ik
        RANGE gbar
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 0.2e-6 
	

	A_anS = 0.00122 (/ms) 
	B_anS = -10.5 (mV)
	C_anS = 23.6 (mV)

	A_bnS = 0.000739 (/ms) 
	B_bnS = 57.1 (mV)
	C_bnS = 21.8 (mV)
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
	g = gbar * n 
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
	if (-Vm-B_anS != 0) {
		alphan=A_anS*(Vm+B_anS)/(1-exp((-Vm-B_anS)/C_anS))
	} else {
		alphan=A_anS*C_anS
	}
}

FUNCTION betan(Vm (mV)) (/ms) {
	if (Vm+B_bnS != 0) {
		betan=A_bnS*(-B_bnS-Vm)/(1-exp((Vm+B_bnS)/C_bnS))
	} else {
		betan=A_bnS*C_bnS
	}
}

FUNCTION rates(Vm (mV)) (/ms) {
	tau_n = 1.0 / (alphan(Vm) + betan(Vm))
	ninf = alphan(Vm) * tau_n
}