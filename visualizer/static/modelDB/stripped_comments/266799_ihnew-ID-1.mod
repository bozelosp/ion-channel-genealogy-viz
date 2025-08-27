NEURON {
	SUFFIX hpkj
	NONSPECIFIC_CURRENT i
	RANGE ghbar, eh, i
	GLOBAL ninf, ntau

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
CONSTANT {
	q10 = 3

}
PARAMETER {
	v	 	(mV)
	celsius (degC)
	ghbar = .0001	(S/cm2)

	eh = -30	(mV)
}

ASSIGNED {
	i (mA/cm2)
	qt
	ninf
	ntau
}

STATE {
	n
}

INITIAL {
    qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	n = ninf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = ghbar*n*(v - eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/ntau
}

PROCEDURE rates(v (mV)) {

	ninf = 1/(1+exp((v+90.3+3+3)/9.67))

	ntau = 1000/(0.62*(exp((v+68)/-22)+exp((v+68)/7.14)))/qt/1.3
}