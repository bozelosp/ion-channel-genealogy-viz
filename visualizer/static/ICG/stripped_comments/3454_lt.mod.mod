NEURON {
	SUFFIX LT
	USEION k READ ek WRITE ik
	RANGE gbar, g, ik
	GLOBAL linf, ltau, rinf, rtau, al, bl, ar, br
}



UNITS {
	(mV) = (millivolt)
	(S) = (mho)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = .002 (S/cm2) 
	gamma = .1

	kal = 1.2 (/ms)
	eal = .03512 (/mV)
	kbl = .2248 (/ms)
	ebl = -.0319 (/mV)

	kar = .0438 (/ms)
	ear = -.0053 (/mV)
	kbr = .0562 (/ms)
	ebr = -.0047 (/mV)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)

	linf
	ltau (ms)
	rinf
	rtau (ms)

	al (/ms)
	bl (/ms)
	ar (/ms)
	br (/ms)
}

STATE {
	l r
}

INITIAL {
	rates(v)
	l = linf
	r = rinf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gbar*l*r*(v - ek) 
}

DERIVATIVE state {
	rates(v)
	l' = (linf - l)/ltau
	r' = (rinf - r)/rtau
}

PROCEDURE rates(v(mV)) {
	al = kal*exp(eal*v)
	bl = kbl*exp(ebl*v)

	ar = kar*exp(ear*v)
	br = kbr*exp(ebr*v)

	linf = al/(al + bl)
	ltau = 1/(al + bl)
	rinf = ar/(ar + br)
	rtau = 1/(ar + br)
}