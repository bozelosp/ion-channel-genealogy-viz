NEURON {
	SUFFIX NaCh
	USEION na READ ena WRITE ina
	RANGE gbar, g, ina
	GLOBAL minf, mtau, hinf, htau, am, bm, ah, bh
}



UNITS {
	(mV) = (millivolt)
	(S) = (mho)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = .05 (S/cm2) 
	gamma = .1

	kam = 76.4 (/ms)
	eam = .037 (/mV)
	
	kbm = 6.930852 (/ms)	
	ebm = -.043 (/mV)

	kah = .00013 (/ms)
	eah = -.1216 (/mV)
	kbh = 1.999 (/ms)
	ebh = .0384 (/mV)
}

ASSIGNED {
	v (mV)
	ena (mV)
	ina (mA/cm2)

	minf
	mtau (ms)
	hinf
	htau (ms)

	am (/ms)
	bm (/ms)
	ah (/ms)
	bh (/ms)
}

STATE {
	m h
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ina = gbar*m^3*h*(v - ena)
}

DERIVATIVE state {
	rates(v)
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
}

PROCEDURE rates(v(mV)) {
	am = kam*exp(eam*v)
	bm = kbm*exp(ebm*v)

	ah = kah*exp(eah*v)
	bh = kbh*exp(ebh*v)

	minf = am/(am + bm)
	mtau = 1/(am + bm)
	hinf = ah/(ah + bh)
	htau = 1/(ah + bh)
}