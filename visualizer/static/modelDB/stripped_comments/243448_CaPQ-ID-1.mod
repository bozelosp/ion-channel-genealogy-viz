NEURON {
	SUFFIX CaPQ
	USEION ca READ cai, cao WRITE ica
	RANGE minf, mtau
	RANGE pmax, ica
}

UNITS {
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	(mM)	= (milli/liter)
	FARADAY = 96489 (coul)
	R       = 8.314 (volt-coul/degC)
}

PARAMETER {
	v				(mV)
	pmax = 8e-6		(cm/s)
	celsius			(degC)
}

STATE {
	m
}

ASSIGNED {
	ica		(mA/cm2)
	mtau		(ms)
	minf
	cai (mM)
	cao	(mM)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	ica = pmax*m*ghk(v,cai,cao,2)
}

DERIVATIVE state {
	rates(v)
	m'= (minf-m) / mtau
}

INITIAL {
	rates(v)
	m = minf
}


FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w
	w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
	if (fabs(w)>1e-4) 
	  { e = w / (exp(w)-1) }
	else 
	  { e = 1-w/2 }
	ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}

UNITSOFF

PROCEDURE rates(v(mV)) {
	minf = 1 / (1+exp((-5.1-v)/3.1))
    mtau = 0.35 + (124.9/(18.14*1.25))*exp(-2*((v+9.73)/18.14)^2) 
}

UNITSON