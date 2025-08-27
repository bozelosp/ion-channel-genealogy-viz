INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaT
	USEION Ca READ Cai, Cao WRITE iCa VALENCE 2
        RANGE minf, mtau, hinf, htau, iCa, p
	GLOBAL pmax
}

UNITS {
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	(mM)	= (milli/liter)
        FARADAY = 96489 (coul)
        R       = 8.314 (volt-coul/degC)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
	Cai		(mM)
	Cao		(mM)
	pmax = 4e-7	(cm/s)		
}

STATE {
	m h
}

ASSIGNED {
	iCa		(mA/cm2)
	mtau		(ms)
	minf
	hinf
	htau		(ms)
	p		(cm/s)
}

BREAKPOINT { 
	
	SOLVE state METHOD euler

	p = pmax*m*m*m*h
	iCa = p*ghk(v,Cai,Cao,2)	
	
}

DERIVATIVE state {
	rates(v)
	m'= (minf-m) / mtau
	h'= (hinf-h) / htau
	}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else 
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}
FUNCTION_TABLE tabmtau(v(mV)) (ms)
FUNCTION_TABLE tabhtau(v(mV)) (ms)
UNITSOFF

PROCEDURE rates(v(mV)) { 
	mtau = tabmtau(v)
	htau = tabhtau(v)
	minf = 1 / (1+exp((v+51.73)/-6.53))
	hinf = 1 / (1+exp((v+80)/6.7))
}