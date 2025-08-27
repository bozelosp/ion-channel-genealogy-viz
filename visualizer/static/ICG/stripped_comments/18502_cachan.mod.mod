UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}

NEURON {
	SUFFIX cachan
	USEION ca READ cai, cao WRITE ica
	RANGE pcabar, ica
}

UNITS {
	
	
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	taufactor=.5	
	celsius		(degC) 
	v		(mV)
	pcabar=.2e-3	(cm/s)	
	cai		(mM)
	cao		(mM)
}

ASSIGNED { ica		(mA/cm2)}

STATE {	oca }		

BREAKPOINT {
	SOLVE castate METHOD cnexp
	ica = pcabar*oca*oca*ghk(v, cai, cao)
}

DERIVATIVE castate {
	oca' = (oca_ss(v) - oca)/oca_tau(v)
}

INITIAL {
	oca = oca_ss(v)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	
	
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION oca_ss(v(mV)) {
	LOCAL a, b
	TABLE FROM -150 TO 150 WITH 200
	
	v = v+65
	a = 1(1/ms)*efun(.1(1/mV)*(25-v))
	b = 4(1/ms)*exp(-v/18(mV))
	oca_ss = a/(a + b)
}

FUNCTION oca_tau(v(mV)) (ms) {
	LOCAL a, b
	TABLE FROM -150 TO 150 WITH 200

	v = v+65
	a = 1(1/ms)*efun(.1(1/mV)*(25-v))
	b = 4(1/ms)*exp(-v/18(mV))
	oca_tau = taufactor/(a + b)
}