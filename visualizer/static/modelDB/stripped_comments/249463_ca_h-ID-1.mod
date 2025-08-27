INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cah
	USEION ca READ cai,cao WRITE ica
	RANGE pbar, minf, taum, hinf, tauh, shift, shifth
	GLOBAL qm, qh
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
	pbar	=.2e-3	(cm/s)	
	shift	= 2 	(mV)	
	shifth	= 0     (mV)	
	cai	  (mM) 
	cao		(mM)
	qm	= 4		
	qh	= 2		

	offm = -14.17
	offh = -22.63
	offmt = -26.31
	offht = 19.73
	slom = 9.76
	sloh = 6.6
	slomt = 31.25
	sloht = 21.2765957
	taummax = 0.97
	tauhmax = 70
	mmax = 1.092
	mmin = 0.75

}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	minf
	taum	(ms)
	hinf
	tauh	(ms)
	phim
	phih
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	ica = pbar * m*m*h * ghk(v, cai, cao)
}

DERIVATIVE castate {
	evaluatefct(v)

	m' = (minf - m) / taum
	h' = (hinf - h) / tauh
}


UNITSOFF
INITIAL {
	phim = qm ^ ((celsius-24)/10)
	phih = qh ^ ((celsius-24)/10)

	evaluatefct(v)

	m = minf
	h = hinf
}

PROCEDURE evaluatefct(v(mV)) {

	minf = mmax/(1+exp((offm-(v+shift))/slom))
	hinf = mmin/(1+exp(-(offh-(v+shifth))/sloh))

	taum = (taummax/(cosh(-(offmt-(v+shift))/slomt)))/phim
	tauh = (tauhmax/(cosh(-(offht-(v+shifth))/sloht)))/phih
	
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
FUNCTION nongat(v,cai,cao) {	
	nongat = pbar * ghk(v, cai, cao)
}
UNITSON