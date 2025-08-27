NEURON {
	SUFFIX xiong
	USEION k READ ko, ki, ek WRITE ik
	USEION na READ nai, nao, ena WRITE ina
	USEION ca READ cao
	GLOBAL tauavg, rmax, ec50, Nh, n
	RANGE ik, ina, g, gpresent, itot
}

UNITS {
	(molar) = 	(1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	FARADAY	= (faraday) (coulomb)
	
	R = (k-mole) (joule/degC)
}



PARAMETER {
	celsius		(degC)
	g=1e-3	(mho/cm2)
	tau_ina=100	(ms)
	tau_act=992	(ms)
	rmax=1
	ec50=.145	(mM) 
	Nh=1.4		
	n=4		
	tauavg=300	(ms) 
}

ASSIGNED { 
	v	(mV)	
	ina	(mA/cm2)
	ik	(mA/cm2)
	ena	(mV)
	ek	(mV)
	cao	(mM)
	ki
	ko
	nai
	nao
	gpresent
	itot
}

STATE { m }

BREAKPOINT {
	SOLVE gatestate METHOD cnexp
	gpresent=g*m^n
	ina = gpresent*(v-ena)
	ik = gpresent*(v-ek)
	itot=ik+ina
}

INITIAL {
	m=hill(cao)
	gpresent=g*m^n
	ina = gpresent*(v-ena)
	ik = gpresent*(v-ek)
	itot=ik+ina
}

DERIVATIVE gatestate {
	m' = ( (hill(cao)-m)/tauavg )
}

FUNCTION hill(co) {
	TABLE DEPEND rmax, ec50, Nh, n FROM 0 TO 15 WITH 150
	hill = (rmax/(1+(co/ec50)^Nh))^(1/n)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*1*FARADAY*v/(R*(celsius+273.11247574))
	eco = co*efun(z)
	eci = ci*efun(-z)
	ghk = (.001)*1*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}