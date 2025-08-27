INCLUDE "Unit.inc"
INCLUDE "Volume.inc"
NEURON {
	SUFFIX ifun
	USEION k READ ek, ko WRITE ik
	USEION na READ ena WRITE ina
	
	RANGE gK, gNa, ifun, ik, ina
	GLOBAL minf, mtau 
}

PARAMETER {
	gK = 3	(uS)
	gNa= 3	(uS)
	Kmf = 45	(mM)
}

STATE { 
	m 
}

ASSIGNED {
	v (mV)
	celsius (degC) 
	ik (mA/cm2)
	ina (mA/cm2)
	ifun (mA/cm2)
	minf 
	mtau (ms)  
	ek (mV)
	ena (mV)
	ko (mM)
}

INITIAL {
	rate(v)
	m = minf
}

BREAKPOINT { LOCAL kc, ifk, ifna
	SOLVE states METHOD derivimplicit

	kc = m * (ko/(ko+Kmf))
	ifna = (1e-06)*kc * (gNa/S*(v-ena))
	ifk  = (1e-06)*kc * (gK/S*(v-ek))
	ik = ifk
	ina = ifna
	ifun = ifk + ifna
}

DERIVATIVE states {
	rate(v)
	m' = (minf - m)/mtau
}

FUNCTION alp(v(mV)) (/ms) { 



	alp = (0.001)* 0.05(/s)*exp(-0.067(/mV)*(v + 52(mV)-10))
}

FUNCTION bet(v(mV)) (/ms) { 
LOCAL Eo 
	Eo= v + 52 - 10
	if (fabs(Eo*1(/mV)) < 1e-5)
	{
		bet = (0.001)* 2.5 (/s)
	} else {



		bet = (0.001)*1(/mV/s)*Eo/(1 - exp(-0.2(/mV)*Eo))
	}
}


PROCEDURE rate(v(mV)) { LOCAL a,b 
TABLE minf, mtau FROM -100 TO 100 WITH 200
	a = alp(v)  b = bet(v) 
	mtau = 1/(a + b)
	minf = a * mtau
}