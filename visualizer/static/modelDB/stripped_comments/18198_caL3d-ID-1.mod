INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX caL
	USEION ca READ cai, cao WRITE ica
	RANGE O, C, I
	RANGE a,b
	GLOBAL Ra, Rb, q, th, p
	GLOBAL q10, temp, tadj
}

UNITS {
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
} 

PARAMETER {
	p    = 0.2e-3  	(cm/s)		
	v 		(mV)

	th   = 5	(mV)		
	q   = 13	(mV)		

	

	Ra   = 1.6	(/ms)		
	Rb   = 0.2	(/ms)		

	celsius		(degC)
	temp = 22	(degC)		
	q10  = 3			
} 


ASSIGNED {
	ica 		(mA/cm2)
	cao		(mM)
	cai		(mM)
	a (/ms)	b (/ms)
	tadj
}
 

STATE { C O }

INITIAL { 
	C = 1 
}


BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse
	ica = O * p * ghk(v,cai,cao)
} 


KINETIC kstates {
	~ C <-> O 	(a,b)	
	CONSERVE C+O = 1
}	
	
PROCEDURE rates(v(mV)) {
	TABLE a, b
	DEPEND Ra, Rb, th, celsius, temp, q10
	FROM -100 TO 100 WITH 200

	tadj = q10 ^ ((celsius - temp)/10 (degC))

	a = Ra / (1 + exp(-(v-th)/q)) * tadj
	b = Rb / (1 + exp((v-th)/q)) * tadj
}






FUNCTION ghk(v(mV), ci(mM), co(mM)) (0.001 coul/cm3) {
	LOCAL z

	z = (0.001)*2*F*v/(R*(celsius+273.15))
	ghk = (.001)*2*F*(ci*efun(-z) - co*efun(z))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}