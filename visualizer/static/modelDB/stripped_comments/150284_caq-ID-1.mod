UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

NEURON {
	SUFFIX caq
	USEION ca READ cai, cao WRITE ica
	RANGE pcaqbar, ica
}

PARAMETER {
	pcaqbar = 6.0e-6(cm/s)		

	mvhalf = -9.0	(mV)		
	mslope = -6.6	(mV)		
	mtau = 1.13	(ms)			
	mshift = 0	(mV)
	
	qfact = 3					
}

ASSIGNED { 
    v 		(mV)
    eca		(mV)
    ica 	(mA/cm2)
    minf

    celsius	(degC)
    cai		(mM)
    cao		(mM)
}

STATE {
    m
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica  = ghk(v,cai,cao) * pcaqbar * m * m    
}						

INITIAL {
    settables(v)
    m = minf
}

DERIVATIVE states {  
    settables(v)
    m' = (minf - m) / (mtau/qfact)
}

PROCEDURE settables( v (mV) ) {
	TABLE minf DEPEND mshift
        FROM -100 TO 100 WITH 201

		minf = 1  /  ( 1 + exp( (v-mvhalf-mshift) / mslope) )
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