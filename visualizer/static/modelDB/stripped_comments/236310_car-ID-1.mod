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
	SUFFIX car
	USEION ca READ cai, cao WRITE ica
	RANGE pcarbar, ica
}

PARAMETER {
	pcarbar = 2.6e-5(cm/s)		

	mvhalf = -10.3	(mV)		
	mslope = -6.6	(mV)		
	mtau = 5.1	(ms)			
	mshift = 9.163	(mV)

	hvhalf = -33.3	(mV)		
	hslope = 17.0	(mV)		
	hshift	= 11.819	(mV)		
	
	qfact = 2.16					
}					

ASSIGNED { 
    v 		(mV)
    ica 	(mA/cm2)
    eca		(mV)
    
    celsius	(degC)
    cai		(mM)
    cao		(mM)

    minf
    hinf
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica  = ghk(v,cai,cao) * pcarbar * m * m * m * h	
}				    

INITIAL {
    settables(v)
	m = minf
	h = hinf
}

DERIVATIVE states {  
	settables(v)
	m' = (minf - m) / (mtau/qfact)
	h' = (hinf - h) / (htau(v)/qfact)
}

FUNCTION_TABLE htau(v(mV))	(ms)	

PROCEDURE settables( v (mV) ) {
	TABLE minf, hinf DEPEND mshift, hshift
        FROM -100 TO 100 WITH 201

		minf = 1  /  ( 1 + exp( (v-mvhalf-mshift) / mslope) )
		hinf = 1  /  ( 1 + exp( (v-hvhalf-hshift) / hslope) ) 
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