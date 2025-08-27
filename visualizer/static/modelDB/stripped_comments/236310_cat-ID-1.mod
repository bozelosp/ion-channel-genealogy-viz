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
	SUFFIX cat
	USEION cal READ cali, calo WRITE ical VALENCE 2
	RANGE pcatbar, ical
}

PARAMETER {
	pcatbar = 7.6e-7(cm/s)		

	mvhalf = -51.73	(mV)		
	mslope = -6.53	(mV)		
	mshift = 9.163	(mV)
	
	hvhalf = -80	(mV)		
	hslope = 6.7	(mV)		
	hshift = 11.819	(mV)
	
	qfact = 2.16				
}					

ASSIGNED { 
    v 		(mV)
    ical 	(mA/cm2)
    ecal		(mV)
    
    celsius	(degC)
    cali		(mM)
    calo		(mM)

    minf
    hinf
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ical  = ghk(v,cali,calo) * pcatbar * m * m * m * h	
}

INITIAL {
    settables(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  
	settables(v)
	m' = (minf - m) / (mtau(v)/qfact)
	h' = (hinf - h) / (htau(v)/qfact)
}

FUNCTION_TABLE mtau(v(mV))	(ms)		
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