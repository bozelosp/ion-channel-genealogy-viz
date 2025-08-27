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
	SUFFIX can
	USEION ca READ cai, cao WRITE ica
	RANGE pbar, ica
}

PARAMETER {
	pbar = 1.0e-5	(cm/s)		

	mvhalf = -8.7	(mV)		
	mslope = -7.4	(mV)		
	mshift = 0	(mV)

	hvhalf = -74.8	(mV)		
	hslope = 6.5	(mV)		
	hshift = 0	(mV)
	htau = 70.0	(ms)		

	vm = -17.19  	(mV)		
	k = 15.22   	(mV)		
	kpr = 23.82   	(mV)		
	c = 0.03856   	(/ms-mV)	
	cpr = 0.3842	(/ms)		

	a = 0.21		
	qfact = 3		
}

ASSIGNED { 
    v		(mV)
    ica 	(mA/cm2)

    celsius	(degC)
    cai		(mM)
    cao		(mM)
    
    minf
    mtau	(ms)

    hinf
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica  = ghk(v,cai,cao) * pbar * m * m * (a*h + (1-a))    
}

INITIAL {
    settables(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  
    settables(v)
    m' = (minf - m) / (mtau/qfact)
    h' = (hinf - h) / (htau/qfact)
}

PROCEDURE settables( v (mV) ) {
	LOCAL malpha, mbeta
	
	TABLE minf, hinf, mtau DEPEND mshift, hshift
        FROM -100 TO 100 WITH 201

		minf = 1  /  ( 1 + exp( (v-mvhalf-mshift) / mslope) )
		hinf = 1  /  ( 1 + exp( (v-hvhalf-hshift) / hslope) )

		malpha = c * (v-vm) / ( exp((v-vm)/k) - 1 )
		mbeta = cpr * exp(v/kpr)		
		mtau = 1 / (malpha + mbeta)
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