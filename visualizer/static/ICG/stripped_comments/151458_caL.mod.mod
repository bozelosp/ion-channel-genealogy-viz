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
	SUFFIX caL
	
	USEION ca READ cai, cao WRITE ica
        RANGE pbar, ical, mshift, hshift, mu
	
}

PARAMETER {
	pbar = 6.7e-6   (cm/s)	

	mvhalf = -8.9	(mV)		
	mslope = -6.7	(mV)		
	mshift = 0	(mV)
	mu = 1
	vm = -8.124  	(mV)		
	k = 9.005   	(mV)		
	kpr = 31.4   	(mV)		
	c = 0.0398   	(/ms-mV)	
	cpr = 0.99	(/ms)			

	hvhalf = -13.4	(mV)		
	hslope = 11.9	(mV)		
	htau = 44.3	(ms)			
	hshift = 0	(mV)
	a = 0.17					
	
	qfact = 3					

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
    mtau	(ms)


}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ica  = mu*ghk(v,cai,cao) * pbar * m * m * (a*h + (1-a))    
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