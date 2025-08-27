UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {           
        dt            (ms)
	v             (mV)
        tBase = 23.5  (degC)
	celsius = 22  (degC)
	gcatbar = 1.0   (mho/cm2)  
	ki = 0.001    (mM)
	cai (mM)       
	cao (mM)       
        tfa = 1                  
        tfi = 0.68               
        
}

NEURON {
	SUFFIX catp
	USEION ca READ cai,cao WRITE ica 
        
        
        
        
        RANGE gcatbar, hinf, minf, taum, tauh
}

STATE {	m h }  

ASSIGNED {     
	ica (mA/cm2)
        gcat  (mho/cm2) 
        minf
        hinf
        taum
        tauh
}

INITIAL {

	rates(v)
        m = minf
        h = hinf
	gcat = gcatbar*m*m*h*h2(cai)
}

BREAKPOINT {
	SOLVE states
	gcat = gcatbar*m*m*h*h2(cai) 
	ica = gcat*ghk(v,cai,cao)    

}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) { LOCAL nu,f
        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {   
        KTF = ((25./293.15)*(celsius + 273.15))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alph(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alph = 1.6e-4*exp(-(v+57)/19)
}

FUNCTION beth(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
	beth = 1/(exp((-v+15)/10)+1.0)
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	betm = 0.046*exp(-v/22.73)
}

UNITSON
LOCAL facm,fach




PROCEDURE states() {     
        rates(v)
        m = m + facm*(minf - m)
        h = h + fach*(hinf - h)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { 
        LOCAL a
        a = alpm(v)
        taum = 1/(tfa*(a + betm(v))) 
        minf =  a/(a+betm(v))        
        facm = (1 - exp(-dt/taum))
        a = alph(v)
        tauh = 1/(tfi*(a + beth(v))) 
        hinf = a/(a+beth(v))         
        fach = (1 - exp(-dt/tauh))
}