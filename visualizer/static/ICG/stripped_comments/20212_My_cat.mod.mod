UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        dt  (ms)
	v (mV)
        tBase = 23.5  (degC)
	celsius = 22	(degC)
	gcatbar = 1.0 (mho/cm2)
	ki=.001 (mM)
	cai=5.e-5 (mM)
	cao = 2  (mM)
        eca = 140
        VhalfH= -85
        VhalfM= -47
        gh = 0.4
        zh = 5
        gm = 0.40
        zm = 5
        aoh = 0.3
        aom =  0.03   
        tfi = 10
        tfa = 1              
        k
        tadj
        taum
        tauh
}


NEURON {
	SUFFIX mycat
	USEION ca READ cai,cao WRITE ica 

        
        
        
        RANGE gcatbar, hinf, minf, taum, tauh, ica

}

STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        gcat  (mho/cm2) 
        minf
        hinf
}

INITIAL {
        k = 1000*8.3134*(celsius+273.15)/96520
        tadj = 3^((celsius-tBase)/10)   
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


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

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
	alph = aoh*exp(-gh*zh*(v-VhalfH)/k)
}

FUNCTION beth(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
	beth = (aoh/10)*exp((1-gh)*zh*(v-VhalfH)/k)
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
        alpm =  aom*exp((1-gm)*zm*(v-VhalfM)/k)
}

FUNCTION betm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	betm = aom*exp(-gm*zm*(v-VhalfM)/k)
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

        taum = 20
        minf =  a/(a+betm(v))
        facm = (1 - exp(-dt/taum))
        a = alph(v)

        tauh = 60
        hinf = a/(a+beth(v))
        fach = (1 - exp(-dt/tauh))
}