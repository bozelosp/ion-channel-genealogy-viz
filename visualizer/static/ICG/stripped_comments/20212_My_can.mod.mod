UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX mycan
	USEION ca READ cai,cao WRITE ica
        RANGE gcanbar, gcan, hinf, minf, taum, tauh, facm, fach, ica
}

PARAMETER {
        dt  (ms)
	v (mV)
	celsius = 6.3	(degC)
	gcanbar = 1.0 (mho/cm2)
	cai=5.e-5 (mM)
	cao = 2  (mM)
        tadjm
        tadjh
        eca = 140
}


STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        gcan  (mho/cm2) 
        minf
        hinf
        taum
        tauh
        facm
        fach
}

INITIAL {
         tadjm= 0.3*3.55^((celsius-23.5)/10)
         tadjh= 0.4*2.8^((celsius-23.5)/10) 
	 rates(v)
         m = minf
         h = hinf
         gcan = gcanbar*m*m*h
}

BREAKPOINT {
	SOLVE states
	gcan = gcanbar*m*m*h

	ica = gcan*ghk(v,cai,cao,2)

}




PROCEDURE states() {     
        rates(v)
        m = m + facm*(minf - m)
        h = h + fach*(hinf - h)
        VERBATIM
        return 0;
        ENDVERBATIM
}

PROCEDURE rates(v (mV)) { 




        taum = 15*(0.055*(25.01 - v)/(exp((25.01-v)/10) - 1)+9.4*exp((-63.01-v)/20))
       minf = 1/(1+exp(-(v+11)/8.3))
        facm = (1 - exp(-dt/taum))






         tauh = 10/(0.01*(10.01 - v)/(exp((10.01-v)/10) - 1)+20*exp((-110.01-v)/20))
        hinf = 1/(1+exp((v+44)/9.2))
        fach = (1 - exp(-dt/tauh))
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) {
        LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        e = w / (exp(w)-1)
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else
        
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}