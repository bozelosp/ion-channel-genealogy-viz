UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX hhsoma
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, q10m, q10n 
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
	THREADSAFE 
}
 
PARAMETER {
        gnabar = .48 (S/cm2)	<0,1e9>
        gkbar = 1.088 (S/cm2)	<0,1e9>
        gl = .0016 (S/cm2)	<0,1e9>
        el = -60.0 (mV)
	q10m = 1
	q10n = 1
	q10h = 1
	}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	mtau (ms) htau (ms) ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 



? rates
PROCEDURE rates(v(mV)) {  
                      
        LOCAL	q10 
        TABLE minf, mtau, hinf, htau, ninf, ntau FROM -100 TO 100 WITH 200

UNITSOFF

                
	minf = 1/(1+exp(-0.4*(36+v)))
	mtau = 2*exp(-0.05*(v+40))
                
	htau = 40*exp(-0.025*(v+55))
        hinf = 1/(1+exp(39.5+v))
                
        ninf = 1/(1+exp(0.125*(-33-v)))
        ntau = 55*exp(-0.015*(v+28))
}
 

 
UNITSON