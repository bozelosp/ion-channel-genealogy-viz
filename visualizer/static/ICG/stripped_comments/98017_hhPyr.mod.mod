UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}

? interface

NEURON {
        SUFFIX hhPyr
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gna, ena, gkbar, gk, ek, gl, el, ina, ik, il
	GLOBAL hinf, ninf, htau, ntau
}
 
PARAMETER {
        gnabar	= 0.0 		(mho/cm2)	<0,1e9>
        gkbar 	= .015 		(mho/cm2)	<0,1e9>
        gl	= 0.0 	(mho/cm2)	<0,1e9>
        el 	= -70 		(mV)
}
 
STATE {
        h n
}
 
ASSIGNED {
        v (mV)
	celsius (degC)
	gna (mho/cm2)
        ina (mA/cm2)
	ena (mV)
	gk (mho/cm2)
        ik (mA/cm2)
	ek (mV) 
        il (mA/cm2)
        minf hinf ninf
	htau (ms) ntau (ms)
}
 
LOCAL mexp, hexp, nexp        
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*minf*minf*minf*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}


INITIAL {
	rates(v)
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {
        rates(v)
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  
                          
		      
        LOCAL  alpha, beta, sum
        TABLE minf, hinf, ninf, htau, ntau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

        
        alpha = .1 * vtrap(-(v+31),10)
        beta =  4 * exp(-(v+56)/18)
        sum = alpha + beta
        minf = alpha/sum

        
        alpha = .07 * exp(-(v+47)/20)
        beta = 1 / (exp(-(v+17)/10) + 1)
        sum = alpha + beta
	htau = 1/(q10*sum)
        hinf = alpha/sum

        
        alpha = .01*vtrap(-(v+34),10) 
        beta = .125*exp(-(v+44)/80)
	sum = alpha + beta
        ntau = 1/(q10*sum)
        ninf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON