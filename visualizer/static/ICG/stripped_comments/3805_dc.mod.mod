UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX dc
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, shift
        GLOBAL minf, hinf, ninf, mtau, htau, ntau, vrest
}
 
PARAMETER {
	
        gnabar = 0.0 (S/cm2)	<0,1e9>
        gkbar = .1 (S/cm2)	<0,1e9>
        gl = 0.0 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
	vrest = 0 (mV)
	shift = 0
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
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
	rates(v - vrest - shift)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v - vrest - shift)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 
? rates
PROCEDURE rates(v(mV)) {  
                      
        LOCAL  alpha, beta, sum

UNITSOFF
                
        alpha = .4 * vtrap(25 - v, 5)


        beta =  .4 * vtrap(v - 45, 5)
        sum = alpha + beta
	mtau = 1/sum
        minf = alpha/sum
                
        alpha = .28 * exp((10 - v)/20)
        beta = 4 / (exp((40 - v)/10) + 1)
        sum = alpha + beta
	htau = 1/sum
        hinf = alpha/sum
                


        alpha = .02*vtrap(20 - v,10) 
        beta = .25*exp((10 - v)/80)
	sum = alpha + beta
        ntau = 1/sum
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