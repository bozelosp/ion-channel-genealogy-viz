UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

? interface

NEURON {
        SUFFIX hh_wbm
        NONSPECIFIC_CURRENT ina,ik,il

        RANGE gnabar,gna,egna,m, gkbar,gk,egk, gl,el
	GLOBAL hinf, ninf, htau, ntau

}
 
PARAMETER {
        gnabar = .035 (mho/cm2)	<0,1e9>
	egna	= 55 (mV)	
        gkbar = .009 (mho/cm2)	<0,1e9>
	egk	= -90 (mV)	
        gl = .0001 (mho/cm2)	<0,1e9>
        el = -65 (mV)
	
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
	celsius (degC)

	gna (mho/cm2)
        ina (mA/cm2)
	gk (mho/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
	htau (ms) ntau (ms)
}
 
LOCAL mexp, hexp, nexp        
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
	m = minf
        gna = gnabar*m*m*m*h
	ina = gna*(v - egna)
        gk = gkbar*n*n*n*n
	ik = gk*(v - egk)      
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
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  
                          
		      
        LOCAL  alpha, beta, sum
        TABLE minf, hinf, htau, ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 200


UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

               
        alpha = .1 * vtrap(-(v+35),10)
        beta =  4 * exp(-(v+60)/18)
        sum = alpha + beta
        minf = alpha/sum

                
        alpha =.35 * exp(-(v+58)/20)
        beta = 5 / (exp(-(v+28)/10) + 1)
        sum = alpha + beta
	htau = 1/(q10*sum)
        hinf = alpha/sum

                
        alpha =.05*vtrap(-(v+34),10) 
        beta = .625*exp(-(v+44)/80)
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