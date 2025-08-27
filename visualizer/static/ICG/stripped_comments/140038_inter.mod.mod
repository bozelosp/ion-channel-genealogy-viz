UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
        SUFFIX inter
        USEION na READ ena WRITE ina
        RANGE gnabar, gna, ina
        GLOBAL minf, mtau, hinf, htau
}
 
PARAMETER {
        gnabar = 1.0 (S/cm2)	<0,1e9>            
}
 
STATE {
        m h
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        
	gna (S/cm2)	
	ina (mA/cm2)
        minf hinf
	mtau (ms) htau (ms)
}
 
LOCAL mexp, hexp
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        m = minf
        gna = gnabar*m*h
	ina = gna*(v - ena)	
}
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

? states
DERIVATIVE states {  
        rates(v)
        h' = (hinf-h)/htau      
}
 
LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  
                      
        LOCAL  alpha, beta, sum
        TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
	q10 = 3^((celsius - 6.3)/10)
                
        alpha = .1 * vtrap(-(v+25),10)
        beta =  4 * exp(-(v+50)/18)
        sum = alpha + beta
        minf = alpha/sum
                
	htau = 0.2218*exp(-0.06883*v)  
        hinf = (1+exp((v+72.5)/8))^-1  
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON