UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

? interface

NEURON {
        SUFFIX INap
        USEION na READ ena WRITE ina
	RANGE gnabar, gna, ena, ina
	GLOBAL htau
}
 
PARAMETER {
        gnabar  = .00015 (mho/cm2)	<0,1e9>
}
 
STATE {
        h
}
 
ASSIGNED {
        v 	(mV)
	celsius (degC)
	gna 	(mho/cm2)
        ina 	(mA/cm2)
        minf
	hinf
	htau 	(ms)
	ena	(mV)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*minf*minf*minf*h
	ina = gna*(v - ena)
}


INITIAL {
	rates(v)
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
        TABLE minf, hinf, htau FROM -100 TO 100 WITH 200

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

        
        minf = 1 / ( 1 + exp(-(v+55.7)/7.7) )

        
        alpha = .001 * exp(-(v+85)/30)
        beta =  0.0034 / ( 1 + exp(-(v+17)/10) )
        sum = alpha + beta
	htau = 1/sum
	hinf = alpha/sum
}
 
 
UNITSON