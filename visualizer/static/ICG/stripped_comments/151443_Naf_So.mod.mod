UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		(S) = (siemens)
}
 
NEURON {
        SUFFIX Naf_So
        USEION na READ ena WRITE ina
        RANGE gnamax, gna
        RANGE minf, hinf, mtau, htau
}
 
PARAMETER {
        gnamax = .1 (S/cm2)   <0,1e9>
}
 
STATE {
        m h
}
 
ASSIGNED {
        v (mV)
        ena (mV)

		gna (S/cm2)
        ina (mA/cm2)
        minf hinf
		mtau (ms) htau (ms)
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnamax*m*m*m*h
		ina = gna*(v - ena)
} 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {
        rates(v)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}
 
PROCEDURE rates(v(mV)) {  
		
        LOCAL  alpha, beta, sum

UNITSOFF
                
        alpha = .4 * vtrap(-66 - v, 5)
        beta =  .4 * vtrap(v + 32, 5)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum
                
		htau = 30/(exp((60+v)/15)+exp(-(60+v)/16))
        hinf = 1/(1+exp((65+v)/7))
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON