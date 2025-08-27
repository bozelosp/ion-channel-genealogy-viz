UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
		(S) = (siemens)
}
 
NEURON {
        SUFFIX Nap_No
        USEION na READ ena WRITE ina
        RANGE gnamax, gna
        RANGE pinf, ptau
}
 
PARAMETER { 
        gnamax = 0.01 (S/cm2)   <0,1e9>
}
 
STATE {
        p
}
 
ASSIGNED {
        v (mV)
        ena (mV)

		gna (S/cm2)
        ina (mA/cm2)
        pinf
		ptau (ms)
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnamax*p*p*p
		ina = gna*(v - ena)
}
 
INITIAL {
	rates(v)
	p = pinf
}

DERIVATIVE states {
        rates(v)
        p' = (pinf-p)/ptau
}
 
PROCEDURE rates(v(mV)) {  
		
        LOCAL  alpha, beta, sum

UNITSOFF
                
        alpha = 0.0353 * vtrap(-17 - v, 10.2)
        beta =  0.000883 * vtrap(v + 24, 10)
        sum = alpha + beta
		ptau = 1/sum
        pinf = alpha/sum
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON