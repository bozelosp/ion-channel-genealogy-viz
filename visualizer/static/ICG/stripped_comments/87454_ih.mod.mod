UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX ih
        NONSPECIFIC_CURRENT i
        RANGE ghbar, gh, ih
        GLOBAL rinf, rtau
	GLOBAL eh
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)
        dt (ms)
        ghbar = 0.00318 (mho/cm2) <0,1e9>
        
}

STATE {
        r
}

ASSIGNED {
        eh (mV)
	gh (mho/cm2)
	i (mA/cm2)
	rinf
    rtau (ms)
}

LOCAL rexp

BREAKPOINT {
	SOLVE states
    
	gh = ghbar*r
    i = gh*(v - eh)
    }

UNITSOFF

INITIAL {
    trates(v)
    r = rinf
}

PROCEDURE states() {  
	trates(v)      
	r = r + rexp*(rinf-r)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10
PROCEDURE rates(v) {  
                      

	q10 = 3^((celsius - 22)/10)
    rinf = 1 / (1+exp((v + 76) / 7))
    rtau = (100000 / (237*exp((v+60) / 12) + 17*exp(-(v+60) / 14))) + 25

}

PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE rinf, rexp
	DEPEND dt, celsius FROM -200 TO 150 WITH 350

    rates(v)    
        
        

	tinc = -dt * q10
	rexp = 1 - exp(tinc/rtau)
}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON