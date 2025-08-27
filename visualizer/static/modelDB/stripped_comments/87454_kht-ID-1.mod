UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX kht

        USEION k WRITE ik
        RANGE gkhtbar, gkht, ik
        GLOBAL ninf, pinf, ntau, ptau
	RANGE ek
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  
        dt (ms)
        ek = -77 (mV)
        gkhtbar = 0.01592 (mho/cm2) <0,1e9>
		nf = 0.85 <0,1> 
}

STATE {
        n p
}

ASSIGNED {
    ik (mA/cm) 
    gkht (mho/cm2)
    pinf ninf
    ptau (ms) ntau (ms)
    }

LOCAL nexp, pexp

BREAKPOINT {
	SOLVE states
    
	gkht = gkhtbar*(nf*(n^2) + (1-nf)*p)
    ik = gkht*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    p = pinf
    n = ninf
}

PROCEDURE states() {  
	trates(v)      
	n = n + nexp*(ninf-n)
	p = p + pexp*(pinf-p)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  
                      

	q10 = 3^((celsius - 22)/10) 

    ninf =   (1 + exp(-(v + 15) / 5))^-0.5
    pinf =  1 / (1 + exp(-(v + 23) / 6))

	ntau =  (100 / (11*exp((v+60) / 24) + 21*exp(-(v+60) / 23))) + 0.7
    ptau = (100 / (4*exp((v+60) / 32) + 5*exp(-(v+60) / 22))) + 5
}

PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE ninf, nexp, pinf, pexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    
        
        

	tinc = -dt * q10
	nexp = 1 - exp(tinc/ntau)
	pexp = 1 - exp(tinc/ptau)
	}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON