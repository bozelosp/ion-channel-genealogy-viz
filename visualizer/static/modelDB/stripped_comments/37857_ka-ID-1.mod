UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gkabar, gka, ik
        GLOBAL ainf, binf, cinf, atau, btau, ctau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 22 (degC)  
        dt (ms)
        ek = -77 (mV)
        gkabar = 0.00477 (mho/cm2) <0,1e9>
}

STATE {
        a b c
}

ASSIGNED {
    ik (mA/cm2) 
    gka (mho/cm2)
    ainf binf cinf
    atau (ms) btau (ms) ctau (ms)
    }

LOCAL aexp, bexp, cexp

BREAKPOINT {
	SOLVE states
    
	gka = gkabar*(a^4)*b*c
    ik = gka*(v - ek)

}

UNITSOFF

INITIAL {
    trates(v)
    a = ainf
    b = binf
    c = cinf
}

PROCEDURE states() {  
	trates(v)      
	a = a + aexp*(ainf-a)
	b = b + bexp*(binf-b)
	c = c + cexp*(cinf-c)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10

PROCEDURE rates(v) {  
                      

	q10 = 3^((celsius - 22)/10) 

    ainf = (1 / (1 + exp(-1*(v + 31) / 6)))^0.25
    binf = 1 / (1 + exp((v + 66) / 7))^0.5
    cinf = 1 / (1 + exp((v + 66) / 7))^0.5

    atau =  (100 / (7*exp((v+60) / 14) + 29*exp(-(v+60) / 24))) + 0.1
    btau =  (1000 / (14*exp((v+60) / 27) + 29*exp(-(v+60) / 24))) + 1
    ctau = (90 / (1 + exp((-66-v) / 17))) + 10
}

PROCEDURE trates(v) {  
                      
	LOCAL tinc
	TABLE ainf, aexp, binf, bexp, cinf, cexp
	DEPEND dt, celsius FROM -150 TO 150 WITH 300

    rates(v)    
        
        

	tinc = -dt * q10
	aexp = 1 - exp(tinc/atau)
	bexp = 1 - exp(tinc/btau)
	cexp = 1 - exp(tinc/ctau)
	}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON