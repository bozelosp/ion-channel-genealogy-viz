UNITS {
        (mA) =		(milliamp)
        (mV) =		(millivolt)
        (uF) = 		(microfarad)
	(molar) = 	(1/liter)
	(nA) = 		(nanoamp)
	(mM) = 		(millimolar)
	(um) = 		(micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
SUFFIX tca

USEION ca READ eca WRITE ica
RANGE gtca
RANGE gcatbar
RANGE ainf, atau, binf, btau, itca
GLOBAL eca
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
	gcatbar = 1.0 (mho/cm2)
}
 
STATE {
	a b
}
 
ASSIGNED {
        gtca (mho/cm2)
	ica (mA/cm2)
	eca (mV)

	ainf binf
	atau (ms) btau (ms) 
	aexp bexp      
} 

BREAKPOINT {
	SOLVE states
        gtca = gcatbar*a*a*b
	ica = gtca*(v-eca)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	a = ainf
	b = binf
}

PROCEDURE states() {	
        trates(v)	
	a = a + aexp*(ainf-a) 
	b = b + bexp*(binf-b)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

PROCEDURE rates(v) {  
                      
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10) 
                
        alpha = -0.2*vtrap(v-19.26,-10)		
	beta = 0.009*exp(-v/22.03)		
	sum = alpha+beta        
	atau = 1/sum      ainf = alpha/sum
                
	alpha = 1e-6*exp(-v/16.26)		
	beta = 1/(exp((29.79-v)/10)+1)		
	sum = alpha+beta        
	btau = 1/sum      binf = alpha/sum
}
 
PROCEDURE trates(v) {  
                      
	LOCAL tinc
        TABLE  ainf, aexp, binf, bexp, atau, btau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	
		
		

	       tinc = -dt * q10
	aexp = 1 - exp(tinc/atau)
	bexp = 1 - exp(tinc/btau)
}
 
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON