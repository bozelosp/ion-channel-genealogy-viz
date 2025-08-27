UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX IA
        USEION k READ ek WRITE ik
        RANGE gkAbar,ik
        GLOBAL ainf, binf, aexp, bexp, tau_b
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        p = 5 (degC)
        dt (ms)
        gkAbar = 0.0165 (mho/cm2)	
        ek = -90 (mV)
	tau_a = 5 (ms)
}
 
STATE {
        a b
}
 
ASSIGNED {
        ik (mA/cm2)
	ainf binf aexp bexp
	tau_b
}
 
BREAKPOINT {
        SOLVE deriv METHOD derivimplicit
        ik = gkAbar*a*b*(v - ek)
}
 
INITIAL {
	rates(v)
	a = ainf
	b = binf
}

DERIVATIVE deriv {  
		
        a' = (ainf - a)/(tau_a)
        b' = (binf - b)/(tau_b)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL alpha_b, beta_b
	TABLE ainf, aexp, binf, bexp, tau_a, tau_b  DEPEND dt, p FROM -200
TO 100 WITH 300
	alpha_b = 0.000009/exp((v-26)/18.5)
	beta_b = 0.014/(exp((v+70)/(-11))+0.2)
        ainf = 1/(1 + exp(-(v + 14)/16.6))
        aexp = 1 - exp(-dt/(tau_a))
	tau_b = 1/(alpha_b + beta_b)
        binf = 1/(1 + exp((v + 71)/7.3))
        bexp = 1 - exp(-dt/(tau_b))
}
 
UNITSON