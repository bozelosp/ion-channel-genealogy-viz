UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX IA
        USEION k READ ek WRITE ik
        RANGE gkAbar,ik, tau_a, tau_b, tau_b_tempInsensitive
        RANGE ainf, binf, aexp, bexp
        GLOBAL qt
}
 
PARAMETER {
        v (mV)
        dt (ms)
        gkAbar = 0.0165 (mho/cm2)	
        ek = -90 (mV)
	tau_a_tempInsensitive = 5 (ms)
        Q10 = 3
        Q10TEMP = 24 (degC)
}
 
STATE {
        a b
}
 
ASSIGNED {
        ik (mA/cm2)
	ainf binf aexp bexp
	tau_b tau_a
        tau_b_tempInsensitive
        celsius (degC)
        qt (1)
}
 
BREAKPOINT {
        SOLVE deriv METHOD cnexp
        ik = gkAbar*a*b*(v - ek)
}
 
INITIAL {
        qt = Q10 ^ ((celsius - Q10TEMP) / 10)
        rates(v)
	a = ainf
	b = binf
}

DERIVATIVE deriv {  
		
        rates(v)
        a' = (ainf - a)/(tau_a)
        b' = (binf - b)/(tau_b)
}
 
PROCEDURE rates(v) {  
                      
        LOCAL alpha_b, beta_b
	alpha_b = 0.000009/exp((v-26)/18.5)
	beta_b = 0.014/(exp((v +70)/(-11))+0.2)
        ainf = 1/(1 + exp(-(v + 14)/16.6))
        tau_a = tau_a_tempInsensitive / qt
        aexp = 1 - exp(-dt/(tau_a))
        tau_b_tempInsensitive = 1/(alpha_b + beta_b)
	tau_b = tau_b_tempInsensitive / qt
        binf = 1/(1 + exp((v + 71)/7.3))
        bexp = 1 - exp(-dt/(tau_b))
}
 
UNITSON