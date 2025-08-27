UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (mS) = (millisiemens)
	(molar) = (1/liter)
	(uM)    = (micromolar)
}

NEURON {
    THREADSAFE
    
    
    SUFFIX gKahp
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gmax, g, i, q_inf, tau_q
    
}

PARAMETER {
    gmax=0.0008 (mho/cm2)  
    beta_q = 0.001
}

STATE {
    q
}

ASSIGNED {
    g		(S/cm2)
    ek		(mV)
    ik		(mA/cm2)
    q_inf	(1)
    tau_q	(ms)
    i		(mA/cm2)
    cai		(uM)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = gmax*q
    i = g*(v-ek)
    ik = i
}

INITIAL {
    	rates(cai)
    	q = q_inf
}

FUNCTION alpha_q(cai) {
    alpha_q = 0.2e-4*cai
    if ( alpha_q > 0.01) {
	alpha_q = 0.01
    }
}

DERIVATIVE states {  
    rates(cai)
    q' = (q_inf - q)/tau_q
}

PROCEDURE rates(cai) { 
    LOCAL alpha, beta
    TABLE q_inf, tau_q  
    FROM 0 TO 500 WITH 2000
    
    alpha = alpha_q(cai)
    beta = beta_q
    tau_q = 1/(alpha + beta)     
    q_inf = alpha*tau_q
}