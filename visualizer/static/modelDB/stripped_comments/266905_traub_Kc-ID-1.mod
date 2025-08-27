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
    
    
    SUFFIX gKc
    USEION ca READ cai
    USEION k READ ek WRITE ik
    RANGE gmax, g, i, c_inf, tau_c, cxi
	
}

PARAMETER {
    gmax=0.01 (mho/cm2)  
}

STATE {
    c
}

ASSIGNED {
    v		(mV)
    ek		(mV)
    ik		(mA/cm2)
    c_inf	(1)
    cxi		(1)
    tau_c	(ms)
    g		(S/cm2)
    i		(mA/cm2)
    cai		(uM)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    cxi = c*xi(cai)
    g = gmax*cxi
    i = g*(v-ek)
    ik = i
}

INITIAL {
    rates(v)
    c = c_inf
}

FUNCTION xi(cai) {
    xi = cai/250
    if ( xi > 1 ) {
		xi = 1
    }
}

FUNCTION alpha_c(v(mV)) {
    if ( v <= 50 ) {
		alpha_c = exp( (v - 10)/11 - (v - 6.5)/27 )/18.975
    } else {
		alpha_c = 2*exp( -(v - 6.5)/27 )
    }
}

FUNCTION beta_c(v(mV)) {
    if ( v<= 50 ) {
		beta_c = 2*exp( -(v - 6.5)/27 ) - alpha_c(v)
    } else {
		beta_c = 0
    }
}

DERIVATIVE states {  
    rates(v)
    c' = (c_inf - c)/tau_c
}

PROCEDURE rates(v (mV)) { 
    LOCAL alpha, beta
    TABLE c_inf, tau_c  
    FROM -20 TO 130 WITH 750
    
    alpha = alpha_c(v)
    beta = beta_c(v)
    tau_c = 1/(alpha + beta)     
    c_inf = alpha*tau_c
}