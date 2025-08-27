UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX HH
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar,gkbar,gl,el,m_inf,h_inf,n_inf,tau_m,tau_n,tau_h,ina,ik 
}
 
PARAMETER {
        v		(mV)
        celsius		(degC)
        dt		(ms)
        gnabar= 0.1	(mho/cm2)
        ena		(mV)
        gkbar= 0.01	(mho/cm2)
        ek		(mV)
        gl= 0		(mho/cm2)
        el		(mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        ina	(mA/cm2)
        ik	(mA/cm2)
        il	(mA/cm2)
        m_inf h_inf n_inf tau_m tau_h tau_n
	tadj
}
 
BREAKPOINT {
        SOLVE states
        ina = gnabar * m*m*m*h * (v - ena)
        ik = gkbar * n*n*n*n * (v - ek)      
        il = gl * (v-el)
}
 
UNITSOFF
INITIAL {
	tadj = 3.0^((celsius-23.5)/10)
	rates(v)
	m = m_inf
	h = h_inf
	n = n_inf
}

PROCEDURE states() { 
        rates(v)   
        m = m + (1-exp(-dt/tau_m)) * (m_inf-m)
        h = h + (1-exp(-dt/tau_h)) * (h_inf-h)
        n = n + (1-exp(-dt/tau_n)) * (n_inf-n)
}
 
PROCEDURE rates(v) { LOCAL alpha, beta, q10, tinc


	
          alpha = .091 * vtrap(v+38,5)
          beta =  .062 * vtrap(-(v+38),5) 
       	  tau_m = 1 / (alpha+beta) / tadj
       	  m_inf = alpha/(alpha+beta)
	
       	  alpha = .016 * exp(-(v+55)/15)
       	  beta = 2.07 / (1+exp((17-v)/21))
       	  tau_h = 1 / (alpha+beta) / tadj
       	  h_inf = alpha/(alpha+beta)
	
       	  alpha = .01*vtrap(v+45,5) 
       	  beta = .17*exp(-(v+50)/40)
       	  tau_n = 1 / (alpha+beta) / tadj
       	  n_inf = alpha/(alpha+beta)
}
 
FUNCTION vtrap( x, b) {
	
	if (fabs(x/b) < 1e-6) {
	  vtrap = b+x/2 }
	else {
	  vtrap = x / (1-exp(-x/b)) }
}
UNITSON