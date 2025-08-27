NEURON {
	SUFFIX calcb
	USEION ca READ cai, eca WRITE ica
        RANGE gcalbar, ica, po
	GLOBAL inf, s_inf, tau_m
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {     
  	ki     = 0.025  (mM)            
	gcalbar = 0   (mho/cm2)  
 	taumin  = 180    (ms)            
        vhalf = -1 (mV)       
	zeta=-4.6
	t0=1.5(ms)
	b = 0.01 	(mM) 
        ba = 0.01	(mM)
	bo = 8
}


ASSIGNED {      
        v               (mV)
 	celsius         (degC)
	cai             (mM)      
	ica             (mA/cm2)
	eca             (mV)

	po
        inf
	s_inf
	tau_m           (ms)
}

STATE {	
	m 
	s 
} 


INITIAL {
	rates(v,cai)
	m = inf    
	s = s_inf
}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	po = m*m*h2(cai)
	ica = gcalbar*(po+s*s*bo)*(v-eca)
}


DERIVATIVE states {
	rates(v,cai)
	m' = (inf-m)/t0
	s' = (s_inf-s)/tau_m
}



FUNCTION alp(v(mV)) {       
UNITSOFF
  alp = exp(1.e-3*zeta*(v-vhalf)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

PROCEDURE rates(v(mV), cai(mM)) {LOCAL a, alpha2
		a = alp(v)
		inf = 1/(1+a)
		alpha2 = (cai/b)^2
		s_inf = alpha2 / (alpha2 + 1)
		tau_m = taumin+ 1(ms)*1(mM)/(cai+ba)
}