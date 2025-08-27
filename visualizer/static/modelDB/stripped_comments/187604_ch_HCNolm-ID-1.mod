UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX ch_HCNolm
        USEION h READ eh WRITE ih VALENCE 1
        RANGE gmax,ih, g
        GLOBAL rinf, rexp, tau_r
        RANGE myi
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gmax = 0.001385 (mho/cm2)
		g (mho/cm2)
        eh = -32.9 (mV)
}
 
STATE {
        r
}
 
ASSIGNED {
        ih (mA/cm2)
	rinf rexp
	tau_r
	myi (mA/cm2)
}
 
BREAKPOINT {
	SOLVE deriv METHOD derivimplicit
	g = gmax*r
	ih = g*(v - eh)
	myi = ih
}
 
INITIAL {
	rates(v)
	r = rinf
}

DERIVATIVE deriv { 
	rates(v)
	r' = (rinf - r)/tau_r
}


PROCEDURE rates(v) {  
                      
        TABLE rinf, tau_r, rexp DEPEND dt FROM -200
TO 100 WITH 300
	rinf = 1/(1 + exp((v+84.1)/10.2))
	tau_r = 100 + 1/(exp(-17.9-0.116*v)+exp(-1.84+0.09*v))
	rexp = 1 - exp(-dt/(tau_r))
}

FUNCTION efun(z) {
  if (fabs(z) < 1e-4) {
    efun = 1 - z/2
  } else {
    efun = z/(exp(z) - 1)
  }
}

 
UNITSON