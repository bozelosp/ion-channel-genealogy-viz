NEURON {
  POINT_PROCESS gradNMDA
  RANGE gmax, g, i, alpha, beta, thetasyn,e, sigma
  GLOBAL mg
  NONSPECIFIC_CURRENT i
  POINTER vpre
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  gmax  = 1e-3 (uS)
  alpha = 0.0163 (/ms)   
  beta  = 0.00292 (/ms)  
  e     = 0	  (mV)    
  thetasyn = 0 (mV)   
  mg       = 1   (mM) 
  sigma    = 2  
}

ASSIGNED { vpre (mV)
           v (mV) 
		   i (nA)
		   g (uS)
	       B       
}

STATE { s }

INITIAL {
  s =  alpha*F(vpre)/(alpha*F(vpre)+beta)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  B = mgblock(v)
  g = gmax*s*B
  i = g*(v - e)
}

DERIVATIVE state {
  s' = alpha*F(vpre)*(1-s) - beta*s
}

FUNCTION F (v1 (mV)) {
  F = 1/(1 + exp(-(v1-thetasyn)/sigma))
}  

FUNCTION mgblock(v(mV)) {
        TABLE 
        DEPEND mg
        FROM -140 TO 80 WITH 1000
        
        mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}