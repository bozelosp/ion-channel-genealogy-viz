INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX iar
  
  NONSPECIFIC_CURRENT i
  USEION ca READ cai
  RANGE ghbar, m, o1, o2, p0, p1, k2, alpha, beta
  GLOBAL cac, k4, Pc, nca, nexp, ginc, qt, origtemp,eh
}

UNITS {
  (molar)	= (1/liter)
  (mM)	= (millimolar)
  (mA) 	= (milliamp)
  (mV) 	= (millivolt)
  (msM)	= (ms mM)
}

PARAMETER {
  
  celsius = 37	(degC)
  ghbar	= 2e-5 (mho/cm2)
  cac = 0.006 (mM)		
  k2 = 0.0001 (1/ms)		
  Pc = 0.01			
  k4 = 0.001	(1/ms)		
  nca = 4			
  nexp= 1			
  ginc= 2			
  q10 = 2.2
  origtemp = 26 
  qt = 1.2668546920110242 
}

STATE {
  c1	
  o1	
  o2	
  p0	
  p1	
}

ASSIGNED {
  eh (mV)
  v	(mV)
  cai	(mM)
  i	(mA/cm2)
  gh	(mho/cm2)
  alpha	(1/ms)
  beta	(1/ms)
  k1ca	(1/ms)
  k3p	(1/ms)
  m
  tadj
}

BREAKPOINT {
  SOLVE ihkin METHOD sparse
  m = o1 + ginc * o2
  i = ghbar * m * (v - eh)
}

KINETIC ihkin {








  evaluate_fct(v,cai)
  ~ c1 <-> o1		(alpha,beta)
  ~ p0 <-> p1		(k1ca,k2)
  ~ o1 <-> o2		(k3p,k4)
  CONSERVE p0 + p1 = 1
  CONSERVE c1 + o1 + o2 = 1
}

INITIAL {




  
  
  
  qt = q10^((celsius-origtemp)/10)
  evaluate_fct(v,cai)
  c1 = 1
  o1 = 0
  o2 = 0
  p0 = 1
  p1 = 0
}

UNITSOFF
PROCEDURE evaluate_fct(v (mV), cai (mM)) {
  alpha = qt / exp(9.63 + 0.0458 * v) 
  beta = qt / exp(1.30  - 0.0447 * v) 
  k1ca = k2 * (cai/cac)^nca
  k3p = k4 * (p1/Pc)^nexp
}



PROCEDURE activation(v (mV), cai (mM)) { LOCAL cc
  evaluate_fct(v,cai)
  cc = 1 / (1 + (cac/cai)^nca ) 		
  m = 1 / ( 1 + beta/alpha + (cc/Pc)^nexp )
  m = ( 1 + ginc * (cc/Pc)^nexp ) * m
}
UNITSON