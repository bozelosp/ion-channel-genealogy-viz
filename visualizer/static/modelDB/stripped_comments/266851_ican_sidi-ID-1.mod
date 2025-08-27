NEURON {
  SUFFIX ican
  NONSPECIFIC_CURRENT i
  USEION ca READ cai
  USEION na WRITE ina
  RANGE gbar, m_inf, tau_m
  GLOBAL beta, cac, taumin
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (molar) = (1/liter)
  (mM) = (millimolar)
}

PARAMETER {
  v               (mV)
  celsius = 36    (degC)
  erev      = -20   (mV)            	
  cai     	     (mM)           	
  gbar    = 0.0001 (mho/cm2)
  beta    = 0.0001 (1/ms) 	 	
  cac     = 0.0004 (mM)
  
  taumin  = 0.1   (ms)            	
}

STATE {
  m
}

ASSIGNED {
  i      (mA/cm2)
  ina     (mA/cm2)
  m_inf
  tau_m   (ms)
  tadj
  g (mho/cm2)
}

PROCEDURE iassign () {
  g = gbar * m * m
  i = g * (v - erev) 
  ina = 0.7 * i
}

BREAKPOINT { 
  SOLVE states METHOD cnexp  
  iassign()
}

DERIVATIVE states { 
  evaluate_fct(v,cai)  
  m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
  
  
  tadj = 3.0 ^ ((celsius-22.0)/10)  
  evaluate_fct(v,cai)
  m = m_inf
  iassign()
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2  
  alpha2 = beta * (cai/cac)^2  
  tau_m = 1 / (alpha2 + beta) / tadj
  m_inf = alpha2 / (alpha2 + beta)  
  if(tau_m < taumin) { tau_m = taumin } 
}
UNITSON