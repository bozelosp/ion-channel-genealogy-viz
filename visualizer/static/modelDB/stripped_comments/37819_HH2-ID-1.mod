INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX hh2ad
  USEION na READ ena WRITE ina
  USEION k READ ek WRITE ik
  RANGE gnabar, gkbar, vtraub, ikhh2, inahh2
  RANGE m_inf, h_inf, n_inf
  RANGE tau_m, tau_h, tau_n
  RANGE m_exp, h_exp, n_exp
}


UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
}

PARAMETER {
  gnabar  = .003  (mho/cm2)
  gkbar   = .005  (mho/cm2)

  ena        (mV)
  ek        (mV)
  celsius    (degC)
  dt              (ms)
  v               (mV)
  vtraub  = -63   (mV)
}

STATE {
  m h n
}

ASSIGNED {
  ina     (mA/cm2)
  ik      (mA/cm2)
  il      (mA/cm2)
  inahh2     (mA/cm2)
  ikhh2      (mA/cm2)
  m_inf
  h_inf
  n_inf
  tau_m
  tau_h
  tau_n
  m_exp
  h_exp
  n_exp
  tadj
}


BREAKPOINT {
  SOLVE state METHOD cnexp
  inahh2 = gnabar * m*m*m*h * (v - ena)
  ikhh2  = gkbar * n*n*n*n * (v - ek)
  ina =   inahh2
  ik  =    ikhh2 
}


DERIVATIVE state {   
  evaluate_fct(v)
  m' = (m_inf - m) / tau_m
  h' = (h_inf - h) / tau_h
  n' = (n_inf - n) / tau_n
}

UNITSOFF
INITIAL {
  
  
  
  tadj = 3.0 ^ ((celsius-36)/ 10 )
  evaluate_fct(v)

  m = m_inf
  h = h_inf
  n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

  v2 = v - vtraub 

  
  a = 0.32 * vtrap(13-v2, 4)
  
  b = 0.28 * vtrap(v2-40, 5)
  tau_m = 1 / (a + b) / tadj
  m_inf = a / (a + b)

  a = 0.128 * Exp((17-v2)/18)
  b = 4 / ( 1 + Exp((40-v2)/5) )
  tau_h = 1 / (a + b) / tadj
  h_inf = a / (a + b)

  
  a = 0.032 * vtrap(15-v2, 5)
  b = 0.5 * Exp((10-v2)/40)
  tau_n = 1 / (a + b) / tadj
  n_inf = a / (a + b)

  m_exp = 1 - Exp(-dt/tau_m)
  h_exp = 1 - Exp(-dt/tau_h)
  n_exp = 1 - Exp(-dt/tau_n)
}
FUNCTION vtrap(x,y) {
  if (fabs(x/y) < 1e-6) {
    vtrap = y*(1 - x/y/2)
  }else{
    vtrap = x/(Exp(x/y)-1)
  }
}

FUNCTION Exp(x) {
  if (x < -100) {
    Exp = 0
  }else{
    Exp = exp(x)
  }
}