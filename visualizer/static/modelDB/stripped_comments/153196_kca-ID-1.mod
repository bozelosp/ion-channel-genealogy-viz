NEURON {
  THREADSAFE
  SUFFIX kca
  USEION k READ ek WRITE ik
  USEION ca READ cai
  RANGE n, gk, gbar
  RANGE ninf, ntau
  GLOBAL Ra, Rb, caix
  GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (um) = (micron)
  (mM) = (milli/liter)
}

PARAMETER {
  gbar = 10 (pS/um2)  


  caix = 1

  Ra   = 0.01  (/ms)    
  Rb   = 0.02  (/ms)    



  temp = 23  (degC)    
  q10  = 2.3      

  vmin = -120  (mV)
  vmax = 100  (mV)
}

ASSIGNED {
  v     (mV)
  celsius (degC)
  cai   (mM)
  a     (/ms)
  b     (/ms)
  ik    (mA/cm2)
  gk    (pS/um2)
  ek    (mV)
  ninf
  ntau  (ms)
  tadj
}

STATE { n }

INITIAL {
  
  
  tadj = q10^((celsius - temp)/(10 (degC)))
  rates(cai)
  n = ninf
}

BREAKPOINT {
  SOLVE states METHOD cnexp

  gk = gbar*n
  ik = (1e-4) * gk * (v - ek)
}



DERIVATIVE states { 
  rates(cai)        
  n' =  (ninf-n)/ntau
  
}

UNITSOFF
PROCEDURE rates(cai(mM)) {
  a = Ra * cai^caix
  b = Rb


  tadj = q10^((celsius - temp)/(10 (degC)))
  ntau = 1/tadj/(a+b)
  ninf = a/(a+b)
  
  
  
  
}
UNITSON