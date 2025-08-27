NEURON {
  SUFFIX bkr
  USEION ip3 WRITE iip3 VALENCE 1
  NONSPECIFIC_CURRENT ix
  RANGE jbar, beta, j, i
  GLOBAL del
}

UNITS {
  (um)    = (micron)
  (molar) = (1/liter)
  (uM)    = (micromolar)
  (mA)	  = (milliamp)
  FARADAY = (faraday)  (coulomb)
}

PARAMETER {
  del (ms)  

  k = 1.188e-3 (/ms) 

  jbar = 20.86 (uM um/s) 
  beta = 1 (1) 
}

ASSIGNED {
  jip3 (micro/um2 ms) 
  iip3 (mA/cm2) 
  ix   (mA/cm2) 
}

INITIAL {
  jip3 = 0
  iip3 = 0
  ix = 0
}

BREAKPOINT {
  at_time(del) 
  at_time(del + (10/k)) 


  if (t > del && t <= del + (10/k)) {
    jip3 = (1e-18)*beta*jbar*exp(-k*(t-del))
  }else{
    jip3 = 0
  }

  iip3 = -(1e8)*jip3*FARADAY
  ix = -iip3
}