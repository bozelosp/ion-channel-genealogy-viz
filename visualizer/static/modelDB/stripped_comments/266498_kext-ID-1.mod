NEURON {
  SUFFIX kext
  USEION k READ ik WRITE ko
  GLOBAL kbath
  RANGE fhspace, txfer
}
UNITS {
  (mV)    = (millivolt)
  (mA)    = (milliamp)
  FARADAY = (faraday) (coulombs)
  (molar) = (1/liter)
  (mM)    = (millimolar)
}
PARAMETER {
  kbath   =  5.4 (mM)        
  fhspace = 1e-4 (cm)  
  txfer   =  50 (ms)  
}
ASSIGNED { ik  (mA/cm2) }
STATE { ko  (mM) }
BREAKPOINT { SOLVE state METHOD cnexp }
DERIVATIVE state {
  ko' = ik/(fhspace*FARADAY) + (kbath - ko)/txfer
}