NEURON {
  SUFFIX acaext
  USEION ca READ ica WRITE cao
  GLOBAL cabath
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
  cabath   =  2 (mM)        
  fhspace = 1e-4 (cm)  
  txfer   =  4511.0 (ms)  
  SA = 2.82743E-05 (cm2) 
  Vol_peri = 1.46136E-09 (cm3)
}
ASSIGNED { ica  (mA/cm2) }
STATE { cao  (mM) }
BREAKPOINT { SOLVE state METHOD cnexp }
DERIVATIVE state {
  cao' = ica*SA/(2*Vol_peri*FARADAY) + (cao - cabath)/txfer
}