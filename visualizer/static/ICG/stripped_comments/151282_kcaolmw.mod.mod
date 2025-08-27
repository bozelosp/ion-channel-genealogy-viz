UNITS {
  (mA)     = (milliamp)
  (mV)     = (millivolt)
  (mS)     = (millisiemens)
  (mollar) = (1/liter)
  (mM)     = (millimollar)
}

NEURON {
  SUFFIX KCaolmw
  USEION k READ ek WRITE ik
  USEION ca READ cai
  RANGE gkca,kd
}
	
PARAMETER {
  gkca =  10 (mS/cm2)
  
  kd   =  30 (mM)
}
    
ASSIGNED {  
  ek (mV)  
  cai (mM) 
  v   (mV)
  ik  (mA/cm2) 
}

PROCEDURE iassign () { ik  = (1e-3) * gkca * cai/(cai+kd) * (v-ek) }

INITIAL {
  iassign()
}

BREAKPOINT { iassign() }