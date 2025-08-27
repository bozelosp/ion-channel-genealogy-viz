UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mS) = (millisiemens)
}

NEURON {
  SUFFIX ICaolmw
  USEION ca READ eca WRITE ica
  RANGE gca
}

PARAMETER {
  gca = 1    (mS/cm2)
  
}
    
ASSIGNED { 
  eca (mV)
  ica (mA/cm2)    
  v   (mV)
}

PROCEDURE iassign () { ica = (1e-3) * gca * mcainf(v)^2 * (v-eca) }

INITIAL {
  iassign()
}

BREAKPOINT { iassign() }

FUNCTION mcainf(v(mV)) { mcainf = fun2(v, -20, 1,  -9)*1(ms) }

INCLUDE "custom_code/inc_files/151282_aux_fun.inc"