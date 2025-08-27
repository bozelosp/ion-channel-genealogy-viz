UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
}

NEURON {
  SUFFIX Iholmw
  NONSPECIFIC_CURRENT i
  RANGE gbar,gh,gfactor
  GLOBAL eh
}
	
PARAMETER {
  gbar = 0.00015 (mho/cm2)
  
  gfactor = 1
}
    
ASSIGNED { 
  eh (mV)
  v (mV)
  i (mA/cm2)
  gh (mho/cm2)
}

STATE { q }

PROCEDURE giassign () { 
  
  gh = gbar * q * gfactor 
  i =  gh * (v-eh)
}

INITIAL { 
  q  = qinf(v) 
  giassign()
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  giassign()
}

DERIVATIVE states { q' = (qinf(v)-q)/qtau(v) }

FUNCTION qinf(v(mV))     { qinf = fun2(v, -80, 1, 10)*1(ms) }
FUNCTION qtau(v(mV))(ms) { qtau = 200(ms)/(exp((v+70(mV))/20(mV))+exp(-(v+70(mV))/20(mV))) + 5(ms) }

INCLUDE "custom_code/inc_files/151282_aux_fun.inc"