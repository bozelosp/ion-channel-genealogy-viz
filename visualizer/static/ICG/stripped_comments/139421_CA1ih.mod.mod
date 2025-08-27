UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
}
 
NEURON {
  SUFFIX hcurrent
  NONSPECIFIC_CURRENT i
  RANGE g, v50, htau, hinf
  RANGE gfactor
  GLOBAL eh
}
 
PARAMETER {
  celsius	(degC)
  g= 0.0001	(mho/cm2)
  
  v50=-82	(mV)
  gfactor = 1
}
 
STATE {
  h
}
 
ASSIGNED {
  eh (mV)
  i	  (mA/cm2) 
  hinf
  htau    (ms)
  v	  (mV)
}

PROCEDURE iassign () { i=g*h*(v-eh)*gfactor }
 
BREAKPOINT {
  SOLVE states METHOD cnexp
  iassign()
}
 
DERIVATIVE states { 
  rates(v)
  h'= (hinf- h)/ htau
}

INITIAL { 
  rates(v)
  h = hinf
  iassign()
}

PROCEDURE rates(v (mV)) {
  UNITSOFF
  
  
  

  
  hinf = 1/(1+exp((v-v50)/10.5))
  htau = (1/(exp(-14.59-0.086*v)+exp(-1.87+0.0701*v))) 
  UNITSON
}