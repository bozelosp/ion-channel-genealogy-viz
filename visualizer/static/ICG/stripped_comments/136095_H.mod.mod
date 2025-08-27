UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX H
    GLOBAL eh, VhlfMaxm, minf
    
    NONSPECIFIC_CURRENT i
    RANGE g,gmax,slopem,taum
    THREADSAFE
}

ASSIGNED { 
  eh (mV)
  i (mA/cm2)
  v (mV)
  g (mho/cm2)
  minf
  iother (mA/cm2)
}

STATE {
  m
}

PARAMETER {
  
  gmax  = 5e-07 (mho/cm2)
  VhlfMaxm = -74
  slopem = -10
  taum = 50 (ms)
}

PROCEDURE iassign () {
  i = g*(v-eh)
  iother = i
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g = m * gmax
    iassign()
}

INITIAL {
    settables(v)
    m = minf
    g = m * gmax
    iassign()
}

DERIVATIVE states {  
    settables(v)      
    m' = ( minf - m ) / taum 
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE minf
          FROM -200 TO 200 WITH 401
    minf = 1 / (1 + exp((VhlfMaxm - v)/ slopem ) )
}

UNITSON