NEURON {
  POINT_PROCESS WC
  RANGE ii0,ii1,ie0,ie1,w,aee,aie,aei,aii,ze,zi,tau,i_e,i_i
}
 
PARAMETER {
  aee = 10
  aie = 8
  aei = 12
  aii = 3
  ze = 0.2
  zi = 4
  tau = 1
  ie0 = 0
  ie1 = 0
  ii0 = 0
  ii1 = 0
  w = 0.25
}
 
STATE { uu vv }
 
ASSIGNED {
  i_e
  i_i  
}
 
BREAKPOINT {
  SOLVE states METHOD derivimplicit
}
 
INITIAL {
  uu = 0.1
  vv = 0.05
}

DERIVATIVE states {  
  i_e = ie0+ie1*sin(w*t)
  i_i = ii0+ii1*sin(w*t)
  uu' = -uu+f(aee*uu-aie*vv-ze+i_e)
  vv' = (-vv+f(aei*uu-aii*vv-zi+i_i))/tau
}

FUNCTION f(x) {
  f = 1/(1+exp(-x))
}