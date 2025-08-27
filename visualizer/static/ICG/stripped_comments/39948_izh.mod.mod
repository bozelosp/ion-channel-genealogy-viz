NEURON {
  POINT_PROCESS IZH
  RANGE a,b,c,d,e,f,g,I,fflag,thresh,erev,gsyn 
}

INITIAL {
  vv=-65
  u=0.2*vv
  gsyn=0
  net_send(0,1)
}

PARAMETER {
  a = 0.02
  b = 0.2
  c = -65
  d = 2
  e = 0.04
  f = 5
  g = 140
  I = 10
  thresh=30
  erev = 0
  fflag = 1
}

STATE { u vv } 

ASSIGNED {
  gsyn
}

BREAKPOINT {
  SOLVE states METHOD derivimplicit
}

DERIVATIVE states {
  vv' = e*vv*vv + f*vv + g - u + I - gsyn*(vv-erev)
  u' = a*(b*vv-u) 
  
  
  
}

NET_RECEIVE (w) {
  if (flag == 1) {
    WATCH (vv>thresh) 2
  } else if (flag == 2) {
    net_event(t)
    vv = c
    u = u+d
  } else { 
    gsyn = gsyn+w
  }
}


PROCEDURE version () {
  printf("$Id
}