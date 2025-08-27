NEURON {
  POINT_PROCESS Izhi2003a
  RANGE a,b,c,d,f,g,Iin,fflag,thresh,erev,taug
}

INITIAL {
  V=-65
  u=0.2*V
  gsyn=0
  net_send(0,1)
}

PARAMETER {
  a = 0.02
  b = 0.2
  c = -65
  d = 2
  f = 5
  g = 140
  Iin = 10
  taug = 1
  thresh=30
  erev = 0
  fflag = 1
}

STATE { u V gsyn } 

ASSIGNED {
}

BREAKPOINT {
  SOLVE states METHOD derivimplicit
}

DERIVATIVE states {
  V' = 0.04*V*V + f*V + g - u + Iin - gsyn*(V-erev)
  u' = a*(b*V-u) 
  gsyn' = -gsyn/taug
}

NET_RECEIVE (w) {
  if (flag == 1) {
    WATCH (V>thresh) 2
  } else if (flag == 2) {
    net_event(t)
    V = c
    u = u+d
  } else { 
    gsyn = gsyn+w
  }
}


PROCEDURE version () {

}