NEURON {
  POINT_PROCESS ZoidSyn
  RANGE trf, tp
  RANGE start, interval, number
  RANGE e, gmax, g, i
  NONSPECIFIC_CURRENT i
}

UNITS {
  (mV) = (millivolt)
  (nS) = (nanosiemens)
  (nA) = (nanoamp)
}

PARAMETER {
  trf (ms) <0, 1e9> 
  tp  (ms) <0, 1e9> 
  start (ms) <0, 1e9> 
  interval (ms) <0, 1e9> 
  number 
  e   (mV) 
  gmax (nS) <0, 1e9> 
}

ASSIGNED {
  v (mV)
  i (nA)
  on
  tally 
  m (1/ms)
  b (1)
  dur (ms) 
  t0 (ms)
  g (nS)
}

INITIAL {
  if (trf <= 0) {
    trf = 1
UNITSOFF
    printf("trf must be longer than 0\n")
    printf("trf has been increased to %g ms\n", trf)
UNITSON
  }
  if (tp < 0) {
    tp = 0
UNITSOFF
    printf("tp must not be negative\n")
    printf("tp has been changed to %g ms\n", tp)
UNITSON
  }
  dur = 2*trf + tp
  if (interval <= dur) {
    interval = dur + 1 (ms)
UNITSOFF
    printf("interval must be longer than trapezoid duration\n")
    printf("interval has been increased to %g ms\n", interval)
UNITSON
  }
  on = 0
  m = 0
  b = 0
  tally = number
  if (tally > 0) {
    net_send(start, 1)
    tally = tally - 1
  }
}

BREAKPOINT {
  g = gmax * (m*(t-t0) + b)
  i = (0.001)*g*(v-e)
}

NET_RECEIVE (w) {
  if ((on == 0) && (flag == 1)) {
    
    t0 = t
    m = 1/trf
    b = 0
    on = 1
    
    net_send(trf, 2)
  }
  if (flag == 2) {
    
    m = 0
    b = 1
    
    net_send(tp, 3)
  }
  if (flag == 3) {
    
    t0 = t
    m = -1/trf
    b = 1
    
    net_send(trf, 4)
  }
  if (flag == 4) {
    
    m = 0
    b = 0
    on = 0
    if (tally > 0) {
      
      net_send(interval - dur, 1)
      tally = tally - 1
    }
  }
}