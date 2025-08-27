NEURON {
  POINT_PROCESS MultIClamp
  ELECTRODE_CURRENT i
  RANGE trf, tp
  RANGE del, per, number
  RANGE amp, i

}

UNITS {
  (mV) = (millivolt)
  (nS) = (nanosiemens)
  (nA) = (nanoamp)
}

PARAMETER {
  trf (ms) <0, 1e9> 
  tp  (ms) <0, 1e9> 
  del (ms) <0, 1e9> 
  per (ms) <0, 1e9> 
  number 
  amp (nA) <0, 1e9> 
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
}

INITIAL {
  if (trf <= 0) {
    trf = 0.025 
UNITSOFF
    printf("time rise fall must be longer than 0
    printf("increased to trf = %g ms\n", trf)
UNITSON
  }
  if (tp < 0) {
    tp = 0
UNITSOFF
    printf("time plateau must not be negative
    printf("changed to tp = %g ms\n", tp)
UNITSON
  }
  dur = 2*trf + tp
  if (per <= dur) {
    per = dur 
UNITSOFF
    printf("period must be longer than trapezoid duration %g
    printf("increased to per = %g ms\n", per)
UNITSON
  }
  on = 0
  m = 0
  b = 0
  tally = number
  if (tally > 0) {
    net_send(del, 1)
    tally = tally - 1
  }
}

BREAKPOINT {
  i = amp * (m*(t-t0) + b)
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
      
      net_send(per - dur, 1)
      tally = tally - 1
    }
  }
}