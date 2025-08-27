NEURON {

  POINT_PROCESS FluctExp2Syn
  RANGE tau_rise, tau_fall, cn, mean_amp, cv, type, e, i, s, g, std, std0, pf, plas, tau_plas
  RANGE seed, flag_print
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (uS) = (microsiemens)
}

PARAMETER {
  tau_rise = 1.0 (ms) <1e-9,1e9>
  tau_fall = 2.0 (ms) <1e-9,1e9>
  cn = 4         
  mean_amp = 0.001
  cv = 0.0
  type = 1
  e=0	(mV)
  std0 = 1.0
  pf = 0.0
  plas = 0.75
  tau_plas = 120.0
  seed = 12
  flag_print = 0
}

ASSIGNED {
  v (mV)
  i (nA)
}

STATE {
  s (uS)
  g (uS)
  std
}

INITIAL {
  s = 0
  g = 0
  std = std0
  set_seed(seed)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  
  if (type == -1) {  
    i = -g               
  } else {
    i = g * (v-e) / 70    
    }
}

DERIVATIVE state {
  s' = -s/tau_rise
  g' = (cn*s-g)/tau_fall
  std' = -(std-1.0)/tau_plas
}

NET_RECEIVE(w (uS)) {
  LOCAL ww, pfail, prob

  
  
  ww = fabs(w*(1+cv*normrand(0,1)))
  
  if (ww > 3*fabs(w)){
    ww = 3*fabs(w)
  }
  
  if (w<0) {
    ww = -ww
  }

  
  ww = std * ww

  
  pfail = pf / std
  
  
  prob = scop_random()
  if (prob >= pfail) {

    s = s + ww         

    std = std * plas   
    
    if (std >= 5) {
      std = 5
    }
    if (std <= 0.4) {
      std = 0.4
    }
  } 
}