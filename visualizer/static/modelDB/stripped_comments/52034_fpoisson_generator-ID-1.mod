INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS fpoisson_generator
  RANGE mean_rate, magnitude, on, off, bg_rate
  POINTER out_stim
}

ASSIGNED {
  rn
  prob1
  prob2
  dt
  out_stim
}

PARAMETER {
  mean_rate (kHz)		
  magnitude (uS)		
  on = 0
  off = 1e11
  bg_rate = 0
}

INITIAL {
  prob1 = bg_rate * dt * 2^31
  prob2 = mean_rate * dt * 2^31
  
}

BREAKPOINT {
SOLVE dum
}

PROCEDURE dum() {

  out_stim = 0
  if (t >= on && t < off) {
    if (scop_random() < prob2)  {out_stim = magnitude}
  } else {
    if (scop_random() < prob1)  {out_stim = magnitude}
  }  
}