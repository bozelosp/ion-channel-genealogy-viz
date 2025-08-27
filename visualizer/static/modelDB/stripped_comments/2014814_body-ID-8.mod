NEURON {
  SUFFIX body 
}

PARAMETER {
  b0 = 1.0  
  w = 0.628 
}

STATE { b }

BREAKPOINT {
  SOLVE states METHOD derivimplicit
}

INITIAL {
  b = 1.0 
}

DERIVATIVE states {
  b' = -b0 * w * sin(w * t)
}