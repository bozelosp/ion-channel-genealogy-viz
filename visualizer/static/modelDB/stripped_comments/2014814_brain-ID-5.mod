NEURON {
  SUFFIX brain 
  POINTER bPointer 
}

ASSIGNED { bPointer } 

STATE { a } 

BREAKPOINT {
  SOLVE states METHOD derivimplicit
}

INITIAL {
  a = 1.0 
}

DERIVATIVE states {
  
  a' = da_dt()
}

FUNCTION da_dt() {
  LOCAL rhs 
  rhs = a * (1 - a) - bPointer
  if (a > 0 || rhs >= 0) {
    da_dt = rhs
  } else {
    da_dt = 0
  }
}