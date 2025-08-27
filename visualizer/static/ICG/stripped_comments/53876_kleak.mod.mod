NEURON {
  SUFFIX kleak
  
  USEION k READ ek WRITE ik
  RANGE g
}

PARAMETER {
  g = 132.8e-6 (siemens/cm2) < 0, 1e9 >
  
}

ASSIGNED {
  ek (millivolt)
  ik (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { ik = g * (v - ek) }