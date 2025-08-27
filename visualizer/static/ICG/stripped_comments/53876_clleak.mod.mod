NEURON {
  SUFFIX clleak
  NONSPECIFIC_CURRENT i
  RANGE i, Erev, g
}

PARAMETER {
  g = 59.6e-6 (siemens/cm2) < 0, 1e9 >
  Erev = -70 (millivolt)
}

ASSIGNED {
  i (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { i = g * (v - Erev) }