NEURON {
  SUFFIX naleak
  
  USEION na READ ena WRITE ina
  RANGE g
}

PARAMETER {
  g = 30.2e-6 (siemens/cm2) < 0, 1e9 >
  
}

ASSIGNED {
  ena (millivolt)
  ina (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { ina = g * (v - ena) }