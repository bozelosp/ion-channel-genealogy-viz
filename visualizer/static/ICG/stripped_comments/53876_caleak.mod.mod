NEURON {
  SUFFIX caleak
  
  USEION ca READ eca WRITE ica
  RANGE  g
}

PARAMETER {
  g = 1.33e-6 (siemens/cm2) < 0, 1e9 >
  
}

ASSIGNED {
  eca (mV)
  
  ica (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { 
  ica = g * (v - eca)
  
}