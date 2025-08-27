NEURON {
  SUFFIX caleak
  
  USEION ca WRITE  ica
  RANGE i, Erev, g
}

PARAMETER {
  g = 1.33e-6 (siemens/cm2) < 0, 1e9 >
  Erev = 80 (millivolt)
}

ASSIGNED {
  i (milliamp/cm2)
  ica (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { 
  ica = g * (v - Erev)
  i = ica	
}