UNITS {
   (mV) = (millivolt)
   (mA) = (milliamp)
}

INDEPENDENT { v FROM -100 TO 50 WITH 50   (mV) }

NEURON {
   SUFFIX PasSA
   NONSPECIFIC_CURRENT i
   RANGE g, erev
}

PARAMETER {
   g = .000677254 (mho/cm2)
   erev = -65     (mV)
}

ASSIGNED { i (mA/cm2) }

BREAKPOINT {
   i = g*(v - erev)
}