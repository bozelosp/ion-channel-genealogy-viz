NEURON {
       SUFFIX I2
       NONSPECIFIC_CURRENT i
       RANGE i,E,g,m,h,Voffa,Vsloa,taua
}

UNITS { 
      (mV) = (millivolt)
      (S) = (siemens) 
}

STATE { m }


PARAMETER {
	  g = 20 (S/cm2)
	  E = -90 (mV)
          Voffa = -50 (mV)
          Vsloa = 8 (mV)
          taua = 0.35 (ms)
}

ASSIGNED {
	 i (mA/cm2)
	 v (mV)
}

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   i = g*(m^4*(v-E))
}

INITIAL {
	m = 1/(1+exp(-(v-Voffa)/Vsloa))
}

DERIVATIVE states {
	   m' = 1/taua*(1/(1+exp(-(v-Voffa)/Vsloa))-m)
}