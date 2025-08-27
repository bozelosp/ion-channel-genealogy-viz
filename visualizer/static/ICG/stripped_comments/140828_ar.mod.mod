UNITS {
         (mV) = (millivolt)
         (mA) = (milliamp)
	   (S)  = (siemens)
}

INDEPENDENT {v FROM -100 TO 50 WITH 50 (mV)}

NEURON {
         SUFFIX ar
         NONSPECIFIC_CURRENT i
         RANGE g, c
	 GLOBAL eh
}

PARAMETER {
         g0 = .0001      (S/cm2)       <0,1e9>
         
         c = 1000000 	 (cm4 ohm2/mV)
}
  
ASSIGNED {
	i (mA/cm2)
	g (S/cm2)
	eh (mV)
}

BREAKPOINT {
	if (c==0) {
		g = g0
		i = g0*(v-eh)
	} else {
		g = 1/sqrt(1/g0^2+4*c*(v-eh)) 
		i = ( -1/g0 + 1/g)/(2*c)
	}
}