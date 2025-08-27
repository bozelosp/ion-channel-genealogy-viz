NEURON {
	SUFFIX VDpas

	NONSPECIFIC_CURRENT i
	
        RANGE i, e, e50, g, gmin, gmax, s, f 
	
}

UNITS {
   (S) = (siemens)
   (mV) = (millivolt)
   (mA) = (milliamp)
}

PARAMETER {
      
	gmax = 0.002800 (S/cm2)
        gmin = 0.000660 (S/cm2)
        e50 = -31 (mV)
        s = -6
        f = 0.051 
        e = -60  
}

ASSIGNED {
	v	(mV)		
	i	(nA)		
	g 	(uS)		
}


BREAKPOINT {
          
        g = f*(gmax/(1+exp((v-e50)/s))+gmin)
        
	i = g * (v-e)
}