UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Leakm
        NONSPECIFIC_CURRENT il
        RANGE  gl, el
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
       
        gl = .000075 (mho/cm2) 
        el = -75 (mV)
}
  
ASSIGNED {
	 v (mV)
        il (mA/cm2)
}
 
BREAKPOINT {

        il = gl*(v - el)
}