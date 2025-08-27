NEURON {
	SUFFIX xtra
	RANGE rx
	RANGE x, y, z
	GLOBAL is
	POINTER im, ex
}

PARAMETER {
	
	rx = 1 (megohm) 
	x = 0 (1) 
	y = 0 (1)
	z = 0 (1)
}

ASSIGNED {
	v (millivolts)
	is (milliamp)
	ex (millivolts)
	im (milliamp/cm2)
	area (micron2)
}

INITIAL {
	ex = is*rx*(1e6)




}

BEFORE BREAKPOINT { 
  ex = is*rx*(1e6)
}