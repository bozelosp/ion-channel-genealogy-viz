NEURON {
	SUFFIX xtra
	RANGE rx, rx2, er
	RANGE x, y, z
	GLOBAL is, is2
	POINTER im, ex
}

PARAMETER {
	
	rx  = 1 (megohm) 
	rx2 = 1 (megohm) 
	x = 0 (1) 
	y = 0 (1)
	z = 0 (1)
}

ASSIGNED {
	v (millivolts)
	is (milliamp)
	is2 (milliamp)
	ex (millivolts)
	im (milliamp/cm2)
	er (microvolts)
	area (micron2)
}

INITIAL {
	ex = is*rx*(1e6) + is2*rx2*(1e6)
	er = (10)*rx*im*area




}
















BEFORE BREAKPOINT { 
  ex = is*rx*(1e6) + is2*rx2*(1e6)
}
AFTER SOLVE { 
  er = (10)*rx*im*area
}