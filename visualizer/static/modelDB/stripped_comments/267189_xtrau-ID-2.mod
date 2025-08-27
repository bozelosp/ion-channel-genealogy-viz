NEURON {
	SUFFIX xtrau
	RANGE rx, er, d
	RANGE x, y, z
	GLOBAL E
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

	E (volt/m) 
	d (micron) 
	ex (millivolts)
	im (milliamp/cm2)
	er (microvolts)
	area (micron2)
}

INITIAL {

	ex = -E*d*(1e-3)
	er = (10)*rx*im*area




}
















BEFORE BREAKPOINT { 

  ex = -E*d*(1e-3)
}
AFTER SOLVE { 
  er = (10)*rx*im*area
}