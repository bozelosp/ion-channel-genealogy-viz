NEURON {
	SUFFIX xtrau_nl
	RANGE rx, er, d
	RANGE x, y, z
	RANGE Vx
	GLOBAL E
	POINTER im, ex
}

PARAMETER {
	
	rx = 1 (megohm) 
	Vx = 1 (millivolts)
	x = 0 (1) 
	y = 0 (1)
	z = 0 (1)
}

ASSIGNED {
	v (millivolts)
	E (1) 
	ex (millivolts)
	im (milliamp/cm2)
	er (microvolts)
	area (micron2)
}

INITIAL {
	ex = E*Vx
	er = (10)*rx*im*area




}

BEFORE BREAKPOINT { 
  ex = E*Vx
}
AFTER SOLVE { 
  er = (10)*rx*im*area
}