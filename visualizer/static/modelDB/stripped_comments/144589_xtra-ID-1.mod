NEURON {
	SUFFIX xtra
	RANGE rx, er, coeff, recx
	RANGE x, y, z
	GLOBAL is
	POINTER im, ex
}

PARAMETER {
	
      coeff=1
recx=1 
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
	er (microvolts)
	area (micron2)
}

INITIAL {
	ex = is*rx*coeff*(1e6)
	er = (10)*recx*im*area




}









BREAKPOINT { 
 

  


}

PROCEDURE f() {
	
	
	ex = is*rx*(1e6)
	er = (10)*recx*im*area
}

BEFORE BREAKPOINT { 
  ex = is*rx*coeff*(1e6)
}
AFTER SOLVE { 
  er = (10)*recx*im*area

}