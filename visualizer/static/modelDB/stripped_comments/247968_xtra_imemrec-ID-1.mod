NEURON {
	SUFFIX xtraimemrec
	RANGE rx, er
	RANGE x, y, z,LFPtemp
	GLOBAL is,LFP
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
	er (microvolts)
	area (micron2)
	LFP (millivolts)
   LFPtemp(millivolts)
}

INITIAL {
	ex = is*rx*(1e6)
	er = (10)*0.1*im*area
	LFP = LFP + er
	LFPtemp=LFP




}
















BEFORE BREAKPOINT { 
  ex = is*rx*(1e6)
  LFP = 0
  LFPtemp=LFP
  
}
AFTER SOLVE { 
  er = (10)*0.1*im*area
  LFP = LFP + er
  LFPtemp=LFP
  
}