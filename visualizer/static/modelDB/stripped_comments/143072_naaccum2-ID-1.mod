NEURON {
	SUFFIX naaccum2
	
	USEION na READ ina, nao WRITE  nao
	RANGE fhspace, k
}

UNITS {
	(molar) = (1/liter)
	(mV) = (millivolt)
	(um) = (micron)
	(mM) = (millimolar)
	(mA) = (milliamp)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}

PARAMETER {
	nabath = 143	(mM)		
	diam		(um)
	ina		(mA/cm2)
	fhspace	= 20000	(angstrom)	
					
	
	k = .05				
	
	nao0 = 143	(mM)		
}

STATE {
	
	nao 	(mM)
}


INITIAL {
	
	
	nao = nao0
	
}

BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state {
	
	

	
	nao' = ina/fhspace/FARADAY*(1e8) + (nabath - nao)/k
	
}