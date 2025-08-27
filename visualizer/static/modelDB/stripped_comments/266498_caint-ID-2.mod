NEURON {
	  SUFFIX caint						
	  USEION ca READ ica WRITE cai		
	  GLOBAL nb							
	  RANGE a, ku, kr, Bi, diffOc		
	}


	UNITS {
	  (mV)    = (millivolt)
	  (mA)    = (milliamp)
	  FARADAY = 96500 (coulombs)
	  (molar) = (1/liter)
	  (mM)    = (millimolar)
	}


	PARAMETER {
	  a = 15e-4 (cm)  			
	  ku = 100 (/mM/ms) 		
	  kr = 0.238 (/ms)		 	
	  nb = 4 					
	  Bi = 0.001 (mM) 			
	  SA = 2.82743E-05 (cm2) 	
	  Vol = 1.27e-8 (cm3) 		
	}


	ASSIGNED { 
		
		
		ica  (mA/cm2)
		
		
		diffOc (/ms)


	}


	STATE { cai  (mM) Oc } 


	BREAKPOINT { SOLVE state METHOD derivimplicit }


	INITIAL {

	Oc=.05

	}


	DERIVATIVE state {
	  LOCAL diffOc											
	  diffOc=ku*cai*(1-Oc)-kr*Oc 							
	  Oc'=diffOc											
	  cai' = -ica*(SA)/Vol/(2*FARADAY) + -nb*Bi*diffOc		
	}