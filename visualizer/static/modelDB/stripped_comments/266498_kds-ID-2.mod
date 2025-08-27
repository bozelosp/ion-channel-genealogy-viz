NEURON {
		   SUFFIX kds						
		   USEION k READ ek WRITE ik		
		   RANGE gbar, ek, ik, shiftkds		

	}


	UNITS {
		  (S) = (siemens)
		  (mV) = (millivolts)
		  (mA) = (milliamp)
	}


PARAMETER {
    gbar =0.000106103  (S/cm2)		
	Q10kds=1.93 					
	Q10TempA = 22.85	(degC)		
	Q10TempB = 10	(degC)

 
	shiftkds=3.0 (mV) 				
	
	
	
		
			V0p5x=-39.59 (mV)		
			S0p5x=14.68(mV)
		
		
			A_taux=5.0	(ms)	
			B_taux=0.022	(/mV)
			C_taux=2.5	(ms)
			Vpx=-65.0		(mV)
	
	
	
		
			V0p5y=-48.0 (mV)
			S0p5y=-7.0 (mV)
		
		
			tau_y22=7500 (ms) 

}

	ASSIGNED {
	
		
		 v	(mV) 
		 ik	(mA/cm2)
		 celsius (degC)
		 ek	(mV)
		 
		 
		 g	(S/cm2)
		 tau_x	(ms)
		 tau_y	(ms)
		 xinf
		 yinf
		 
			 
	}


	STATE { x y1 } 


	BREAKPOINT {
		   SOLVE states METHOD cnexp
		   g = gbar * x^3 * y1
		   ik = g * (v-ek)
	}


	INITIAL {
		rates(v) 
		
		

		x = xinf
		y1 = yinf
	}


	DERIVATIVE states {
		   rates(v)
		   x' = (xinf - x)/tau_x
		   y1' = (yinf - y1)/tau_y
	}



	
		FUNCTION rates(Vm (mV)) (/ms) {
			tau_x = A_taux*exp(-(B_taux)^2*(Vm-Vpx)^2)+C_taux
				xinf = 1.0/(1.0+exp((Vm-V0p5x+shiftkds)/(-S0p5x)))
				 
			tau_y=tau_y22
				yinf = 1.0/(1.0+exp((Vm-V0p5y+shiftkds)/(-S0p5y)))
			
			
			tau_x=tau_x*Q10kds^((Q10TempA-celsius)/Q10TempB)
			tau_y=tau_y*Q10kds^((Q10TempA-celsius)/Q10TempB)
		}