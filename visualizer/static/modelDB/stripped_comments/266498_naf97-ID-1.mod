NEURON {
		   SUFFIX naf97						
		   USEION na READ ena WRITE ina		
		   RANGE gbar, ena, ina, shiftnaf	

	}


	UNITS {
		  (S) = (siemens)
		  (mV) = (millivolts)
		  (mA) = (milliamp)
	}


	PARAMETER {
		gbar =0.0012	(S/cm2)	
		Q10nafm=2.30				
		Q10nafh=1.50				
		Q10TempA = 22	(degC)		
		Q10TempB = 10	(degC)
		
		
		
		
			
			V0p5m=-37.75 (mV)	
			S0p5m=6.98 (mV)
			
			
			A_taum=1.2575	(ms)	
			B_taum=0.0625	(/mV)
			C_taum=0.175	(ms)
			Vpm=-39.0		(mV)
		
		
		
		
			
				V0p5h=-65.99 (mV)
				S0p5h=-5.97 (mV)
			
			
				A_tauh=25.5	(ms)
				B_tauh=0.035	(/mV)
				C_tauh=1.05	(ms)
				Vph=-72.5	(mV)
	}


	ASSIGNED {
		
		
		 v	(mV) 
		 celsius (degC)
		 ina	(mA/cm2)
		 ena	(mV)
		 
		 
		 g	(S/cm2)
		 tau_h	(ms)
		 tau_m	(ms)
		 minf
		 hinf
		 
			 
	}


	STATE { m h } 


	BREAKPOINT {
		   SOLVE states METHOD cnexp
		   g = gbar * m^3 * h
		   ina = g * (v-ena)
	}


	INITIAL {
		rates(v) 
		
		

		m = minf
		h = hinf
	}


	DERIVATIVE states {
		   rates(v)
		   m' = (minf - m)/tau_m
		   h' = (hinf - h)/tau_h
	}



	
	
		FUNCTION rates(Vm (mV)) (/ms) {
			 tau_m = A_taum*exp(-(B_taum)^2*(Vm-Vpm)^2)+C_taum
				 minf = 1.0/(1.0+exp((V0p5m-Vm)/S0p5m))

			 tau_h = A_tauh*exp(-(B_tauh)^2*(Vm-Vph)^2)+C_tauh
				 hinf = 1.0/(1.0+exp((V0p5h-Vm)/S0p5h))

			
				tau_m=tau_m*Q10nafm^((Q10TempA-celsius)/Q10TempB)
				tau_h=tau_h*Q10nafh^((Q10TempA-celsius)/Q10TempB)
		}