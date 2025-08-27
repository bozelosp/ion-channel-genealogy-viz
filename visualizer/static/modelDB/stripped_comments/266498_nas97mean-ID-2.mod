NEURON {
		   SUFFIX nas97mean						
		   USEION na READ ena WRITE ina		
		   RANGE gbar, ena, ina, shiftnas	

	}


	UNITS {
		  (S) = (siemens)
		  (mV) = (millivolts)
		  (mA) = (milliamp)
	}


	PARAMETER {
		gbar =0.001043349 (S/cm2)	
		Q10nasm=2.30				
		Q10nash=1.50				
		Q10TempA = 22	(degC)		
		Q10TempB = 10	(degC)

		
		
		
			
				V0p5m=-11.29 (mV)
				S0p5m=5.54 (mV)
				
			
				A_taum=1.45	(ms)	
				B_taum=0.058	(/mV)
				C_taum=0.26	(ms)
				Vpm=-14.5		(mV)
			
		
			
			
				V0p5h=-31.00 (mV)
				S0p5h=-5.20 (mV)

			
				A_tauh=10.75	(ms)
				B_tauh=0.067	(/mV)
				C_tauh=3.15	(ms)
				Vph=-13.5	(mV)

	}


	ASSIGNED {
		
		
		 v	(mV) 
		 ina	(mA/cm2)
		 celsius (degC)
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
			
			tau_m=tau_m*Q10nasm^((Q10TempA-celsius)/Q10TempB)
			tau_h=tau_h*Q10nash^((Q10TempA-celsius)/Q10TempB)
		}