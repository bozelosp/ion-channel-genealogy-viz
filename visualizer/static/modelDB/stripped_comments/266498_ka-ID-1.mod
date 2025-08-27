NEURON {
       SUFFIX ka						
       USEION k READ ek WRITE ik		
       RANGE gbar, ek, ik, shiftka		

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}


	PARAMETER {
		gbar =0.000141471 (S/cm2)	
		Q10ka=1.93 					
		Q10TempA = 22.85	(degC)			
		Q10TempB = 10	(degC)

	 
		shiftka=3.0 (mV) 			
		
		
			
			
				V0p5p=-28.0 (mV)
				S0p5p=28.0 (mV)
			
			
				A_taup=5.0	(ms)	
				B_taup=0.022	(/mV)
				C_taup=2.5	(ms)
				Vpp=-65.0		(mV)
		
		
		
			
				V0p5q=-58.0 (mV)
				S0p5q=-7.0 (mV)
			
			
				A_tauq=100.0	(ms)
				B_tauq=0.035	(/mV)
				C_tauq=10.5	(ms)
				Vpq=-30.0	(mV)

	}


	ASSIGNED {

		
		 v	(mV) 
		 ik	(mA/cm2)
		 celsius  (degC)
		 ek	(mV)
		 
		 
		 g	(S/cm2)
		 tau_q	(ms)
		 tau_p	(ms)
		 pinf
		 qinf
		 
			 
	}


	STATE { p q } 


	BREAKPOINT {
		   SOLVE states METHOD cnexp
		   g = gbar * p^3 * q
		   ik = g * (v-ek)
	}


	INITIAL {
		rates(v) 
		
		

		p = pinf
		q = qinf
	}


	DERIVATIVE states {
		   rates(v)
		   p' = (pinf - p)/tau_p
		   q' = (qinf - q)/tau_q
	}



	
		FUNCTION rates(Vm (mV)) (/ms) {
			 tau_p = A_taup*exp(-(B_taup)^2*(Vm-Vpp)^2)+C_taup
				 pinf = 1.0/(1.0+exp((Vm-V0p5p+shiftka)/(-S0p5p)))

			 tau_q = A_tauq*exp(-(B_tauq)^2*(Vm-Vpq)^2)+C_tauq
				 qinf = 1.0/(1.0+exp((Vm-V0p5q+shiftka)/(-S0p5q)))
			
			
			tau_p=tau_p*Q10ka^((Q10TempA-celsius)/Q10TempB)
			tau_q=tau_q*Q10ka^((Q10TempA-celsius)/Q10TempB)
		}