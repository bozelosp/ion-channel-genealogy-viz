NEURON {
       SUFFIX akd
       USEION k READ ek WRITE ik
       RANGE gbar, ek, ik

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
    gbar =0.001061033 (S/cm2)
	Q10kdn=1.40

	
	V0p5n=-14.62 (mV)
	S0p5n=18.38 (mV)
	
	A_alphan=.001265 (/ms-mV) 
	B_alphan=14.273 (mV)
	C_alphan=-10.0	(mV)
	
	A_betan=0.125 (/ms)	
	B_betan=55.0 (mV)
	C_betan=-2.5 (mV)
}


ASSIGNED {
	 v	(mV) 
	 ik	(mA/cm2)
	 celsius (degC)
	 g	(S/cm2)
	 tau_n	(ms)
	 ninf
	 alphan (/ms)
	 betan (/ms)
     ek	(mV)
         
}

STATE { n }

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * n
	   ik = g * (v-ek)
}

INITIAL {
	rates(v) 
	
	
    n = ninf
}

DERIVATIVE states {
	   rates(v)
	   n' = (ninf - n)/tau_n
}

FUNCTION alpha(Vm (mV)) (/ms) {
	alphan=(A_alphan*(Vm+B_alphan))/(1.0-exp((Vm+B_alphan)/C_alphan))
}

FUNCTION beta(Vm (mV)) (/ms) {
	betan=A_betan*exp((Vm+B_betan)/C_betan)
}

FUNCTION rates(Vm (mV)) (/ms) {
	alpha(Vm)
	beta(Vm)
	
	tau_n = 1/(alphan+betan)+1.0
	ninf = 1.0/(1.0+exp((Vm-V0p5n)/(-S0p5n)))
	
	if (celsius >= 37) {
	tau_n=tau_n*(1/Q10kdn)
	}
}