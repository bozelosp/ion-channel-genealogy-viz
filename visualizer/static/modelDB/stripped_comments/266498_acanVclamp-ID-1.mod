NEURON {
       SUFFIX acan
       USEION ca READ cao, cai WRITE ica
       RANGE gbar, ecan, ica

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
	  F = (faraday) (coulombs)
	  (molar) = (1/liter)
	  (mM)    = (millimolar)
}

PARAMETER {
    gbar =0.000795775 (S/cm2)
	Q10can=4.30 
	
	V0p5d=-20.0 (mV)
	S0p5d=4.5 (mV)
	
	V0p5f1=-20.0 (mV)
	S0p5f1=-25.0 (mV)
	
	V0p5f2=-40.0 (mV)
	S0p5f2=-10.0 (mV)
	
	A_taud=3.25	(ms)	
	B_taud=0.042	(/mV)
	C_taud=0.395	(ms)
	Vpd=-31.0		(mV)
	
	A_tauf1=33.5	(ms)
	B_tauf1=.0395	(/mV)
	C_tauf1=5.0	(ms)
	Vpf1=-30.0	(mV)
	
	A_tauf2=225.0	(ms)
	B_tauf2=0.0275	(/mV)
	C_tauf2=75.00	(ms)
	Vpf2=-40.0		(mV)
	
	A_rn=5.0 (mV)
	B_rn=-10.0 (mV)
	
	R=8.314 (joule/degC)
	z=2 
	ecaoffset=78.7 (mV)

}


ASSIGNED {
	 v	(mV) 
	 ica	(mA/cm2)
	 celsius (degC)
	 g	(S/cm2)
	 tau_f1	(ms)
	 tau_d	(ms)
	 tau_f2	(ms)
	 dinf
	 f1inf
	 f2inf
	 rn
     ecan	(mV)
	 cao (mM)
	 cai (mM)
         
}

STATE { d f1 f2 } 


BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * d * (0.55*f1+0.45*f2)
	   ica = g * (v-ecan)
}

INITIAL {
	rates(v) 
	
	

    d = dinf
	f1 = f1inf
	f2 = f2inf
}

DERIVATIVE states {
	   rates(v)
	   d' = (dinf - d)/tau_d
	   f1' = (f1inf - f1)/tau_f1
	   f2' = (f2inf - f2)/tau_f2
}



FUNCTION rates(Vm (mV)) (/ms) {

	rn=0.2/(1.0+exp((Vm +A_rn)/B_rn))
	
	tau_d = A_taud*exp(-(B_taud)^2*(Vm-Vpd)^2)+C_taud
        dinf = 1.0/(1.0+exp((Vm-V0p5d)/(-S0p5d)))

	tau_f1 = A_tauf1*exp(-(B_tauf1)^2*(Vm-Vpf1)^2)+C_tauf1
        f1inf = 1.0/(1.0+exp((Vm-V0p5f1)/(-S0p5f1)))

	tau_f2 = A_tauf2*exp(-(B_tauf2)^2*(Vm-Vpf2)^2)+C_tauf2
        f2inf =rn+(1.0/(1.0+exp((Vm-V0p5f2)/(-S0p5f2))))
		
	ecan=(1000)*(R*(celsius+273.15)/z/F*log(cao/cai))-ecaoffset 
		
	if (celsius >= 37) {
		tau_d=tau_d*(1/Q10can)
		tau_f1=tau_f1*(1/Q10can)
		tau_f2=tau_f2*(1/Q10can)
	}
}