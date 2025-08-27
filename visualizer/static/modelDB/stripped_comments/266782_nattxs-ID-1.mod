NEURON {
       SUFFIX nattxs
       USEION na READ ena WRITE ina
       RANGE gbar, ena, ina, celsiusT, Tshift, tau_m, tau_h, tau_s

}

UNITS {
      (S) = (siemens)
      (mV) = (millivolts)
      (mA) = (milliamp)
}

PARAMETER {
	  gbar = 0 (S/cm2)
          enainit (mV)
          kvot_qt
          celsiusT


  A_am = 15.5 (/ms)  
  B_am = -5 (mV)    
  C_am = -12.08 (mV)   

  A_ah = 0.38685 (/ms) 
  B_ah = 122.35 (mV)     
  C_ah = 15.29 (mV)   

  A_as = 0.00092 (/ms) 
  B_as = 93.9 (mV)     
  C_as = 16.6 (mV)   

  A_bm = 35.2 (/ms)   
  B_bm = 72.7 (mV)    
  C_bm = 16.7 (mV)    

  A_bh = 2.00283 (/ms)    
  B_bh = 5.5266 (mV)    
  C_bh = -12.70195 (mV) 

  A_bs = -132.05 (/ms)    
  B_bs = -384.9 (mV)    
  C_bs = 28.5 (mV) 

  shift=0 (mV) 
  Tshift=0 (mV)

}

ASSIGNED {
	 v	(mV) 
	 ina	(mA/cm2)
	 g	(S/cm2)
	 tau_h	(ms)
	 tau_m	(ms)
	 tau_s	(ms)
	 minf
	 hinf
	 sinf
         ena	(mV)
         
}

STATE { m h s }

BREAKPOINT {
	   SOLVE states METHOD cnexp
	   g = gbar * m^3 * h *s
	   ina = g * (v-ena)
}

INITIAL {
	rates(v) 
	
	

  m = minf
	h = hinf
	s = sinf
}

DERIVATIVE states {
	   rates(v)
	   m' = (minf - m)/tau_m
	   h' = (hinf - h)/tau_h
	   s' = (sinf - s)/tau_s
}

FUNCTION alpham(Vm (mV)) (/ms) {
	 alpham=A_am/(1+exp((Vm+shift+B_am)/C_am))
}

FUNCTION alphah(Vm (mV)) (/ms) {
	 alphah=A_ah/(1+exp((Vm+shift+B_ah)/C_ah))
}

FUNCTION alphas(Vm (mV)) (/ms) {
	 alphas=0.00003+A_as/(1+exp((Vm+shift+B_as+Tshift)/C_as))
}

FUNCTION betam(Vm (mV)) (/ms) {
	 betam=A_bm/(1+exp((Vm+shift+B_bm)/C_bm))
}

FUNCTION betah(Vm (mV)) (/ms) {
	 betah=-0.00283+A_bh/(1+exp((Vm+shift+B_bh)/C_bh))
}

FUNCTION betas(Vm (mV)) (/ms) {
	 betas=132.05+A_bs/(1+exp((Vm+shift+B_bs+Tshift)/C_bs))
}

FUNCTION rates(Vm (mV)) (/ms) {
	 tau_m = 1.0 / (alpham(Vm) + betam(Vm))
         minf = alpham(Vm) * tau_m

	 tau_h = 1.0 / (alphah(Vm) + betah(Vm))
         hinf = alphah(Vm) * tau_h

	 tau_s = 1.0 / (alphas(Vm) + betas(Vm))
         sinf = alphas(Vm) * tau_s

         kvot_qt=1/((2.5^((celsiusT-21)/10)))
         tau_m=tau_m*kvot_qt
         tau_h=tau_h*kvot_qt
         tau_s=tau_s*kvot_qt

}