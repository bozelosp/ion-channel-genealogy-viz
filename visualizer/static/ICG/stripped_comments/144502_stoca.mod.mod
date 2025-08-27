NEURON {
       SUFFIX stoca
       USEION ca READ eca WRITE ica
       
       RANGE mMidV,gbar,g,i,minf,hinf,tauh,m,h 
}

UNITS {
      (S) = (siemens)
      (mS) = (millisiemens)
      (mV) = (millivolt)
      (mA) = (milliamp)
}

PARAMETER {
	  
	  gbar = 0.4 (mS/cm2)
      mMidV=-61 (mV) 
}

ASSIGNED {
         eca (mV)
	 v (mV)
	 ica (mA/cm2)
	 g (mS/cm2)
	 minf 
	 hinf 
	 tauh (ms)
}

STATE {
      m 
      h
}

INITIAL {
	rates(v)
	h = hinf
	m = minf
}

BREAKPOINT {
	   rates(v)
	   SOLVE states METHOD cnexp
	   g = gbar *minf*h
	   ica = g * (v - eca)*(0.001)
	   
}

DERIVATIVE states {
	   h' = (hinf -h)/tauh
}

PROCEDURE rates(v (mV)) {
	  
	  
	  
	  

	  
	  UNITSOFF
	  hinf =1/( 1+exp( (v+85.5)/8.6 ) ) 
	  tauh=40+30*(1/( 1+exp((v+84.0)/7.3) ))*exp((v+160.0)/30.0)
 	  minf = 1/(  (1+exp((mMidV-v)/4.2)) *(1+exp((mMidV-v)/4.2))* (1+exp((mMidV-v)/4.2)) )
	  m = minf
	  UNITSON
}