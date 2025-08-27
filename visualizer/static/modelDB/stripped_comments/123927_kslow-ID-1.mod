NEURON {
  SUFFIX kslow
  USEION k READ ek WRITE ik
  RANGE gbar, g, i
}

UNITS {
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)	
}

PARAMETER {
  gbar = 6.59090909e-5 (S/cm2) 
  eK = -95 (mV)	
}

ASSIGNED {
  v	(mV)
  ek	(mV)
  ik 	(mA/cm2)
  i 	(mA/cm2)
  g	(S/cm2)
  
}

STATE {m h}

BREAKPOINT {
  SOLVE states METHOD cnexp
  g = gbar*h*m*m 
  i = g*(v-eK)
  ik = i
}

INITIAL {
  m = minf(v)
  h = hinf(v)
}

DERIVATIVE states {
 m'= (minf(v)-m)/mtau(v)
 h'= (hinf(v)-h)/htau(v)
}

FUNCTION minf (Vm (mV)) () {

  UNITSOFF
    minf = 1/(1+exp(-(Vm+14.3)/14.6)) 

}

FUNCTION mtau (Vm (mV)) (ms) {

  UNITSOFF
	
    if (Vm < 50) {
    	mtau = 1.25+1.15*exp(-0.026*Vm) 
	} else {
    	mtau = 1.25+13*exp(-0.026*Vm) 
	}
  UNITSON

}


FUNCTION hinf (Vm (mV)) () {

  UNITSOFF
    hinf = 1/(1+exp((Vm+54)/11)) 
  UNITSON

}

FUNCTION htau (Vm (mV)) (ms) {

  UNITSOFF
    htau = 360+(1010+24*(Vm+55))*exp(-((Vm+75)/48)*((Vm+75)/48)) 
  UNITSON

}