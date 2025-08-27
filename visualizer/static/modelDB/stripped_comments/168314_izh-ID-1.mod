NEURON {
  POINT_PROCESS IZH
  NONSPECIFIC_CURRENT vv
  RANGE a,b,c,d,e,f,I,vv,thresh, vr, vt, vpeak, aACH, cACH, dACH,alphaShutdown, bACH, ACHshutdown, aMin, aMax, g, vrACH, k, Cap, uinit
}
UNITS {
	(mV) = (millivolt)
    (nA) = (nanoamp)
	(pA) = (picoamp)
	(uS) = (microsiemens)
	(nS) = (nanosiemens)
}

INITIAL {
u = uinit
net_send(0,1)
}

PARAMETER {  

  k = 0.0011 (nA/mV2) 
  a = 0.01 (1/ms)
  b = 0.0002 (uS)
  c = -65 (mV)
  d = .001 (nA)

 
  vpeak= 30 (mV)
  vv = 0 (mV)
  vr = - 70 (mV)
  vt = - 45 (mV)

  a_OLM = 0.002
  ACH = 1 	 		
  dACH = 0 (nA)		
  cACH = 0 (mV)		
  vrACH = 0 (mV)	
  bACH = 1.25
  
  ACHshutdown = 0	
  
  uinit = 0 (nA)
  
}

STATE { u }

ASSIGNED {
  }

BREAKPOINT {
  SOLVE states METHOD derivimplicit
  vv = -(k*(v - (vr + vrACH * (-1+ ACH) ))*(v - vt) - u)
  
  a_OLM = 0.023*ACH^2 - 0.022*ACH + 0.002 
}

DERIVATIVE states {
    u' = (a_OLM * ACHshutdown + a )*(b*(v - (vr + vrACH * (-1+ ACH) ))-u)
}

NET_RECEIVE (w) {
  if (flag == 1) {
    WATCH (v>vpeak) 2
  } else if (flag == 2) {
    net_event(t)
    v = c + cACH * (-1+ACH) 
    u = u+d + dACH * (1-ACH)
  }
}