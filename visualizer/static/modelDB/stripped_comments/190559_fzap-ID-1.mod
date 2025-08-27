NEURON {
  POINT_PROCESS Fzap
  RANGE del, dur, f0, f1, amp, f
  POINTER x
}

UNITS {
  PI = (pi) (1)
}

PARAMETER {
  del (ms)
  dur (ms)
  f0 (1/s)  
  f1 (1/s)
  amp (1)
}

ASSIGNED {
  f (1/s)
  x (1)
  on (1)
}

INITIAL {
  f = 0
  x = 0
  on = 0

  if (del<0) { del=0 }
  if (dur<0) { dur=0 }
  if (f0<=0) { f0=0 (1/s) }
  if (f1<=0) { f1=0 (1/s) }

  
  if (dur>0) {
    net_send(del, 1)  
    net_send(del+dur, 1)  
  }
}



BEFORE BREAKPOINT {
  if (on==0) {
    f = 0
    x = 0
  } else {
    f = f0 + (f1 - f0)*(t-del)/dur
    x = amp * sin( 2*PI * (t-del) * (f0 + (f1 - f0)*(t-del)/(2*dur)) * (0.001) )
  }
}

NET_RECEIVE (w) {
  
  if (flag == 1) {
VERBATIM
	on = (double)(on == 0.0);
ENDVERBATIM
  }
}