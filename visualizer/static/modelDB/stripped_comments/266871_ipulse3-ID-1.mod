NEURON {
  POINT_PROCESS Ipulse3
  RANGE dur, amp, i
  ELECTRODE_CURRENT i
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  dur (ms) <0, 1e9> 
  amp (nA) 
}

ASSIGNED {
  ival (nA)
  i (nA)
  on
}

INITIAL {
  on = 0
  i = 0
  ival = 0
}

BREAKPOINT {
  i = ival
}

NET_RECEIVE (w) {
  if (flag == 0) { 
    if (on == 0) {
      
      ival = amp
      on = 1
      
      net_send(dur, 1)
    }
  } else { 
    
    ival = 0
    on = 0
  }
}