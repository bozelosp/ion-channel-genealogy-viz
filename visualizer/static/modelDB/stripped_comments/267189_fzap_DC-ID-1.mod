NEURON {
  POINT_PROCESS Fzap_DC
  RANGE del, dur, amp
  POINTER x
}

PARAMETER {
  del (ms)
  dur (ms)
  amp (1)
}

ASSIGNED {
  x (1)
  on (1)
}

INITIAL {
  x = 0
  on = 0

  if (del<0) { del=0 }
  if (dur<0) { dur=0 }

  
  if (dur>0) {
    net_send(del, 1)  
    net_send(del+dur, 1)  
  }
}

BEFORE BREAKPOINT {
  if (on==0) {
    x = 0
  } else {
    x = amp
  }
}

NET_RECEIVE (w) {
  
  if (flag == 1) {
VERBATIM
	on = (double)(on == 0.0);
ENDVERBATIM
  }
}