NEURON {
	POINT_PROCESS DAsyn
	RANGE tau, tau0, tau1, dur, msg, msginf, spkcnt, ip3i
}

PARAMETER {
	tau0 = 100 (ms)
	tau1 = 30 (ms)
	dur = 60
}

ASSIGNED {
	tau (ms)
	msginf (1)
	
spkcnt
	ip3i (mM)

}

STATE {
	msg (1)
}

INITIAL {
	tau = tau0
	msginf = 1
	msg = 1

	spkcnt=0
		ip3i =0
}

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	msg' = (msginf - msg)/tau
}

NET_RECEIVE (w(1)) {
  if (flag == 0) { 
	if (w>0) { 
		if (msginf==1) {
			msginf = 1 + w
			tau = tau1
			net_send(dur, 1)
		} else {
			net_move(t + dur)			
		}
	}
  } else { 
	msginf = 1
	tau = tau0
  }
  
  	spkcnt = spkcnt +1
}