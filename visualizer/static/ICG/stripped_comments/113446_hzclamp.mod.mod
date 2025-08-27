NEURON {
	POINT_PROCESS Hzclamp
	ELECTRODE_CURRENT i

	RANGE del, dur, amp, freq, width, i, telpulse, stand
}

UNITS {
	(nA) = (nanoamp)
	(uS) = (micromho)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (nA)
	dt (ms)
	freq	(1/s)
	width	(ms)
}

ASSIGNED {
	i (nA)
	i_amp (nA)
	telpulse
	stand
	notify
}

INITIAL {
	i_amp = 0
	stand = 0
	telpulse = 0
	notify=del
}

BREAKPOINT {
	SOLVE state METHOD after_cvode
	at_time(notify)
	i = i_amp
}

PROCEDURE state() {
	if (t <= del + dur && t >= del) {
	  if (telpulse/(freq/1000) < t-del && t-del <= telpulse/(freq/1000)+width ) {
	    notify = del + telpulse/(freq/1000)+width
	    i_amp=amp
	    stand=1 
	  } else if (stand == 1 ) {
	    stand = 0
	    telpulse = telpulse + 1
	    i_amp = 0
	    notify = del + telpulse/(freq/1000)
	  } else {
	    i_amp = 0
	  }
	}else{
	  i_amp = 0
	}
}