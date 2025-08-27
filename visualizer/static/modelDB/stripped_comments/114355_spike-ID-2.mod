INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS Spike
	RANGE thresh
	RANGE spike
	RANGE spike_on
	RANGE spike_time
	RANGE spike_count
	RANGE spike_freq_isi	
	RANGE spike_freq_count	
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	thresh = 0 (mV)	
}

ASSIGNED {
	spike
	spike_on
	spike_time
	spike_freq_isi
	spike_freq_count
	spike_count
}

INITIAL {
	spike=0
	spike_time=0
	spike_freq_isi=0
	spike_on = 0
	spike_count=0
	spike_freq_count=0
	if (v > thresh) {
		spike_on = 1
	}
}

BREAKPOINT {
	SOLVE check
}

PROCEDURE check() {
	if (spike) {
		spike=0 
	}
	if (spike_on && v < thresh) {
		spike_on = 0
	}
	if (!spike_on && v > thresh) {
		spike=1
		spike_on = 1
		spike_count=spike_count+1
		if (t) {
			spike_freq_isi=1000/(t-spike_time)
			spike_freq_count=spike_count/t*1000
		}
		spike_time = t
	}
}