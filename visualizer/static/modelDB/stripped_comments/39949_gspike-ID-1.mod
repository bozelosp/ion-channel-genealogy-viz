NEURON {
	POINT_PROCESS GSpike
	RANGE vt, vh, vc, spike, minisi
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	vt = -58 (mV)
	vh = -55 (mV)
	vc = 2.5 (mV)
	minisi = 20 (ms)
	spikewidth = 1 (ms)
}

ASSIGNED {
	v (mV)
	tp (ms)
	spike (1)
}

INITIAL {
	tp = -1e9 (ms)
	spike = 0
}

AFTER SOLVE {
	if (at_time(tp+spikewidth)) {
		spike = 0
	}
	if (v>vt) {
		if (t>tp+minisi*(1+exp(-(v-vh)/vc))) {
			spike = 1
			tp = t
		}
	}
}