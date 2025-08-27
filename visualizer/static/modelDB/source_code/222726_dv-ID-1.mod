COMMENT
	deltav, deltavmax and vrest
ENDCOMMENT

NEURON {
	SUFFIX dv
	RANGE vrest, deltav, vmax, vmaxt, vmin, vmint
	RANGE dvmax, dvmaxt, dvmin, dvmint
}

ASSIGNED {
	v (millivolt)
	vrest (millivolt)
	deltav (millivolt)
	vmax (millivolt)
	vmaxt (ms)
	vmin (millivolt)
	vmint (ms)
	dvmax (millivolt)
	dvmaxt (ms)
	dvmin (millivolt)
	dvmint (ms)
}

INITIAL {
	vrest = v
	deltav = 0
	vmax = v
	vmaxt = 0
	vmin = v
	vmint = 0
	dvmax = 0
	dvmaxt = 0
	dvmin = 0
	dvmint = 0
}

BREAKPOINT {
if  (t<950) {
		vrest=v
	}
	deltav=v-vrest
	if (v>vmax) {
		vmax=v
		vmaxt=t
	}
	if (v<vmin) {
		vmin=v
		vmint=t
	}
	if (deltav>dvmax) {
		dvmax=deltav
		dvmaxt=t
	}
	if (deltav<dvmin) {
		dvmin=deltav
		dvmint=t
	}
}
