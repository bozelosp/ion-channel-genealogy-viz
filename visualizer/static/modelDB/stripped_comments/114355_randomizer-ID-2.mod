NEURON {
	POINT_PROCESS Randomizer
	RANGE Eintrvl			
	RANGE NEproc			
	POINTER proc_num		
	RANGE special
}

PARAMETER {
	Eintrvl
	NEproc
	special = 0	
}

ASSIGNED {
	proc_num
	last_event
	tshift
	skip
	count
}

INITIAL {
	LOCAL j
	skip = 0
	count = 0
	if (special) {
		Eintrvl = Eintrvl / 2.5
	}
	last_event = 0
	tshift = 0
	j = 0
	while (j < NEproc) {
		net_send(exprand(Eintrvl),j)
		j = j + 1
	}
}

NET_RECEIVE(w) {
	if (t == last_event) {
		tshift = tshift + dt/10
		net_send(tshift,flag)

	} else {
		proc_num = flag
		net_send(exprand(Eintrvl),flag)
		
		
		if (special) {
			count = count + 1
			if (count == 2 || count == 4) {
				skip = 0
			} else {
				skip = 1
			}
			if (count == 5) {
				count = 0
			}
		}
		if (!skip) {
			net_event(t)
		}
		tshift = 0
	}
	last_event = t
}

PROCEDURE seed(x) {
	set_seed(x)
}