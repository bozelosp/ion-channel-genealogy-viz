TITLE all current pointer


NEURON {
	POINT_PROCESS FlagALL
	RANGE f,a,del
	
	
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	
}

PARAMETER {
	a
    del (ms)
}


ASSIGNED {
	f				: flag
		
	}


INITIAL {
	f=0
}

BREAKPOINT {

    at_time(del)
	
	if (t >= del) {
		f = a
	}else{
		f = 0
	}

}

