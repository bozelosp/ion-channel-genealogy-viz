NEURON {
	SUFFIX brain
	POINTER bPointer 
}

ASSIGNED { bPointer } 

STATE { a } 

BREAKPOINT {
	SOLVE states METHOD derivimplicit 
	if ((a < 0) && (a * (1 - a) - bPointer < 0)) {
		a = 0
	}
}

INITIAL {
	a = 1.0 
}

DERIVATIVE states {
	a' = a * (1 - a) - bPointer
}