NEURON {
	SUFFIX prey
	POINTER bPointer
}

PARAMETER {}

ASSIGNED { bPointer }

STATE { a }

BREAKPOINT {
	SOLVE states METHOD derivimplicit
}

INITIAL { a = 10 }

DERIVATIVE states {
	a' = 1.1 * a - 0.4 * a * bPointer
}