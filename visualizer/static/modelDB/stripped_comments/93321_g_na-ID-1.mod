NEURON {
	SUFFIX gbarna
	POINTER F, S, D
	RANGE gbarna, tau, Fbar, Sbar, Dbar, A, B, C






	RANGE gbarna_init
}

UNITS {
	(mA) = (milliamp)
}

PARAMETER {
	Fbar = 0.1 (1) 
	Sbar = 0.1 (1)
	Dbar = 0.1 (1)
	tau = 5000 (ms) 
	A = 1 
	B = 0
	C = 0
}

ASSIGNED {
	F (1)
	S (1)
	D (1)
	gbarna_init (mA/cm2)
}

INITIAL {
	gbarna = gbarna_init
}

STATE { gbarna (mA/cm2) }

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	gbarna' = ( A*(Fbar-F) + B*(Sbar-S) + C*(Dbar-D) ) * gbarna / tau
}