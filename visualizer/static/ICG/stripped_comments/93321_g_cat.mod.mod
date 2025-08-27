NEURON {
	SUFFIX gbarcat
	POINTER F, S, D
	RANGE gbarcat, tau, Fbar, Sbar, Dbar, A, B, C






	RANGE gbarcat_init
}

UNITS {
	(mA) = (milliamp)
}

PARAMETER {
	Fbar = 0.1 (1) 
	Sbar = 0.1 (1)
	Dbar = 0.1 (1)
	tau = 5000 (ms) 
	A = 0 
	B = 1
	C = 0
}

ASSIGNED {
	F (1)
	S (1)
	D (1)
	gbarcat_init (mA/cm2)
}

INITIAL {
	gbarcat = gbarcat_init
}

STATE { gbarcat (mA/cm2) }

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	gbarcat' = ( A*(Fbar-F) + B*(Sbar-S) + C*(Dbar-D) ) * gbarcat / tau
}