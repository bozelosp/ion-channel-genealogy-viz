NEURON {
	SUFFIX gbarkd
	POINTER F, S, D
	RANGE gbarkd, tau, Fbar, Sbar, Dbar, A, B, C






	RANGE gbarkd_init
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
	B = -1
	C = 0
}

ASSIGNED {
	F (1)
	S (1)
	D (1)
	gbarkd_init (mA/cm2)
}

INITIAL {
	gbarkd = gbarkd_init
}

STATE { gbarkd (mA/cm2) }

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	gbarkd' = ( A*(Fbar-F) + B*(Sbar-S) + C*(Dbar-D) ) * gbarkd / tau
}