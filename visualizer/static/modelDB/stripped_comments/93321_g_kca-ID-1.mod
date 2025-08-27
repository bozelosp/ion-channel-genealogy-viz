NEURON {
	SUFFIX gbarkca
	POINTER F, S, D
	RANGE gbarkca, tau, Fbar, Sbar, Dbar, A, B, C






	RANGE gbarkca_init
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
	B = -1
	C = -1
}

ASSIGNED {
	F (1)
	S (1)
	D (1)
	gbarkca_init (mA/cm2)
}

INITIAL {
	gbarkca = gbarkca_init
}

STATE { gbarkca (mA/cm2) }

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	gbarkca' = ( A*(Fbar-F) + B*(Sbar-S) + C*(Dbar-D) ) * gbarkca / tau
}