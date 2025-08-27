NEURON {
	SUFFIX gbarcas
	POINTER F, S, D
	RANGE gbarcas, tau, Fbar, Sbar, Dbar, A, B, C






	RANGE gbarcas_init
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
	gbarcas_init (mA/cm2)
}

INITIAL {
	gbarcas = gbarcas_init
}

STATE { gbarcas (mA/cm2) }

BREAKPOINT {
	SOLVE state METHOD cnexp
}

DERIVATIVE state {
	gbarcas' = ( A*(Fbar-F) + B*(Sbar-S) + C*(Dbar-D) ) * gbarcas / tau
}