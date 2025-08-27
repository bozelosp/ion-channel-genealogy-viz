NEURON {
	SUFFIX xtrau_nl_const
	RANGE rx, er, d
	RANGE x, y, z
	RANGE Vx
	POINTER im, ex
}

PARAMETER {
	: default transfer resistance between recording electrode and node
	rx = 1 (megohm) : mV/nA
	Vx = 1 (millivolts)
	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)
}

ASSIGNED {
	v (millivolts)
	ex (millivolts)
	im (milliamp/cm2)
	er (microvolts)
	area (micron2)
}

INITIAL {
	ex = Vx
	er = (10)*rx*im*area
: this demonstrates that area is known
: UNITSOFF
: printf("area = %f\n", area)
: UNITSON
}

BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  ex = Vx
}
AFTER SOLVE { : after each solution step
  er = (10)*rx*im*area
}
