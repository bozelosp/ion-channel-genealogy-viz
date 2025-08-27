INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS two_exp_syn1
	RANGE Ar, Af, As, onset, risetime, tslow, tfast, gmax, e, i, g
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset = 0 (ms)
        risetime = 0.5
	tfast = 7.5 (ms)
        tslow = 33
	gmax  = 0.01 (umho)
	e = 0	(mV)
	v	(mV)
        Af = 0.76
        Ar = 1.00
        As = 0.24
        cst
        T
}

INITIAL {
  T = -log(Af*risetime/(Ar*tfast))/(1/risetime - 1/tfast)
  cst = 1/(Af*exp(-T/tfast)-Ar*exp(-T/risetime)+(As*exp(-T/tslow)))
}

ASSIGNED { i (nA)  g (umho)}

BREAKPOINT {
        g = get_g()
	i = g*(v - e)
}

FUNCTION get_g() {
   if ((t>onset) && ((((t-onset)/risetime)<10) && (((t-onset)/tfast)<10) && (((t-onset)/tslow)<10))) {
     get_g = gmax*cst*(-Ar*exp(-(t-onset)/risetime) + Af*exp(-(t-onset)/tfast) + As*exp(-(t-onset)/tslow))
   }
   else
   {
   get_g = 0
   }
}