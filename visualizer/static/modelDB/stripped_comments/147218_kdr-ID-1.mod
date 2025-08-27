NEURON {
     SUFFIX kdr
     USEION k READ ek WRITE ik
     RANGE gk, ik
}

UNITS {
     (S)  = (siemens)
     (mV) = (millivolt)
     (mA) = (milliamp)
}

PARAMETER { gk = 42e-4 (S/cm2) }

ASSIGNED {
     v       (mV)
     ek      (mV)
     ik      (mA/cm2)
}


STATE { m }

BREAKPOINT {
     SOLVE states METHOD cnexp
     ik = gk * m^4 * (v - ek)
}


INITIAL {
     
     m = minf(v)
}


DERIVATIVE states {
     
     m' = (minf(v)-m) / taum(v)

}


FUNCTION minf (Vm (mV)) {
     UNITSOFF
     minf = 1 / (1 + exp(-(Vm+25)/11.5)) 
     UNITSON
}

FUNCTION taum (Vm (mV)) (ms) {
     LOCAL x
     UNITSOFF
     x = exp((Vm+22.5)/17)
     taum = 0.2 + 4.15 / (x + 0.6/x)
     UNITSON
}