NEURON {
	POINT_PROCESS FDSExp2Syn
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE g
	GLOBAL total
        RANGE f, tau_F, d1, tau_D1, d2, tau_D2
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau1 = 0.1 (ms) < 1e-9, 1e9 >
	tau2 = 10 (ms) < 1e-9, 1e9 >
	e = 0	(mV)
        
	
        f = 0.917 (1) < 0, 1e9 >    
        tau_F = 94 (ms) < 1e-9, 1e9 >
        d1 = 0.416 (1) < 0, 1 >     
        tau_D1 = 380 (ms) < 1e-9, 1e9 >
        d2 = 0.975 (1) < 0, 1 >     
        tau_D2 = 9200 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	factor
	total (umho)
}

STATE {
	A (umho)
	B (umho)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (umho), F, D1, D2, tsyn (ms)) {
INITIAL {

        F = 1
        D1 = 1
        D2 = 1
        tsyn = t


}

        F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)
        D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)
        D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)

        tsyn = t

	state_discontinuity(A, A + weight*factor*F*D1*D2)
	state_discontinuity(B, B + weight*factor*F*D1*D2)
	total = total+weight*F*D1*D2

        F = F + f
        D1 = D1 * d1
        D2 = D2 * d2

}