NEURON {

	ARTIFICIAL_CELL FDSIntFire
	RANGE tau, y, spikedur, refrac, f, tau_F, d1, tau_D1, d2, tau_D2
	
}

PARAMETER {
	tau = 10 (ms) < 1e-9, 1e9 >
	spikedur = 1 (ms) < 0, 1e9 >
	refrac = 5 (ms) < 0, 1e9 >
	
	
	f = 0.917 (1) < 0, 1e9 >	
	tau_F = 94 (ms) < 1e-9, 1e9 >
	d1 = 0.416 (1) < 0, 1 >	
	tau_D1 = 380 (ms) < 1e-9, 1e9 >
	d2 = 0.975 (1) < 0, 1 >	
	tau_D2 = 9200 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	y
	t0 (ms)
	refractory
}

INITIAL {
	y = 0
	t0 = t
	refractory = 0 
}

FUNCTION yy() {
	yy = y*exp(-(t - t0)/tau)
}

FUNCTION Fval(F, time) {
	Fval = 1 + (F-1)*exp(-time/tau_F)
}

FUNCTION D1val(D1, time) {
	D1val = 1 - (1-D1)*exp(-time/tau_D1)
}

FUNCTION D2val(D2, time) {
	D2val = 1 - (1-D2)*exp(-time/tau_D2)
}

NET_RECEIVE (w, F, D1, D2, tsyn (ms)) {
INITIAL {

	F = 1
	D1 = 1
	D2 = 1
	tsyn = t


}

	if (flag == 0) {
		
		



		F = Fval(F, t-tsyn)
		D1 = D1val(D1, t-tsyn)
		D2 = D2val(D2, t-tsyn)

		tsyn = t
	}
	if (refractory == 0) { 

		y = yy()
		y = y + w*F*D1*D2
		t0 = t
		if (y > 1) {
			refractory = 1
			y = 2
			net_send(spikedur, refractory)
			net_event(t)
		}
	} else if (flag == 1) { 
		
		refractory = 2
		y = -1
		net_send(refrac, refractory)
	} else if (flag == 2) { 
		refractory = 0
		y = 0
	}
	if (flag == 0) {
		
		F = F + f
		D1 = D1 * d1
		D2 = D2 * d2

	}
}