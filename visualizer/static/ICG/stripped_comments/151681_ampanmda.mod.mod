NEURON {
	POINT_PROCESS AmpaNmda
	RANGE R, g, mg, inmda, iampa, gnmda, gampa
	RANGE x, mgid, ggid, srcgid, gmax
	NONSPECIFIC_CURRENT i
	GLOBAL Cdur, Alpha, Beta, E, Rinf, Rtau, ampatau
	GLOBAL gampafactor, nmdafactor
	GLOBAL ltdinvl, ltpinvl, sighalf, sigslope
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cdur	= 1		(ms)	
	Alpha	= 0.35		(/ms)	
	Beta	= 0.035		(/ms)	
	E	= 0	(mV)		
	mg	= 1    (mM)		
	gmax = 2 (umho)		
	gampafactor = 0.001 (1)
	nmdafactor = 0.0035 (1)
	ltdinvl = 250 (ms)		
	ltpinvl = 33.33 (ms)		
	sighalf = 50 (1)
	sigslope = 10 (1)
	ampatau = 3 (ms)
	x = 0 (um) 
	mgid = -1 
	ggid = -1 
	srcgid = -1 
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	inmda 		(nA)		
	iampa 		(nA)		
	gnmda 		(umho)		
	Rinf				
	Rtau		(ms)		
	synon
}

STATE {Ron Roff
	gampa 		(umho)
}

INITIAL {
	PROTECT Rinf = Alpha / (Alpha + Beta)
	PROTECT Rtau = 1 / (Alpha + Beta)
	synon = 0
	gampa = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	gnmda = mgblock(v)*(Ron + Roff)*gmax*nmdafactor
	inmda = gnmda*(v - E)
	iampa = gampa*(v - E)
	i = iampa + inmda
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
	gampa' = -gampa/ampatau
}








FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}

FUNCTION plast(step(1))(1) {
	plast = 1 - 1/(1 + exp((step - sighalf)/sigslope))
}

NET_RECEIVE(weight, s, w, tlast (ms), r0, t0 (ms)) {
	INITIAL {
		
		w = weight*plast(s)
		tlast = -1e9 (ms)
		r0 = 0
		t0 = -1e9 (ms)
	}
	
        if (flag == 0) { 
		
		
		
		if (t - tlast < ltpinvl) { 
			s = s + 1
			if (s > 2*sighalf) { s = 2*sighalf }
		}else if (t - tlast > ltdinvl) { 
		}else{ 
			s = s - 1
			if (s < 0) { s = 0 }
		}
		tlast = t

		w = weight*plast(s)
		gampa = gampa + w*gmax*gampafactor
		r0 = r0*exp(-Beta*(t - t0))
		t0 = t
		synon = synon + w
		Ron = Ron + r0
		Roff = Roff - r0
		
		net_send(Cdur, w + 1)
        }else{ 
		r0 = (flag-1)*Rinf + (r0 - (flag-1)*Rinf)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - (flag-1)
		Ron = Ron - r0
		Roff = Roff + r0
	}
}