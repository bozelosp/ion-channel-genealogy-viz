NEURON {
	POINT_PROCESS ltpltd
	RANGE e, i, g2, peso, U, tau_facil, tau_rec, tau_1, Rin, Ase
	RANGE   u0, Pini, Nini, deltap, deltad, f, nip, nid, gamma, eta, VVini
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	
	
	e = 0	(mV)
	
	
	tau_1 = 3 (ms) < 1e-9, 1e9 >
	
	
	tau_rec = 50 (ms) < 1e-9, 1e9 >
	
	
	tau_facil = 200 (ms) < 0, 1e9 >
	
	
	
	
	U = 0.36(1) < 0, 1 >
	
	u0 = 0 (1) < 0, 1 > 
	
	f = 0.05e-3 
	deltap = 400
	deltad = 400
	gamma = 0.2
	eta = 2e-3
	nip = 0.0987
	nid = 0.07
	lambdap = 1e-3
	lambdad = 2e-3
	mp = 3e-3
	md = 3e-3
	ap = 2
	ad = 0.5
	taum = 40
	Rin=10e7
	Ase=2.5e-7
	Pini=0
	Nini=0
	VVini=0
	g2=43 (umho)
	peso=0.5e-6 
	}

ASSIGNED {
	v (mV)
	i (nA)
	x
	y
	}

STATE {
	g (umho)
	C
	Np
	Nd
	VV (mV)
	}

INITIAL {
	VV=VVini
	g=0
	C =0
	Np=Pini
	Nd=Nini 
	}

	
FUNCTION v2(v) { if (v<-65) { v2 = 0	} else {v2 = v+65}
	}
	
BREAKPOINT {
	
	SOLVE state METHOD derivimplicit
	i = -g2*VV
}

DERIVATIVE state {

	g' = -g/tau_1
	C' = gamma*(VV)-eta*C+peso*(v2(v))
	Np' = nip*C-(lambdap+g*deltap)*Np+((mp*Np*Np)/(ap+Np*Np))
	Nd' = nid*C-(lambdad+g*deltad)*Nd+((md*Nd*Nd)/(ad+Nd*Nd))
	VV' = -(VV/taum)+Rin*Ase*g*((1/taum)+f*(deltap*Np-deltad*Nd))
		}



		

NET_RECEIVE(weight (umho), y, z, u, tsyn (ms)) 
{
INITIAL {

	y = 0
	z = 0

	u = u0
	tsyn = t


}
	
	
	z = z*exp(-(t - tsyn)/tau_rec)
	z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )
	
	y = y*exp(-(t - tsyn)/tau_1)

	x = 1-y-z

	
	if (tau_facil > 0) {
		u = u*exp(-(t - tsyn)/tau_facil)
	} else {
		u = U
	}

 

	if (tau_facil > 0) {
		state_discontinuity(u, u + U*(1-u))
	}



	state_discontinuity(g, g + weight*x*u)
	state_discontinuity(y, y + x*u)
	
	tsyn = t



}