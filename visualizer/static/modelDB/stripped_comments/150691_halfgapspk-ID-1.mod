NEURON {
	POINT_PROCESS HalfGapSpk
	
	
	
	ELECTRODE_CURRENT i
	RANGE g, i, vgap, meang, meant, rg, rt, drift
	THREADSAFE 
	POINTER donotuse 
	RANGE id 
		 
		 
		 
		 
	RANGE gmax, gmin, vhalf 
		
		
		
		
}

PARAMETER {
	gmax = 1 (nanosiemens)
	gmin = 1 (nanosiemens)
	vhalf = 0 (millivolt)
	slope4 = 10 (/millivolt)
	meang = 30 (nanosiemens)
	meant = 1000000 (ms)
	drift = 0
	rg=0
	rt=0
	id = 0
}

ASSIGNED {
	g (nanosiemens)
	v (millivolt)
	vgap (millivolt)
	i (nanoamp)
	donotuse
}

INITIAL {
}






FUNCTION gv(x(millivolt))(nanosiemens) {
	
	gv = (gmax - gmin)/(1 + exp(slope4*(vhalf - x))) + gmin
}

BREAKPOINT {
	LOCAL x
	if (gmax == gmin) { 
		g = gmax
		i = g * (vgap - v) * (.001)
	}else{
		
		if (id > 0 ) {
			x = v - vgap 
		}else if (id < 0){
			x = vgap - v 
		}else{
VERBATIM
			assert(0);
ENDVERBATIM

		}
		g = gv(x)
		i = g * (vgap - v) * (.001)
	}
}

PROCEDURE getpar() {
	gmax=mynormrand(meang/1(nanosiemens),rg)*1(nanosiemens)
	if (gmax<0) {gmax=0}
	if (gmin != 0) {
		gmin = gmax
	}
	meang=meang+drift*meang
	rg=rg+drift*rg
}

NET_RECEIVE (w) {
	getpar()		
}













VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
double nrn_random_pick(void* r); 
void* nrn_random_arg(int argpos);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif
ENDVERBATIM


FUNCTION mynormrand(mean, var) {
VERBATIM
	if (_p_donotuse) {
		double x = nrn_random_pick(RANDCAST _p_donotuse);
		_lmynormrand = x*_lvar + _lmean;
	}else{
		_lmynormrand = _lmean;
	}
ENDVERBATIM
}

PROCEDURE setRandom() {
VERBATIM
 {
        void** pv = (void**)(&_p_donotuse);
        if (ifarg(1)) {
                *pv = nrn_random_arg(1);
        }else{
                *pv = (void*)0;
        }
 }
ENDVERBATIM
}