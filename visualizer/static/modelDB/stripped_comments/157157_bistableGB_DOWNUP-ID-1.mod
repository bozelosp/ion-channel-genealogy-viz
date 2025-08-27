INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX BistableGBdownup
	USEION ca READ cai
	RANGE pDOWN,pUP,pDOWN0,pUP0, C1, tau, ps, gp, gd, thp, thd, tp, td, beta, w0, w1, b, weight
	
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

PARAMETER {
	pDOWN0	= 0 (1)	
        pUP0	= 1 (1)	
	tau	= 700000 (ms)	
	ps	= 0.5 (1)	
	gp	= 1600 (1)	
	gd	= 300 (1)	
	thp	= 1.3 (uM)	
	thd	= 1.0 (uM)	
	tp	= 20000 (ms)	
	td	= 20000 (ms)	
        beta    = 0.5 (1)	
	w0	= 0.5 (1)	
	w1	= 1.5 (1)	

}

ASSIGNED {
	cai (mM)	
	C1 (uM)	        
        b (1)		
	weight (1)      
}

STATE {
	pDOWN	
        pUP	

}

INITIAL {
	pDOWN=pDOWN0
        pUP=pUP0
        b=w1/w0
        weight=(((1-wpDOWN(pDOWN))*beta+wpUP(pUP)*(1-beta))+b*(wpDOWN(pDOWN)*beta+(1-wpUP(pUP))*(1-beta)))/(beta+(1-beta)*b)

}

BREAKPOINT {
	SOLVE states METHOD euler
        
        weight=(((1-wpDOWN(pDOWN))*beta+wpUP(pUP)*(1-beta))+b*(wpDOWN(pDOWN)*beta+(1-wpUP(pUP))*(1-beta)))/(beta+(1-beta)*b)

}

DERIVATIVE states {
	
	C1 = 1000*cai

	pDOWN' = (-pDOWN*(1-pDOWN)*(ps-pDOWN)+gp*(1-pDOWN)*heavi(C1-thp)-gd*pDOWN*heavi(C1-thd))/tau

        pUP' = (-pUP*(1-pUP)*(ps-pUP)+gp*(1-pUP)*heavi(C1-thp)-gd*pUP*heavi(C1-thd))/tau

 
}

FUNCTION heavi(x (uM)) (1) {
	heavi = 0
	if (x>0) {
	  heavi = 1
	}
}


FUNCTION wpUP(x) (1) {
	wpUP=0
	if (x<=0.5) {
	   wpUP = 1
	}
}


FUNCTION wpDOWN(x) (1) {
	wpDOWN=0
	if (x>=0.5) {
	   wpDOWN = 1
	}
}