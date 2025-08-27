NEURON {
	POINT_PROCESS GRC_GABA
	NONSPECIFIC_CURRENT i
	RANGE g,Cdur,Erev,Open,OpenScaled,ScaleFactor

	
	RANGE kon,koff,d3,r3,d1d2,r1r2,a1,b1,a2,b2,r1,r2,d1,d2

	RANGE Tmax,gmax,onSET	
	
	RANGE tau_1,tau_rec,tau_facil,U,T 
	RANGE diff_flag,M,Rd,Diff,lamd
	RANGE diffusione,nd	
}

UNITS {
	(nA) 	= (nanoamp)
	(mV) 	= (millivolt)
	(umho)  = (micromho)
	(mM) 	= (milli/liter)
	(pS) 	= (picosiemens)
	PI   	= (pi)(1)
}

PARAMETER {
	
	gmax	=  756.35	(pS)	
	Cdur	= 0.3		(ms)	

	kon	= 20		(/ms/mM)	
	koff	= 2		(/ms) 		 
	d3	= 15		(/ms) 		 
	r3	= 3.75		(/ms) 		

	d1d2	= 15		(/ms/mM)	
	r1r2	= 0.007		(/ms)

	a1	= 0.06		(/ms)
	b1	= 0.03		(/ms)
	a2	= 0.4		(/ms)
	b2	= 10		(/ms)
	
	r1	= 7e-4		(/ms)
	r2	= 6e-3		(/ms)
	d1	= 3.3e-4	(/ms)
	d2	= 1.2		(/ms)

	Erev	= -65	(mV)

	
	tau_1 		= 0.1 (ms) 	< 1e-9, 1e9 >
	tau_rec 	= 43.4 (ms) 	< 1e-9, 1e9 > 	
	tau_facil 	= 6.22 (ms) 	< 0, 1e9 >    	
	U 		= 0.35		< 0, 1 >	
	
	Tmax	= 1  (mM)	
	onSET	= 1
		
	
	

	M		= 52.76		

	Rd		= 4.79 (um)	
	Diff		= 0.223 (um2/ms)
	lamd		= 20 	(nm)
	diff_flag	= 1			
	nd		= 1			

	ScaleFactor	= 1 			
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(pS)		
	Open
	OpenScaled	
	
	T		(mM)	
	Trelease	(mM)
	Mres		(mM)	
	tpre		(ms)

	tspike[50]	(ms)	
	PRE[50]
	numpulses
	tzero
}

STATE {	
	C
	CA1
	CA2
	DA1
	DA2
	DA2f
	OA1
	OA2	
}

INITIAL {
	C=1
	CA1=0
	CA2=0
	DA1=0
	DA2=0
	DA2f=0
	OA1=0  	
	OA2=0
	CA1=0
	CA2=0
	Open=0
	T=0 		(mM)
	
		
	numpulses=0
	Mres=1e3* (1e3 * 1e15 / 6.022e23 * M)     
	FROM i=1 TO 50{ PRE[i-1]=0 tspike[i-1]=0}
	tspike[0]=1e12	(ms)
	if(tau_1>=tau_rec){ 
		printf("Warning
		tau_rec=tau_1+1e-5
		
	} 
}

FUNCTION diffusione(){
	LOCAL DifWave,i	
	DifWave=0
	FROM i=1 TO numpulses{
		tzero=tspike[i-1]
		if(t>tzero){
			DifWave=DifWave+PRE[i-1]*Mres*exp(-Rd*Rd/(4*Diff*(t-tzero)))/((4*PI*Diff*(1e-3)*lamd)*(t-tzero))^nd
		}
	}	
	diffusione=DifWave 
}


BREAKPOINT {
	SOLVE kstates METHOD sparse
	Open = OA1 + OA2
	OpenScaled=Open*ScaleFactor
	g = gmax * Open
	i = (1e-6) * g * (v - Erev)
}

KINETIC kstates {
	if ( diff_flag ) { Trelease = T + diff_flag * diffusione() } else { Trelease = T }
	
	~	C  	<-> 	CA1	(2*kon*Trelease,koff)
	~	CA1 	<-> 	CA2	(kon*Trelease,2*koff)
	~	CA2	<->	DA2f	(d3,r3)
	
	~ 	DA1  	<-> 	DA2	(d1d2*Trelease,r1r2)
	
	~ 	OA1  	<-> 	CA1	(a1,b1)
	~ 	OA2  	<-> 	CA2	(a2,b2)
	
	~	DA1	<->	CA1	(r1,d1)
	~	DA2	<->	CA2	(r2,d2)
	CONSERVE C+CA1+CA2+DA1+DA2+DA2f+OA1+OA2 = 1
}


NET_RECEIVE(weight, on, nspike, tzero (ms),x,y, z, u, tsyn (ms)) {

	INITIAL {
		x = 0
		y = 0
		z = 0
		u = 0 
		tsyn = t
		nspike = 1
	}

	if(onSET){on=0 onSET=0}
        if (flag == 0) { 
		
		nspike = nspike + 1
		if (!on) {
			tzero = t	
			tpre=t	
			on = 1				
			z = z*exp(-(t - tsyn)/tau_rec)		
			z = z + ( y*(exp(-(t - tsyn)/tau_1) - exp(-(t - tsyn)/tau_rec)) / ((tau_1/tau_rec)-1) )
			y = y*exp(-(t - tsyn)/tau_1)			
			x = 1-y-z
				
			if (tau_facil > 0) { 
				u = u*exp(-(t - tsyn)/tau_facil)
				u = u + U * ( 1 - u )							
			} else { u = U }
			
			y = y + x * u
			T=Tmax*y
				
			PRE[numpulses]=y
			tspike[numpulses]=t
			numpulses=numpulses+1
			tsyn = t			
		}
		net_send(Cdur, nspike)
        }
	if (flag == nspike) { 
			T = 0
			on = 0
	}
}