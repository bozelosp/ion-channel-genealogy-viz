NEURON {
	POINT_PROCESS GRANULE_Gaba_stoch_vi
	 
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,Q10_channel
	RANGE g, ic ,Cdur ,Tmax, Erev
	RANGE Open,Open_a6	
	

	RANGE kon,koff,d3,r3,d1d2,r1r2,a1,b1,a2,b2,r1,r2,d1,d2
	RANGE kon_a6,koff_a6,d3_a6,r3_a6,d1d2_a6,r1r2_a6,a1_a6,b1_a6,a2_a6,b2_a6,r1_a6,r2_a6,d1_a6,d2_a6

	RANGE gmaxA1,gmaxA6,onSET	
	RANGE gA1,gA6

	RANGE xout,yout,zout,uout	
	
	
    
    RANGE T,U_2 
	RANGE diff_flag,diff_flag2,M,Rd,Diff,lamd
	RANGE nd,diffus ,Trelease, gmax_factor, syntype 
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
	syntype
	gmax_factor = 1
	
	gmaxA1	= 918.807	(pS)	
	gmaxA6	= 132.842	(pS)	
	Q10_diff	= 1.5
	Q10_channel	= 2.4
	Cdur	= 0.3		(ms)	
	Cdur_2 	= 0.3 		(ms)

	
	kon	= 17.7		(/ms/mM)	
	koff	= 2		(/ms) 		 
	d3	= 21.703	(/ms) 	
	r3	= 0.97325	(/ms) 	
	d1d2	= 0		(/ms/mM) 
	r1r2	= 0		(/ms)	
	a1	= 0.4144	(/ms)	
	b1	= 0.03		(/ms)
	a2	= 1.0002	(/ms)	
	b2	= 10		(/ms)	
	r1	= 7e-4		(/ms)
	r2	= 0.14208	(/ms)	
	d1	= 3.3e-4	(/ms)
	d2	= 3.4898	(/ms)	
	
	
	kon_a6	= 54.8		(/ms/mM)	
	koff_a6	= 0.31061	(/ms) 		
	d3_a6	= 15		(/ms) 		 
	r3_a6	= 7.4132	(/ms) 		
	d1d2_a6	= 24.2		(/ms/mM)	
	r1r2_a6	= 0.091668	(/ms)		
	a1_a6	= 0.06		(/ms)
	b1_a6	= 0.03		(/ms)
	a2_a6	= 0.4		(/ms)
	b2_a6	= 10		(/ms)	
	r1_a6	= 0.040001	(/ms)		
	r2_a6	= 0.4316	(/ms)		
	d1_a6	= 0.86042	(/ms)		
	d2_a6	= 2.7012	(/ms)		
	
	Erev	= -65		(mV)

	
	
	
	
	
	 
	T_factor = 0.5
	Tmax	= 1  		(mM)	
	onSET	= 1
	 
		
	
	

	M	= 7.506	
	Rd	= 0.978 (um)	
	Diff	= 0.223 (um2/ms)

	lamd	= 20 		(nm)
	diff_flag	= 1			
	 
	nd		= 1			

	 
	celsius (degC)
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	ic 		(nA)		
	g 		(pS)		
	gA1		(nA)
	gA6		(nA)
	
	Open
	Open_a6
 
	diffus	
 
	T		(mM)	
	 
	Trelease	(mM)
	 	
	Mres		(mM)	
	tpre		(ms)

	xout
	yout
	zout
	uout
	
	tspike[100]	(ms)	
	PRE[100]
	
	numpulses
	tzero
	gbar_Q10 (mho/cm2)
	Q10 (1)
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

	
	C_a6
	CA1_a6
	CA2_a6
	DA1_a6
	DA2_a6
	DA2f_a6
	OA1_a6
	OA2_a6
 
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

	
	C_a6=1
	CA1_a6=0
	CA2_a6=0
	DA1_a6=0
	DA2_a6=0
	DA2f_a6=0
	OA1_a6=0  	
	OA2_a6=0
	CA1_a6=0
	CA2_a6=0

 
	Open=0
	  
	T=0 		(mM)
	
	gbar_Q10 = Q10_diff^((celsius-30)/10)
	Q10 = Q10_channel^((celsius-30)/10)
 
	numpulses=0
	Mres=1e3* (1e3 * 1e15 / 6.022e23 * M)     
	FROM i=1 TO 100 { PRE[i-1]=0 tspike[i-1]=0 } 
	tspike[0]=1e12	(ms)
    
	
	
	
	
	
    
	onSET	= 1
	 
}
	FUNCTION imax(a,b) {
	    if (a>b) { imax=a }
	    else { imax=b }
	}
	
FUNCTION diffusione(){	 
	LOCAL DifWave,i,cntc,fi
	DifWave=0
	cntc=imax(numpulses-100,0) 
	FROM i=cntc  TO numpulses{
	    fi=fmod(i,100)
	    tzero=tspike[fi]
		if(t>tzero){
			DifWave=DifWave+PRE[fi]*Mres*exp(-Rd*Rd/(4*Diff*(t-tzero)))/((4*PI*Diff*(1e-3)*lamd)*(t-tzero))^nd
		}
	}	
	diffusione=DifWave
}



BREAKPOINT {
	SOLVE kstates METHOD sparse
	
	Open=OA1+OA2
	Open_a6=OA1_a6+OA2_a6
	
	gA1 = gmaxA1 * gbar_Q10 * Open
	gA6 = gmaxA6 * gbar_Q10 * Open_a6
	
	g = (gA1 + gA6) * gmax_factor
	
	i = (1e-6) * g * ( v - Erev )
	ic = i
}

KINETIC kstates {

	
	
	
	

	diffus=diffusione() 	
	Trelease = T + diffus
	
	
	
	
	~	C  	<-> 	CA1	(2*kon*Trelease*Q10,koff*Q10)
	~	CA1 	<-> 	CA2	(kon*Trelease*Q10,2*koff*Q10)
	~	CA2	<->	DA2f	(d3*Q10,r3*Q10)
	
	~ 	DA1  	<-> 	DA2	(d1d2*Trelease*Q10,r1r2*Q10)
	
	~ 	OA1  	<-> 	CA1	(a1*Q10,b1*Q10)
	~ 	OA2  	<-> 	CA2	(a2*Q10,b2*Q10)
	
	~	DA1	<->	CA1	(r1*Q10,d1*Q10)
	~	DA2	<->	CA2	(r2*Q10,d2*Q10)	
	CONSERVE C+CA1+CA2+DA1+DA2+DA2f+OA1+OA2=1
	
	
	
	
	~	C_a6  	<-> 	CA1_a6	(2*kon_a6*diffus*Q10,koff_a6*Q10)
	~	CA1_a6 	<-> 	CA2_a6	(kon_a6*diffus*Q10,2*koff_a6*Q10)
	~	CA2_a6	<->	DA2f_a6	(d3_a6*Q10,r3_a6*Q10)
	
	~ 	DA1_a6  <-> 	DA2_a6	(d1d2_a6*diffus*Q10,r1r2_a6*Q10)
	
	~ 	OA1_a6  <-> 	CA1_a6	(a1_a6*Q10,b1_a6*Q10)
	~ 	OA2_a6  <-> 	CA2_a6	(a2_a6*Q10,b2_a6*Q10)
	
	~	DA1_a6	<->	CA1_a6	(r1_a6*Q10,d1_a6*Q10)
	~	DA2_a6	<->	CA2_a6	(r2_a6*Q10,d2_a6*Q10)	
	CONSERVE C_a6+CA1_a6+CA2_a6+DA1_a6+DA2_a6+DA2f_a6+OA1_a6+OA2_a6=1	
	
}


NET_RECEIVE(weight,on, nspike,tzero (ms),x,y,z,u, tsyn (ms)) {LOCAL fi
	INITIAL {
		
		
		
		
		 		
		
		
		
		
        
		
        
		nspike = 1
	}

	if(onSET){
		on=0 
		onSET=0
		 
	}
	
	if (flag == 0) { 
		
		nspike = nspike + 1
		
		if (!on) {
			tzero = t	
			tpre=t	
			on = 1		
					
			
			
			
			
			
			
			
			
			
            
			T=Tmax*T_factor 

			fi=fmod(numpulses,100)
			PRE[fi]=T_factor 
			 
			tspike[fi] = t
			numpulses=numpulses+1

			
		}
		net_send(Cdur, nspike)
		 
	}

	if (flag == nspike) { 
		T = 0
		on = 0
	}

	 


}