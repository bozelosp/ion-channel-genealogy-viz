NEURON {
	SUFFIX na_stoch
	USEION na READ ena WRITE ina
	RANGE gbar,m0h0,m1h0,m2h0,m3h0,m0h1,m1h1,m2h1,m3h1,g,N
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 1000   	(pS/um2)	
	vshift = -5	(mV)		

	tha  = -43	(mV)		
	qa   = 7	(mV)		
	Ra   = 0.182	(/ms)		
	Rb   = 0.124	(/ms)		

	thi1  = -50	(mV)		
	thi2  = -75	(mV)		
	qi   = 5	(mV)	        
	thinf  = -72	(mV)		
	qinf  = 6.2	(mV)		
	Rg   = 0.0091	(/ms)		
	Rd   = 0.024	(/ms)		

	temp = 23	(degC)		
	q10  = 2.3			

	v 		(mV)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)

	gamma=10e-12		(S)		
	seed
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)

	minf 		hinf		
	mtau (ms)	htau (ms)	
	tadj				

	dt			(ms)
	area			(um2)
   	N			
	m0h0			
	m1h0			
	m2h0			
	m3h0			
	m0h1			
	m1h1			
	m2h1			
	m3h1			
}
 

STATE { m h }		

FUNCTION trap0(v,th,a,q) {
	if (fabs((v-th)/th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}


INITIAL { 
	m = 0
	h = 0
	N=abs(((1e-8*area*1e-4*gbar)/gamma)+0.5)				

	m0h0=0
	m1h0=0
	m2h0=0
	m3h0=0
	m0h1=N		
	m1h1=0
	m2h1=0
	m3h1=0

	set_seed(seed)
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ina = tadj*((m3h1*gamma)/(1e-8*area))*(v-ena)		
} 


DERIVATIVE states {   
        m' =  m
        h' =  h
	noise()
}

PROCEDURE noise() {
	LOCAL a,b,vm,p,m0h0_merk,m1h0_merk,m2h0_merk,m3h0_merk,m0h1_merk,m1h1_merk,m2h1_merk,m3h1_merk,am,bm,ah,bh

	vm=v+vshift
	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
        tadj = q10^((celsius - temp)/10)
	mtau = 1/tadj/(a+b)
	minf = a/(a+b)
	am=tadj*a
	bm=tadj*b


		
	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/tadj/(a+b)
	hinf = 1/(1+exp((vm-thinf)/qinf))	




	ah=hinf/htau
	bh=(1-hinf)/htau

	m0h0_merk=m0h0
	m1h0_merk=m1h0
	m2h0_merk=m2h0
	m3h0_merk=m3h0
	m0h1_merk=m0h1
	m1h1_merk=m1h1
	m2h1_merk=m2h1
	m3h1_merk=m3h1



	
	p=1-exp(-dt*(ah+3*am))				
	FROM ii=1 TO m0h0_merk {
		if (scop_random()<= p)	{					
			if (scop_random()<= ah/(ah+3*am))	{		

				m0h0=m0h0-1
				m0h1=m0h1+1
			}
			else {							
				m0h0=m0h0-1
				m1h0=m1h0+1
			}
		}
	}

	
	p=1-exp(-dt*(ah+2*am+bm))
	FROM ii=1 TO m1h0_merk {
		if (scop_random()<= p) {	
			if (scop_random()<= ah/(ah+2*am+bm)) {	
				m1h0=m1h0-1
				m1h1=m1h1+1
			}
			else {							
				if (scop_random()<= 2*am/(2*am+bm)) {	
					m1h0=m1h0-1
					m2h0=m2h0+1
				}
				else {
					m1h0=m1h0-1
					m0h0=m0h0+1
				}
			}
		}
	}

	
	p=1-exp(-dt*(ah+am+2*bm)) 
	FROM ii=1 TO m2h0_merk {
		if (scop_random()<= p){	
			if (scop_random()<= ah/(ah+am+2*bm))	{
				m2h0=m2h0-1
				m2h1=m2h1+1
			}
			else {
				if (scop_random()<= am/(am+2*bm))	{
					m2h0=m2h0-1
					m3h0=m3h0+1
				}
				else {
					m2h0=m2h0-1
					m1h0=m1h0+1
				}
			}
		}
	}

	
	p=1-exp(-dt*(ah+3*bm)) 
	FROM ii=1 TO m3h0_merk {
		if (scop_random()<= p){	
			if (scop_random()<= ah/(ah+3*bm))	{	
				m3h0=m3h0-1
				m3h1=m3h1+1
			}
			else {	
				m3h0=m3h0-1
				m2h0=m2h0+1
			}
		}
	}


	
	p=1-exp(-dt*(bh+3*am)) 
	FROM ii=1 TO m0h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+3*am))	{	
				m0h1=m0h1-1
				m0h0=m0h0+1
			}
			else {
				m0h1=m0h1-1
				m1h1=m1h1+1
			}
		}
	}

	
	p=1-exp(-dt*(bh+2*am+bm)) 
	FROM ii=1 TO m1h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+2*am+bm))	{	
				m1h1=m1h1-1
				m1h0=m1h0+1
			}
			else {
				if (scop_random()<= 2*am/(2*am+bm))	{	
					m1h1=m1h1-1
					m2h1=m2h1+1
				}
				else {
					m1h1=m1h1-1
					m0h1=m0h1+1
				}
			}
		}
	}

	
	p=1-exp(-dt*(bh+am+2*bm)) 
	FROM ii=1 TO m2h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+am+2*bm))	{	
				m2h1=m2h1-1
				m2h0=m2h0+1
			}
			else {
				if (scop_random()<= am/(am+2*bm))	{	
					m2h1=m2h1-1
					m3h1=m3h1+1
				}
				else {
					m2h1=m2h1-1
					m1h1=m1h1+1
				}
			}
		}
	}

	
	p=1-exp(-dt*(bh+3*bm)) 
	FROM ii=1 TO m3h1_merk {
		if (scop_random()<= p){	
			if (scop_random()<= bh/(bh+3*bm))	{	
				m3h1=m3h1-1
				m3h0=m3h0+1
			}	
			else {
				m3h1=m3h1-1
				m2h1=m2h1+1
			}
		}
	}
}