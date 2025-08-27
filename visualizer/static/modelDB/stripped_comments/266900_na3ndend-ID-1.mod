NEURON {
	SUFFIX na3d
	POINTER stim_i
	USEION na READ ena WRITE ina
	RANGE flag, curr,  gbar, sh,thegna,ena2,alpha,vrun,count,vvrun,vrun2,delta2, stim_moltNa
	GLOBAL  minf, hinf, mtau, htau, sinf, taus,qinf, thinf,ar,sh2
}

PARAMETER { curr
   
	sh   = 0	(mV)
	gbar = 0.010   	(mho/cm2)	
								
	tha  =  -25	(mV)		
	qa   = 7.2	(mV)		
	Ra   = 0.4	(/ms)		
	Rb   = 0.124 	(/ms)		

	thi1  = -45	(mV)		
	thi2  = -45 	(mV)		
	qd   = 1.5	(mV)	        
	qg   = 1.5      (mV)
	mmin=0.02	
	hmin=0.5			
	q10=2
	Rg   = 0.01 	(/ms)		
	Rd   = .03 	(/ms)		
	qq   = 10        (mV)
	tq   = -55      (mV)

	thinf  = -50 	(mV)		
	qinf  = 2 	(mV)		

        vhalfs=-60	(mV)		
        a0s=0.0003	(ms)		
        zetas=12	(1)
        gms=0.2		(1)
        smax=10		(ms)
        vvh=-58		(mV) 
        vvs=2		(mV)
        ar=1	(1)		
	ena		(mV)            
	ena2		(mV)  
	celsius
	v 		(mV)
	vrun=0  (mV)
	vinit=-76.2  (mV)
	delta=0
    count=1
	alphaena=0.25
	sh2
	alphash0=0
	alphash1=0.13
	 timestep=1000
	vrun2
	v0
	dv0
	ddv
	flag=0
	FNa = 2
	PNa = 1
	BNa = 2.6
	CNa = 60
	stim_moltNa=1
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	thegna		(mho/cm2)
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
	sinf (ms)	taus (ms)
	vvrun
	stim_i
}
 

STATE { m h s}

BREAKPOINT {
        SOLVE states METHOD cnexp
        thegna = gbar*m*m*m*h*s
	    ena2=ena-alphaena*vvrun
	    ina = thegna * (v - ena2)
} 

INITIAL {
	trates(v,ar,sh2)
	m=minf  
	h=hinf
	s=sinf
	vrun=0
	vvrun=vrun
}


FUNCTION alpv(v(mV)) {
         alpv = 1/(1+exp((v-vvh)/vvs))
}
        
FUNCTION alps(v(mV)) {  
  alps = exp(1.e-3*zetas*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bets(v(mV)) {
  bets = exp(1.e-3*zetas*gms*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
}

LOCAL mexp, hexp, sexp

DERIVATIVE states {   
        trates(v,ar,sh2)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        s' = (sinf - s)/taus
}

BEFORE STEP { LOCAL i
       
	    if(stim_i==0 && flag==0){ 
		  vrun=0
		  vvrun=0
		  
	    }else{
		 flag=1
		             		
		delta=v-vinit
	    if (count<timestep+1){
		  vrun= (delta-vrun)*(FNa/(count+1))+vrun
          vrun2=vrun 
		 }else{

		vrun2= (delta)*(FNa/(timestep+1))+vrun2*pow((1-FNa/(timestep+1)),PNa)
			
			}
	    
	   vvrun=(BNa*vrun2/(1+vrun2/CNa))
	   
		count=count+1   
        }					
		 sh2=sh+alphash1*vvrun
	
}
  
   
PROCEDURE trates(vm,a2,sh2) {  
        LOCAL  a, b, c, qt,i
        qt=q10^((celsius-24)/10)
	   	
	a = trap0(vm,tha+sh2,Ra,qa)
	b = trap0(-vm,-tha-sh2,Rb,qa)
	mtau = (1/(a+b))/qt
        if (mtau<mmin) {mtau=mmin}
	minf = a/(a+b)

	a = trap0(vm,thi1+sh2,Rd,qd)
	b = trap0(-vm,-thi2-sh2,Rg,qg)
	htau =  (1/(a+b))/qt
        if (htau<hmin) {htau=hmin}
	hinf = 1/(1+exp((vm-thinf-sh2)/qinf))
	c=alpv(vm-sh2)
        sinf = c+a2*(1-c)
        taus = bets(vm-sh2)/(a0s*(1+alps(vm-sh2)))
        if (taus<smax) {taus=smax}
		
	
  
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}