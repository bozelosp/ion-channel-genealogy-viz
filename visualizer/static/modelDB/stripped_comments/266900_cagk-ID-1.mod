UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


NEURON {
	SUFFIX mykca
	POINTER stim_i
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE flag, curr, gbar,ik,sh,taua,taub,ggbar,ek2,vrun,count,vvrun,taun2,vrun2,delta,delta2, stim_moltCa
	GLOBAL  oinf, tau,alpha
}

UNITS {
	FARADAY = (faraday)  (kilocoulombs)
	R = 8.313424 (joule/degC)
}

PARAMETER { curr
      v		(mV)
      sh=0		(mV)

	dt		(ms)
	ek		(mV)
	ek2=-80 (mV)
	celsius = 20	(degC)
	gbar = 0.01	(mho/cm2)	
	cai = 1e-3	(mM)

      d1 =0.92
     	d2 = 1
		d11 =0.92
     	d22 = 1
	k1 = 0.01 (mM)   
	k2 = 0.000003	(mM)  
	bbar = 0.48	(/ms) 
	abar = 0.28	(/ms)  
	k1b = 0.01	(mM)   
	k2b = 0.000003	(mM)  
	bbar2 = 0.48	(/ms) 
	abar2 = 0.28	(/ms)  
	count=1
	vrun=0 (mV)
	vvrun=0
	delta=0
	vinit=-76.2  
	alpha=1.1
	alphash1=0.15
	sh2
	delta1=0
	 timestep=1000
	vrun2
	v0
	dv0
	ddv
	FK = 2
	PK = 1
	BK= 2.11
	CK = 48
	stim_moltK=1
	flag=0



        st=1            (1)
}

ASSIGNED {
	ik		(mA/cm2)
	oinf
	tau		(ms)
	taub   (ms)
	taua      (ms)
	ggbar  (mho/cm2)
	taun2
	stim_i
	
      
}

INITIAL {
        vrun=0
	   vvrun=vrun
	    rate(v,cai,sh2)
        o=oinf
				
}

STATE {	o }		

BREAKPOINT {
	SOLVE state METHOD cnexp
	ggbar=gbar*o^st
	ek2=ek+vvrun*alpha
	ik = gbar*o^st*(v - ek2)
}

DERIVATIVE state {	
	rate(v, cai,sh2)
	if (o<oinf) {tau=taua} else {tau=taub}
	
	oinf=oinf*tau
	o' = (oinf - o)/tau
   
}

FUNCTION alp(v (mV), c (mM)) (1/ms) { 
	alp = c*abar/(c + exp1(k1,d1,v))
}

FUNCTION bet(v (mV), c (mM)) (1/ms) { 
	bet = bbar/(1 + c/exp1(k2,d2,v))
}

FUNCTION alp2(v (mV), c (mM)) (1/ms) { 
	alp2 = c*abar2/(c + exp1(k1b,d11,v))
}

FUNCTION bet2(v (mV), c (mM)) (1/ms) { 
	bet2 = bbar2/(1 + c/exp1(k2b,d22,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { 
	exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius))
}

BEFORE STEP { LOCAL i

         if(stim_i==0 && flag==0){ 
		  vrun=0
		  vvrun=0
		  
	    }else{
		 flag=1
		            		  
		delta=v-vinit
		if (count<timestep+1){
			vrun= (delta-vrun)*(FK/(count+1))+vrun
			vrun2=vrun 
		 }else{
		
		vrun2= (delta)*(FK/(timestep+1))+vrun2*pow((1-FK/(timestep+1)),PK)
			
			}
		 
	   vvrun=(BK*vrun2/(1+vrun2/CK))
	    
		count=count+1   
        }
							
		 sh2=sh+alphash1*vvrun
	
}
 
   

PROCEDURE rate(v (mV), c (mM),sh2) { 
	LOCAL a,i,ema
	
	    
	 
	a = alp(v-sh2,c)
	taua = 1/(a + bet(v-sh2, c))
	taub=1/(alp2(v-sh2,c) + bet2(v-sh2, c))
	
	
	oinf=a
	
	
	

	
}