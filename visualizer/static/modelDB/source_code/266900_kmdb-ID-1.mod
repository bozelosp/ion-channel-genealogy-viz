TITLE CA1 KM channel from Mala Shah
: M. Migliore June 2006

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER { 
  curr

	v 		(mV)
	ek
	ek2=-80 (mV)
      sh=0 (mV)
	celsius 	(degC)
	gbar=.0001 	(mho/cm2)
       vhalfl=-42  	(mV)
   	kl=-4
      : kl=-2
       vhalft=-42 (mV)
        a0t=0.04  	(/ms)
      : a0t=0.009      	(/ms)

      : zetat=7    	(1)
      zetat=4
        gmt=.7	(1)
       	q10=5
	b0=60
	st=1
	count=1
	vrun (mV)
	delta=0
	vinit=-76.2 (mV)
	alpha=1.06
	 zetat2=4
     gmt2=.7	(1)
	 a0t2=0.04
	 b02=60
	 sh2
	alphash0=0
	alphash1=0.15
	vvrun=0
	 timestep=1000
	vrun2
	v0
	dv0
	ddv
	flag=0
	FK = 2
	PK = 1
	BK = 2.11
	CK = 48	
	stim_moltK=1
}


NEURON {
	SUFFIX    kmdb
	POINTER stim_i
	USEION k READ ek WRITE ik
      RANGE flag, curr,  gbar,ik,sh,taua,taub,ggbar,ek2,vrun,count,vvrun,vrun2,delta2, stim_moltK
      GLOBAL  inf, tau,alpha
}

STATE {
        m
}

ASSIGNED {
	ik (mA/cm2)
        inf
	tau
        taua
	taub
	ggbar (mho/cm2)
	
	stim_i

}

INITIAL {
    vrun=0
	vvrun=vrun
	rate(v,sh2)
	m=inf
	
}


BREAKPOINT {
	SOLVE state METHOD cnexp
	ggbar=gbar*m^st
    ek2=ek+vvrun*alpha
	ik = gbar*m^st*(v-ek2)
}


FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

FUNCTION alpt2(v(mV)) {
  alpt2 = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett2(v(mV)) {
  bett2 = exp(0.0378*zetat2*gmt2*(v-vhalft)) 
}
DERIVATIVE state {
        rate(v,sh2)
        if (m<inf) {tau=taua} else {tau=taub}
	m' = (inf - m)/tau
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
   
   
PROCEDURE rate(v (mV),sh2) { :callable from hoc
        LOCAL a,qt,i
        qt=q10^((celsius-35)/10)
		
		
        inf = (1/(1 + exp((v-vhalfl-sh2)/kl)))
        a = alpt(v-sh2)
        taua = b0 + bett(v-sh2)/(a0t*(1+a))
        taub = b02 + bett2(v-sh2)/(a0t2*(1+alpt2(v-sh2)))

}















