TITLE I-h channel from Magee 1998 for distal dendrites

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER { 
   
   curr
	v 		(mV)
	ehd  = -30		(mV)        
	celsius 	(degC)
	ghdbar=.0001 	(mho/cm2)
        vhalfl=-81   	(mV)
	kl=-8
        vhalft=-75   	(mV)
        a0t=0.011      	(/ms)
        zetat=2.2    	(1)
        gmt=.4   	(1)
	q10=4.5
	qtl=1
    sh=0 (mV)
	time1=600
	time0=100
	alphash0=0
	alphash1=0.13
	alphahd=0.3
	ehd2
	sh2
	count=1
	vrun (mV)
	delta=0
	vinit=-76.2  
	vvrun=0
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


NEURON {
	SUFFIX hd
	POINTER stim_i
	NONSPECIFIC_CURRENT i
        RANGE flag, curr, ghdbar, vhalfl, sh,count,delta2,vrun2, stim_moltNa, ehd2
        GLOBAL  linf,taul
}

STATE {
        l
}

ASSIGNED {
	i (mA/cm2)
        linf      
        taul
        ghd
	stim_i	
}

INITIAL {
	rate(v,sh2)
	l=linf
	vrun=0
	vvrun=vrun
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l
	ehd2=ehd-alphahd*vvrun
	i = ghd*(v-ehd2)

}


FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rate(v,sh2)
        l' =  (linf - l)/taul
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
   
   
PROCEDURE rate(v (mV),sh2) { :callable from hoc
        LOCAL a,qt,i
        qt=q10^((celsius-33)/10)
		
			   
        a = alpt(v-sh2)
        linf = 1/(1 + exp(-(v-vhalfl-sh2)/kl))
:       linf = 1/(1+ alpl(v-sh2))
        taul = bett(v-sh2)/(qtl*qt*a0t*(1+a))
}














