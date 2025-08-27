TITLE K-DR channel
: from Klee Ficker and Heinemann
: modified to account for Dax et al.
: M.Migliore 1997

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
    curr

	v (mV)
        ek (mV)		: must be explicitely def. in hoc
	 ek2=-80 (mV)		
	celsius		(degC)
	gkdrbar=.003 (mho/cm2)
        vhalfn=13   (mV)
		 vhalfn2=-60  (mV)
		
		
        a0n=0.02      (/ms)
        zetan=-3    (1)
		
		 a0n2=0.02    (/ms)
        zetan2=-3  (1)
		 gmn2=0.68 (1) :0.7
		 
        gmn=0.68  (1) :0.7
	nmax=2  (1)
	q10=1
	sh= 0	(mV)
	vrun=0  (mV)
	vinit=-76.2  (mV)
	alpha=1.06
	delta=0
    count=1
	sh2=0
	alphash0=0
	alphash1=0.15
	vvrun=0
	 timestep=1000
	vrun1
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
	SUFFIX kdr
	POINTER stim_i
	USEION k READ ek WRITE ik
    RANGE flag, curr, gkdr,gkdrbar,sh,taua,taub,ek2,vrun,count,delta,vvrun,vrun2,delta2,vrun1, stim_moltK
	GLOBAL  ninf,taun,zetan,vhalfn,vhalfn2,zetan2,alpha
}

STATE {
	n
	wrun
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkdr
        taun
		taua
		taub
		taun2
		
	stim_i	
		
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gkdrbar*n	
	ek2=ek+vvrun*alpha
	ik = gkdr*(v-ek2)

}

INITIAL {
    vrun=0
	vrun1=vrun
	vvrun=vrun
	rates(v,sh2)
	n=ninf
	 v0=vinit

	
	
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}



FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpn2(v(mV)) {
  alpn2 = exp(1.e-3*zetan2*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}



FUNCTION betn2(v(mV)) {
  betn2 = exp(1.e-3*zetan2*gmn2*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v,sh2)
		if (n<ninf) {taun=taua} else {taun=taub}
        n' = (ninf - n)/taun
	
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
		 :printf("stim %f %\n",stim_i) 
		 
	
}
   

   
PROCEDURE rates(v (mV),sh2) { :callable from hoc
        LOCAL a,qt,aa,i
        qt=q10^((celsius-24)/10)
		
	    a = alpn(v-sh2)
        ninf = 1/(1+a)
        taua= betn(v-sh2)/(qt*a0n*(1+a))
		 
	if (taua<nmax) {taua=nmax}
	    aa = alpn2(v-sh2)
	    taub= betn2(v-sh2)/(qt*a0n2*(1+aa))
		    

}














