TITLE  K-D channel
: M.Migliore jun 2006

UNITS {
	(mA) = (milliamp) 
	(mV) = (millivolt)

}

PARAMETER {
    flag=0
      curr
	v (mV)
        ek (mV)		: must be explicitely def. in hoc
	celsius		(degC)
	gkdbar=.0 (mho/cm2)
        vhalfn=-33   (mV)
        a0n=0.01      (/ms)
        zetan=3    (1)
        gmn=0.7  (1)
	nmax=2  (1)
	q10=1
	sh = 0
	ek2 (mV)
	count=1
	vrun (mV)
	vvrun=0
	vrun2
	delta=0
	vinit=0 (mV)
	alpha=1.06
	sh2=0
	alphash0=0
	alphash1=0.15
	FK = 2
	PK = 1
	BK = 2.11
	CK = 48
	timestep=1000
	stim_moltK=1
}


NEURON {
	SUFFIX kd
	POINTER stim_i
	USEION k READ ek WRITE ik
        RANGE flag, curr,gkd,gkdbar, sh,ek2,vvrun, vrun2, stim_moltK
	GLOBAL ninf,taun,alpha
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkd
        taun
		stim_i
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkd = gkdbar*n
	ek2=ek+vvrun*alpha
	ik = gkd*(v-ek2)

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
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
	
}
   

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
		
			
	   
        a = alpn(v-sh2)
        ninf = 1/(1+a)
        taun = betn(v-sh2)/(qt*a0n*(1+a))
	if (taun<nmax) {taun=nmax/qt}
	
}














