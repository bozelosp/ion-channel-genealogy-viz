TITLE K-A channel from Klee Ficker and Heinemann
: modified to account for Dax A Current ----------
: M.Migliore Jun 1997

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER { 
     
     curr
	celsius
        v (mV)
        gkabar=.008 (mho/cm2)
        vhalfn=-1   (mV)
        vhalfl=-56   (mV)
        a0l=0.05      (/ms)
        a0n=.1    (/ms)
        zetan=-1.8    (1)
        zetal=3    (1)
        gmn=0.39   (1)
        gml=1   (1)
        lmin=2  (mS)
        nmin=0.2  (mS)
        pw=-1    (1)
        tq=-40
        qq=5
        q10=5
        qtl=1
	ek
	ek2=-80 (mV)
	  sh=0
	  count=1
	delta=0
	vinit=-76.2
	alpha=1.06
	sh2
	alphash0=0
	alphash1=0.15
	vvrun=0
	vrun=0
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
        SUFFIX kad
		POINTER stim_i
        USEION k READ ek WRITE ik
        RANGE flag, curr, gkabar,gka,sh,ek2,vrun,count,vvrun,vrun2,delta2,stim_moltK
        GLOBAL  ninf,linf,taul,taun,lmin,gmn,zetal,alpha
}

STATE {
        n
        l
}

ASSIGNED {
        ik (mA/cm2)
        ninf
        linf      
        taul
        taun
        gka
		stim_i

		
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gka = gkabar*n*l
		ek2=ek+vvrun*alpha
	    ik = gka*(v-ek2)

}

INITIAL {
     vrun=0
	vvrun=vrun
	rates(v,sh2)
	n=ninf
	l=linf
}


FUNCTION alpn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  alpn = exp(1.e-3*zeta*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
LOCAL zeta
  zeta=zetan+pw/(1+exp((v-tq)/qq))
  betn = exp(1.e-3*zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alpl(v(mV)) {
  alpl = exp(1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betl(v(mV)) {
  betl = exp(1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
 
}

DERIVATIVE states {  
        rates(v,sh2)
        n' = (ninf - n)/taun
        l' = (linf - l)/taul
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
    
PROCEDURE rates(v (mV),sh2) { :callable from hoc
        LOCAL a,qt,i
        qt=q10^((celsius-24)/10)
		
				
        a = alpn(v-sh2)
        ninf = 1/(1 + a)
        taun = betn(v-sh2)/(qt*a0n*(1+a))
        if (taun<nmin) {taun=nmin}
        a = alpl(v-sh2)
        linf = 1/(1+ a)
        taul = 0.26*(v+50-sh2)/qtl
        if (taul<lmin/qtl) {taul=lmin/qtl}
		

}

