DEFINE SYN_MAX_NUM 6000 

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(S) = (siemens)
}

NEURON{ 
	POINT_PROCESS VectorSynNS

	GLOBAL eI,eE 
	GLOBAL gmax1,gmax2,gmax3,gmaxI 
	GLOBAL factor1,factor2,factor3,factorI
	GLOBAL CL3f 
	
	GLOBAL CL1tau1,CL2tau1,CL3tau1,CLItau1,CL1tau2,CL2tau2,CL3tau2,CLItau2
	GLOBAL tau_FCL3,tau_FCLI 
	GLOBAL tau_DCL1,tau_DCL3,tau_DCLI 
	GLOBAL FCL3,FCLI,DCL1,DCL3,DCLI 
	
	
	
	RANGE intrains_CL1,intrains_CL2,intrains_CL3,intrains_CLI 
	
	RANGE F,D,time
	RANGE syntype,syncount 
	
	RANGE tsyn
	RANGE count
	RANGE k
	
	RANGE g1,g2,g3,gi
	
	NONSPECIFIC_CURRENT i
	
}
	
PARAMETER {

	syncount	 
	syntype[SYN_MAX_NUM]	
	tsyn[SYN_MAX_NUM]		
	count[SYN_MAX_NUM]	
  	F[SYN_MAX_NUM]		
  	D[SYN_MAX_NUM]		
  	
  	eI=-75
  	eE=0
  	gmax1=0.0015 (umho)
  	CL1tau1=0.4 (ms)
	CL1tau2=0.5 (ms)
  	DCL1=0.17946757
  	tau_DCL1=121.34669 (ms)
  	
  	gmax2=0.0008 (umho)	
  	CL2tau1=0.589999 (ms)
  	CL2tau2=0.59 (ms)
  	
  	gmax3=0.0009 (umho)	
  	CL3tau1=0.4 (ms)
  	CL3tau2=0.6 (ms)
  	FCL3=2.6705809
  	tau_FCL3=21.896552 (ms)
  	DCL3=0.0017831198
  	tau_DCL3=11.319578 (ms)
	
  	gmaxI=0.00053 (umho)
  	CLItau1=0.55 (ms)
  	CLItau2=6.5 (ms)
  	FCLI=2.6370474
  	DCLI=0.56560689
  	tau_DCLI=45.241561 (ms)
	tau_FCLI=6.6086312 (ms)
	
	factor1
	factor2
	factor3
	factorI
	time
	CL3f
	k
}

VERBATIM
	double *hoc_pgetarg();
	double *Syntimes[150];  
ENDVERBATIM

STATE {
	A1 (umho)
	B1 (umho)
	A2 (umho)
	B2 (umho)
	A3 (umho)
	B3 (umho)
	AI (umho)
	BI (umho)
}

ASSIGNED {
	v (mV)
	i (nA)
	g1
	g2
	g3
	gi
}

PROCEDURE addspiketrain() {
VERBATIM
	int x=(int)*getarg(1);
	Syntimes[x]=hoc_pgetarg(2);
ENDVERBATIM
}

INITIAL {LOCAL ari,tp1,tp2,tp3,tpI
	A1 =0
	B1 =0
	A2 =0
	B2 =0
	A3 =0
	B3 =0
	AI =0
	BI =0
	g1=0
	g2=0
	g3=0
	gi=0
	
	
	if (CL1tau1/CL1tau2 > 0.9999) {CL1tau1 = 0.9999*CL1tau2}
	if (CL2tau1/CL2tau2 > 0.9999) {CL2tau1 = 0.9999*CL2tau2}
	if (CL3tau1/CL3tau2 > 0.9999) {CL3tau1 = 0.9999*CL3tau2}
	if (CLItau1/CLItau2 > 0.9999) {CLItau1 = 0.9999*CLItau2}
	
	tp1 = (CL1tau1*CL1tau2)/(CL1tau2 - CL1tau1) * log(CL1tau2/CL1tau1)
	factor1 = -exp(-tp1/CL1tau1) + exp(-tp1/CL1tau2)
	factor1 = 1/factor1
	
	tp2 = (CL2tau1*CL2tau2)/(CL2tau2 - CL2tau1) * log(CL2tau2/CL2tau1)
	factor2 = -exp(-tp2/CL2tau1) + exp(-tp2/CL2tau2)
	factor2 = 1/factor2
	
	tp3 = (CL3tau1*CL3tau2)/(CL3tau2 - CL3tau1) * log(CL3tau2/CL3tau1)
	factor3 = -exp(-tp3/CL3tau1) + exp(-tp3/CL3tau2)
	factor3 = 1/factor3
	
	tpI = (CLItau1*CLItau2)/(CLItau2 - CLItau1) * log(CLItau2/CLItau1)
	factorI = -exp(-tpI/CLItau1) + exp(-tpI/CLItau2)
	factorI = 1/factorI
		
	FROM ari=0 TO syncount-1{
	    F[ari]=1
		D[ari]=1
		tsyn[ari]=-1000
		time=0
		
	}
	
	
	
	
	
	
	
	
	
	
    
}

BREAKPOINT{
	SOLVE state METHOD cnexp
	g1=B1 - A1
	g2=B2 - A2
	g3=B3 - A3
	gi=BI - AI
	i=((g1+g2+g3)*(v-eE)) + gi*(v-eI)
	
}

DERIVATIVE state{
	
	A1'=-A1/CL1tau1
	B1'=-B1/CL1tau2
	
	A2'=-A2/CL2tau1
	B2'=-B2/CL2tau2
	
	A3'=-A3/CL3tau1
	B3'=-B3/CL3tau2
	
	AI'=-AI/CLItau1
	BI'=-BI/CLItau2
}

NET_RECEIVE(w){ LOCAL x,j
        k = w
		
		if (syntype[k] == 1){				
			
			
        	D[k] = 1 - (1-D[k])*exp(-(t - tsyn[k])/tau_DCL1)
	        tsyn[k] = t
	        
			state_discontinuity(A1, A1 + gmax1*factor1*D[k])
			state_discontinuity(B1, B1 + gmax1*factor1*D[k])
	        
	        D[k] = D[k] * DCL1
	        
	        
	        
	        
	        
	        
	        
	        
		} else if (syntype[k] == 2){			
			
       		
        	
			state_discontinuity(A2, A2 + gmax2*factor2)
			state_discontinuity(B2, B2 + gmax2*factor2)
			
        	
        	
        	
	    } else if (syntype[k] == 3){			
			F[k] = 1 + (F[k]-1)*exp(-(t - tsyn[k])/tau_FCL3)
       		D[k] = 1 - (1-D[k])*exp(-(t - tsyn[k])/tau_DCL3)
        	tsyn[k] = t
			
			state_discontinuity(A3, A3 + gmax3*factor3*F[k]*D[k]) 
			state_discontinuity(B3, B3 + gmax3*factor3*F[k]*D[k])
        	F[k] = F[k] + FCL3
        	D[k] = D[k] * DCL3
        	
			
	    } else if (syntype[k] == 4){				
			F[k] = 1 + (F[k]-1)*exp(-(t - tsyn[k])/tau_FCLI)
        	D[k] = 1 - (1-D[k])*exp(-(t - tsyn[k])/tau_DCLI)
	       	tsyn[k] = t
			state_discontinuity(AI, AI + gmaxI*factorI*F[k]*D[k])
			state_discontinuity(BI, BI + gmaxI*factorI*F[k]*D[k])
	       	F[k] = F[k] + FCLI
	       	D[k] = D[k] * DCLI
	       	
	       	
	       	
	    }
            
    count[k]=count[k]+1	
}

PROCEDURE add(x){
	if (syncount<SYN_MAX_NUM-1) {
		syntype[syncount]=x
		syncount=syncount+1
	} else {
		printf("\n\n\n EXCEEDING ALLOCATED # OF SYNAPSES\n\n\n")
	}
}