NEURON {
POINT_PROCESS bipNMDA
	RANGE Vpre,Vdel,Vdur,Vamp,Vbase,locx,locy,local_v,i,g,release,numves 
	RANGE gAMPA,gNMDA,s_inf,t1,A,B,Vinf	
	GLOBAL maxves,newves,gAMPAsingle,gNMDAsingle,Vtau,Voff,Vset,VampK	
	GLOBAL icaconst ,gama,n,e,tauAMPA,tau1NMDA ,tau2NMDA
	NONSPECIFIC_CURRENT iAMPA,iNMDA
	
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}
PARAMETER {
							
	maxves=10				
	newves=0.01				
	
	Vdel=50 	(ms)		
	Vdur=100	(ms)		
	Vamp=10		(mV)		
	Vbase=0		(mV)		
	VampK=2					
	
	
	Vtau=30	(/ms)		
	
							
	gAMPAsingle=0.2	(nS)	
	gNMDAsingle=0.2	(nS)	
	tau1NMDA=50	(ms)		
	tau2NMDA=2	(ms)		
	tauAMPA=2	(ms)		
	n=0.25 		(/mM)		
	gama=0.08 	(/mV)		
	e = 0 		(mV)		
	locx=0					
	locy=0					
	icaconst =0.1			
	Voff=0					
	Vset=-60				
}

ASSIGNED {
	
	Vinf 		(mV)
	
	s_inf
	t1	
	numves
	release

	
	v 			(mV)
	i 			(nA)
	g           (nS)
	iNMDA		(nA)
	iAMPA		(nA)
	gNMDA		(nS)
	local_v		(mV)
	
}

STATE {
	A
	B
	gAMPA 		(nS)
	
	
	Vpre		(mV)
}

BREAKPOINT {
	SOLVE state METHOD euler
	if (t>t1){										
		
		
		
		



		releasefunc(Vpre)
		t1=t1+1
	}
	
	
	
	
	
	
	
	

	local_v=v*(1-Voff)+Vset*Voff					
	gNMDA=(A-B)/(1+n*exp(-gama*local_v) )
	iAMPA = (1e-3)*gAMPA * (v - e)
	iNMDA = (1e-3)*gNMDA * (v - e)
	i= iAMPA+iNMDA									
	g=gNMDA+gAMPA									

	
	
}

INITIAL {
	
	s_inf=0
	release=0
	numves=maxves
	t1=0
	Vinf=0	
	Vpre=Vinf
	
	
	

	
	gAMPA=0
	gNMDA=0
	A=0
	B=0
}
 
FUNCTION releasefunc(vpre){
	LOCAL rand,addves
	s_inf=vpre/100
	release=0	
	FROM rand=0 TO numves-1 {			
		if (scop_random()<s_inf){
			release=release+1
		}
	}
	if (release>0){						
		numves=numves-release
		if (numves<0){numves=0}
		state_discontinuity( gAMPA, gAMPA+ release*gAMPAsingle)
		state_discontinuity( A, A+ release*gNMDAsingle)
		state_discontinuity( B, B+ release*gNMDAsingle)
	}
	addves=0							
	FROM rand=0 TO maxves-numves-1 {
		if (scop_random()<newves){addves=addves+1}
	}
	numves=numves+addves
}
DERIVATIVE state {
	A'=-A/tau1NMDA
	B'=-B/tau2NMDA
	gAMPA'=-gAMPA/tauAMPA
	Vpre'=(-Vpre+Vinf)/Vtau
	

}