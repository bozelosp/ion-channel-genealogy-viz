NEURON {
POINT_PROCESS SACexc
	RANGE Vpre ,Vinf
	GLOBAL tau ,Vtau ,e,maxves,gsingle ,newves
	RANGE release,numves,g,s_inf,t1,i,g
	RANGE locx,locy,local_v
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}
PARAMETER {
								
	maxves=10				
	newves=0.01				
	Vtau=30	(/ms)		
							
	gsingle=0.2	(nS)
	tau=3		(ms)
	e = 0 	(mV)
	locx=0		
	locy=0		

}

ASSIGNED {
	
	Vinf 		(mV)
	s_inf
	t1
	numves
	release
	
	v 			(mV)
	i 			(nA)
	local_v		(mV)
}

STATE {
	g	 		(nS)
	Vpre 		(mV)
}
 
BREAKPOINT {
	SOLVE state METHOD euler
	if (t>t1){										
		releasefunc(Vpre)
		t1=t1+1
	}
	i = (1e-3)*g * (v - e)
	local_v=v
}

INITIAL {
	
	s_inf=0
	release=0
	numves=maxves
	t1=0
	Vinf=0
	Vpre=Vinf
	
	g =0
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
		state_discontinuity( g, g+ release*gsingle)
	}
	addves=0							
	FROM rand=0 TO maxves-numves-1 {
		if (scop_random()<newves){addves=addves+1}
	}
	numves=numves+addves
}
DERIVATIVE state {
	g'=-g/tau
	Vpre'=(-Vpre+Vinf)/Vtau
}