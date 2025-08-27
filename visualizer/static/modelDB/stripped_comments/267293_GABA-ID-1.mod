NEURON {

        POINT_PROCESS GABA
        RANGE gGABAmax
	RANGE E_Cl, tau_sGABA
        RANGE i, i_GABA, g_GABA, sGABA
        NONSPECIFIC_CURRENT i, i_GABA

}

PARAMETER {

	gGABAmax = 0.01  (uS)

	E_Cl = -80       (mV)
	tau_sGABA = 2   (ms)

}

ASSIGNED {

        v (mV)
        i (nA)
	i_GABA (nA)
        g_GABA (uS)
}

STATE {

        sGABA       
}

INITIAL{

	sGABA = 0
        
}

BREAKPOINT {

        SOLVE state METHOD cnexp
        g_GABA = gGABAmax*sGABA          
        i_GABA = g_GABA*(v-E_Cl) 
	i = i_GABA
}

DERIVATIVE state{

        sGABA' = -sGABA/tau_sGABA
}



NET_RECEIVE (weight){
	
        sGABA = sGABA + 1
        if (sGABA > 1) { 
          sGABA = 1
        }
}