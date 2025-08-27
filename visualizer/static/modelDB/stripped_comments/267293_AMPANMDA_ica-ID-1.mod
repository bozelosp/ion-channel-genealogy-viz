NEURON {

        POINT_PROCESS AMPANMDA_ica
        RANGE gAMPAmax, gNMDAmax, MgCon
	RANGE E_Glu, tau_sAMPA, tau_sNMDA, tau_xNMDA, alphas
        RANGE i, i_AMPA, i_NMDA, g_AMPA, g_NMDA, sAMPA, sNMDA, xNMDA
        NONSPECIFIC_CURRENT i, i_AMPA,i_NMDA
        USEION ca WRITE ica

}

PARAMETER {

	gAMPAmax = 0.01  (uS)
	gNMDAmax = 0.007 (uS)
	MgCon = 0.69     
	mggate

	E_Glu = 0       (mV)
	tau_sAMPA = 2   (ms)
	tau_sNMDA = 100 (ms)
	tau_xNMDA = 2   (ms)
	alphas = 0.5    (kHz)

}

ASSIGNED {

        v (mV)
        i (nA)
	i_AMPA (nA)
	i_NMDA (nA)
        g_AMPA (uS)
	g_NMDA (uS)
	ica (nA)
}

STATE {

        sAMPA       
	sNMDA       
        xNMDA       
}

INITIAL{

	sAMPA = 0
	sNMDA = 0
	xNMDA = 0
        
}

BREAKPOINT {

        SOLVE state METHOD cnexp
        if (sNMDA > 1) { 
          sNMDA = 1
        }
	mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (MgCon / 3.57 (mM))) 
        g_AMPA = gAMPAmax*sAMPA          
	g_NMDA = gNMDAmax*sNMDA * mggate 
        i_AMPA = g_AMPA*(v-E_Glu) 
	i_NMDA = g_NMDA*(v-E_Glu) 
        i = i_AMPA + i_NMDA
        ica = 0.1*i																		   
}

DERIVATIVE state{

        sAMPA' = -sAMPA/tau_sAMPA
	sNMDA' = -sNMDA/tau_sNMDA + alphas*xNMDA*(1-sNMDA)
        xNMDA' = -xNMDA/tau_xNMDA
}



NET_RECEIVE (weight){
	
        sAMPA = sAMPA + 1
	xNMDA = xNMDA + 1
        if (sAMPA > 1) { 
          sAMPA = 1
        }
}