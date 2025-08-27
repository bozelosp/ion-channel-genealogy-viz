NEURON {
	THREADSAFE
    SUFFIX cadp
    USEION ca READ cao, cai, ica WRITE cai, ica
    RANGE ica_pmp, diam_factor, TotalPump, ca, CaEndBuffer, EndBuffer, pump, pumpca
    GLOBAL vrat, TotalEndBuffer, k1, k2, k3, k4, k1bufend, k2bufend


    NONSPECIFIC_CURRENT icont                                                       
}



UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    FARADAY = (faraday) (10000 coulomb)
    PI = (pi) (1)
    (mol) = (1)

}

PARAMETER {

    DCa = 0.6 (um2/ms)
    k1bufend = 100 (/mM-ms) 
    k2bufend = 0.1 (/ms)
    TotalEndBuffer = 0.003 (mM)

    k1bufex = 100 (/mM-ms) 
    k2bufex = 0.017 (/ms)       
    TotalExBuffer = 0.04 (mM)   

    k1 = 1 (/mM-ms)
    
    k3 = 1 (/ms)
                            
    TotalPump = 1e-14 (mol/cm2)

    fl_ratio=14 (1)

    diam_factor=1 (1)



}

ASSIGNED {
    diam (um)
    L (um)
    ica (mA/cm2)
    cai (mM)

                 
                 
    

    k2          (/ms)
    k4          (/mM-ms)

    ka_end (/mM)


    B0end (mM)


    cao (mM)
    ica_pmp (mA/cm2)
    icont (mA/cm2)
    parea (um)


    
    diamf   (um)


}

CONSTANT { volo = 1e10 (um2) }

STATE {
    
    
    
    ca (mM) <1e-10>
    CaEndBuffer (mM)
    EndBuffer (mM)
    
    

    pump (mol/cm2)
    pumpca (mol/cm2)

}


BREAKPOINT {


    SOLVE state METHOD sparse

    ica = ica_pmp
    icont = -ica_pmp
    
    
    
    
    

    
    
    
}

LOCAL factors_done

INITIAL {
    k2=sqrt(cai/cao)    
    k4=sqrt(cai/cao)
    diamf=diam*diam_factor

    parea = PI*diamf
    vrat = PI*0.25
    pump = TotalPump/(1 + (cai*k1/k2))
    pumpca = TotalPump - pump





    ka_end = k1bufend/k2bufend
    

    B0end = TotalEndBuffer/(1 + ka_end*cai)
    


    EndBuffer = B0end
    CaEndBuffer = TotalEndBuffer - B0end
    ca=cai









    

}




















LOCAL dsq, dsqvol   
                    

KINETIC state {
    COMPARTMENT diamf*diamf*vrat {ca CaIndBuffer IndBuffer}
    COMPARTMENT (1e10)*parea {pump pumpca}
    COMPARTMENT volo {cao}

    

    
    ~ ca + pump <-> pumpca (k1*parea*(1e10), k2*parea*(1e10))
    ~ pumpca <-> pump + cao (k3*parea*(1e10), k4*parea*(1e10))

    CONSERVE pump + pumpca = TotalPump * parea * (1e10)
    ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

    
    ~ ca << (-(ica-ica_pmp)*PI*diamf/(2*FARADAY)) 

    
    
    
    dsqvol = diamf*diamf*vrat
    ~ ca + EndBuffer <-> CaEndBuffer (k1bufend*dsqvol, k2bufend*dsqvol)
    cai=ca







}