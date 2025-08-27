NEURON {
    SUFFIX cadp
    USEION ca READ cao, cai, ica WRITE cai, ica
    RANGE ica_pmp, ex_buffer_ratio, f, AvgCaExBuffer, diam_factor, TotalPump, ca, CaEndBuffer, EndBuffer, CaExBuffer, ExBuffer, pump, pumpca
    GLOBAL vrat, TotalEndBuffer, TotalExBuffer, k1, k2, k3, k4, k1bufend, k2bufend, k1bufex, k2bufex, DCa, fl_ratio, dep_factor
    USEION dep READ depi VALENCE 1

    NONSPECIFIC_CURRENT icont                                                       
}

DEFINE Nannuli 4

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
    dep_factor=1 (1)


}

ASSIGNED {
    diam (um)
    L (um)
    ica (mA/cm2)
    cai (mM)
    vrat[Nannuli] 
                 
                 
    
    depi
    k2          (/ms)
    k4          (/mM-ms)

    ka_end (/mM)
    ka_ex (/mM)

    B0end (mM)
    B0ex (mM)

    cao (mM)
    ica_pmp (mA/cm2)
    icont (mA/cm2)
    parea (um)
    AvgCaExBuffer (mM)  

    f    (1)
    diamf   (um)


}

CONSTANT { volo = 1e10 (um2) }

STATE {
    
    
    ca[Nannuli] (mM) <1e-10>
    CaEndBuffer[Nannuli] (mM)
    EndBuffer[Nannuli] (mM)
    CaExBuffer[Nannuli] (mM)
    ExBuffer[Nannuli] (mM)

    pump (mol/cm2)
    pumpca (mol/cm2)

}


BREAKPOINT {


    SOLVE state METHOD sparse

    ica = ica_pmp
    icont = -ica_pmp
    AvgCaExBuffer = 0.0
    
    FROM i=0 TO Nannuli-1 {
        AvgCaExBuffer = AvgCaExBuffer + (CaExBuffer[i] * vrat[i])
    }

    AvgCaExBuffer = AvgCaExBuffer * (4/PI)
    f = TotalExBuffer + (fl_ratio - 1) * AvgCaExBuffer
    
}

LOCAL factors_done

INITIAL {

    k2=sqrt(cai/cao)    
    k4=sqrt(cai/cao)
    diamf=diam*diam_factor

    parea = PI*diamf
    pump = TotalPump/(1 + (cai*k1/k2))
    pumpca = TotalPump - pump
    if (factors_done == 0) {    
        factors_done = 1        
        factors()               
    }

    ka_end = k1bufend/k2bufend
    ka_ex = k1bufex/k2bufex

    B0end = TotalEndBuffer/(1 + ka_end*cai)
    B0ex = TotalExBuffer/(1 + ka_ex*cai)


    FROM i=0 TO Nannuli-1 {
        ca[i] = cai
        EndBuffer[i] = B0end
        CaEndBuffer[i] = TotalEndBuffer - B0end
        ExBuffer[i] = B0ex
        CaExBuffer[i] = TotalExBuffer - B0ex
        

    }
    


}

LOCAL frat[Nannuli]     

PROCEDURE factors() {
    LOCAL r, dr2
    r = 1/2                 
    dr2 = r/(Nannuli-1)/2   
                            
    vrat[0] = 0
    frat[0] = 2*r
    FROM i=0 TO Nannuli-2 {
        vrat[i] = vrat[i] + PI*(r-dr2/2)*2*dr2  
        r = r - dr2
        frat[i+1] = 2*PI*r/(2*dr2)              
                                                
        r = r - dr2
        vrat[i+1] = PI*(r+dr2/2)*2*dr2 
    }
}

LOCAL dsq, dsqvol   
                    

KINETIC state {
    COMPARTMENT i, diamf*diamf*vrat[i] {ca CaIndBuffer IndBuffer CaExBuffer ExBuffer}
    COMPARTMENT (1e10)*parea {pump pumpca}
    COMPARTMENT volo {cao}

    

    
    ~ ca[0] + pump <-> pumpca (k1*(1-depi*dep_factor)*parea*(1e10), k2*(1-depi*dep_factor)*parea*(1e10))
    ~ pumpca <-> pump + cao (k3*(1-depi*dep_factor)*parea*(1e10), k4*(1-depi*dep_factor)*parea*(1e10))

    CONSERVE pump + pumpca = TotalPump * parea * (1e10)
    ica_pmp = 2*FARADAY*(f_flux - b_flux)/parea

    
    ~ ca[0] << (-(ica-ica_pmp)*PI*diamf/(2*FARADAY)) 

    FROM i=0 TO Nannuli-2 {
        ~ ca[i] <-> ca[i+1] (DCa*frat[i+1], DCa*frat[i+1])
    }
    dsq = diamf*diamf
    FROM i=0 TO Nannuli-1 {
        dsqvol = dsq*vrat[i]
        ~ ca[i] + EndBuffer[i] <-> CaEndBuffer[i] (k1bufend*dsqvol, k2bufend*dsqvol)
        ~ ca[i] + ExBuffer[i] <-> CaExBuffer[i] (k1bufex*dsqvol, k2bufex*dsqvol)

    }
    cai = ca[0]
}