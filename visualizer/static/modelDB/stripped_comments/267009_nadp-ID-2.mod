NEURON {
    SUFFIX nadp
    USEION na READ nai, ina WRITE nai, ina, nao
    NONSPECIFIC_CURRENT ik_pump
    RANGE ina_pmp, TotalPump, ik_ratio, na, pump, pumpna, na3,DNa, k1, k2, k3, k4, pcu,nai_i
    GLOBAL vrat                                                        
}

DEFINE Nannuli 4

UNITS {
    (mol) = (1)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    (mV) = (millivolt)
    FARADAY = (faraday) (10000 coulomb)
    R = (k-mole) (joule/degC)    
    PI = (pi) (1)

}

PARAMETER {
	
	nai_i = 10
	nao_i = 140
    DNa = 0.6 (um2/ms)
    k1 = 2.0 (/mM3-ms)

    k2 = 0.001 (/ms)
    
    k3 = 0.6 (/ms)
                            
    TotalPump = 1e-14 (mol/cm2)
    
    ik_ratio = -0.66666666 (1)

}

ASSIGNED {
    diam (um)
    L (um)
    ina (mA/cm2)
    nai (mM)
    vrat[Nannuli] 
                 
                 
    
    k4          (/mM3-ms)

    ina_pmp (mA/cm2)
    parea (um)

    ik_pump (mA/cm2)
	pcu (mA/cm2)
}

CONSTANT { volo = 1e14 (um2) }

STATE {
    na[Nannuli] (mM) <1e-3>
    nao (mM)
    pump (mol/cm2)
    pumpna (mol/cm2)
}

BREAKPOINT {

    SOLVE state METHOD sparse

    ina = ina_pmp
    ik_pump = ik_ratio*ina_pmp
    pcu = ik_pump + ina_pmp
}

LOCAL factors_done

INITIAL {
	nai = nai_i
	nao = nao_i
    k4=(((10/nao)^3)*k1*k3)/k2    
    parea = PI*diam
    pump = TotalPump/(1 + (nai*k1/k2))
    pumpna = TotalPump - pump
    if (factors_done == 0) {    
        factors_done = 1        
        factors()               
    }

    FROM i=0 TO Nannuli-1 {
        na[i] = nai
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


KINETIC state {
    COMPARTMENT i, diam*diam*vrat[i] {na}
    COMPARTMENT (1e10)*parea {pump pumpna}
    COMPARTMENT volo {nao}

    LONGITUDINAL_DIFFUSION i, DNa*diam*diam*vrat[i] {na}

    
    ~ 3 na[0] + pump <-> pumpna (k1*parea*(1e10), k2*parea*(1e10))
    ~ pumpna <-> pump + 3 nao (k3*parea*(1e10), k4*parea*(1e10))
    ina_pmp = FARADAY*(f_flux - b_flux)/parea
    CONSERVE pump + pumpna = TotalPump * parea * (1e10)

    
    ~ na[0] << (-(ina-ina_pmp)*PI*diam/(FARADAY)) 
	
    FROM i=0 TO Nannuli-2 {
        ~ na[i] <-> na[i+1] (DNa*frat[i+1], DNa*frat[i+1])
    }

    nai = na[0]
}