NEURON {
    THREADSAFE
    SUFFIX nadp
    USEION na READ nao, nai, ina WRITE nai, ina, ena
    NONSPECIFIC_CURRENT ik_pump
    RANGE ina_pmp, TotalPump, ik_ratio, na, pump, pumpna, k4_coeff, DNa_coeff, initial_na
    GLOBAL vrat, DNa, k1, k2, k3, k4, change_ena, fix_na
}



UNITS {
    (molar) = (1/liter)
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    (mV) = (millivolt)
    FARADAY = (faraday) (10000 coulomb)
    R = (k-mole) (joule/degC)    
    PI = (pi) (1)
    (mol) = (1)

}

PARAMETER {

    DNa = 0.6 (um2/ms)

    k1 = 1.0 (/mM3-ms)

    k2 = 0.001 (/ms)
    
    k3 = 0.3 (/ms)
                            
    TotalPump = 1e-14 (mol/cm2)
    
    ik_ratio = -0.66666666 (1)

    change_ena = 1 (1)

    k4_coeff = 1.0 (1)
    DNa_coeff = 1.0 (1)

    fix_na = 0 (1)

}

ASSIGNED {
    diam (um)
    L (um)
    ina (mA/cm2)
    nai (mM)

                 
                 
    
    k4          (/mM3-ms)

    nao (mM)
    ena (mV)

    ina_pmp (mA/cm2)
    parea (um)

    ik_pump (mA/cm2)
	celsius

    initial_na (mM)

    



}

CONSTANT { volo = 1e10 (um2) }

STATE {
    
    na (mM) <1e-3>

    pump (mol/cm2)
    pumpna (mol/cm2)

}


BREAKPOINT {


    SOLVE state METHOD sparse

    ina = ina_pmp
    ik_pump = ik_ratio*ina_pmp

}

LOCAL factors_done

INITIAL {
	MUTEXLOCK
    k4=(((nai/nao)^3)*k1*k3)/k2    
    parea = PI*diam
    vrat = PI*0.25
    pump = TotalPump/(1 + (nai*k1/k2))
    pumpna = TotalPump - pump





    na = nai
    initial_na = nai




	MUTEXUNLOCK
}







                            






                                                






KINETIC state {
    COMPARTMENT diam*diam*vrat {na}
    COMPARTMENT (1e10)*parea {pump pumpna}
    COMPARTMENT volo {nao}

    LONGITUDINAL_DIFFUSION DNa*DNa_coeff*diam*diam*vrat {na}

    
    ~ 3 na + pump <-> pumpna (k1*parea*(1e10), k2*parea*(1e10))
    ~ pumpna <-> pump + 3 nao (k3*parea*(1e10), k4*k4_coeff*parea*(1e10))

    CONSERVE pump + pumpna = TotalPump * parea * (1e10)
    ina_pmp = FARADAY*(f_flux - b_flux)/parea

    
    ~ na << (-(ina-ina_pmp)*PI*diam/(FARADAY)) 

    
    
    
    if (fix_na) {
        na = initial_na
    }
    nai = na
    if (change_ena == 1) {
        ena = ((R*(273.15+celsius))/(FARADAY*10))*log(nao/nai)
    }
    else {
        ena = 60.
    }
}