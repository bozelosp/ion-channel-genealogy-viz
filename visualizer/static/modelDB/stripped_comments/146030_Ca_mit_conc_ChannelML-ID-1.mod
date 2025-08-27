?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file



    

? Creating ion concentration






UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (um) = (micrometer)
    (l) = (liter)
    (molar) = (1/liter)
    (mM) = (millimolar)
}

    
NEURON {
    SUFFIX Ca_mit_conc_ChannelML
    USEION ca READ ica WRITE cai VALENCE 2
    
    RANGE cai
    
    RANGE rest_conc
    
    
    RANGE tau
    
    
    GLOBAL total_current
    
    
    RANGE thickness, F
    
    GLOBAL volume, surf_area
    
    
}

ASSIGNED {

    ica (mA/cm2)
    diam (um)
}

INITIAL {
    
        
    LOCAL shell_inner_diam

    shell_inner_diam = diam - (2*thickness)
    
    volume = (diam*diam*diam)*3.14159/6 - (shell_inner_diam*shell_inner_diam*shell_inner_diam)*3.14159/6
    
    surf_area = (diam*diam)*3.14159
    
    cai = rest_conc

}

PARAMETER {

    total_current
    rest_conc = 0.0000052 (mM)
          
    
    tau = 10 (ms)
   
    F = 96494 (C)
    
    thickness = 10 (um)   
                
    volume
    surf_area
    
    
}

STATE {

    cai (mM)

}

BREAKPOINT {

    SOLVE conc METHOD derivimplicit
    

}

DERIVATIVE conc {
    
    LOCAL thickness_cm, surf_area_cm2, volume_cm3 ? Note, normally dimensions are in um, but curr dens is in mA/cm2, etc
    
    thickness_cm = thickness *(1e-4)
    surf_area_cm2 = surf_area * 1e-8
    volume_cm3 = volume * 1e-12
    
    total_current = ica * surf_area_cm2


    cai' =  ((-1 * total_current)/(2 * F * volume_cm3)) - ((cai - rest_conc)/tau)
    

}