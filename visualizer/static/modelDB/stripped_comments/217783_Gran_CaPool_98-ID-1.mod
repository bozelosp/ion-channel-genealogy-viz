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
    SUFFIX Gran_CaPool_98
    USEION ca READ ica WRITE cai VALENCE 2
    
    RANGE cai
    
    RANGE rest_conc
    
    
    RANGE tau
    
    RANGE thickness, F


    RANGE total_current
    RANGE volume_pool
    
    
    
}

ASSIGNED {

    ica (mA/cm2)
    diam (um)
    area (um)
}

INITIAL {
    
        
    LOCAL pi, shell_inner_diam, cylinderLen, circumference, circumference_shell, volumeOuter, volumeInner, volumeSph, volumeCyl

    pi = 3.14159265

    shell_inner_diam = diam - (2*thickness)


    ?  Volume of the pool if it is a shell inside a sphere of diameter diam

    volumeSph = (diam*diam*diam) * pi / 6 - (shell_inner_diam*shell_inner_diam*shell_inner_diam)* pi / 6


    ? Volume of the pool if it is a cylinder

    circumference = diam * pi
    circumference_shell = shell_inner_diam * pi

    cylinderLen = area/circumference

    volumeOuter = (diam * diam/4) * pi * cylinderLen
    volumeInner = (shell_inner_diam * shell_inner_diam/4) * pi * cylinderLen
    volumeCyl = volumeOuter - volumeInner

    if ((area - (pi * diam * diam)) < 1e-3 && (area - (pi * diam * diam)) > -1e-3 ) {

        ? Assume the segment is a sphere
        

        volume_pool = volumeSph
        
    } else {

        ? assume segment is a cylinder
        

        volume_pool = volumeCyl
    }

    cai = rest_conc

}

PARAMETER {

    total_current
    rest_conc = 0.0000755 (mM)
          
    
    tau = 10 (ms)
   
    F = 96494 (C)
    
    thickness = 0.08399999999999999 (um)   
                
    volume_pool
    
    
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
    surf_area_cm2 = area * 1e-8
    volume_cm3 = volume_pool * 1e-12
    
    total_current = ica * surf_area_cm2


    cai' =  ((-1 * total_current)/(2 * F * volume_cm3)) - ((cai - rest_conc)/tau)
    

}