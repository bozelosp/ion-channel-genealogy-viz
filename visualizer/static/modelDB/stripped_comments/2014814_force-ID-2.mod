NEURON {
    POINT_PROCESS force
    POINTER aPointer, xmPointer

    RANGE a0, b0, c0, d0
    RANGE p0, g1, g2, g3
    RANGE Kse, A, Fc, F
    RANGE xm_init, xce_init, xce
    RANGE spk_index, t_axon
}

PARAMETER {
    a0 = 2.35        
    b0 = 24.35        
    c0 = -7.4        
    d0 = 30.3        
    p0 = 23            
    g1 = -8            
    g2 = 21.4        
    xm_init = -8    
    xce_init = -8    
    Kse = 0.4        
}

STATE {
    xce
}

ASSIGNED {
    aPointer
    xmPointer
    F
    Fc
}

BREAKPOINT { LOCAL i
    SOLVE state METHOD cnexp
    F = p0*Kse*xse(xmPointer, xce)
}

DERIVATIVE state {
    Fc = p0*g(xmPointer)*aPointer
    xce' = dxdt (F, Fc)
}

FUNCTION xse (x, y) { LOCAL d_xm, d_xce, d_se
    d_xm = xmPointer - xm_init
    d_xce = xce - xce_init
    d_se = d_xm - d_xce

    if (d_se <= 0) {xse = 0}
    else {xse = d_se}
}

FUNCTION g (x) {
    g = exp(-((x-g1)/g2)^2)
}

FUNCTION dxdt (x, xc) {LOCAL gain_length
    if (x <= xc) {
        dxdt = (10^-3)*(-b0*(xc-x))/(x+a0*xc/p0)
    } else {
        gain_length = (-d0*(xc-x))/(2*xc-x+c0*xc/p0)
        if (gain_length <= 0) {dxdt = (10^-3)*1e5}
        else {dxdt = (10^-3)*gain_length}
    }
}

INITIAL {
    xce = xce_init
    F = 1e-5
}