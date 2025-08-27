NEURON {
    SUFFIX kir
    USEION k WRITE ik
    RANGE gbar,Offset,Slope,ik,ek,g
}

PARAMETER {
    gbar (siemens/cm2)
    ek=-95 (mV)
    Offset=15 (mV)
    Slope=-10 (mV)
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
}

STATE {
    g
}

BREAKPOINT {
    LOCAL Arg
    Arg=-(v-ek+Offset)/Slope

    if (Arg<-50) {g=gbar}
    else if (Arg>50) {g=0}
    else {g=gbar/(1.0+exp(Arg))}

    ik=(v-ek)*g
}