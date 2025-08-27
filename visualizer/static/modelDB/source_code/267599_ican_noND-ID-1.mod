TITLE Slow Ca-dependent cation current
:
:   We've moved to the Model described by Nillus in 2004, while keeping the description of the nanodomain
: modified by Canavier too include separate pool for ICan calcium microdomain
: at this point ican doesn't activate other pools of calcium need to declare a new ion species
: because at this point the code requires pools for SK+BK+T inactivation and ICAN that decay at different rates

INDEPENDENT {
    t FROM 0 TO 1 WITH 1 (ms)
}

NEURON {
    SUFFIX icannoND
    USEION ca READ cai 
    RANGE erev
    RANGE gbar,itrpm4
    RANGE beta,taumin
    NONSPECIFIC_CURRENT itrpm4
}

UNITS {
    (mA)=(milliamp)
    (mV)=(millivolt)
    (molar)=(1/liter)
    (mM)=(millimolar)
    (um)=(micron)
    (msM)=(ms mM)
    FARADAY=(faraday) (coulomb)
}

PARAMETER {
    v (mV)
    erev=0 (mV)             : reversal potential
    cai (mM)                : will now decay to bulk cai
    gbar=0.0001 (mho/cm2)
    : middle point of activation fct, for ip3 as somacar, for current injection
    taumin=0.1 (ms)         : minimal value of time constant
    Kd = 87e-3 (mM)		:20e-3
}

STATE {
    can (mM) 
    Po
}

ASSIGNED {
	jip3p (mM/ms)
    ican (mA/cm2)
    itrpm4 (mA/cm2)
    Po_inf
    Tau (ms)
    :cai (mM) 
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    itrpm4=gbar*Po*(v-erev)
}

DERIVATIVE states {
    evaluate_fct(v,cai)
    Po'=(Po_inf-Po)/Tau
}

FUNCTION MyExp(x) {
    if (x<-50) {MyExp=0}
    else if (x>50) {MyExp=exp(50)}
    else {MyExp=exp(x)}
}

UNITSOFF

INITIAL {
    : activation kinetics are assumed to be at 22 deg. C
    : Q10 is assumed to be 3
    evaluate_fct(v,cai)
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {
    LOCAL alpha, alpha2, beta
    alpha=0.0057*MyExp(0.0060*v)
    beta=0.033*MyExp(-0.019*v)
    alpha2=alpha/(1+(Kd/cai))
    Po_inf=alpha2/(alpha2+beta)
    Tau=1/(alpha2+beta)
    if (Tau<taumin) {Tau=taumin}                        : min value of time cst
}

UNITSON
