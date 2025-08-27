TITLE Slow Ca-dependent cation current
:
:   Ca++ dependent nonspecific cation current ICAN
:   Differential equations
:
:   Model based on a first order kinetic scheme
:
:       + n cai <->     (alpha,beta)
:
:   Following this model, the activation fct will be half-activated at 
:   a concentration of Cai = (beta/alpha)^(1/n) = cac (parameter)
:
:   The mod file is here written for the case n=2 (2 binding sites)
:   ---------------------------------------------
:
:   Kinetics based on: Partridge & Swandulla, TINS 11: 69-72, 1988.
:
:   This current has the following properties:
:      - inward current (non specific for cations Na, K, Ca, ...)
:      - activated by intracellular calcium
:      - NOT voltage dependent
:
:   A minimal value for the time constant has been added
:
:   Ref: Destexhe et al., J. Neurophysiology 72: 803-818, 1994.
:   See also:  http://www.cnl.salk.edu/~alain , http://cns.fmed.ulaval.ca
:
: Updated by Kiki Sidiropoulou (2010) so that dADP has slow inactivation kinetics and it is activated after 5 spikes
: modified by Canavier too include separate pool for ICan calcium microdomain
: at this point ican doesn't activate other pools of calcium need to declare a new ion species
: because at this point the code requires pools for SK+BK+T inactivation and ICAN that decay at different rates

INDEPENDENT {
    t FROM 0 TO 1 WITH 1 (ms)
}

NEURON {
    SUFFIX ican
    USEION ca READ cai 
    USEION na WRITE ina
    RANGE depth,taur,frac,ica,ican
    RANGE gbar,m_inf,tau_m,in,mystart
    RANGE beta,cac,taumin
    POINTER cip3p
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
    depth=0.08 (um)         : depth of shell 0.00027
    frac=0
    taur=100 (ms)           : rate of calcium removal
    :cainf=100e-6 (mM)
    en=-20 (mV)             : reversal potential
    eca=140 (mV)            : reversal potential
    cai (mM)                : will now decay to bulk cai
    gbar=0.0001 (mho/cm2)
    beta=0.0001 (1/ms)      : backward rate constant
    cac=0.15 (mM)	:0.02
    : middle point of activation fct, for ip3 as somacar, for current injection
    taumin=0.1 (ms)         : minimal value of time constant
}

STATE {
    can (mM) 
    m
}

ASSIGNED {
	cip3p (mA/cm2)
    ica (mA/cm2)
    ican (mA/cm2)
    drive_channel (mM/ms)
    in (mA/cm2)
    ina (mA/cm2)
    m_inf
    tau_m (ms)
    tadj
    :cai (mM) 
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    in=gbar*m*m*(v-en)
    ican=gbar*frac*m*m*(v-eca)
    ina=0.7*in
}

DERIVATIVE states {
    evaluate_fct(v,can)
    m'=(m_inf-m)/tau_m
    drive_channel=-(10000)*(ican-100*cip3p)/(2*FARADAY*depth) :
    if (drive_channel<=0.0) {drive_channel=0.0}             : cannot pump inward 
    can'=drive_channel/18+(cai-can)/(taur)
}

UNITSOFF

INITIAL {
    : activation kinetics are assumed to be at 22 deg. C
    : Q10 is assumed to be 3
    can=cai
    tadj=3.0^((celsius-22.0)/10)
    evaluate_fct(v,can)
    m=m_inf
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {
    LOCAL alpha2
    alpha2=beta*(cai/cac)^2
    if (cai<1.09e-4) {alpha2=0}
    tau_m=1/(alpha2+beta)/tadj
    m_inf=alpha2/(alpha2+beta)
    if (tau_m<taumin) {tau_m=taumin}                        : min value of time cst
}

UNITSON
