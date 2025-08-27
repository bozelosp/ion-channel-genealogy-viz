:---------------------------------------------------------------
:  Translation of the cooperatve Na+ current model 
:   used for Fig 5 in the publication 
:   Unique features of action potential initiation in cortical neurons
:   Bjoern Naundorf, Fred Wolf, Maxim Volgushev
:   Nature 440:1060-3 (2006)
:
:   K+ currents are added to their model 
:   Parameters optimized to fit the recorded APs of immature hippocampal CA3 neurons 
:   last changed: 27.07.2021
:   Implemented by W. Kilb
:
:----------------------------------------------------------------------

UNITS {
 	(mA) = (milliamp)
 	(mV) = (millivolt)
	(S) = (siemens)		
}

NEURON {
	SUFFIX ndrfAP
	NONSPECIFIC_CURRENT i
	RANGE iNa,iK
	RANGE eNa, eK
	RANGE gNaMax, gKaMax
	RANGE AlphaActNa, AlphaInaNa, BetaActNa, BetaInaNa
	RANGE AlphaActKa, BetaInaKa 
        RANGE CoopNa, CoopKa  
        RANGE TauNa, TauKa 
        RANGE TauActNa, VActNa, kActNa 
        RANGE TauInaNa, VInaNa, kInaNa
        RANGE TauActKa, VActKa, kActKa 
        RANGE TauInaKa, VInaKa, kInaKa
 }


	
PARAMETER {
        gNaMax = 0.006 (S/cm2)	
        gKaMax = 0.0013 (S/cm2) 	
        eK = -90 (mV)
        eNa = 66 (mV)
        TauActNa = 0.51
        VActNa = -32.8
        kActNa = 3
        TauInaNa = 1.1
        VInaNa = -35.7
        kInaNa = 5
        TauNa = 999
        CoopNa = 14
        TauActKa = 1.5
        VActKa = -30
        kActKa = 4
        TauInaKa = 1.66 
        VInaKa = -10
        kInaKa = 1
        TauKa = 500
        CoopKa = 0
}

 
STATE {
      NaO
      NaH
      KaO        
      KaH
}
 
ASSIGNED {
        v (mV)
        i (mA/cm2)
        cm (uF)
        celsius
        iNa iK(mA/cm2)
        gNa gK (S/cm2)

        : rate Variables for the transition between channel states
        AlphaActNa 
        AlphaInaNa
        BetaActNa
        BetaInaNa
        AlphaActKa 
        BetaInaKa

        Q10                 : Temperatur Factor



}

INITIAL {
	Q10 = 1
:3^((celsius-31(degC))/10 (degC)) 
      NaO = 0
      NaH = 0
      KaO = 0       
      KaH = 0
      if (TauNa == 0) {TauNa = 10e-10}
      if (TauKa == 0) {TauKa = 10e-10}
}

BREAKPOINT {
        SOLVE states METHOD cnexp 
         
                : ---- calculate conductance and currents for Na+ and K+
        gNa = gNaMax * NaO
	iNa = gNa*(v-eNa)
		
        gK = gKaMax * KaO 
	iK = gK*(v-eK)
	
        i = iK + iNa      
}
 

DERIVATIVE states {  
        : Just for debugging  - can be deleted in the final version
	AlphaActNa = Alpha((v+CoopNa*NaO),TauActNa,VActNa,kActNa)
        BetaActNa = Beta((v+CoopNa*NaO),TauActNa,VActNa,kActNa)
        AlphaInaNa = Alpha(v,TauInaNa,VInaNa,kInaNa)
        BetaInaNa = Beta(v,TauInaNa,VInaNa,kInaNa)
	AlphaActKa = Alpha((v+CoopKa*KaO),TauActKa,VActKa,kActKa)
        BetaInaKa = Alpha(v,TauInaKa,VInaKa,kInaKa)


        NaO' = Alpha((v+CoopNa*NaO),TauActNa,VActNa,kActNa)*(1-NaH-NaO)-Beta((v+CoopNa*NaO),TauActNa,VActNa,kActNa)*NaO-NaO/TauNa
        NaH' = Alpha(v,TauInaNa,VInaNa,kInaNa)*(1-NaH)-Beta(v,TauInaNa,VInaNa,kInaNa)*NaH-NaH/TauNa
        KaO' = Alpha((v+CoopKa*KaO),TauActKa,VActKa,kActKa)*(1-KaO)-Beta(v,TauInaKa,VInaKa,kInaKa)*KaO-KaO/TauKa
}

: Definition of the alpha and beta function
: Inputs V   = voltage
:        Tau = tau value (resp. 1/gain value)
:        Vh  = v(1/2) halfmaximum activation voltage
:        k   = slope of the sigmoidal curve
:        note that alpha and beta functions just differes in the direction of the voltage terms
FUNCTION Alpha(V,Tau,Vh,k) {
     Alpha = Q10/Tau*(1/(1+exp(Vh-V)/k))
}

FUNCTION Beta(V,Tau,Vh,k) {
     Beta = Q10/Tau*(1/(1+exp(V-Vh)/k))
}

