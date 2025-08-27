INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS modulator
	RANGE  Alpha_Max, Alpha_Delay, Alpha_tau, DC_Level, DC_Delay, DC_Off, Ramp_Max, Ramp_Delay, Ramp_Off, Slope_UP, Slope_DOWN
	GLOBAL mod
}

UNITS {
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
        Alpha_Max=0     (umho)
	Alpha_Delay=0 (ms)
	Alpha_tau=.1 (ms)
	e=0	(mV)
	v	(mV)
	mod=0	(mM)
        DC_Level (mM)
        DC_Delay (ms)
        DC_Off (ms)
        Ramp_Max=0.5 (mM)
        Ramp_Delay=0 (ms)
        Ramp_Off=0 (ms)
        Slope_UP=1 (ms)
        Slope_DOWN=-.0001 (ms)
}


BREAKPOINT {
        mod = mramp(t) + (Alpha_Max * alpha( (t - Alpha_Delay)/Alpha_tau ))
}

FUNCTION alpha(x) {
	if (x < 0 || x > 10) {
		alpha = 0
	}else{
		alpha = x * exp(1 - x)
	}
}



FUNCTION mramp(x)
{
VERBATIM
double timeramp, x, Dc;
x = _lx;

if (x >= DC_Delay & x < DC_Off)
   Dc = DC_Level;
else
   Dc = 0;

if (x < Ramp_Delay)
   timeramp = Dc;
else
  {
  if (x < Ramp_Off)
   {
    timeramp = Dc + (Slope_UP * (x - Ramp_Delay));
    if (timeramp >= Ramp_Max)
       timeramp = Ramp_Max;
   }
  else
   {
    timeramp = Ramp_Max + (Slope_DOWN * (x - Ramp_Off));
    if (timeramp <= Dc) 
       timeramp = Dc;
   }
  }
return (timeramp);
ENDVERBATIM
}