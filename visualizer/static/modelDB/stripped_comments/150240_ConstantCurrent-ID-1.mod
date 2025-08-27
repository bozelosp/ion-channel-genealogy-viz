NEURON {
	POINT_PROCESS ConstantCurrent
	RANGE Enabled
	RANGE Amplitude
	ELECTRODE_CURRENT I
}

UNITS
{
}

PARAMETER
{
	Enabled = 0
	Amplitude = 0 (nA)
}

ASSIGNED
{
	I (nA)
}

INITIAL
{
}

PROCEDURE UpdateParameters()
{
}

BREAKPOINT
{
	if (Enabled == 1)
	{
		I = Amplitude
	}
	else
	{
		I = 0
	}
}