UNITS { }

NEURON {
	SUFFIX CaSP
	
	
	RANGE k1, k2, k3, k4, k5, k6, k, k5i, k6i 
	RANGE Umax, Rmax, t1, t2, R, vth
	RANGE phi0, phi1, phi2, phi3, phi4

	
	RANGE c1, c2, c3, c4, c5 
	RANGE AMinf, AMtau, SF_AM
	RANGE acm, alpha, alpha1, alpha2, alpha3, beta, gamma

	
	RANGE t_axon 
	USEION mg WRITE mgi VALENCE 2
	USEION cl READ cli
}

PARAMETER {
	
	k1 = 3000		
	k2 = 3			
	k3 = 400		
	k4 = 1			
	k5i = 4e5		
	k6i = 150		
	k = 850			
	SF_AM = 5
	Rmax = 10		
	Umax = 2000		
	t1 = 3			
	t2 = 25			
	phi1 = 0.03
	phi2 = 1.23
	phi3 = 0.01		
	phi4 = 1.08		
	CS0 = 0.03     	
	B0 = 0.00043	
	T0 = 0.00007 	
	
	
	c1 = 0.128 		 
	c2 = 0.093		 
	c3 = 61.206	 
	c4 = -13.116	 
	c5 = 5.095		
	alpha = 2
	alpha1 = 4.77
	alpha2 = 400
	alpha3 = 160
	beta = 0.47
	gamma = 0.001

	
	vth = -40		
	t_axon = 0.01	
}

STATE {
	CaSR			
	CaSRCS			
	Ca				
	CaB				
	CaT				
	AM				
	mgi
	R_On
	spk_index
	Spike_On				
}

ASSIGNED {
	dt
	v 	(mV)
	R
	t_shift 		
	k5
	k6
	AMinf
	AMtau
	cli				
	spk[1000] 		
	xm[2] 			
	vm			    
	acm	
}

BREAKPOINT { LOCAL i, temp_R

	SPK_DETECT (v, t) 
	CaR (CaSR, t)

	SOLVE state METHOD cnexp
	
	xm[0]=xm[1]
	xm[1]=cli
	
	vm = (xm[1]-xm[0])/(dt*10^-3)
	
	
	mgi = AM^alpha

	
	VERBATIM
	FILE *outfile_cli2;
	outfile_cli2=fopen("cli_output_results.txt","a");
	fprintf(outfile_cli2,"%g \n", cli);
	fclose(outfile_cli2);
	ENDVERBATIM
	VERBATIM
	FILE *outfile_mgi2;
	outfile_mgi2=fopen("mgi_output_results.txt","a");
	fprintf(outfile_mgi2,"%g \n", mgi);
	fclose(outfile_mgi2);
	ENDVERBATIM
}

DERIVATIVE state {
	rate (cli, CaT, AM, t)
	
	CaSR' = -k1*CS0*CaSR + (k1*CaSR+k2)*CaSRCS - R + U(Ca)
	CaSRCS' = k1*CS0*CaSR - (k1*CaSR+k2)*CaSRCS
	
	Ca' = - k5*T0*Ca + (k5*Ca+k6)*CaT - k3*B0*Ca + (k3*Ca+k4)*CaB + R - U(Ca)
	CaB' = k3*B0*Ca - (k3*Ca+k4)*CaB
	CaT' = k5*T0*Ca - (k5*Ca+k6)*CaT
	
	AM' = (AMinf -AM)/AMtau
	mgi' = 0
}

PROCEDURE SPK_DETECT (v (mv), t (ms)) {LOCAL spk_output_value
	if (Spike_On == 0 && v > vth) {
	Spike_On = 1
	spk[spk_index] = t + t_axon

	
	spk_output_value = spk[spk_index]
	VERBATIM
	FILE *outfile_spk;
	outfile_spk=fopen("spk_output.txt","a");
	fprintf(outfile_spk,"%g \n", _lspk_output_value);
	fclose(outfile_spk);
	ENDVERBATIM

	spk_index = spk_index + 1

	
	VERBATIM
	FILE *outfile_spk_index;
	outfile_spk_index=fopen("spk_index_output.txt","w");
	fprintf(outfile_spk_index,"%g", spk_index);
	fclose(outfile_spk_index);
	ENDVERBATIM

	R_On = 1
	} else if (v < vth) {
	Spike_On = 0
	}
}

FUNCTION U (x) {
	if (x >= 0) {U = Umax*(x^2*k^2/(1+x*k+x^2*k^2))^2}
	else {U = 0}
}

FUNCTION phi (x) {
	if (x <= -8) {phi = phi1*x + phi2}
	else {phi = phi3*x + phi4}
}

PROCEDURE CaR (CaSR (M), t (ms)) { LOCAL i, temp_R  
	if (R_On == 1) {
		FROM i=0 TO spk_index-1 {
			temp_R = temp_R + CaSR*Rmax*(1-exp(-(t-spk[i])/t1))*exp(-(t-spk[i])/t2)
		}
		R = temp_R
		temp_R = 0
	}
	else {R = 0}
}

PROCEDURE rate (cli (M), CaT (M), AM (M), t(ms)) {
	k5 = phi(cli)*k5i
	k6 = k6i/(1 + SF_AM*AM)
	AMinf = 0.5*(1+tanh(((CaT/T0)-c1)/c2))
	AMtau = c3/(cosh(((CaT/T0)-c4)/(2*c5)))
}

INITIAL {LOCAL i, j
	CaSR = 0.0025  		
	CaSRCS = 0			
	Ca = 1e-10			
	CaB = 0				
	CaT = 0				
	AM = 0				
	mgi = 0
	R_On = 0
	
	
	
	xm[0] = 0
	
	VERBATIM
	float xm_value;
	FILE *infile_xm;
	infile_xm=fopen("xm_tracker.txt","r+");
	fscanf(infile_xm,"%g",&xm_value);
	fclose(infile_xm);
	xm[1] = xm_value;
	ENDVERBATIM

	
	
	
	VERBATIM
	float spk_index_value;
	FILE *infile_spk_index;
	infile_spk_index=fopen("spk_index_output.txt","r+");
	fscanf(infile_spk_index,"%g",&spk_index_value);
	fclose(infile_spk_index);
	ENDVERBATIM

	if (spk_index_value > 0) {
		VERBATIM
		float spk_input_value;
		int j;
		FILE *infile_spk;
		infile_spk=fopen("spk_output.txt","r+");
		for (j = 0 ; j <= spk_index_value - 1 ; j ++) {
			fscanf(infile_spk,"%g",&spk_input_value);
			spk[j] = spk_input_value;
		}
		fclose(infile_spk);
		ENDVERBATIM

		FROM i = spk_index_value TO 999 {
			spk[i] = 0
		}
	}
	else {
		FROM i = 0 TO 999 {
			spk[i] = 0
		}
	}
}