NEURON {
	POINT_PROCESS Syn4P
	
	RANGE tau_a, tau_b, e, i

	NONSPECIFIC_CURRENT i
	
	RANGE A_LTD_pre, A_LTP_pre, A_LTD_post, A_LTP_post
	RANGE tau_G_a, tau_G_b, m_G
	RANGE w_pre, w_post, w_pre_init, w_post_init, w
	RANGE s_ampa, s_nmda
	
	RANGE tau_u_T, theta_u_T, m_T
	RANGE theta_u_N, tau_Z_a, tau_Z_b, m_Z, tau_N_alpha, tau_N_beta, m_N_alpha, m_N_beta, theta_N_X
	RANGE theta_u_C, theta_C_minus, theta_C_plus, tau_K_alpha, tau_K_gamma, m_K_alpha, m_K_beta, s_K_beta
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	
	tau_a = 0.2 (ms) <1e-9,1e9>			
	tau_b = 2 (ms) <1e-9,1e9>			
	e = 0 (mV)							
	
	w_pre_init = 0.5					
	w_post_init = 2.0					
	
	s_ampa = 0.5						
	s_nmda = 0.5						

	
	tau_G_a = 2 (ms) <1e-9,1e9>			
	tau_G_b = 50 (ms) <1e-9,1e9>		
	m_G = 10							
	
	A_LTD_pre = 8.5e-7					
	A_LTP_pre = 8.5e-7					
	A_LTD_post = 3.6e-7					
	A_LTP_post = 5.5e-5					
	
	tau_u_T = 10 (ms) <1e-9,1e9>		
	theta_u_T = -60						
	m_T = 1.7							
	
	theta_u_N = -30						
	tau_Z_a = 1	(ms) <1e-9,1e9>			
	tau_Z_b = 15 (ms) <1e-9,1e9>		
	m_Z = 6								
	tau_N_alpha = 7.5 (ms) <1e-9,1e9>	
	tau_N_beta = 30	(ms) <1e-9,1e9>		
	m_N_alpha = 2						
	m_N_beta = 10						
	theta_N_X = 0.2						
	
	theta_u_C = -68						
	theta_C_minus = 15					
	theta_C_plus = 35					
	tau_K_alpha = 15 (ms) <1e-9,1e9>	
	tau_K_gamma = 20 (ms) <1e-9,1e9>	
	m_K_alpha = 1.5						
	m_K_beta = 1.7						
	s_K_beta = 100						
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	g_ampa (uS)
	g_nmda (uS)
	epsilon
	epsilon_G
	epsilon_Z
	C
	T
	G
	E
	P
	K_alpha
	Rho
	K_beta
	K
	N
	X
	Z
	flag_D
	w_pre
	w_post
	w
	A_n
	B_n
	last_weight
}

STATE {
	A (uS)
	B (uS)
	G_a
	G_b
	u_bar
	K_alpha_bar
	K_gamma
	Z_a
	Z_b
	N_alpha_bar
	N_beta_bar
}

INITIAL {
	LOCAL omega, omega_G, omega_Z
	g = 0
	g_ampa = 0
	g_nmda = 0
	u_bar = 0
	K_alpha_bar = 0
	K_gamma = 0
	flag_D = -1
	N_alpha_bar = 0
	N_beta_bar = 0
	w_pre = w_pre_init
	w_post = w_post_init
	w = w_pre * w_post
	last_weight = 0
	
	
	
	if (tau_a/tau_b > .9999) {
		tau_a = .9999*tau_b
	}
	A = 0
	B = 0
	omega = (tau_a*tau_b)/(tau_b - tau_a) * log(tau_b/tau_a)
	epsilon = -exp(-omega/tau_a) + exp(-omega/tau_b)
	epsilon = 1/epsilon
	
	
	if (tau_G_a/tau_G_b > .9999) {
		tau_G_a = .9999*tau_G_b
	}
	G_a = 0
	G_b = 0
	omega_G = (tau_G_a*tau_G_b)/(tau_G_b - tau_G_a) * log(tau_G_b/tau_G_a)
	epsilon_G = -exp(-omega_G/tau_G_a) + exp(-omega_G/tau_G_b)
	epsilon_G = 1/epsilon_G
	
	
	if (tau_Z_a/tau_Z_b > .9999) {
		tau_Z_a = .9999*tau_Z_b
	}
	Z_a = 0
	Z_b = 0
	omega_Z = (tau_Z_a*tau_Z_b)/(tau_Z_b - tau_Z_a) * log(tau_Z_b/tau_Z_a)
	epsilon_Z = -exp(-omega_Z/tau_Z_a) + exp(-omega_Z/tau_Z_b)
	epsilon_Z = 1/epsilon_Z
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	G = sigmoid_sat(m_G, G_b - G_a)			
	Z = sigmoid_sat(m_Z, Z_b - Z_a)			
	
	g_ampa = s_ampa * w_post * (B - A)							
	g_nmda = s_nmda * w_post_init * (B_n - A_n) * mgblock(v)	
	
	g = w_pre * (g_ampa + g_nmda)			
	i = g * (v - e)
}

DERIVATIVE state {
	LOCAL D, u, Eta, LTD_pre, LTP_pre, LTD_post, LTP_post, g_update, N_alpha, N_beta, N
	
	if(flag_D == 1) {	
		D = 1
	    flag_D = -1
	}
	else {
	    D = 0
	}
	
	UNITSOFF
	u = v		
	Eta = dt	
	UNITSON
	
	
	u_bar' = (- u_bar + positive(u - theta_u_T)) / tau_u_T
	T = sigmoid_sat(m_T, u_bar)									
	E = D * T													
	
	
	N_alpha_bar' = (- N_alpha_bar + positive(u - theta_u_N)) / tau_N_alpha
	N_beta_bar'	= (- N_beta_bar + N_alpha_bar) / tau_N_beta
	N_alpha = sigmoid_sat(m_N_alpha, N_alpha_bar)
	N_beta = sigmoid_sat(m_N_beta, N_beta_bar)
	N = positive(N_alpha * N_beta - theta_N_X)
	X = Z * N
	
	
	C = G * positive(u - theta_u_C)
	P = positive(C - theta_C_minus) * positive(theta_C_plus - C) / pow((theta_C_plus - theta_C_minus) / 2, 2)
	
	
	K_alpha = sigmoid_sat(m_K_alpha, positive(C - theta_C_plus)) * Rho
	K_alpha_bar' = (- K_alpha_bar + K_alpha) / tau_K_alpha
	K_beta = sigmoid_sat(m_K_beta, (K_alpha_bar * s_K_beta))
	Rho = 1.0 - K_beta
	K_gamma' = (- K_gamma + K_beta) / tau_K_gamma
	K = K_alpha * K_beta * K_gamma
	
	
	LTD_pre = - A_LTD_pre * E
	LTP_pre	= A_LTP_pre * X * Eta			
	LTD_post = - A_LTD_post * P * Eta
	LTP_post = A_LTP_post * K * Eta
	
	
	w_pre = w_pre + LTD_pre + LTP_pre
	if(w_pre > 1.0) {
		w_pre = 1.0
	}
	if(w_pre < 0.0) {
		w_pre = 0.0
	}
	w_post = w_post + LTD_post + LTP_post
	if(w_post > 5.0) {
		w_post = 5.0
	}
	if(w_post < 0.0) {
		w_post = 0.0
	}
	w = w_pre * w_post
	
	A' = -A / tau_a
	B' = -B / tau_b
	G_a' = -G_a / tau_G_a
	G_b' = -G_b / tau_G_b
	A_n = G_a * last_weight	
	B_n = G_b * last_weight
	Z_a' = -Z_a / tau_Z_a
	Z_b' = -Z_b / tau_Z_b
}

NET_RECEIVE(weight (uS)) {
	
	A = A + weight * epsilon	
	B = B + weight * epsilon
	G_a = G_a + epsilon_G
	G_b = G_b + epsilon_G
	Z_a = Z_a + epsilon_Z
	Z_b = Z_b + epsilon_Z
	flag_D = 1
	last_weight = weight	
}

FUNCTION positive(value) {	
	if(value < 0) {
		positive = 0
	}
	else {
		positive = value
	}
}

FUNCTION sigmoid_sat(slope, value) {	
	sigmoid_sat = 2.0 / (1.0 + pow(slope, -value)) - 1.0
}

FUNCTION mgblock(v(mV)) {	
	LOCAL u
	UNITSOFF
	u = v
	UNITSON
	
	mgblock = 1 / (1 + exp(0.080 * -u) * (1.0 / 3.57))
}