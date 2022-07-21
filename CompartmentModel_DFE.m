syms rs ri phi qs qi m_e m_l muL muV b c_l kl p_mh gl ga km1 km2 cV d_l
syms p_hm1 p_hm2 p_hm3 omega1 omega2 omega3 p_hh1 p_hh2 p_hh3 g1 g2 g3 
syms gamma1 gamma2 gamma3 mu_h1 mu_h2 mu_h3 c_h1 c_h2 c_h3 d_h1 d_h2 d_h3 Lambda1 Lambda2 Lambda3
syms Hs1 Hs2 Hs3 Hi1 Hi2 Hi3 Hr1 Hr2 Hr3 NH1 NH2 NH3 rec1 rec2 rec3 Vs Es Ls

rec1 = Hs1 + Hr1
rec2 = Hs2 + Hr2
rec3 = Hs3 + Hr3
NH1 = Hs1 + Hi1 + Hr1
NH2 = Hs2 + Hi2 + Hr2
NH3 = Hs3 + Hi3 + Hr3

dHs1 = Lambda1*Hs1 - d_h1*Hs1 - mu_h1*Hs1 == 0
solve(dHs1, Hs1)
% Hs1 = c_h1
% Hs2 = c_h2
% Hs3 = c_h3
% Ls = c_l
% Vs = c_l*m_l/muV
% Es = c_l*m_l*rs/(muV*m_e)

