syms Hi1 Hi2 Hi3 Ei Li Ve Vi n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 j1 j2 j3 j4 j5 j6 j7 j8
syms rs ri phi qs qi m_e m_l muL muV b c_l kl p_mh gl ga km1 km2 cV d_l Vs
syms p_hm1 p_hm2 p_hm3 omega1 omega2 omega3 p_hh1 p_hh2 p_hh3 g1 g2 g3 
syms gamma1 gamma2 gamma3 mu_h1 mu_h2 mu_h3 c_h1 c_h2 c_h3 d_h1 d_h2 d_h3 Lambda1 Lambda2 Lambda3
syms Hs1 Hs2 Hs3 Hr1 Hr2 Hr3 Es Ls Vs Ul Ua x_int

NH1 = Hs1 + Hi1 + Hr1; % total hosts group 1
NH2 = Hs2 + Hi2 + Hr2; % total hosts group 2
NH3 = Hs3 + Hi3 + Hr3; % total hosts group 3

% Vector ODEs % updated 7/20/2022

dEs = rs*(Vs+Ve)-m_e*Es;
dEi = ri*(Vi)-m_e*Ei;
dLs = m_e*qs*Es + m_e*qi*(1-phi)*Ei - muL*Ls - m_l*Ls - d_l*Ls*(Ls+Li)/c_l - km1*Ls*Ul;
dLi = m_e*qi*phi*Ei - muL*Li - m_l*Li - d_l*Li*(Li+Ls)/c_l - km1*Li*Ul;
dVs = m_l*Ls - b*p_hm1*Vs*Hi1/NH1 - b*p_hm2*Vs*Hi2/NH2 - b*p_hm3*Vs*Hi3/NH3 - muV*Vs - km2*Vs*Ua;
dVe = b*p_hm1*Vs*Hi1/NH1 + b*p_hm2*Vs*Hi2/NH2 + b*p_hm3*Vs*Hi3/NH3 - kl*Ve - muV*Ve - km2*Ve*Ua;
dVi = m_l*Li + kl*Ve - muV*Vi - km2*Vi*Ua;

% Host ODEs

dHs1 = Lambda1*NH1 - b*p_mh*Vi*Hs1/NH1 - omega1*p_hh1*Hi1*Hs1/NH1 - d_h1*NH1*Hs1/c_h1 - mu_h1*Hs1; % Updated 07/20/2022
dHs2 = Lambda2*NH2 - b*p_mh*Vi*Hs2/NH2 - omega2*p_hh2*Hi2*Hs2/NH2 - d_h2*NH2*Hs2/c_h2 - mu_h2*Hs2;
dHs3 = Lambda3*NH3 - b*p_mh*Vi*Hs3/NH3 - omega3*p_hh3*Hi3*Hs3/NH3 - d_h3*NH3*Hs3/c_h3 - mu_h3*Hs3;
dHi1 = b*p_mh*Vi*Hs1/NH1 + omega1*p_hh1*Hi1*Hs1/NH1 - (gamma1+g1)*Hi1 - d_h1*NH1*Hi1/c_h1 - mu_h1*Hi1; % Updated 07/20/2022
dHi2 = b*p_mh*Vi*Hs2/NH2 + omega2*p_hh2*Hi2*Hs2/NH2 - (gamma2+g2)*Hi2 - d_h2*NH2*Hi2/c_h2 - mu_h2*Hi2;
dHi3 = b*p_mh*Vi*Hs3/NH3 + omega3*p_hh3*Hi3*Hs3/NH3 - (gamma3+g3)*Hi3 - d_h3*NH3*Hi3/c_h3 - mu_h3*Hi3;
dHr1 = g1*Hi1 - d_h1*NH1*Hr1 - mu_h1*Hr1; % Updated 07/20/2022
dHr2 = g2*Hi2 - d_h2*NH2*Hr2 - mu_h2*Hr2;
dHr3 = g3*Hi3 - d_h3*NH3*Hr3 - mu_h3*Hr3;

%Chemical ODEs

dUl = -gl*Ul;
dUa = -ga*Ua;
    
    
% Integral State ODE

% dx_int = cV*(Vi+Vs+Ve); 
dx_int=(Vi+Hi1 + Hi2 + Hi3);


% NEEDS UPDATE

f=[dEs; dEi; dLs; dLi; dVs; dVe; dVi; dHs1; dHs2; dHs3; dHi1; dHi2; dHi3; dHr1; dHr2; dHr3; dUl; dUa; dx_int]; 

dfdx = jacobian(f, [Es Ei Ls Li Vs Ve Vi Hs1 Hs2 Hs3 Hi1 Hi2 Hi3 Hr1 Hr2 Hr3 Ul Ua x_int])

sz = size(dfdx)
