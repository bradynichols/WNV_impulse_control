syms Hs Hi Hr Es Ei Ls Li Vs Ve Vi Ua Ul x_int
syms rs ri phi qs qi m_e m_l mu_h muL muV b c_l kl p_hm p_mh d_h g gl ga km1 km2 cV d_l gamma omega p_hh c_h Lambda NH 
NH = Hs + Hi + Hr

    %Vector ODEs
    dEs=  rs*(Vs+Ve)-m_e*Es;
    dEi=  ri*(Vi)-m_e*Ei;
    dLs = m_e*qs*Es+m_e*qi*(1-phi)*Ei-muL*Ls-m_l*Ls-d_l*Ls*(Ls+Li)-km1*Ls*Ul;
    dLi = m_e*qi*phi*Ei-muL*Li-m_l*Li-d_l*Li*(Li+Ls)-km1*Li*Ul;
    dVs = m_l*Ls-b*p_hm*Vs*Hi/NH-muV*Vs-km2*Vs*Ua;
    dVe = b*p_hm*Vs*Hi/NH-kl*Ve-muV*Ve-km2*Ve*Ua;
    dVi = m_l*Li+kl*Ve-muV*Vi-km2*Vi*Ua;

    %Host ODEs
    dHs = Lambda*(Hs+Hr) - b*p_mh*Vi*Hs/NH - omega*p_hh*Hi*Hs/NH - d_h*NH*Hs - mu_h*Hs;                    
    dHi = b*p_mh*Vi*Hs/NH + omega*p_hh*Hi*Hs/NH - gamma*Hi - g*Hi - d_h*NH*Hi - mu_h*Hi;
    dHr = g*Hi - d_h*NH*Hr - mu_h*Hr;
    
    %Chemical ODEs
    dUl=-gl*Ul;
    dUa=-ga*Ua;
    
    %Integral State ODE
    dx_int=(Vi+Hi);

f=[dEs; dEi; dLs; dLi; dVs; dVe; dVi; dHs; dHi; dHr; dUl; dUa; dx_int];

dfdx = jacobian(f, [Es Ei Ls Li Vs Ve Vi Hs Hi Hr Ul Ua x_int])

sz = size(dfdx)
