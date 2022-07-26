syms Hs Hi Hr Es Ei Ls Li Vs Ve Vi Ua Ul
syms rs ri phi qs qi m_e m_l mu_h muL muV b c_l kl p_hm p_mh d_h g gl ga km1 km2 cV d_l gamma omega p_hh c_h Lambda NH
NH = Hs + Hi + Hr

f=[dEs; dEi; dLs; dLi; dVs; dVe; dVi; dHs; dHi; dHr; dUl; dUa; dx_int];

jacobian(f, [Es Ei Ls Li Vs Ve Vi Hs Hi Hr Ul Ua x_int])