syms Lambda1 Hs1 b p_mh Vi Hi1 Hr1 omega1 p_hh1 d_h1 mu_h1

% (Hs1 + Hi1 + Hr1)

y = Lambda1*Hs1 - b*p_mh*Vi*Hs1/(Hs1 + Hi1 + Hr1) - omega1*p_hh1*Hi1*Hs1/(Hs1 + Hi1 + Hr1) - d_h1*(Hs1 + Hi1 + Hr1)*Hs1 - mu_h1*Hs1;

diff(y,Hr1)