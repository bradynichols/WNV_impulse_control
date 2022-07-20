syms Lambda1 Hs1 b p_mh Vi Hi1 Hr1 omega1 p_hh1 d_h1 mu_h1 Ls p_hm1 Vs NH1 p_hm2 Hi2 NH2 p_hm3 Hi3 NH3 Ua

% (Hs1 + Hi1 + Hr1)

y = m_l*Ls - b*p_hm1*Vs*Hi1/NH1 - b*p_hm2*Vs*Hi2/NH2 - b*p_hm3*Vs*Hi3/NH3 - muV*Vs - km2*Vs*Ua;

diff(y,Vs)