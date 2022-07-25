% Updated 07-14-2022

function [dxdt] = West_Nile_ModelRL(t,x,p)
% This file gives the derivatives of the continuous state variables.
x1=x(1:19);
% x2=dx/dx0 is the derivative of the state with respect to its initial value.
% x2=[dx1/dx1(0), dx1/dx2(0),. . .dx1/dx19(0);
% dx2/dx1(0), dx2/dx2(0), . . .dx2/dx19(0);
%    .
%    .
%    .
% dx19/dx1(0), dx19/dx2(0), . . .,dx19/dx19(0)];
x2 = x(length(x1)+1:length(x1)+length(x1)^2);
x2 = reshape(x2,length(x1),length(x1));

% Model Parameters

rs = p(1); % egg laying rate of S and E mosquitoes
ri = p(2); % egg laying rate of I mosquitoes

phi = p(3); % fraction of eggs born to infected mothers that are infected (fraction of eggs infected)
qs = p(4); % fraction of eggs from uninfected mosquitoes that hatch
qi = p(5); % fraction of eggs laid to infected mosquitoes that hatch

m_e = p(6); % hatch rate 
m_l = p(7); % larval maturation rate

muL = p(8); % larval death rate
muV = p(9); % adult death rate

b = p(10); % mosquito biting rate
c_l = p(11); % mosquito carrying capacity (larval)

kl = p(12); % disease progression in mosquitoes (1/latency period)
p_mh = p(13); % mosquito-to-host transmission

p_hm1 = p(14); % host-to-mosquito transmission host group 1
p_hm2 = p(15); % host-to-mosquito transmission host group 2
p_hm3 = p(16); % host-to-mosquito transmission host group 3

omega1 = p(17); % direct transmission rate host group 1
omega2 = p(18); % direct transmission rate host group 2
omega3 = p(19); % direct transmission rate host group 3

p_hh1 = p(20); % contact rate host group 1
p_hh2 = p(21); % contact rate host group 2
p_hh3 = p(22); % contact rate host group 3 

g1 = p(23); % WNV recovery host group 1
g2 = p(24); % WNV recovery host group 2
g3 = p(25); % WNV recovery host group 3

gamma1 = p(26); % WNV death host group 1
gamma2 = p(27); % WNV death host group 2
gamma3 = p(28); % WNV death host group 3

Lambda1 = p(29); % recruitment rate host group 1
Lambda2 = p(30); % recruitment rate host group 2
Lambda3 = p(31); % recruitment rate host group 3

mu_h1 = p(32); % natural death rate host group 1
mu_h2 = p(33); % natural death rate host group 2
mu_h3 = p(34); % natural death rate host group 3

c_h1 = p(35); % carrying capacity host group 1
c_h2 = p(36); % carrying capacity host group 2
c_h3 = p(37); % carrying capacity host group 3

km1 = p(38); % max rate at which larvicide kills larvae
km2 = p(40); % max rate at which adulticide kills adult vectors

gl = p(39); % larvicide decay rate
ga = p(41); % adulticide decay rate

cV = p(42); % weight of cost of vectors in objective functional

d_l=((rs*m_l*qs/muV)-muL-m_l); % density-dependent death rate for larvae
d_h1 = (Lambda1 - mu_h1); % density-dependent death rate for host group 1
d_h2 = (Lambda2 - mu_h2); % density-dependent death rate for host group 2
d_h3 = (Lambda3 - mu_h3); % density-dependent death rate for host group 3

% State variables

% Vector

Es = x(1); % eggs laid by susceptible and exposed mothers
Ei = x(2); % eggs laid by infected mothers
Ls = x(3); % susceptible larvae
Li = x(4); % infected larvae
Vs = x(5); % susceptible vectors
Ve = x(6); % exposed vectors
Vi = x(7); % infected vectors

% Host

Hs1 = x(8); % susceptible host group 1
Hs2 = x(9); % susceptible host group 2
Hs3 = x(10); % susceptible host group 3

Hi1 = x(11); % infected host group 1
Hi2 = x(12); % infected host group 2
Hi3 = x(13); % infected host group 3

Hr1 = x(14); % recovered host group 1
Hr2 = x(15); % recovered host group 2
Hr3 = x(16); % recovered host group 3
    
% Control 

Ul = x(17); % larvacide
Ua = x(18); % adultacide

NH1 = Hs1 + Hi1 + Hr1; % total hosts group 1
NH2 = Hs2 + Hi2 + Hr2; % total hosts group 2
NH3 = Hs3 + Hi3 + Hr3; % total hosts group 3

Hi1 = 0;
Hi2 = 0;
Hi3 = 0;
Hr1 = 0;
Hr2 = 0;
Hr3 = 0;
Ei = 0;
Li = 0;
Ve = 0;
Vi = 0;

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

dx_int = cV*(Vi+Vs+Ve); 

% NEEDS UPDATE

f=[dEs; dEi; dLs; dLi; dVs; dVe; dVi; dHs1; dHs2; dHs3; dHi1; dHi2; dHi3; dHr1; dHr2; dHr3; dUl; dUa; dx_int]; 

% dfdx = jacobian(f, [Es Ei Ls Li Vs Ve Vi Hs1 Hs2 Hs3 Hi1 Hi2 Hi3 Hr1 Hr2 Hr3 Ul Ua])

dfdx = 
[  -m_e,                 0,                                             0,                                             0,                                                                                                                   rs,                  rs,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0,       0,       0, 0;
     0,              -m_e,                                             0,                                             0,                                                                                                                    0,                   0,                              ri,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0,       0,       0, 0;
m_e*qs, -m_e*qi*(phi - 1), - m_l - muL - Ls*d_l/c_l - Ul*km1 - d_l*(Li + Ls)/c_l,                                       -Ls*d_l/c_l,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0, -Ls*km1,       0, 0;
     0,        m_e*phi*qi,                                       -Li*d_l/c_l, - m_l - muL - Li*d_l/c_l - Ul*km1 - d_l*(Li + Ls)/c_l,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0, -Li*km1,       0, 0;
     0,                 0,                                           m_l,                                             0, - muV - Ua*km2 - (Hi1*b*p_hm1)/(Hi1 + Hr1 + Hs1) - (Hi2*b*p_hm2)/(Hi2 + Hr2 + Hs2) - (Hi3*b*p_hm3)/(Hi3 + Hr3 + Hs3),                   0,                               0,                                                                                                                                                                          (Hi1*Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1)^2,                                                                                                                                                                          (Hi2*Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2)^2,                                                                                                                                                                          (Hi3*Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3)^2,                                                                                                             (Hi1*Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1)^2 - (Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1),                                                                                                             (Hi2*Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2)^2 - (Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2),                                                                                                             (Hi3*Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3)^2 - (Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3),                                                          (Hi1*Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1)^2,                                                          (Hi2*Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2)^2,                                                          (Hi3*Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3)^2,       0, -Vs*km2, 0;
     0,                 0,                                             0,                                             0,                  (Hi1*b*p_hm1)/(Hi1 + Hr1 + Hs1) + (Hi2*b*p_hm2)/(Hi2 + Hr2 + Hs2) + (Hi3*b*p_hm3)/(Hi3 + Hr3 + Hs3), - kl - muV - Ua*km2,                               0,                                                                                                                                                                         -(Hi1*Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1)^2,                                                                                                                                                                         -(Hi2*Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2)^2,                                                                                                                                                                         -(Hi3*Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3)^2,                                                                                                             (Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1) - (Hi1*Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1)^2,                                                                                                             (Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2) - (Hi2*Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2)^2,                                                                                                             (Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3) - (Hi3*Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3)^2,                                                         -(Hi1*Vs*b*p_hm1)/(Hi1 + Hr1 + Hs1)^2,                                                         -(Hi2*Vs*b*p_hm2)/(Hi2 + Hr2 + Hs2)^2,                                                         -(Hi3*Vs*b*p_hm3)/(Hi3 + Hr3 + Hs3)^2,       0, -Ve*km2, 0;
     0,                 0,                                             0,                                           m_l,                                                                                                                    0,                  kl,                  - muV - Ua*km2,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0,       0, -Vi*km2, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0, -(Hs1*b*p_mh)/(Hi1 + Hr1 + Hs1), Lambda1 - mu_h1 - Hs1*d_h1/c_h1 - d_h1*(Hi1 + Hr1 + Hs1)/c_h1 - (Hi1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1) - (Vi*b*p_mh)/(Hi1 + Hr1 + Hs1) + (Hi1*Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1)^2 + (Hs1*Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)^2,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                (Hi1*Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1)^2 - (Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1) - Hs1*d_h1/c_h1 + (Hs1*Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)^2,                                                                                                                                                                                 0,                                                                                                                                                                                 0,   (Hi1*Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1)^2 - Hs1*d_h1/c_h1 + (Hs1*Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)^2,                                                                                             0,                                                                                             0,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0, -(Hs2*b*p_mh)/(Hi2 + Hr2 + Hs2),                                                                                                                                                                                                             0, Lambda2 - mu_h2 - Hs2*d_h2/c_h2 - d_h2*(Hi2 + Hr2 + Hs2)/c_h2 - (Hi2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2) - (Vi*b*p_mh)/(Hi2 + Hr2 + Hs2) + (Hi2*Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2)^2 + (Hs2*Vi*b*p_mh)/(Hi2 + Hr2 + Hs2)^2,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                (Hi2*Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2)^2 - (Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2) - Hs2*d_h2/c_h2 + (Hs2*Vi*b*p_mh)/(Hi2 + Hr2 + Hs2)^2,                                                                                                                                                                                 0,                                                                                             0,   (Hi2*Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2)^2 - Hs2*d_h2/c_h2 + (Hs2*Vi*b*p_mh)/(Hi2 + Hr2 + Hs2)^2,                                                                                             0,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0, -(Hs3*b*p_mh)/(Hi3 + Hr3 + Hs3),                                                                                                                                                                                                             0,                                                                                                                                                                                                             0, Lambda3 - mu_h3 - Hs3*d_h3/c_h3 - d_h3*(Hi3 + Hr3 + Hs3)/c_h3 - (Hi3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3) - (Vi*b*p_mh)/(Hi3 + Hr3 + Hs3) + (Hi3*Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3)^2 + (Hs3*Vi*b*p_mh)/(Hi3 + Hr3 + Hs3)^2,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                (Hi3*Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3)^2 - (Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3) - Hs3*d_h3/c_h3 + (Hs3*Vi*b*p_mh)/(Hi3 + Hr3 + Hs3)^2,                                                                                             0,                                                                                             0,   (Hi3*Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3)^2 - Hs3*d_h3/c_h3 + (Hs3*Vi*b*p_mh)/(Hi3 + Hr3 + Hs3)^2,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,  (Hs1*b*p_mh)/(Hi1 + Hr1 + Hs1),                                            (Hi1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1) - Hi1*d_h1 + (Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)*c_h1 - (Hi1*Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1)^2 - (Hs1*Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)^2,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0, (Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1) - gamma1 - mu_h1 - Hi1*/c_h1 - d_h1*(Hi1 + Hr1 + Hs1)/c_h1 - g1 - (Hi1*Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1)^2 - (Hs1*Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)^2,                                                                                                                                                                                 0,                                                                                                                                                                                 0, - Hi1*d_h1/c_h1 - (Hi1*Hs1*omega1*p_hh1)/(Hi1 + Hr1 + Hs1)^2 - (Hs1*Vi*b*p_mh)/(Hi1 + Hr1 + Hs1)^2,                                                                                             0,                                                                                             0,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,  (Hs2*b*p_mh)/(Hi2 + Hr2 + Hs2),                                                                                                                                                                                                             0,                                            (Hi2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2) - Hi2*d_h2/c_h2 + (Vi*b*p_mh)/(Hi2 + Hr2 + Hs2) - (Hi2*Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2)^2 - (Hs2*Vi*b*p_mh)/(Hi2 + Hr2 + Hs2)^2,                                                                                                                                                                                                             0,                                                                                                                                                                                 0, (Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2) - gamma2 - mu_h2 - Hi2*d_h2/c_h2 - d_h2*(Hi2 + Hr2 + Hs2)/c_h2 - g2 - (Hi2*Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2)^2 - (Hs2*Vi*b*p_mh)/(Hi2 + Hr2 + Hs2)^2,                                                                                                                                                                                 0,                                                                                             0, - Hi2*d_h2/c_h2 - (Hi2*Hs2*omega2*p_hh2)/(Hi2 + Hr2 + Hs2)^2 - (Hs2*Vi*b*p_mh)/(Hi2 + Hr2 + Hs2)^2,                                                                                             0,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,  (Hs3*b*p_mh)/(Hi3 + Hr3 + Hs3),                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                            (Hi3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3) - Hi3*d_h3/c_h3 + (Vi*b*p_mh)/(Hi3 + Hr3 + Hs3) - (Hi3*Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3)^2 - (Hs3*Vi*b*p_mh)/(Hi3 + Hr3 + Hs3)^2,                                                                                                                                                                                 0,                                                                                                                                                                                 0, (Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3) - gamma3 - mu_h3 - Hi3*d_h3/c_h3 - d_h3*(Hi3 + Hr3 + Hs3)/c_h3 - g3 - (Hi3*Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3)^2 - (Hs3*Vi*b*p_mh)/(Hi3 + Hr3 + Hs3)^2,                                                                                             0,                                                                                             0, - Hi3*d_h3/c_h3 - (Hi3*Hs3*omega3*p_hh3)/(Hi3 + Hr3 + Hs3)^2 - (Hs3*Vi*b*p_mh)/(Hi3 + Hr3 + Hs3)^2,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                     -Hr1*d_h1/c_h1,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                     g1 - Hr1*d_h1/c_h1,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                   - mu_h1 - Hr1*d_h1/c_h1 - d_h1*(Hi1 + Hr1 + Hs1)/c_h1,                                                                                             0,                                                                                             0,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                     -Hr2*d_h2/c_h2,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                     g2 - Hr2*d_h2/c_h2,                                                                                                                                                                                 0,                                                                                             0,                                                   - mu_h2 - Hr2*d_h2/c_h2 - d_h2*(Hi2 + Hr2 + Hs2),                                                                                             0,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                     -Hr3*d_h3/c_h3,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                     g3 - Hr3*d_h3/c_h3,                                                                                             0,                                                                                             0,                                                   - mu_h3 - Hr3*d_h3/c_h3 - d_h3*(Hi3 + Hr3 + Hs3)/c_h3,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0,     -gl,       0, 0;
     0,                 0,                                             0,                                             0,                                                                                                                    0,                   0,                               0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0,       0,     -ga, 0;
     0,                 0,                                             0,                                             0,                                                                                                                   cV,                  cV,                              cV,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                                             0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                                                                                                                 0,                                                                                             0,                                                                                             0,                                                                                             0,       0,       0, 0];



dx2dt=dfdx*x2;

dx2dt=reshape(dx2dt,361,1);

dxdt=[f; dx2dt];
end