% Updated 07-14-2022

function [dxdt] = West_Nile_ModelRL(t,x,p)
%This file gives the derivatives of the continuous state variables.
%x=[Es,Ei,Ls,Li,Vs,Ve,Vi,Hs,Hi,Hr,Ul,Ua,arti_state,reshape(dx/dx0,13^2,1)]
%  Es = x(1);  Ei = x(2);  Ls = x(3);   Li = x(4); Vs=x(5); Ve=x(6); Vi=x(7)
%  Hs = x(8);  Hi = x(9);  Hr = x(10); Ul=x(11), Ua=x(12), int_0^t{cV*(Vi(s)+Vs(s)+Ve(s))ds}=x(13)
x1=x(1:13);
%x2=dx/dx0 is the derivative of the state with respect to its initial value.
%x2=[dx1/dx1(0), dx1/dx2(0),. . .dx1/dx13(0);
%    dx2/dx1(0), dx2/dx2(0), . . .dx2/dx13(0);
%    .
%    .
%    .
%     dx13/dx1(0), dx13/dx2(0), . . .,dx13/dx13(0)];
x2=x(length(x1)+1:length(x1)+length(x1)^2);
x2=reshape(x2,length(x1),length(x1));

%Model Parameters

%Vector
rs = p(1); %intrinsic rate of increase of uninfected mosquitoes
ri = p(2); %intrinsic rate of increase of infected mosquitoes
phi=p(3); %fraction of eggs born to infected mothers that are infected
qs=p(4); %fraction of eggs from uninfected mosquitoes that hatch
qi=p(5); %fraction of eggs laid to infected mosquitoes that hatch
m_e = p(6); %egg maturatation rate
m_l = p(7); %larval maturation rate

muL = p(8); %larval death rate
muV=p(9); %adult death rate
b = p(10); %mosquito biting rate
c_l = p(11); %mosquito carrying capacity

%Disease
kl = p(12);  %disease progression (1/latency period)
p_mh = p(13); %mosquito-to-host transmission
p_hm = p(14); %host-to-mosquito transmission

gamma = p(15); % WNV- induced host mortality
g = p(16); %host recovery rate

gl=p(17);
ga=p(18);

km1=p(19);
km2=p(20);

cV=p(21); %weight of cost of vectors in objective functional

Lambda = p(30); % Host recruitment rate
mu_h = p(31); % Host natural death rate

p_hh = p(32); %Host-Host contact rate
omega = p(33); % Host-Host transmission probability

c_h = p(34); % Host carrying capacity

d_l=((rs*m_l*qs/muV)-muL-m_l)/c_l; % density-dependent death rate for larvae
d_h = (Lambda - mu_h)/c_h; % density-dependent death rate for host

    %State variables

    %Vector
    Es = x(1);              %eggs laid by susceptible and exposed mothers
    Ei = x(2);              %eggs laid by infected mothers
    Ls = x(3);              %susceptible larvae
    Li = x(4);              %infected larvae
    Vs = x(5);              %susceptible vectors
    Ve = x(6);              %exposed vectors
    Vi = x(7);              %infected vectors

    %Host
    Hs = x(8);              %susceptible hosts
    Hi = x(9);              %infected hosts
    Hr = x(10);              %recovered hosts
    
    %Control 
    Ul = x(11);              %larvacide
    Ua = x(12);             %adultacide

    NH = Hs+Hi+Hr;           %total hosts

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
    dx_int=cV*(Vi+Vs+Ve); 

f=[dEs; dEi; dLs; dLi; dVs; dVe; dVi; dHs; dHi; dHr; dUl; dUa; dx_int];

% Updated to include host demographics

% dfdx=[ -m_e, 0, 0, 0, rs, rs, 0, 0, 0, 0, 0, 0, 0;
%         0, -m_e, 0, 0, 0, 0, ri, 0, 0, 0, 0, 0, 0;
%       m_e*qs, m_e*qi*(1-phi), -muL-m_l-km1*Ul-2*d_l*Ls-d_l*Li, -d_l*Ls, 0, 0, 0, 0, 0, 0, -km1*Ls, 0, 0;
%  0, m_e*qi*phi, -d_l*Li, -muL-m_l-km1*Ul-2*d_l*Li-d_l*Ls, 0, 0, 0, 0, 0, 0, -km1*Li, 0, 0;
%  0, 0, m_l, 0, -b*p_hm*Hi/NH-muV-km2*Ua, 0, 0, b*p_hm*Vs*Hi/(NH)^2, -b*p_hm*Vs/(NH)+b*p_hm*Vs*Hi/(NH)^2, b*p_hm*Vs*Hi/(NH)^2, 0, -km2*Vs, 0;
%  0, 0, 0, 0, b*p_hm*Hi/NH, -kl-muV-km2*Ua, 0, -b*p_hm*Vs*Hi/(NH)^2, b*p_hm*Vs/(NH)-b*p_hm*Hi*Vs/(NH)^2, -b*p_hm*Vs*Hi/(NH)^2, 0, -km2*Ve, 0;
%  0, 0, 0, m_l, 0, kl, -muV-km2*Ua, 0, 0, 0, 0, -km2*Vi, 0;
%  0, 0, 0, 0, 0, 0, -b*p_mh*Hs/NH, Lambda - mu_h - Hs*d_h - d_h*(NH) - (Hi*omega*p_hh)/(NH) - (Vi*b*p_mh)/(NH) + (Hi*Hs*omega*p_hh)/(NH)^2 + (Hs*Vi*b*p_mh)/(NH)^2, (Hi*Hs*omega*p_hh)/(NH)^2 - (Hs*omega*p_hh)/(NH) - Hs*d_h + (Hs*Vi*b*p_mh)/(NH)^2, (Hi*Hs*omega*p_hh)/(NH)^2 - Hs*d_h + (Hs*Vi*b*p_mh)/(NH)^2, 0, 0, 0; % Hs updated 07-19-2022
%  0, 0, 0, 0, 0, 0, b*p_mh*Hs/NH, (Hi*omega*p_hh)/(NH) - d_h*(NH) - Hs*d_h + (Vi*b*p_mh)/(NH) - (Hi*Hs*omega*p_hh)/(NH)^2 - (Hs*Vi*b*p_mh)/(NH)^2, (Hs*omega*p_hh)/(NH) - g - mu_h - Hs*d_h - mu_h - (Hi*Hs*omega*p_hh)/(NH)^2 - (Hs*Vi*b*p_mh)/(NH)^2, - Hs*d_h - (Hi*Hs*omega*p_hh)/(NH)^2 - (Hs*Vi*b*p_mh)/(NH)^2, 0, 0, 0; % Hi updated 07-19-2022
%  0, 0, 0,  0,  0,  0,  0,  - Hs*d_h - d_h*(NH),  g - Hs*d_h,  - mu_h - Hs*d_h,  0,  0,  0; % Hr updated 07/19/2022.
%  0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  -gl, 0, 0;
%  0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  -ga,  0; 
%  0, 0, 0,  0,  cV,  cV,  cV,  0,  0,  0,  0,  0,  0];


dfdx = [  -m_e,                 0,                                             0,                                             0,                                          rs,                  rs,                           0,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0,       0,       0, 0;
     0,              -m_e,                                             0,                                             0,                                           0,                   0,                          ri,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0,       0,       0, 0;
m_e*qs, -m_e*qi*(phi - 1), - m_l - muL - Ls*d_l - Ul*km1 - d_l*(Li + Ls),                                       -Ls*d_l,                                           0,                   0,                           0,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0, -Ls*km1,       0, 0;
     0,        m_e*phi*qi,                                       -Li*d_l, - m_l - muL - Li*d_l - Ul*km1 - d_l*(Li + Ls),                                           0,                   0,                           0,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0, -Li*km1,       0, 0;
     0,                 0,                                           m_l,                                             0, - muV - Ua*km2 - (Hi*b*p_hm)/(Hi + Hr + Hs),                   0,                           0,                                                                                                                                                   (Hi*Vs*b*p_hm)/(Hi + Hr + Hs)^2,                                                                                            (Hi*Vs*b*p_hm)/(Hi + Hr + Hs)^2 - (Vs*b*p_hm)/(Hi + Hr + Hs),                                                         (Hi*Vs*b*p_hm)/(Hi + Hr + Hs)^2,       0, -Vs*km2, 0;
     0,                 0,                                             0,                                             0,                  (Hi*b*p_hm)/(Hi + Hr + Hs), - kl - muV - Ua*km2,                           0,                                                                                                                                                  -(Hi*Vs*b*p_hm)/(Hi + Hr + Hs)^2,                                                                                            (Vs*b*p_hm)/(Hi + Hr + Hs) - (Hi*Vs*b*p_hm)/(Hi + Hr + Hs)^2,                                                        -(Hi*Vs*b*p_hm)/(Hi + Hr + Hs)^2,       0, -Ve*km2, 0;
     0,                 0,                                             0,                                           m_l,                                           0,                  kl,              - muV - Ua*km2,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0,       0, -Vi*km2, 0;
     0,                 0,                                             0,                                             0,                                           0,                   0, -(Hs*b*p_mh)/(Hi + Hr + Hs), Lambda - mu_h - Hs*d_h - d_h*(Hi + Hr + Hs) - (Hi*omega*p_hh)/(Hi + Hr + Hs) - (Vi*b*p_mh)/(Hi + Hr + Hs) + (Hi*Hs*omega*p_hh)/(Hi + Hr + Hs)^2 + (Hs*Vi*b*p_mh)/(Hi + Hr + Hs)^2,                                         (Hi*Hs*omega*p_hh)/(Hi + Hr + Hs)^2 - (Hs*omega*p_hh)/(Hi + Hr + Hs) - Hs*d_h + (Hs*Vi*b*p_mh)/(Hi + Hr + Hs)^2, Lambda - Hs*d_h + (Hi*Hs*omega*p_hh)/(Hi + Hr + Hs)^2 + (Hs*Vi*b*p_mh)/(Hi + Hr + Hs)^2,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                           0,                   0,  (Hs*b*p_mh)/(Hi + Hr + Hs),                                      (Hi*omega*p_hh)/(Hi + Hr + Hs) - Hi*d_h + (Vi*b*p_mh)/(Hi + Hr + Hs) - (Hi*Hs*omega*p_hh)/(Hi + Hr + Hs)^2 - (Hs*Vi*b*p_mh)/(Hi + Hr + Hs)^2, (Hs*omega*p_hh)/(Hi + Hr + Hs) - gamma - mu_h - Hi*d_h - d_h*(Hi + Hr + Hs) - g - (Hi*Hs*omega*p_hh)/(Hi + Hr + Hs)^2 - (Hs*Vi*b*p_mh)/(Hi + Hr + Hs)^2,        - Hi*d_h - (Hi*Hs*omega*p_hh)/(Hi + Hr + Hs)^2 - (Hs*Vi*b*p_mh)/(Hi + Hr + Hs)^2,       0,       0, 0;
     0,                 0,                                             0,                                             0,                                           0,                   0,                           0,                                                                                                                                                                           -Hr*d_h,                                                                                                                                              g - Hr*d_h,                                                    - mu_h - Hr*d_h - d_h*(Hi + Hr + Hs),       0,       0, 0;
     0,                 0,                                             0,                                             0,                                           0,                   0,                           0,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0,     -gl,       0, 0;
     0,                 0,                                             0,                                             0,                                           0,                   0,                           0,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0,       0,     -ga, 0;
     0,                 0,                                             0,                                             0,                                          cV,                  cV,                          cV,                                                                                                                                                                                 0,                                                                                                                                                       0,                                                                                       0,       0,       0, 0];


dx2dt=dfdx*x2;

dx2dt=reshape(dx2dt,169,1);

dxdt=[f; dx2dt];
end