% Updated 08/27/2022
% Single-Compartment Reproductive Number

% Hi=infected hosts, Ei=infected eggs,
%Li=infected larva, Ve=exposed mosquitoes, Vi=infected mosquitoes

clear
syms Hi Ei Li Ve Vi
syms rs ri phi qs qi m_e m_l mu_h muL muV b c_l kl p_hm p_mh d_h g gl ga km1 km2 cV d_l gamma omega p_hh c_h % Included gamma, omega, p_hh

% DFE
Hs = c_h;
Ls = c_l;

% Compute Jacobian
% F=new infections, V=transfer between compartments 
% Only need to focus on infection compartments: [Hi1 Hi2 Hi3 Ei Li Ve Vi]

Ffun=[p_mh*b*Vi+omega*p_hh*Hi, ri*Vi, 0, b*p_hm*m_l*Hi*c_l/(c_h*muV), 0]; 
Vfun=[-gamma*Hi-g*Hi-d_h*Hi*Hs-mu_h*Hi, -m_e*Ei, m_e*qi*phi*Ei-muL*Li-m_l*Li-d_l*Li*Ls, -kl*Ve-muV*Ve, m_l*Li+kl*Ve-muV*Vi]; 

% Compute the jacobian with respect to infection compartments: [Hi Ei Li Ve Vi]

FF=jacobian(Ffun, [Hi Ei Li Ve Vi]);
VV=jacobian(Vfun, [Hi Ei Li Ve Vi]);
FF
VV

% Find matrix F and V
% Evaluate FF and VV at disease free equilibrium
% Only need to set infection compartments [I, As, Is, F, X, Ms, V] as zeros

MatrixF = subs(FF, [Hi Ei Li Ve Vi], [0, 0, 0, 0, 0])
MatrixV = subs(VV, [Hi Ei Li Ve Vi], [0, 0, 0, 0, 0])

%Compute F*V^{-1}
RR = -MatrixF*inv(MatrixV)

% Find eigenvalue of RR, largest eigenvalue = R0
syms lambda
pp=det(RR-lambda*eye(5));
pp_factors=factor(pp);
num_factors=length(pp_factors);
eigen_values=[];
for i=1:num_factors
eigen_values=[eigen_values; solve(pp_factors(i)==0, lambda,'MaxDegree', 5)];
end

p = ES_SC_Parameters(1,90);

rs = p(1); % egg laying rate of S and E mosquitoes
ri = p(2); % egg laying rate of I mosquitoes

phi = p(3); % fraction of eggs infected

qs = p(4); % fraction of eggs laid by uninfected mosquitoes that hatch
qi = p(5); % fraction of eggs laid by infected mosquitoes that hatch

m_e = p(6); % hatch rate

m_l = p(7); % larval maturation rate
muL = p(8); %larval death rate
muV = p(9); % adult death rate 

b = p(10); % mosquito biting rate
c_l = p(11); % mosquito larval carrying capacity

kl = p(12); % disease progression in mosquitoes (1/latency period)

p_mh = p(13); % mosquito-to-host transmission
p_hm = p(14); % host-to-mosquito transmission

omega = p(15) % host-host direct transmission rate
p_hh = p(16) % host-host contact rate

g = p(17); % host recovery rate
gamma = p(18); % WNV induced host mortality

Lambda = p(19) % host recruitment rate
mu_h = p(20) % host natural death rate

c_h = p(21) % host carrying capacity

km1 = p(22); % max rate at which larvicide kills larvae
gl = p(23); % larvicide decay rate

km2 = p(24); % max rate at which adulticide kills adult vectors
ga = p(25); % adulticide decay rate

cV = p(26);  % weight of cost of vectors in objective functional

d_l = ((rs*m_l*qs/muV)-muL-m_l)/c_l; % density-dependent death rate for larvae
d_h = (Lambda - mu_h)/c_h % density-dependent death rate for host

% eigenvalues are copied from  eig=solve(p, lambda)
sol1=double(subs(eigen_values(1)));
sol2=double(subs(eigen_values(2)));
sol3=double(subs(eigen_values(3)));
sol4=double(subs(eigen_values(4)));
sol5=double(subs(eigen_values(5)));
sol=[sol1 sol2 sol3 sol4 sol5];
 %%%%largest eigenvalue is R0
[subR0,max_index]=max(sol)
R0=eigen_values(max_index)