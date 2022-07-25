tic;

%%%% Write on 6/9/2020, 
%%%% To compute the basic reproduction number, R0, of a West Nile Virus
%%%% Model
%%%% Define variables Hi=infected hosts, Ei=infected eggs,
%%%% Li=infected larva, Ve=exposed mosquitoes, Vi=infected mosquitoes
%%%%%
clear
syms Hi1 Hi2 Hi3 Ei Li Ve Vi n1 n2 n3 n4 n5 n6 n7 n8 j1 j2 j3 j4
syms rs ri phi qs qi m_e m_l muL muV b c_l kl p_mh gl ga km1 km2 cV d_l Vs
syms p_hm1 p_hm2 p_hm3 omega1 omega2 omega3 p_hh1 p_hh2 p_hh3 g1 g2 g3 
syms gamma1 gamma2 gamma3 mu_h1 mu_h2 mu_h3 c_h1 c_h2 c_h3 d_h1 d_h2 d_h3 Lambda1 Lambda2 Lambda3

% DFE
% Hs1 = c_h1
% Hs2 = c_h2
% Hs3 = c_h3
% Ls = c_l
% Vs = c_l*m_l/muV
% Es = c_l*m_l*rs/(muV*m_e)

Vs = c_l*m_l/muV;

%%% Compute Jacobian
%%%%% F=new infections, V=transfer between compartments 
%%%% Only need to focus on infection compartments: [Hi1 Hi2 Hi3 Ei Li Ve Vi]

Ffun=[p_mh*b*Vi+p_hh2*omega2*Hi2, ri*Vi, 0, b*p_hm2*Vs*Hi2/c_h2, 0]

Vfun=[-gamma2*Hi2-g2*Hi2-d_h2*Hi2-mu_h2*Hi2, -m_e*Ei, m_e*qi*phi*Ei-muL*Li-m_l*Li-d_l*Li, -kl*Ve-muV*Ve, m_l*Li+kl*Ve-muV*Vi]

%%%% Compute the jacobian with respect to infection compartments: [Hi1 Ei Li Ve Vi]

FF=jacobian(Ffun, [Hi2 Ei Li Ve Vi]);
VV=jacobian(Vfun, [Hi2 Ei Li Ve Vi]);
FF
VV

%%% Find matrix F and V
%%%% Evaluate FF and VV at disease free equilibrium
%%%% Only need to set infection compartments [I, As, Is, F, X, Ms, V] as zeros

MatrixF=subs(FF, [Hi2 Ei Li Ve Vi], [0, 0, 0, 0, 0])
MatrixV=subs(VV, [Hi2 Ei Li Ve Vi], [0, 0, 0, 0, 0])

MatrixF
MatrixV

MatrixF = subs(MatrixF, [omega2*p_hh2 (b*c_l*m_l*p_hm2)/(c_h2*muV) b*p_mh ri], [j1 j2 j3 j4]);
MatrixV = subs(MatrixV, [-d_h2-g2-gamma2-mu_h2 -m_e m_e*phi*qi -d_l-m_l-muL m_l -kl-muV kl -muV], [n1 n2 n3 n4 n5 n6 n7 n8]);

%%%%% Compute F*V^{-1}
RR=-MatrixF*inv(MatrixV)

%%% Find eigenvalue of RR, largest eigenvalue=R0
syms lambda
pp=det(RR-lambda*eye(5));
pp_factors=factor(pp);
num_factors=length(pp_factors);
eigen_values=[];
for i=1:num_factors
eigen_values=[eigen_values; solve(pp_factors(i)==0, lambda,'MaxDegree', 5)];
end
%%%%% eigenvalues are given, the largest one is R0
%%%%% Comment out the following if you just need a symbolic expression of R0
%%%%% Numerical example, Largest eigenvalue is spectral radius (R0)

%%% set values for parameters. Note the argument is for plotting decay in
%%% larvicidal effect and does not impact this computation.
p = System_parametersRL(1,90); % Not quite sure why this goes to 90. 

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

Lambda1 = p(29);
Lambda2 = p(30);
Lambda3 = p(31);

mu_h1 = p(32); % natural death rate host group 1
mu_h2 = p(33); % natural death rate host group 2
mu_h3 = p(34); % natural death rate host group 3

% c_h1 = p(35); % carrying capacity host group 1
% c_h2 = p(36); % carrying capacity host group 2
% c_h3 = p(37); % carrying capacity host group 3
c_h1 = 0.0005
c_h2 = 0.0005
c_h3 = 0.0005

km1 = p(38); % max rate at which larvicide kills larvae
km2 = p(40); % max rate at which adulticide kills adult vectors

gl = p(39); % larvicide decay rate
ga = p(41); % adulticide decay rate

cV = p(42); % weight of cost of vectors in objective functional

d_l=((rs*m_l*qs/muV)-muL-m_l); % density-dependent death rate for larvae
d_h1 = (Lambda1 - mu_h1); % density-dependent death rate for host group 1
d_h2 = (Lambda2 - mu_h2); % density-dependent death rate for host group 2
d_h3 = (Lambda3 - mu_h3); % density-dependent death rate for host group 3

j1 = omega2*p_hh2
j2 = (b*c_l*m_l*p_hm2)/(c_h2*muV)
j3 = b*p_mh
j4 = ri

n1 = -d_h2-g2-gamma2-mu_h2
n2 = -m_e
n3 = m_e*phi*qi
n4 = d_l-m_l-muL
n5 = m_l
n6 = -kl-muV
n7 = kl
n8 = -muV

%eivenvalues are copied from  eig=solve(p, lambda)
sol1=double(subs(eigen_values(1)));
sol2=double(subs(eigen_values(2)));
sol3=double(subs(eigen_values(3)));
sol4=double(subs(eigen_values(4)));
sol5=double(subs(eigen_values(5)));

sol=[sol1 sol2 sol3 sol4 sol5]

% %%%% largest eigenvalue is R0
[subR0,max_index]=max(sol)
R0=eigen_values(max_index)

toc;