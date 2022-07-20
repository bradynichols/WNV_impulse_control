function [tt,x] = West_Nile_Model_run(larvicide_type)
% This function simulates a West Nile virus epidemic or mosquito population with no control. 
% Disease-free or disease-free-equilibrium conditionns are selected by the user prior to running on lines 56-60.

% Vector

%Es = x(1); % eggs laid by susceptible and exposed mothers
%Ei = x(2); % eggs laid by infected mothers
%Ls = x(3); % susceptible larvae
%Li = x(4); % infected larvae
%Vs = x(5); % susceptible vectors
%Ve = x(6); % exposed vectors
%Vi = x(7); % infected vectors

% Host

%Hs1 = x(8); % susceptible host group 1
%Hs2 = x(9); % susceptible host group 2
%Hs3 = x(10); % susceptible host group 3

%Hi1 = x(11); % infected host group 1
%Hi2 = x(12); % infected host group 2
%Hi3 = x(13); % infected host group 3

%Hr1 = x(14); % recovered host group 1
%Hr2 = x(15); % recovered host group 2
%Hr3 = x(16); % recovered host group 3

% Control

%Ul = x(17); % larvacide
%Ua = x(18); % adultacide

% Artifical State

%int_0^t{cV*cV*(Vi(s)+Vs(s)+Ve(s)}ds=x(19)

%NH1 = Hs1 + Hi1 + Hr1; % total hosts group 1
%NH2 = Hs2 + Hi2 + Hr2 % total hosts group 2
%NH3 = Hs3 + Hi3 + Hr3 % total hosts group 3

% The continuous state has 19 components. 
% These components track the value of the continous variables
% [Es,Ei,Ls,Li,Vs,Ve,Vi,Hs1,Hs2,Hs3,Hi1,Hi2,Hi3,Hr1,Hr2,Hr3,Ul,Ua,Int] % post treatment

% Model Parameters
% generate parameters
% Duration of simulation
Tf=2000;
p = System_parametersRL(larvicide_type,Tf);

% Initial conditions for discrete/continuous state variables
% ic(2)=infected egges, % ic(4)=infected larva, % ic(7)=infected mosquitoes

%Vector parameters

rs = p(1); % egg laying rate of S and E mosquitoes
m_e = p(6); % hatch rate 
m_l = p(7); % larval maturation rate
muV = p(9); % adult death rate
c_l = p(11); % mosquito carrying capacity (larval)
c_h1 = p(35); % carrying capacity host group 1
c_h2 = p(36); % carrying capacity host group 2
c_h3 = p(37); % carrying capacity host group 3

ic_V = m_l*c_l/muV;
ic_E = rs*m_l*c_l/(m_e*muV);

c_h1 = p(35);
c_h2 = p(36);
c_h3 = p(37);

% Initial conditions for discrete/continuous state variables
% healthy, summer ic. starts from DFE%
%ic = [ic_E;0;C;0;ic_V;0;0;NH;0;0;0;0;0];
%diseased
ic = [ic_E;0;c_l;0;ic_V;0;.01*ic_V;c_h1;c_h2;c_h3;0;0;0;0;0;0;0;0;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X01=ic;
%the derivative of the state with respect to its initial condition is intially one,
%so we have an identity matrix
X02=reshape(eye(length(X01)),length(X01)^2,1);
X0=[X01; X02];

f=@(t,x)West_Nile_ModelRL(t,x,p);

%solve the state equations
[tt,x]=ode23(f,[0,Tf],X0);

West_Nile_Model_plots(tt,x,0,[],[],Tf,[],[],larvicide_type);
end
