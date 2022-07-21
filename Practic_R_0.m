%%%% Write on 6/9/2020, 
%%%% To compute the basic reproduction number, R0, of a West Nile Virus
%%%% Model
%%%% Define variables Hi=infected hosts, Ei=infected eggs,
%%%% Li=infected larva, Ve=exposed mosquitoes, Vi=infected mosquitoes
%%%%%
clear
syms a b c d e
syms aa bb cc dd ee


%%% Compute Jacobian
%%%%% F=new infections, V=transfer between compartments 
%%%% Only need to focus on infection compartments: [Hi1 Hi2 Hi3 Ei Li Ve Vi]

% UPDATE F AND V VECTORS!!!

Ffun=[a/ee, 0, b*ee, c*d, 0]; % Updated to inlcude new interaction term % + omega*p_hh*Hi % 
Vfun=[aa*b, e*cc, 0, 0, 0]; % Updated to include gamma %-gamma*Hi% 

%%%% Compute the jacobian with respect to infection compartments: [Hi Ei Li Ve Vi]
% FF=jacobian(Ffun, [Hi1 Hi2 Hi3 Ei Li Ve Vi]);
% VV=jacobian(Vfun, [Hi1 Hi2 Hi3 Ei Li Ve Vi]);

FF=jacobian(Ffun, [a b c d e]);
VV=jacobian(Vfun, [a b c d e]);
FF
VV

%%% Find matrix F and V
%%%% Evaluate FF and VV at disease free equilibrium
%%%% Only need to set infection compartments [I, As, Is, F, X, Ms, V] as zeros

% MatrixF=subs(FF, [Hi1 Hi2 Hi3 Ei Li Ve Vi], [0, 0, 0, 0, 0, 0, 0])
% MatrixV=subs(VV, [Hi1 Hi2 Hi3 Ei Li Ve Vi], [0, 0, 0, 0, 0, 0, 0])

MatrixF=subs(FF, [a b c d e], [0, 0, 0, 0, 0])
MatrixV=subs(VV, [a b c d e], [0, 0, 0, 0, 0])

%%%%% Compute F*V^{-1}
RR=-MatrixF*inv(MatrixV)

%%% Find eigenvalue of RR, largest eigenvalue=R0
syms lambda
pp=det(RR-lambda*eye(5)); % Change this from 7 to 5 if something isn't working correctly.
pp_factors=factor(pp);
num_factors=length(pp_factors);
eigen_values=[];
for i=1:num_factors
eigen_values=[eigen_values; solve(pp_factors(i)==0, lambda,'MaxDegree', 5)]; % Change this from 7 to 5 if something isn't working correctly.
end
%%%%% eigenvalues are given, the largest one is R0
%%%%% Comment out the following if you just need a symbolic expression of R0
%%%%% Numerical example, Largest eigenvalue is spectral radius (R0)

%%% set values for parameters. Note the argument is for plotting decay in
%%% larvicidal effect and does not impact this computation.
p = System_parametersRL(1,90); % Not quite sure why this goes to 90. 

% Model Parameters
   
aa = p(14); % host-to-mosquito transmission host group 1
bb = p(15); % host-to-mosquito transmission host group 2
cc = p(16); % host-to-mosquito transmission host group 3

dd = p(17); % direct transmission rate host group 1
ee = p(18); % direct transmission rate host group 2


%% eivenvalues are copied from  eig=solve(p, lambda)
sol1=double(subs(eigen_values(1)));
sol2=double(subs(eigen_values(2)));
sol3=double(subs(eigen_values(3)));
sol4=double(subs(eigen_values(4)));
sol5=double(subs(eigen_values(5)));

% sol1=(subs(eigen_values(1)));
% sol2=(subs(eigen_values(2)));
% sol3=(subs(eigen_values(3)));
% sol4=(subs(eigen_values(4)));
% sol5=(subs(eigen_values(5)));
% sol6=(subs(eigen_values(6))); 
% sol7=(subs(eigen_values(7)));

% sol=[sol1 sol2 sol3 sol4 sol5 sol6 sol7]

sol=[sol1 sol2 sol3 sol4 sol5]

% %%%% largest eigenvalue is R0
[subR0,max_index]=max(sol)
R0=eigen_values(max_index)