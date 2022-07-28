% Updated 07-14-2022

%%%%Write on 6/9/2020, 
%%%%To compute the basic reproduction number, R0, of a West Nile Virus
%%%%Model
%%%%Define variables Hi=infected hosts, Ei=infected eggs,
%%%%Li=infected larva, Ve=exposed mosquitoes, Vi=infected mosquitoes
%%%%%
clear
syms Hi Ei Li Ve Vi
syms rs ri phi qs qi m_e m_l mu_h muL muV b c_l kl p_hm p_mh d_h g gl ga km1 km2 cV d_l gamma omega p_hh c_h % Included gamma, omega, p_hh
%%%Compute Jacobian
%%%%%F=new infections, V=transfer between compartments
%%%%Only need to focus on infection compartments: [Hi Ei Li Ve Vi]

Hs = c_h;
Ls = c_l;

Ffun=[p_mh*b*Vi+omega*p_hh*Hi, ri*Vi, 0, b*p_hm*m_l*Hi*c_l/(c_h*muV), 0]; 
Vfun=[-gamma*Hi-g*Hi-d_h*Hi*Hs-mu_h*Hi, -m_e*Ei, m_e*qi*phi*Ei-muL*Li-m_l*Li-d_l*Li*Ls, -kl*Ve-muV*Ve, m_l*Li+kl*Ve-muV*Vi]; 
%%%%Compute the jacobian with respect to infection compartments: [Hi Ei Li Ve Vi]
FF=jacobian(Ffun, [Hi Ei Li Ve Vi]);
VV=jacobian(Vfun, [Hi Ei Li Ve Vi]);
FF;
VV;
%%%Find matrix F and V
%%%%Evaluate FF and VV at disease free equilibrium
%%%%Only need to set infection compartments [I, As, Is, F, X, Ms, V] as zeros


MatrixF=subs(FF, [Hi Ei Li Ve Vi], [0, 0, 0, 0, 0]);
MatrixV=subs(VV, [Hi Ei Li Ve Vi], [0, 0, 0, 0, 0]);
%%%%%Compute F*V^{-1}
RR=-MatrixF*inv(MatrixV);
%%%Find eigenvalue of RR, largest eigenvalue=R0
syms lambda
pp=det(RR-lambda*eye(5));
pp_factors=factor(pp);
num_factors=length(pp_factors);
eigen_values=[];
for i=1:num_factors
eigen_values=[eigen_values; solve(pp_factors(i)==0, lambda,'MaxDegree', 5)];
end

p = ES_SC_Parameters(1,90);
Nh=p(22);
%Vector parameters
rs = p(1);            %intrinsic rate of increase of uninfected mosquitoes
ri = p(2);           %intrinsic rate of increase of infected mosquitoes
phi=p(3);               %fraction of eggs born to infected mothers that are infected
qs=p(4);                %fraction of eggs from uninfected mosquitoes that hatch
qi=p(5);                %fraction of eggs laid to infected mosquitoes that hatch
m_e = p(6);              %egg maturatation rate
m_l = p(7);             %larval maturation rate
muL = p(8);           %larval death rate
muV=p(9);             %adult death rate
b = p(10);             %mosquito biting rate
%Disease parameters
kl = p(12);           %disease progression (1/latency period)
p_mh = p(13);     %mosquito-to-host transmission
p_hm = p(14);     %host-to-mosquito transmission
gamma = p(15);           % WNV induced host mortality
g = p(16);            %host recovery rate
gl=p(17);
ga=p(18);
km1=p(19);
km2=p(20);
cV=p(21);                   %weight of cost of vectors in objective functional

Lambda = p(30); % Host birth rate
mu_h = p(31); % Host natural death rate
omega = p(32); %Host-to-Host Contact Rate
p_hh = p(33); %Host-to-Host Transmission Probability


d_l=((rs*m_l*qs/muV)-muL-m_l)/c_l; % density-dependent death rate for larvae
d_h = (Lambda - mu_h)/c_h; % density-dependent death rate for host


c_l = 0.01;             %mosquito carrying capacity
c_h = 0.001; % Host carrying capacity
% ending: 0.0265

reps = 50;

R0s = zeros([1 reps]);
xgraph = zeros([1 reps]);

xvalues = {};
Rvalues = {};
for x = 1:reps
    d_l=((rs*m_l*qs/muV)-muL-m_l)/c_l; % density-dependent death rate for larvae
    d_h = (Lambda - mu_h)/c_h; % density-dependent death rate for host
    %%%%eivenvalues are copied from  eig=solve(p, lambda)
    sol1=double(subs(eigen_values(1)));
    sol2=double(subs(eigen_values(2)));
    sol3=double(subs(eigen_values(3)));
    sol4=double(subs(eigen_values(4)));
    sol5=double(subs(eigen_values(5)));
    sol=[sol1 sol2 sol3 sol4 sol5];
     %%%%largest eigenvalue is R0
    [subR0,max_index]=max(sol);
%     subR0 = sol4;
%     subR0;
%     xvalues(end+1) = mat2cell(x);
    R0s(x) = sol4;
    xgraph(x) = c_l / c_h;
    c_h = c_h + 0.01;
    x
    %     R0=eigen_values(max_index);
end
% xvalues
R0s;
plot(xgraph,R0s,'-','LineWidth',4)
% set(gcf,'position',[0,0,1100,450])
% legend('H_S','H_I','H_R','FontSize', 12);
xlabel('c_L/c_H', 'FontSize', 12);
ylabel('Basic Reproductive Number R_0', 'FontSize', 12);
title('R_0 and Host/Vector Carrying Capacities');
% file_name=sprintf('R0_Graph_NewModel.png');
% exportgraphics(gcf,file_name);

% hold off
% figure
% hold on
%     
% 
% c_l = 0.01;             %mosquito carrying capacity
% c_h = 0.00005; % Host carrying capacity
% % ending: 0.0265
% 
% R0s = zeros([1 reps]);
% xgraph = zeros([1 reps]);
% 
% xvalues = {};
% Rvalues = {};
% for x = 1:reps
%     d_l=((rs*m_l*qs/muV)-muL-m_l)/c_l; % density-dependent death rate for larvae
%     d_h = (Lambda - mu_h)/c_h; % density-dependent death rate for host
%     %%%%eivenvalues are copied from  eig=solve(p, lambda)
%     sol1=double(subs(eigen_values(1)));
%     sol2=double(subs(eigen_values(2)));
%     sol3=double(subs(eigen_values(3)));
%     sol4=double(subs(eigen_values(4)));
%     sol5=double(subs(eigen_values(5)));
%     sol=[sol1 sol2 sol3 sol4 sol5];
%      %%%%largest eigenvalue is R0
%     [subR0,max_index]=max(sol);
% %     subR0;
% %     xvalues(end+1) = mat2cell(x);
%     R0s(x) = subR0;
%     xgraph(x) = c_h / c_l;
%     c_h = c_h + 0.0001;
%     x
%     %     R0=eigen_values(max_index);
% end
% % xvalues
% R0s;
% plot(xgraph,R0s,'-','LineWidth',4)
% % set(gcf,'position',[0,0,1100,450])
% % legend('H_S','H_I','H_R','FontSize', 12);
% xlabel('c_H/c_L', 'FontSize', 12);
% ylabel('Basic Reproductive Number R_0', 'FontSize', 12);
% title('R_0 and Host/Vector Carrying Capacities');
% % file_name=sprintf('Single_Compartment_Host_Density_4_Long.eps',Tf);
% % exportgraphics(gcf,file_name);
%     
%     
% 
