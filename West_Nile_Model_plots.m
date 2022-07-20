function [] = West_Nile_Model_plots(tt,x,control_type,Obj_type,N,Tf,J,J_comp,larvicide_type)

hold off
figure
hold on
plot(tt,x(:,1:2),'LineWidth',4)
legend('e_s', 'e_i','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('density (ind/m^2)', 'FontSize', 12);
if control_type==1
    file_name=sprintf('vector_control_eggs_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('eggs with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=sprintf('eggs_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('eggs with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
else 
    file_name=sprintf('eggs_no_control_T=%.2f.eps',Tf);
    figure_title=sprintf('egg density without control');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off
hold off
figure
hold on
plot(tt,x(:,3:7),'LineWidth',4)
legend('l_s', 'l_i', 'v_s', 'v_e', 'v_i','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('density (ind/m^2)', 'FontSize', 12);
if control_type==1
    file_name=sprintf('vector_control_vectors_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('vectors with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=si1 + Hr1; % total hosts group 1
NH2 = Hs2 + Hi2 + Hr2; % total hosts group 2
NH3 = Hs3 + Hi3 + Hr3; % total hosts group 3

% Vector ODEs % updated 7/20/2022

dEs = rs*(Vs+Ve)-m_e*Es;
dEi = ri*(Vi)-m_e*Ei;
dLs = m_e*qs*Es + m_e*qi*(1-phi)*Ei - muL*Ls - m_l*Ls - d_l*Ls*(Ls+Li) - km1*Ls*Ul;
dLi = m_e*qi*phi*Ei - muL*Li - m_l*Li - d_l*Li*(Li+Ls) - km1*Li*Ul;
dVs = m_l*Ls - b*p_hm1*Vs*Hi1/NH1 - b*p_hm2*Vs*Hi2/NH2 - b*p_hm3*Vs*Hi3/NH3 - muV*Vs - km2*Vs*Ua;
dVe = b*p_hm1*Vs*Hi1/NH1 + b*p_hm2*Vs*Hi2/NH2 + b*p_hm3*Vs*Hi3/NH3 - kl*Ve - muV*Ve - km2*Ve*Ua;
dVi = m_l*Li + kl*Ve - muV*Vi - km2*Vi*Ua;

% Host ODEs

dHs1 = Lambda1*Hs1 - b*p_mh*Vi*Hs1/NH1 - omega1*p_hh1*Hi1*Hs1/NH1 - d_h1*NH1*Hs1 - mu_h1*Hs1; % Updated 07/20/2022
dHs2 = Lambda2*Hs2 - b*p_mh*Vi*Hs2/NH2 - omega2*p_hh2*Hi2*Hs2/NH2 - d_h2*NH2*Hs2 - mu_h2*Hs2;
dHs3 = Lambda3*Hs3 - b*p_mh*Vi*Hs3/NH3 - omega3*p_hh3*Hi3*Hs3/NH3 - d_h3*NH3*Hs3 - mu_h3*Hs3;
dHi1 = b*p_mh*Vi*Hs1/NH1 + omega1*p_hh1*Hi1*Hs1/NH1 - (gamma1+g1)*Hi1 - d_h1*NH1*Hi1 - mu_h1*Hi1; % Updated 07/20/2022
dHi2 = b*p_mh*Vi*Hs2/NH2 + omega2*p_hh2*Hi2*Hs2/NH2 - (gamma2+g2)*Hi2 - d_h2*NH2*Hi2 - mu_h2*Hi2;
dHi3 = b*p_mh*Vi*Hs3/NH3 + omega3*p_hh3*Hi3*Hs3/NH3 - (gamma3+g3)*Hi3 - d_h3*NH3*Hi3 - mu_h3*Hi3;
dHr1 = g1*Hi1 - d_h1*NH1*Hr1 - mu_h1*Hr1; % Updated 07/20/2022
dHr2 = g2*Hi2 - d_h2*NH2*Hr2 - mu_h2*Hr2;
dHr3 = g3*Hi3 - d_h3*NH3*Hr3 - mu_h3*Hr3;
    
%Chemical ODEs

dUl = -gl*Ul;
printf('vectors_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('vectors with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
else 
    file_name=sprintf('vectors_no_control_T=%.2f.eps',Tf);
    figure_title=sprintf('vector density without control');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off
figure
hold on
plot(tt,x(:,8:16),'LineWidth',4)
legend('h_1s','h_2s','h_3s','h_1i','h_2i','h_3i','h_1r','h_2r','h_3r','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('density (ind/m^2)', 'FontSize', 12);
if control_type==1
    file_name=sprintf('hosts_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=sprintf('hosts_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
else 
    file_name=sprintf('hosts_no_control_T=%.2f.eps',Tf);
    figure_title=sprintf('host density without control');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off 
figure
hold on
plot(tt,x(:,8), tt, x(:,9), tt, x(:,10),'LineWidth',4)
legend('h_1s','h_2s','h_3s','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('density (ind/m^2)', 'FontSize', 12);
if control_type==1
    file_name=sprintf('hosts_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=sprintf('hosts_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
else 
    file_name=sprintf('hosts_no_control_T=%.2f.eps',Tf);
    figure_title=sprintf('host density without control');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off 
figure
hold on
%
plot(tt,x(:,11), tt, x(:,12), tt, x(:,13),'LineWidth',4)
legend('h_1i','h_2i','h_3i','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('density (ind/m^2)', 'FontSize', 12);
if control_type==1
    file_name=sprintf('hosts_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=sprintf('hosts_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
else 
    file_name=sprintf('hosts_no_control_T=%.2f.eps',Tf);
    figure_title=sprintf('host density without control');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off 
figure
hold on
%
plot(tt,x(:,14), tt, x(:,15), tt, x(:,16),'LineWidth',4)
legend('h_1r','h_2r','h_3r','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('density (ind/m^2)', 'FontSize', 12);
if control_type==1
    file_name=sprintf('hosts_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=sprintf('hosts_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('hosts with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
else 
    file_name=sprintf('hosts_no_control_T=%.2f.eps',Tf);
    figure_title=sprintf('host density without control');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off 
figure
hold on
%
plot(tt,x(:,11:12),'LineWidth',4)
legend('u_l', 'u_a','FontSize', 12);
xlabel('time (days)', 'FontSize', 12);
ylabel('normalized level', 'FontSize', 12);
if control_type==1
    file_name=sprintf('fixed-time_controls_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('fixed-time control (N = %.2f, T = %.2f, J = %.2f, J_c = %.2f, Obj fun = %.0f)',N,Tf,J,J_comp, Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
elseif control_type==2
    file_name=sprintf('variable_times_controls_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
    figure_title=sprintf('variable-time control (N = %.2f, T = %.2f, J = %.2f, J_c = %.2f, Obj fun = %.0f)',N,Tf,J,J_comp, Obj_type);
    title(figure_title)
    exportgraphics(gcf,file_name)
end
hold off

if control_type==0
p=System_parametersRL(larvicide_type,Tf);
%plot percent effectiveness of larvicide through time, according to
%p(7)/[ul(day)*p(19)+p(7)+p(8)] = (1-perc_ef)* p(7)/[p(7)+p(8)]
%perc_ef=1-[p(7)+p(8)]/[ul(day)*p(19)+p(7)+p(8)]
days=0:.1:150;
perc_ef=1-(p(7)+p(8))./(exp(-days*p(17))*p(19)+p(7)+p(8));
plot(days,perc_ef,'LineWidth',4)
xlabel('time (days)', 'FontSize', 12);
ylabel('instantaneous percent control of adult emergence', 'FontSize', 12);
file_name=sprintf('percent_control_adult_emergence.eps');
    figure_title=sprintf('Instantaneous percent control of adult emergence through time');
    title(figure_title)
    exportgraphics(gcf,file_name)
end
end


