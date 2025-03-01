function [] = ES_MC_Plots(tt,x,control_type,Obj_type,N,Tf,J,J_comp,larvicide_type)

% Plots Adult Vector Populations

hold off
figure
hold on
plot(tt,x(:,5:7),'Linewidth', 4)
legend('V_S', 'V_E', 'V_I', 'FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
title('Multiple Compartment Vector Density');
file_name=sprintf('Multiple_Compartment_Model_Vector_Density.eps',Tf);
exportgraphics(gcf,file_name);

% Plots Host 1 Populations

hold off
figure
hold on
plot(tt,x(:,8), tt,x(:,11), tt,x(:,14), 'LineWidth', 4)
legend('H_{S1}', 'H_{I1}', 'H_{R1}', 'FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
title('Multiple Compartment Host 1 Density');
file_name=sprintf('Multiple_Compartment_Model_Host1_Density.eps',Tf);
exportgraphics(gcf,file_name);

% Plots Host 2 Populations

hold off
figure
hold on
plot(tt,x(:,9), tt,x(:,12), tt,x(:,15), 'LineWidth', 4)
legend('H_{S2}', 'H_{I2}', 'H_{R2}', 'FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
title('Multiple Compartment Host 2 Density');
file_name=sprintf('Multiple_Compartment_Model_Host2_Density.eps',Tf);
exportgraphics(gcf,file_name);

% Plots Host 3 Populations

hold off
figure
hold on
plot(tt,x(:,10), tt,x(:,13), tt,x(:,16), 'LineWidth', 4)
legend('H_{S3}', 'H_{I3}', 'H_{R3}', 'FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
title('Multiple Compartment Host 3 Density');
file_name=sprintf('Multiple_Compartment_Model_Host3_Density.eps',Tf);
exportgraphics(gcf,file_name);

% Control, Not Needed for Initial Bird Demographic Model

% hold off 
% figure
% hold on
% %
% plot(tt,x(:,11:12),'LineWidth',4)
% legend('u_l', 'u_a','FontSize', 12);
% xlabel('time (days)', 'FontSize', 12);
% ylabel('normalized level', 'FontSize', 12);
% if control_type==1
%     file_name=sprintf('fixed-time_controls_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('fixed-time control (N = %.2f, T = %.2f, J = %.2f, J_c = %.2f, Obj fun = %.0f)',N,Tf,J,J_comp, Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% elseif control_type==2
%     file_name=sprintf('variable_times_controls_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('variable-time control (N = %.2f, T = %.2f, J = %.2f, J_c = %.2f, Obj fun = %.0f)',N,Tf,J,J_comp, Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% end

% if control_type==0
% p=System_parametersRL(larvicide_type,Tf);
% %plot percent effectiveness of larvicide through time, according to
% %p(7)/[ul(day)*p(19)+p(7)+p(8)] = (1-perc_ef)* p(7)/[p(7)+p(8)]
% %perc_ef=1-[p(7)+p(8)]/[ul(day)*p(19)+p(7)+p(8)]
% days=0:.1:150;
% perc_ef=1-(p(7)+p(8))./(exp(-days*p(17))*p(19)+p(7)+p(8));
% plot(days,perc_ef,'LineWidth',4)
% xlabel('time (days)', 'FontSize', 12);
% ylabel('instantaneous percent control of adult emergence', 'FontSize', 12);
% file_name=sprintf('percent_control_adult_emergence.eps');
%     figure_title=sprintf('Instantaneous percent control of adult emergence through time');
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% end
% end

% hold off
% figure
% hold on
% plot(tt,x(:,1:2),'LineWidth',4)
% legend('e_s', 'e_i','FontSize', 12);
% xlabel('time (days)', 'FontSize', 12);
% ylabel('density (ind/m^2)', 'FontSize', 12);
% if control_type==1
%     file_name=sprintf('vector_control_eggs_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('eggs with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% elseif control_type==2
%     file_name=sprintf('eggs_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('eggs with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% else 
%     file_name=sprintf('eggs_no_control_T=%.2f.eps',Tf);
%     figure_title=sprintf('egg density without control');
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% end
% 
% 
