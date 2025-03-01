function [] = ES_SC_Plots(tt,x,control_type,Obj_type,N,Tf,J,J_comp,larvicide_type)

hold off
figure
hold on
plot(tt,x(:,1:2),'LineWidth',4)
legend('E_S', 'E_I','FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
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

% Adult Vector Density
hold off
hold off
figure
hold on
plot(tt,x(:,5:7),'LineWidth',4)
legend('V_S', 'V_E', 'V_I','FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
title('Single-Compartment Vector Density');
file_name=sprintf('Single_Compartment_Vector_Density_4.eps',Tf);
exportgraphics(gcf,file_name);

% if control_type==1
%     file_name=sprintf('vector_control_vectors_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('vectors with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% elseif control_type==2
%     file_name=sprintf('vectors_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('vectors with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% else 
%     file_name=sprintf('vectors_no_control_T=%.2f.eps',Tf);
%     figure_title=sprintf('Single-Compartment Vector Density');
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% end

% Host Density

hold off
figure
hold on
plot(tt,x(:,8:10),'-','LineWidth',4)
% set(gcf,'position',[0,0,1100,450])
legend('H_S','H_I','H_R','FontSize', 12);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Density (Individuals/m^2)', 'FontSize', 12);
title('Single-Compartment Host Density');
file_name=sprintf('Single_Compartment_Host_Density_4_Long.eps',Tf);
exportgraphics(gcf,file_name);
% 
% if control_type==1
%     file_name=sprintf('hosts_fixed_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('hosts with fixed-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% elseif control_type==2
%     file_name=sprintf('hosts_variable_times_T=%.2f_N=%.2f_Obj_fun=%.0f.eps',Tf,N,Obj_type);
%     figure_title=sprintf('hosts with variable-time control (N = %.2f, T = %.2f, Obj fun = %.0f)',N,Tf,Obj_type);
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% else 
%     file_name=sprintf('hosts_no_control_T=%.2f.eps',Tf);
%     figure_title=sprintf('Single-Compartment Host Density');
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% end
% hold off 
% figure
% hold on
% plot(tt,x(:,11:12),'LineWidth',4)
% legend('u_l', 'u_a','FontSize', 12);
% xlabel('Time (Days)', 'FontSize', 12);
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
% hold off
% 
% if control_type==0
% p=System_parametersRL(larvicide_type,Tf);
% %plot percent effectiveness of larvicide through time, according to
% %p(7)/[ul(day)*p(19)+p(7)+p(8)] = (1-perc_ef)* p(7)/[p(7)+p(8)]
% %perc_ef=1-[p(7)+p(8)]/[ul(day)*p(19)+p(7)+p(8)]
% days=0:.1:150;
% perc_ef=1-(p(7)+p(8))./(exp(-days*p(17))*p(19)+p(7)+p(8));
% plot(days,perc_ef,'LineWidth',4)
% xlabel('Time (Days)', 'FontSize', 12);
% ylabel('instantaneous percent control of adult emergence', 'FontSize', 12);
% file_name=sprintf('percent_control_adult_emergence.eps');
%     figure_title=sprintf('Instantaneous percent control of adult emergence through time');
%     title(figure_title)
%     exportgraphics(gcf,file_name)
% end
% end
% 
% 
