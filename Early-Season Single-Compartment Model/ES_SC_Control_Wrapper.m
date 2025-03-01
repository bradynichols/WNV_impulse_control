function [tau,ul,ua,X,J,J_comp,final_treatment_time,K] = ES_SC_Control_Wrapper(N,Ka,Tf,Obj_type,larvicide_type,ul,ua)
%This code optimizes the timing of the control and the control level
%in order to minimize a user selected objective functional. 
%Obj_type sets the objective functional. Vector control: 1, Disease
%control:2, Host preservation: 3. 
%larvicide_type: 1=long-lasting s-methorpene briquet, 2=VectoBac


p = ES_SC_Parameters(larvicide_type,Tf);
if Obj_type==1
f=@(t,x)ES_SC_Model_Definition(t,x,p);
end

if Obj_type==2
f=@(t,x)ES_SC_Model_Definition_Disease(t,x,p);
end

% This one hasn't existed in a while
if Obj_type==3
f=@(t,x)ES_SC_ModelRL_Hosts(t,x,p);
end

%This part of the code adjusts the constant K until final_treatment_time=Tf.
final_treatment_time_error_b=1;
%[tau,ul,ua,X,J,J_comp,final_treatment_time_a,X0,T] = West_Nile_Control(N,Ka,Tf);
if Obj_type==1
    [tau,ul,ua,X,J,J_comp,final_treatment_time_a,X0,T] = ES_SC_Control_Vector_Control(N,Ka,Tf,larvicide_type);
end
if Obj_type==2
    [tau,ul,ua,X,J,J_comp,final_treatment_time_a,X0,T] = ES_SC_Control_Disease_Control(N,Ka,Tf,larvicide_type,ul,ua);
end
if Obj_type==3
    [tau,ul,ua,X,J,J_comp,final_treatment_time_a,X0,T] = ES_SC_Control_Host_Preservation(N,Ka,Tf,larvicide_type);
end
final_treatment_time_error_a=Tf-final_treatment_time_a;
Kb=Ka+.5;
%Linear approximation of treatment_time as a function of K
%f(x)-f(xa)=(x-xa)*(f(xb)-f(xa))/(xb-xa))
%setting f(x)=0 yields x=xa-f(xa)*(xb-xa)/(f(xb)-f(xa));
while abs(final_treatment_time_error_b)>10^(-3)
    %[tau,ul,ua,X,J,J_comp,final_treatment_time_b,X0,T] = West_Nile_Control(N,Kb,Tf);
    if Obj_type==1
    [tau,ul,ua,X,J,J_comp,final_treatment_time_b,X0,T] = ES_SC_Control_Vector_Control(N,Kb,Tf,larvicide_type);
    end
    if Obj_type==2
    [tau,ul,ua,X,J,J_comp,final_treatment_time_b,X0,T] = ES_SC_Control_Disease_Control(N,Kb,Tf,larvicide_type,ul,ua);
    end
    if Obj_type==3
    [tau,ul,ua,X,J,J_comp,final_treatment_time_b,X0,T] = ES_SC_Control_Host_Preservation(N,Kb,Tf,larvicide_type);
    end
    final_treatment_time_error_b=Tf-final_treatment_time_b;
    Kc=Ka-final_treatment_time_error_a*(Ka-Kb)/(final_treatment_time_error_a-final_treatment_time_error_b);
    Ka=Kb;
    final_treatment_time_a=final_treatment_time_b;
    final_treatment_time_error_a=final_treatment_time_error_b;
    Kb=Kc;
end
K=Ka;
final_treatment_time=final_treatment_time_a;

%This part of the code plots the solution
x=[];
t=[];
for i=1:N
        if i==1
        x0=X0;
        else
            x0=X(:,i-1);
        end
        %if not time passes, we do not need to add any points to the
        %solution graph
    if tau(i)~=0
        %solve the state equations forward in time
        [tt,xx]=ode45(f,[0,tau(i)],x0);
        if i>1
        tt=tt+T(i-1);
        end
        t=cat(1,t,tt);
        x=cat(1,x,xx);  
    end
    
end
%control with optimal schedule
control_type=2;
ES_SC_Plots(t,x,control_type,Obj_type,N,Tf,J,J_comp, larvicide_type)
end



    
    
        
    
    
