function [tau,ul,ua,X,J,J_comp,final_treatment_time,X0,T] = ES_MC_Control_Vector_Control(N,K,Tf,larvicide_type)
%This code optimizes the timing of the control and the control level
%in order to minimize the density of vectors plus control cost.
%Initial conditions are set on line 93 as the disease-free equilibrium. 
% N is the maximal number of treatments
% K is the unknown value of the constant adjoint variable which determines sum(tau(i)).
% A wrapper file adjusts K to meet the constraint sum(tau(i))=Tf
%In the future we will write a wrapper function to find the value of K that satisfies the constraint.
% format short
% optimal control of WNV model
% State variables

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

%Artifical State
%int_0^t{cV*(Vi(s)+Vs(s)+Ve(s)}ds=x(19)

%NH = Hs+Hi+Hr;           %total hosts

%%% Impulse equation
%%%%% lower case u is for the impulse of adulticide and larvacide
% Ul(Ti^+) = Ul(Ti^+)+ul(i);
% Ua(Ti^+) = Ua(Ti^+)+ua(i);

%The discrete state has 19 components. These components track the value of the continous variables
%[Es,Ei,Ls,Li,Vs,Ve,Vi,Hs1,Hs2,Hs3,Hi1,Hi2,Hi3,Hr1,Hr2,Hr3,Ul,Ua,Int] post treatment, and the current time, T(i)=sum_{j=1}^i(tau(i)).
%This final variable is denoted X(i,14);
%The discrete state is denoted by X(i,:)=[X(i,1), . . ., X(i,19)]

%J=X(19,N)+c_e*(X(1,N)+X(2,N))+cl*sum(ul.^2)+ca*sum(ua.^2)+cT*sum(tau.^2);

%The Hamiltonian is -J+sum_{i=1}^{N}{<Y(i),(G(X(i-1),ul(i),ua(i),tau(i))-X(i))>} + <Y(0),X0-X(0)>
% + <K(N),Phi(X(N))>
% + sum_{i=1}^N{mu_t^+(i)*(tau_i-tau_max)-mu_t^-(i)*tau_i-mu_a^-(i)*ua-mu_l^-(i)*ul}
%Here the vector valued function G determines how the discrete state at the next time is determined by
%the discrete state at the previous time. Y(i) has one component for each state component.
% <Y,G> denotes the dot product.
%The second line of the Hamiltonian contains the inequality constraints. The coefficients here are <=0.
%In addition sum_{i=1}^N{mu_t^+(i)*(tau_i-tau_max)-mu_t^-(i)*tau_i-mu_a^-(i)*ua-mu_l^-(i)*ul}=0
%Phi is a vecotor of equality contraints on the final state. We have only one such constraint, X(N,14)-Tf=0.
                                              
%let dx/dt=g(x(t,x0)), where g is vector-valued and x depends on the initial state x0.
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(x(t,x0)) dt} for i=1, . . . 10
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(t,x0) dt} + ul(i) for i=11
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(t,x0) dt} + ua(i) for i=12
%G_i(x0,ul,ua,tau)=x0_i+int_0^tau{g_i(x(t,x0)) dt} for i=13
%G_i(x0,ul,ua,tau)=x0_i+tau(i) for i=19

test=-1;

delta=10^(-7);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Model Parameters
%generate parameters
p = ES_MC_Parameters(larvicide_type,Tf);

c_h1 = p(35); % carrying capacity host group 1
c_h2 = p(36); % carrying capacity host group 2
c_h3 = p(37); % carrying capacity host group 3

%Vector

rs = p(1); % egg laying rate of S and E mosquitoes
m_e = p(6); % hatch rate 
m_l = p(7); % larval maturation rate
muV = p(9); % adult death rate
c_l = p(11); % mosquito carrying capacity (larval)

ic_V=m_l*c_l/muV;
ic_E=rs*m_l*c_l/(m_e*muV);

% Steady state initial conditions for discrete/continuous state variables
%Disease-free equilibirum initial conditions.
ic = [ic_E;0;c_l;0;ic_V;0;0;c_h1+c_h2+c_h3;0;0;0;0;0;0]; % DFE; might not be right.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Duration of treatment
%Tf]
cl = p(43); % weight of cost of larvacide
ca = p(44); % weight of cost of adulticide
            % adulticide treatment is more expensive. Cite Mosquitoes and disease Illinois dept of public health
cT = p(45); % weight of cost of time
ce=p(46); % cost of eggs at the final time

%maximum time between controls, is set to Tf
Maxt=p(28);
%minimum time between controls
mint=p(29);

%state variables with underlying continuous dynamics
X01=ic(1:13);
%the derivative of the state with respect to its initial condition is intially one,
%so we have an identity matrix
X02=reshape(eye(length(X01)),length(X01)^2,1);
%Discrete states at each time will be columns
X0=[X01; X02];

%initial guess
%waiting times
%tau=(Tf/(N-1))*ones(N,1);
%tau(1)=0;
tau=ones(1,N);

T=zeros(size(tau));

ua=zeros(1,N);
ul=zeros(1,N);

%vector of discrete states after control and derivatives before treatment
X=zeros(length(X0),N);
%discrete states just before control is added
XX=zeros((length(X0)),N);
%derivative of discrete state with respect to previous state value (initial value).
DX=zeros(length(X01),length(X01),N);
%discrete adjoint variables
Y=zeros(N,length(ic));


f=@(t,x)ES_MC_Model_Definition(t,x,p);

final_times=[];
J_values=[];

count = 0;

while (test<0)
    
    oldul=ul;
    oldua=ua;   %%% control variable
    oldtau=tau;
    oldX=X;
    
%This loop finds the discrete states just before (XX) and just after (X) each treatment
for i=1:N
        if i==1
        x0=X0;
        T(i)=tau(i);
        else
            T(i)=T(i-1)+tau(i);
            x01=X(1:13,i-1);%The IC in the ODE solver must be a column
        %the derivative of the state with respect to its initial condition is intially one,
        %so we have an identity matrix
        x02=reshape(eye(length(x01)),length(x01)^2,1);
        x0=[x01; x02];
     
        end   
    %solve the state equations forward in time
    if tau(i)==0
        x=x0;  %states just before the next dose is added
    
        %state just before treatment
        XX(:,i) = x(:);
        
        X(:,i) = x(:);
    
        %states after the addition of the next dose
        X(11,i)=X(11,i)+ul(i);
        X(12,i)=X(12,i)+ua(i);
    else
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        sol=ode45(f,[0,tau(i)],x0,options);
        ll=length(sol.x);
        x=sol.y(:,ll);
        %states just before the next dose is added
    
        %state just before treatment
        XX(:,i) = x(:);
        
        X(:,i) = x(:);
    
        %states after the addition of the next dose
        X(11,i)=X(11,i)+ul(i);
        X(12,i)=X(12,i)+ua(i);
    end 
    
    DX(:,:,i)=reshape(x(length(X01)+1:length(X01)^2+length(X01)),length(X01),length(X01));
    %DX=[dx1/dx1(0), dx1/dx2(0),. . .dx1/dx13(0);
       % dx2/dx1(0), dx2/dx2(0), . . .dx2/dx13(0);
       %    .
       %    .
       %    .
       % dx13/dx1(0), dx13/dx2(0), . . .,dx13/dx13(0)];
end

%get the adjoint variables at steps 1,. . .,N.
Y(N,:)= [-ce -ce 0 0 0 0 0 0 0 0 0 0 -1 K];
for i=1:N-1
    j=N-i;
    Y(j,1:length(X01))=Y(j+1,1:length(X01))*DX(:,:,j);
end
%Note the 14th adjoint variable, which corresponds to the aritificial state
%variable is constant.
Y(:,length(ic))=K*ones(N,1);

for i=1:N
        
    ul(i) = Y(i,11)/(2*cl);
    ua(i) =Y(i,12)/(2*ca);

    ul(i)=max(0,min(1,ul(i)));
    ua(i)=max(0,min(1,ua(i)));
    
 %Here there is (was) an error, we should have the derivative of the state with
 %respect to time as well.
 %f=@(t,x)West_Nile_ModelRL2(t,x,p); gives derivatives of states at given
 %time and state.
 
 %dH/dtau(i)=(Y(i,1:13)*f(T(i),XX(i),p))+Y(i,14)-cT*2*tau(i)
 dxdt=f(T(i),XX(:,i));
 tau(i)=((Y(i,1:13)*dxdt(1:13))+Y(i,14))/(2*cT);

 %at least mint days between treatments
 if i>1
    tau(i)=min(Maxt,max(mint,tau(i)));
 end
 %the first treatment can happen at t=0
 if i==1
     tau(i)=min(Maxt,max(0,tau(i)));
 end
 
end
    
    J=X(13,N)+ce*(X(1,N)+X(2,N))+cl*sum(oldul.^2)+ca*sum(oldua.^2)+cT*sum(oldtau.^2);
    J_values=[J_values J];
    final_times=[final_times T(N)];

    %evaluate total relative error
    testua=delta*sum(abs(ua))-sum(abs(ua-oldua));
    testul=delta*sum(abs(ul))-sum(abs(ul-oldul));
    testtau=delta*sum(abs(tau))-sum(abs(tau-oldtau));
    testX=delta*sum(abs(X(1:13,:)))-sum(abs(X(1:13,:)-oldX(1:13,:)));
    
    
    %update control
    ul=(.9*oldul+.1*ul);
    ua=(.9*oldua+.1*ua);
    tau=(.9*oldtau+.1*tau);
        
    test=min([testua testul testtau testX]);
    count=count+1;
    
    final_treatment_time=sum(tau);
end

J_comp=X(13,N)+cl*sum(oldul.^2)+ca*sum(oldua.^2)+ce*(X(1,N)+X(2,N));

M=length(J_values);

figure
yyaxis left
plot(1:M,J_values,'*')

yyaxis right 
plot(1:M,final_times,'o')

yyaxis left 
ylabel('Objective functional value','FontSize', 20)

yyaxis right
ylabel('Final time','FontSize', 20)

xlabel('iteration','FontSize', 20)

set(gca,'fontsize',16)

file_name=sprintf('J_and_final_time_K=%.2f_T=%.2f_N=%.2f.eps',K,Tf,N);

figure_title=sprintf('K=%.2f T=%.2f N=%.2f',K,Tf,N);

title(figure_title)

%legend({'Objective functional value','Final time'},'Location','best','FontSize', 20)
exportgraphics(gcf,file_name)

hold off
    
end