% Updated 07/28/2022
% For citations regarding vector parameters, please refer to comments in West Nile code written by Dr. Wandi Ding and Dr. Rachel Leander.

function [p] = ES_SC_Parameters(larvicide_type, Tf)
%parameter values for simulations
p=zeros(1,33);

% VECTOR PARAMETERS

p(1)=150/8; % egg laying rate of S and E mosquitoes
p(2)=100/8; % egg laying rate of I mosquitoes
p(3) = 0.003; % fraction of eggs infected

p(4)=.56; % fraction of eggs laid by uninfected mosquitoes that hatch
p(5)=.43; % fraction of eggs laid by infected mosquitoes that hatch

p(6) = 1/2; % egg maturation rate  
p(7) = 1/7; % larval maturation rate (1/larval lifespan)

p(8) = 0.16; % Daily death rate of 1-3% of larvae!                
             % Larval stage survival varied between 70-90%. 
             % This would mean p(7)/(p(7)+p(8))=.7-.9.
             % For p(8) and p(7) as above this ratio is .82.                 
                  
p(9) = 1/10.4; % adult death rate (1/adult lifespan)

p(10) = 1/5; % mosquito biting rate.

p(11)=.01; % mosquito larval carrying capacity.  
   
p(12) = 1/10;  %disease progression in mosquitoes (1/latency period)

p(13) = 0.5;   %mosquito-to-host transmission

% HOST PARAMETERS

p(14) = 0.3111312377; % host-to-mosquito transmission CHECK EXCEL

p(15) = 0.02960806276; % host-host direct transmission CHECK EXCEL

p(16) = 2.5768/7; % host-host contact rate CHECK EXCEL
                  % number in typical flock of those species that can transmit WNV horizontally/typical viremic period.

p(17) = 0.1121256665; % host recovery parameter CHECK EXCEL
                      % (1/infection period)*proportion recovered 

p(18) = 0.01672790387; % WNV-induced death rate

p(19) = 0.003952915781; % host recruitment rate

p(20) = 0.001993512447; % host natural death rate

p(21) = 0.0015; % host carrying capacity

% INSECTICIDE PARAMETERS

% 1 CORRESPONDS TO METHOPRENE
if larvicide_type == 1
%%%%%%%Computation of p(22) and p(23) for s-methoprene briquets%%%%%%%%%%

min_ef=.03; % Rough approximation of min_ef and min_ef days

% The product assessment has the product lasting up to 150 days and 69.5% effective at 120 days
% It does not last this long in the field. 
min_ef_day=150;
half_ef_day=100;
max_ef=1/(1+(min_ef/(1-min_ef))^(half_ef_day/(min_ef_day-half_ef_day)));

p(22) = max_ef*(p(8)+p(7))/(1-max_ef); % max rate at which larvicide kills larvae
p(23)= -log((1-max_ef)/max_ef)/half_ef_day; % rate at which larvicide decays

end

% 2 CORRESPONDS TO VECTOBAC
if larvicide_type == 2
%%%%%%%Computation of p(17) and p(19) for vectobac%%%%%%%%%%

min_ef = .22;
min_ef_day = 42;
mid_ef_day = 34;
mid_ef = .57;

max_ef=1/(1+((min_ef/(1-min_ef))^(mid_ef_day/(min_ef_day-mid_ef_day)))*((1-mid_ef)/mid_ef)^(min_ef_day/(min_ef_day-mid_ef_day)));
p(22) = max_ef*(p(8)+p(7))/(1-max_ef); % max rate at which larvicide kills larvae
p(23) = -log(mid_ef*(1-max_ef)/((1-mid_ef)*max_ef))/mid_ef_day; % rate at which larvicide decays

end

p(24) = -log(0.1)*p(18)*2; % max rate at which adulticide kills adult vectors
p(25) = 24; % adulticide decay rate
per_remain_one_hour = exp(p(24)*0.5*exp(-p(25)/24)/p(25)-p(24)*0.5/p(25)); % Percent of adulticide remaining after one hour

% CONTROL PARAMETERS

p(26) = 5000; % the weight of the cost of the infected vectors.    
p(27) = 1; % weight of cost of larvacide
p(28) = 10; % weight of cost of adulticide
         
p(29) = 0.05; % weight of cost of time
p(30) = 5000; % cost of eggs at the final time
p(31) = -100000; % cost of hosts at the final time; value used in paper simulations

p(32) = Tf; % maximum time between controls
p(33) = 1; % minimum time between controls

end