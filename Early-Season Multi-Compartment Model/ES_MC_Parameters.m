% Updated 07/28/2022
% For citations regarding vector parameters, please refer to comments in West Nile code written by Dr. Wandi Ding and Dr. Rachel Leander.

function [p] = ES_MC_Parameters(larvicide_type, Tf)

% Parameter values for simulations:
p = zeros(1,49);

% VECTOR PARAMETERS

p(1) = 150/8; % egg laying rate of S and E mosquitoes
p(2) = 100/8; % egg laying rate of I mosquitoes

p(3) = 0.003; % fraction of eggs infected

p(4)=.56; % fraction of eggs laid by uninfected mosquitoes that hatch
p(5)=.43; % fraction of eggs laid by infected mosquitoes that hatch

p(6) = 1/2; % hatch rate  

p(7) = 1/7; % larval maturation rate (1/larval lifespan) cite MosquitoLifecycleFINAL.
            % Here we combine the larval and pupal stages. 

p(8) = 0.16; % Daily death (of larvae) rate of 1 - 3%.                   
             % Larval stage survival varied between 70-90%. 
             % This would mean p(7)/(p(7)+p(8))=.7-.9.
             % For p(8) and p(7) as above this ratio is .82.                  
                  
p(9) = 1/10.4; % adult death rate (1/adult lifespan)
               % female mosquitoes have a life expectancy of 3-7 days in the wild

p(10) = 1/5; % mosquito biting rate. It seems reasonable to assume a female mosquito will take a blood meal every 5 days.
                  
p(11)=.01; % mosquito larval carrying capacity. 

p(12) = 1/10; % disease progression in mosquitoes (1/latency period)

p(13) = 0.5; % mosquito-to-host transmission

% HOST PARAMETERS
% Host group 1: Birds with zero mortality from WNV.
% Host group 2: Birds with 0.01-0.49 mortality from WNV.
% Host group 3: Birds with 0.50-1.00 mortality from WNV.

p(14) = 0.2146393254; % host-to-mosquito transmission host group 1
p(15) = 0.376731036; % host-to-mosquito transmission host group 2
p(16) = 0.4066812182; % host-to-mosquito transmission host group 3

p(17) = 0; % direct transmission host group 1
p(18) = 0; % direct transmission host group 2
p(19) = 0.1590100415; % direct transmission host group 3

p(20) = 0; % contact rate host group 1
p(21) = 0; % contact rate host group 2
p(22) = 9.202857143/7; % contact rate host group 3

p(23) = 0.1429; % WNV recovery host group 1
p(24) = 0.1018258535; % WNV recovery host group 2
p(25) = 0.0602776217; % WNV recovery host group 3

p(26) = 0; % WNV death host group 1
p(27) = 0.03519053473; % WNV death host group 2
p(28) = 0.09875217396; % WNV death host group 3

p(29) = 0.004396413222; % recruitment rate host group 1
p(30) = 0.003783533594; % recruitment rate host group 2
p(31) = 0.003304357023; % recruitment rate host group 3

p(32) = 0.001989171602; % natural death rate host group 1
p(33) = 0.001870952684; % natural death rate host group 2
p(34) = 0.002269399628; % natural death rate host group 3

p(35) = 0.0005; % carrying capacity host group 1
p(36) = 0.0005; % carrying capacity host group 2
p(37) = 0.0005; % carrying capacity host group 3

% INSECTICIDE PARAMETERS

% 1 CORRESPONDS TO METHOPRENE
if larvicide_type == 1
%%%%%%%Computation of p(38) and p(39) for s-methoprene briquets%%%%%%%%%%%%

min_ef=.03; % Rough approximation of min_ef and min_ef days

% The product assessment has the product lasting up to 150 days and 69.5% effective at 120 days
% It does not last this long in the field. 
min_ef_day=150;
half_ef_day=100;
max_ef=1/(1+(min_ef/(1-min_ef))^(half_ef_day/(min_ef_day-half_ef_day)));

p(38) = max_ef*(p(8)+p(7))/(1-max_ef); % max rate at which larvicide kills larvae
p(39)= -log((1-max_ef)/max_ef)/half_ef_day; % rate at which larvicide decays

end

% 2 CORRESPONDS TO VECTOBAC
if larvicide_type == 2
%%%%%%%Computation of p(38) and p(39) for vectobac%%%%%%%

min_ef = .22;
min_ef_day = 42;
mid_ef_day = 34;
mid_ef = .57;

max_ef=1/(1+((min_ef/(1-min_ef))^(mid_ef_day/(min_ef_day-mid_ef_day)))*((1-mid_ef)/mid_ef)^(min_ef_day/(min_ef_day-mid_ef_day)));
p(38) = max_ef*(p(8)+p(7))/(1-max_ef); % max rate at which larvicide kills larvae
p(39)= -log(mid_ef*(1-max_ef)/((1-mid_ef)*max_ef))/mid_ef_day; % rate at which larvicide decays

end

p(40) = -log(0.1)*p(18)*2; % max rate at which adulticide kills adult vectors
p(41) = 24; % adulticide decay rate
per_remain_one_hour = exp(p(40)*0.5*exp(-p(41)/24)/p(41)-p(40)*0.5/p(41)); % Percent of adulticide remaining after one hour

% CONTROL PARAMETERS

p(42) = 5000; % the weight of the cost of the infected vectors.    
p(43) = 1; % weight of cost of larvacide
p(44) = 10; % weight of cost of adulticide

p(45) = 0.05; % weight of cost of time
p(46) = 5000; % cost of eggs at the final time
p(47) = -100000; % cost of hosts at the final time; value used in paper simulations

p(48) = Tf; % maximum time between controls
p(49) = 1; % minimum time between controls

end