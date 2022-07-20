function [p] = System_parametersRL(larvicide_type, Tf)

% Parameter values for simulations:
p = zeros(1,49);

% VECTOR PARAMETERS

p(1) = 150/8; % egg laying rate of S and E mosquitoes
p(2) = 100/8; % egg laying rate of I mosquitoes
              % infected mosquitoes lay few eggs per raft (83) than nonexposed mosquitoes (127).
              % the first raft contains about 50 fewer eggs.  If we include only the first gonotrophic cycle,
              % we might estimate the rate of egg-laying of unifected females as 150/8 and that of infected mosquitoes as 100/8.
              % Cite West Nile Virus Infection Decreases Fecundity of Culex tarsalis Females


p(3) = 0.003; % fraction of eggs infected about 3/1000
              % Cite importance of Vertical and Horizontal Transmission of West Nile Virus by Culex pipiens in the Northeastern United States.

p(4)=.56; % fraction of eggs laid by uninfected mosquitoes that hatch
p(5)=.43; % fraction of eggs laid by infected mosquitoes that hatch
          % 43% of eggs laid by infected mosquitoes hatch. 56% of eggs laid by nonexposed mosquitoes hatch.

p(6) = 1/2; % hatch rate  
            % eggs hatch in as little as a few days, as long as several months.
            % They require moisture. Cite MosquitoLifecycleFINALfor Culex Mosquitoes.

p(7) = 1/7; % larval maturation rate (1/larval lifespan) cite MosquitoLifecycleFINAL.
            % Here we combine the larval and pupal stages. 

p(8) = 0.16; % Daily death (of larvae) rate of 1 - 3%.                   

             % Larval stage survival varied between 70-90%. 
             % This would mean p(7)/(p(7)+p(8))=.7-.9.
             % For p(8) and p(7) as above this ratio is .82.
             % Cite Ruybal, Geographic variation in the response of 
             % Culex pipiens life history traits to temperature                   
                  
p(9) = 1/10.4; % adult death rate (1/adult lifespan)
               % female mosquitoes have a life expectancy of 3-7 days in the wild
               % Cite Maciel-de-Freitas (Calculating the survival rate and estimated population density of gravid Aedes aegypti
               % (Diptera, Culicidae) in Rio de Janeiro, Brazil)
               % the type that transmit West Nile lives about 10 days in the wild. 
               % Cite Jones, Rainfall Influences Survival of Culex pipiens (Diptera: Culicidae) in a Residential Neighborhood in the mid-Atlantic USA

p(10) = 1/5; % mosquito biting rate. It seems reasonable to assume a female mosquito will take a blood meal every 5 days.
             % Cite Ruybal, Geographic variation in the response of Culex 
             % pipiens life history traits to temperature
             %p(11) = 22; % per m^2 (in a larval habitat) see Fillinger
               
p(11)=.01; % mosquito larval carrying capacity. 
           %residential area; see paper for references MOSQUITO (LARVAE) CARRYING CAPACITY
           % Note that culex mosquitoes may be primarily breeding in drainage basins and ditches.

        
p(12) = 1/10; % disease progression in mosquitoes (1/latency period)
              % Cite West Nile Virus Infection Decreases Fecundity of Culex tarsalis Females

p(13) = 0.5; % mosquito-to-host transmission

% HOST PARAMETERS
% Host group 1: Birds with zero mortality from WNV.
% Host group 2: Birds with 0.01-0.49 mortality from WNV.
% Host group 3: Birds with 0.50-1.00 mortality from WNV.

p(14) = 0; % host-to-mosquito transmission host group 1
p(15) = 0; % host-to-mosquito transmission host group 2
p(16) = 0; % host-to-mosquito transmission host group 3

p(17) = 0; % direct transmission host group 1
p(18) = 0; % direct transmission host group 2
p(19) = 0; % direct transmission host group 3

p(20) = 0; % contact rate host group 1
p(21) = 0; % contact rate host group 2
p(22) = 0; % contact rate host group 3

p(23) = 0; % WNV recovery host group 1
p(24) = 0; % WNV recovery host group 2
p(25) = 0; % WNV recovery host group 3

p(26) = 0; % WNV death host group 1
p(27) = 0; % WNV death host group 2
p(28) = 0; % WNV death host group 3

p(29) = 0; % recruitment rate host group 1
p(30) = 0; % recruitment rate host group 2
p(31) = 0; % recruitment rate host group 3

p(32) = 0; % natural death rate host group 1
p(33) = 0; % natural death rate host group 2
p(34) = 0; % natural death rate host group 3

p(35) = 0; % carrying capacity host group 1
p(36) = 0; % carrying capacity host group 2
p(37) = 0; % carrying capacity host group 3

% urban bird densities varied between 
% 281 birds / 10 ha = 281 birds / 100,000 m^2 = 0.0028 birds/m^2 and 48 birds / 10 ha = .00048 birds / m^2 
% in Clergeau "Bird Abundance and Diversity Along an Urban Rural Gradient," a comparative
% study between two cities on different continents. In a residential area bird density was estimated as 15.3 birds / ha = .00153 birds/m^2. See Jones.
% In a manmade wetland, the density of individual species was as great as 7.15 / ha, although it is not clear how 
%  density was computed. Summing the numbers of individual terestrial birds and dividing by the land area gives
% a density of 200 / ha or .002. This is similar to the density in a residential area.
% "Population density of bird specides in a man-made wetland of peninsular Malaysia."
% The density of birds in low-lanf Hawaii was estimated at 65/ha = .0065 / m^2. 
% "Land Use and Larval Habitat Increase Aedes albopictus
% (Diptera: Culicidae) and Culex quinquefasciatus (Diptera:
% Culicidae) Abundance in Lowland Hawaii"

% INSECTICIDE PARAMETERS

% There are several types of larvicide. Decay rates vary.
% p(17)=rate at which larvicide decays
% p(18)=rate at which adulticide decays
% p(19)= max mate at which larvicide kills larva
% p(20)=max rate at which adulticide kills adults

% 1 CORRESPONDS TO METHOPRENE
if larvicide_type == 1
%%%%%%%Computation of p(17) and p(19) for s-methoprene briquets%%%%%%%%%%%%

%S-methoprene briquets reduce emergence of adult mosquitoes by almost 100% immediately following application. 
% Cite: Effectiveness of S-Methoprene Briquets and Application Method for Mosquito
% Control in Urban Road Gullies/Catch Basins/Gully Pots in a Mediterranean Climate: Implications for Ross R..
% If we assume max_ef effective at application, this gives the maximal larvicide-induced mortality p(19) satisfies 
% p(7)/(p(7)+p(8)+p(19))= (1-max_ef) p(7)/(p(7)+p(8)) 
% (1-max_ef)(p(7)+p(8)+p(19)) = (p(7)+p(8))
% (1-max_ef) p(19) = max_ef(p(7)+p(8))
% p(19) = max_ef*(p(7)+p(8))/(1-max_ef)


% After half_ef_day, the larvicide induced death rate ul*p(19) should satisfy ul(half_ef_day)*p(19)=p(8)+p(7).
% (That is, the number of adults emerging is reduced by 50%.) 
% So ul(half_ef_day)*max_ef/(1-max_ef) = 1, or ul(half_ef_day) = (1-max_ef)/(max_ef) = exp(-half_ef_day*p(17))
% p(17)=-log[(1-max_ef)/(max_ef)]/half_ef_day

% If we require the control be only min_ef=3% effective after min_ef_day days 
% p(7)/[ul(min_ef_day)*p(19)+p(7)+p(8)] = (1-min_ef)* p(7)/[p(7)+p(8)]
% [ul(min_ef_day)*p(19)]=min_ef*(p(7)+p(8))/(1-min_ef)
% exp(-min_ef_day*p(17))*max_ef*(p(8)+p(7))/(1-max_ef) = min_ef*(p(7)+p(8))/(1-min_ef)

% ((1-max_ef)/max_ef))^((min_ef_day-half_ef_day)/half_ef_day) = min_ef/(1-min_ef)
% (min_ef/(1-min_ef))^(half_ef_day/(min_ef_day-half_ef_day)) = (1-max_ef)/max_ef
% max_ef*(1+(min_ef/(1-min_ef))^(half_ef_day/(min_ef_day-half_ef_day))) = 1

%this is just a rough approximation of min_ef and min_ef days
min_ef=.03;

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
%%%%%%%Computation of p(17) and p(19) for vectobac%%%%%%%%%%%%

% Vectobac reduces emergence of adult mosquitoes by ~100% immediately following application. 
% Cite: Slow release formulations of Bacillus thuringiensis israelensis (AM 65-52) and spinosyns: effectiveness
% against the West Nile vector Culex pipiens in Saudi Arabia
% If we assume max_ef effective at application, this gives the maximal larvicide-induced mortality.
% p(19) satisfies p(7)/(p(7)+p(8)+p(19))= (1-max_ef) p(7)/(p(7)+p(8)) 
% (1-max_ef)(p(7)+p(8)+p(19)) = (p(7)+p(8))
% (1-max_ef) p(19) = max_ef(p(7)+p(8))
% p(19) = max_ef*(p(7)+p(8))/(1-max_ef)

% After mid_ef_day, the larvicide induced death rate ul*p(19) should satisfy (1-mid_eff)*[ul(mid_ef_day)*p(19)+p(8)+p(7)]=(p(8)+p(7)).
% (That is, the number of adults emerging is reduced by mid_eff*100%.) 
% So ul(mid_ef_day)*p(19)=[(mid_eff)/(1-mid_eff)](p(8)+p(7))
% ul(mid_ef_day)*max_ef*(p(7)+p(8))/(1-max_ef)=[(mid_eff)/(1-mid_eff)](p(8)+p(7))
% ul(mid_ef_day)=[(mid_eff)(1-max_ef)]/[(1-mid_eff)max_ef], 
% [(mid_eff)(1-max_ef)]/[(1-mid_eff)max_ef] = exp(-mid_ef_day*p(17))
% p(17)=-log([(mid_eff)(1-max_ef)]/[(1-mid_eff)max_ef])/mid_ef_day

% If we require the control be only min_ef=64% effective after min_ef_day days 
% p(7)/[ul(min_ef_day)*p(19)+p(7)+p(8)] = (1-min_ef)* p(7)/[p(7)+p(8)]
% [ul(min_ef_day)*p(19)]=min_ef*(p(7)+p(8))/(1-min_ef)
% exp(-min_ef_day*p(17))*max_ef*(p(8)+p(7))/(1-max_ef) = min_ef*(p(7)+p(8))/(1-min_ef)
% [(1-max_ef)/max_ef]^{(min_ef_day-mid_ef_day)/mid_ef_day}*[mid_eff/(1-mid_eff)]^(min_ef_day/mid_ef_day) = min_ef/(1-min_ef)
% [(1-max_ef)/max_ef] = [min_ef/(1-min_ef)]^{mid_ef_day/(min_ef_day-mid_ef_day)}*[(1-mid_eff)/mid_ef]^{min_ef_day/(min_ef_day-mid_ef_day)} 
% max_ef*(1+((min_ef/(1-min_ef))^(mid_ef_day/(min_ef_day-mid_ef_day)))*((1-mid_eff)/mid_ef)^(min_ef_day/(min_ef_day-mid_ef_day))) = 1

% The following values come from "Slow release formulations of Bacillus thuringiensis israelensis (AM 65-52) and spinosyns: effectiveness
% against the West Nile vector Culex pipiens in Saudi Arabia"
% min_ef=.64;
% min_ef_day=10;
% mid_ef_day=8;
% mid_ef=.95;

% The following values come from "Field Efficacy of Vectobac GR as a Mosquito Larvicide for
% the Control of Anopheline and Culicine Mosquitoes in Natural Habitats in Benin, West Africa"
min_ef=.22;
min_ef_day=42;
mid_ef_day=34;
mid_ef=.57;

max_ef=1/(1+((min_ef/(1-min_ef))^(mid_ef_day/(min_ef_day-mid_ef_day)))*((1-mid_ef)/mid_ef)^(min_ef_day/(min_ef_day-mid_ef_day)));
p(38) = max_ef*(p(8)+p(7))/(1-max_ef); % max rate at which larvicide kills larvae
p(39)= -log(mid_ef*(1-max_ef)/((1-mid_ef)*max_ef))/mid_ef_day; % rate at which larvicide decays

end

p(40)=-log(0.1)*p(18)*2; % max rate at which adulticide kills adult vectors
p(41)=24; % adulticide decay rate
per_remain_one_hour = exp(p(20)*0.5*exp(-p(18)/24)/p(18)-p(20)*0.5/p(18)); % Percent of adulticide remaining after one hour

% Adulticide max kill rate estimated from mortality of caged mosquitoes in a fallow agricultural field.
% Cite: Field trial efficiency of Anvil 10+10 and Biomist 31:66 in against . . . 
% When applied at 50% maximal level Biomist yields 88.6% mortality after 1 hour, and 95% mortality after 12 hours, when 
% mortality in the control group is 4%. Since we do not include a "treated"mosquito class, we use parameterize the model 
% so that adulticide-induced mortality is 90% after treatment. 
% Thus "successfully treated" mosquitoes are removed.
% let A = adulticide leve
% A=A_0*exp(-p(18)*t)
% Let V=mosquitoes, ignorming natural mortality, we want 90% of mosquitoes to eventually perish.
% dV/dt = -p(20)*A_0*exp(-p(18)*t)*V
% log(V(t)/V_0)=p(20)*A_0*exp(-p(18)*t)/p(18)-p(20)*A_0/p(18)
% taking limit as t goes to inifinity
% log(0.1)=-p(20)*A_0/p(18)
% p(20)=-log(0.1)*p(18)/A_0
% A_0=1/2.

% CONTROL PARAMETERS

p(42) = 5000; %the weight of the cost of the infected vectors.    
p(43) = 1; %weight of cost of larvacide
p(44) = 10; %weight of cost of adulticide
            %adulticide treatment is more expensive. Cite Mosquitoes and disease Illinois dept of public health
p(45) = 0.05; % weight of cost of time
p(46) = 5000; %cost of eggs at the final time
p(47) = -100000; %cost of hosts at the final time; value used in paper simulations

p(48)=Tf; % maximum time between controls
p(49)=1; % minimum time between controls

end