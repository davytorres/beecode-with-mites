% Copyright (c) 2020, David Torres and Nicholas Torres
% This MATLAB code simulates a honey bee with mites

clear all
close all
eps = 1.e-15;
iplot = 1;

% Transient bee model
% Egg(i,j) represents the egg population that is j days old at time i(dt)
% B(i,j) represents the brood population that is j days old at time i(dt)

% P(i,j) represents the pupae population that is j days old at time i(dt)
% P_infected(i,j) represents the infected pupae population that is j days
% old at time i(dt)

% Drone(i,j) represents the pupae population that is j days old at time i(dt)
% Drone_infected(i,j) represents the infected pupae population that is j days
% old at time i(dt)

% H(i,j) represents the hive population that is j days old at time i(dt)
% H_infected(i,j) represents the infected hive population that is j days
% old at time i(dt)

% F(i,j) represents the foraging population that is j days old at time i(dt)
% F_infected(i,j) represents the infected forager population that is j days
% old at time i(dt)

% L is the daily egg laying rate of the queen bee, Lmax is the maximum egg-laying rate 

Lmax = 1600; % Maximum Egg-laying rate

% Number of days to run simulation 
ndays = 2*365;


dt = .1; % Time step in days
% n: Number of cycles
n = ceil(ndays/dt);

time = zeros(n+1,1);
% Day of year to start simulation where 
% 1: January 1st and 365: December 31st
time(1) = 32.d0; % February 1st

summer_begin = 65.d0; % Day when summer begins
summer_end = 260.d0; % Day when summer ends

initial_mite = 200.d0; % Initial mite population
%initial_mite = 0.0;

% Initial total number of bees
total_bees = 8000.;

% Grooming rate
grooming_rate = 0.05;
transmission_rate = 0.25;

% Parameter which regulates survival rate of larvae in the absence of
% sufficient hive bees
expl = 5.d0;

% ndays_transfer: Number of days in older hive population which can be age accelerated or
% deaccelerated due to pheromones
ndays_transfer = 6;

ijump = floor(1.0/dt);

% Grams of food gathered in one day by one forager (and processed by hive)
forager_rate = .1d0;

ratio_increase = 1.d0; % Increase in number of nurse bees to tend to brood if infected

% Increase mortality of infected hive, drone, and forager castes
increase_mortality_infected = 2.d0;

% Increase mortality of infected pupae. Survival is further reduced if the
% capped cells are multiply infested.
increase_mortality_infected_pupae = 2.d0;


% Number of days spent in each class 
negg_levels= 3;
nlarvae_levels = 5; % Worker larvae spend 5 days uncapped
nlarvae_drone_levels = 7; % Drone larvae spend 7 days uncapped
npupae_levels = 12; % Worker pupae spend 12 days capped
npupae_drone_levels = 14; % Drone pupae spend 14 days capped
nnurse_levels = 10; % Number of days spent as nurse bee/hive bee
nhive_levels = 21; % Number of bees spent as a hive bee
nforager_levels = 14; % Schmickl has a mean of 11 % Khoury uses mean of 7
ndrone_levels = 11; % Lifespan of drone bees

nmite_levels = 27; % Life expectancy of mites during summer season
nfertile_mites = 14; % Maximum age at which mites can reproduce by invading cells
mite_death_rate = .006d0;  % Summer mite daily death rate
mite_death_rate_winter = .002d0; % Winter mite daily death rate
%mite_death_rate = .008 + .01*double(jjj-1)/double(numv-1);

% Survival values of multiply infested mite cells come from S. Martin/Ecological Modelling 109 (1998)
% p. 274
reprod_drone = 2.91;   % Reproductive rate in capped drone cells
reprod_worker = 1.01;  % Reproductive rate in capped worker cells
xi = 6.0; % Maximum number of mites per infected hive or drone bee
% DeGrandi-Hoffman values
reprod_drone = 2.6;   % Reproductive rate in capped drone cells
%reprod_drone = 1.5;
reprod_worker = 1.5;  % Reproductive rate in capped worker cells
% Rosenkranz 2010 Biology and control of Varroa destructor
%rep_rate_worker = .5d0*(1.3+1.45);
%rep_rate_drone = .5d0*(2.2+2.6);
facd = 1.0;
facw = 1.0;

% Survival rate of mites is reduced if more than one mite invades a capped
% drone cell
max_levels_drone = 4;
survival_multiple_drones(1) = reprod_drone*facd;
survival_multiple_drones(2) = reprod_drone*.84*facd;
survival_multiple_drones(3) = reprod_drone*.65*facd;
survival_multiple_drones(4) = reprod_drone*.66*facd;

% Survival rate of mites is reduced if more than one mite invades a capped
% worker cell
max_levels_worker = 4;
survival_multiple_workers(1) = reprod_worker*facw;
survival_multiple_workers(2) = reprod_worker*.91*facw;
survival_multiple_workers(3) = reprod_worker*.86*facw;
survival_multiple_workers(4) = reprod_worker*.60*facw;

% Survival rates at different levels of worker honey bees
% Degrandi-Hoffman, 2004, International Journal of Acarology
% A mathematical model of Varroa mite ...
reduce_level(1) = 1.0;
reduce_level(2) = .991258;  % .90 = .991258^12 reduction overall reduction over 12 days spent as pupue
reduce_level(3) = .9815765; % .80 reduction overall reduction over 12 days spent as pupue
reduce_level(4) = .9583245; % .60 reduction overall reduction over 12 days spent as pupue
reduce_level(5) = .8744853; % .20 reduction overall reduction over 12 days spent as pupue
reduce_level(6) = .8744853; % .20 reduction overall reduction over 12 days spent as pupue

reduce_level_drone(1) = 1.0;
reduce_level_drone(2) = .9925025; % .90 reduction overall reduction over 12 days spent as pupue
reduce_level_drone(3) = .9841875; % .80 reduction overall reduction over 12 days spent as pupue
reduce_level_drone(4) = .9641701; % .60 reduction overall reduction over 12 days spent as pupue
reduce_level_drone(5) = .8914019; % .20 reduction overall reduction over 12 days spent as pupue
reduce_level_drone(6) = .8914019; % .20 reduction overall reduction over 12 days spent as pupue

% Set reduction values to 1 for sensitivity studies
for i = 1:6
    reduce_level(i) = 1.0;
end
%%
for i = 1:6
    reduce_level_drone(i) = 1.0;
end

% Food consumption
gamma_hive = .007; % grams/(bee*day)
gamma_forager = gamma_hive; % grams/(bee*day)
gamma_drone = 2.0*gamma_hive; % grams/(bee*day)
gamma_larvae = .018; % grams/day % .0326 Russell
gamma_larvae_death = 0.005; % grams/day for cannibalization

% No food survival rates
nofood_surv_forager = 2.d0/3.d0;
nofood_surv_hive = 1.d0/2.d0;
nofood_surv_brood = 1.d0/5.d0;
    
% Daily Death rates (Schmickl) based on data from Sakagami and Fukuda (1968)
egg_death_rate = .03d0;
larvae_death_rate = .01d0;
pupae_death_rate = .001d0;
hive_death_rate = .015d0;
% According to Sumpter and Martin (2004) the daily death rate in the winter
% is 1/190  = .0053
hive_death_rate_winter = 1./190.;
forager_death_rate = .045d0;

pupae_death_rate_infected = increase_mortality_infected_pupae*pupae_death_rate;
hive_death_rate_infected = increase_mortality_infected*hive_death_rate;
hive_death_rate_winter_infected = increase_mortality_infected*hive_death_rate_winter;
forager_death_rate_infected = increase_mortality_infected*forager_death_rate;

% Weight of larvae at different ages - used for cannibalism
weight(1) = .001d0;
weight(2) = .006d0;
weight(3) =  .020d0;
weight(4) = .080d0;
weight(5) = .150d0;


% Death rate of foragers and drones after their last day
forager_death_rate_end = 1.d0;
drone_death_rate_end = 1.d0;

% Daily survival rates
% Egg survival rate
s_e = 1.d0 - egg_death_rate;
% Brood survival rate
s_b = 1.d0 - larvae_death_rate;
% Pupae survival rate
s_p = 1.d0 - pupae_death_rate;
% Infected pupae survival rate
s_pi = 1.d0 - pupae_death_rate_infected;
% s_h calculated in time loop
s_f  = 1.d0 - forager_death_rate;
% Infected forager survival rate
s_fi = 1.d0 - forager_death_rate_infected;

s_mite = 1.d0 - mite_death_rate;  % Mite daily survival rate

%s_pi = .8744853; % Survival rate of pupae with Deformed Wing Virus
%s_pi = 0.0;
% Healthy Hive to Brood ratio
ratio_H_L_eq = 2.d0;

% Healthy Nurse (average 10 days) to Brood (average 21 days) ratio
%ratio_N_L_eq = .5d0*ratio_H_L_eq;


%ratio_N_L_eq_infected = ratio_increase*ratio_N_L_eq;
% There need to be ratio_increase as many infected hive bees to be effective in caring
% for brood
ratio_H_L_eq_infected = ratio_increase*ratio_H_L_eq;

% Healthy Hive to Forager ratio (phermone impact parameter)
ratio_H_F_eq = 2.3d0;

% Establishing initial populations of egg, larvae, pupae, drone, hive and
% foragers
initial_egg = .05*double(total_bees)/double(negg_levels);
initial_larvae = .10*double(total_bees)/double(nlarvae_levels);
initial_pupae = .98*.21*double(total_bees)/double(npupae_levels);
initial_drone = .02*.21*double(total_bees)/double(ndrone_levels);
initial_hive = .45*double(total_bees)/double(nhive_levels);
initial_forager = .19*double(total_bees)/double(nforager_levels);
% .05+.10+.21+.45+.19 = 1.0
fraction_hive_infected = double(initial_mite)/double(total_bees);
fraction_hive_healthy = 1. - fraction_hive_infected;
initial_egg = .0*double(total_bees)/double(negg_levels);
initial_larvae = .0*double(total_bees)/double(nlarvae_levels);
initial_pupae = .0*.21*double(total_bees)/double(npupae_levels);
initial_drone = .0*.21*double(total_bees)/double(ndrone_levels);
initial_hive = fraction_hive_healthy*double(total_bees)/double(nhive_levels);
initial_hive_infected = fraction_hive_infected*double(total_bees)/double(nhive_levels);
initial_forager = 0.*double(total_bees)/double(nforager_levels);

food = zeros(n+1,1);
% Initial food (in grams) 
food(1) = 2000000.0; % 

Egg = zeros(n+1,negg_levels);
Etot = zeros(n+1,1);
Egg_drone = zeros(n+1,negg_levels);
Etot_drone = zeros(n+1,1);

B = zeros(n+1,nlarvae_levels);
Btot = zeros(n+1,1);
B_drone = zeros(n+1,nlarvae_drone_levels);
Btot_drone = zeros(n+1,1);

P = zeros(n+1,npupae_levels);
P_infected = zeros(n+1,npupae_levels);
Ptot = zeros(n+1,1);
Ptot_infected = zeros(n+1,1);

P_drone = zeros(n+1,npupae_drone_levels);
P_infected_drone = zeros(n+1,npupae_drone_levels);
Ptot_drone = zeros(n+1,1);
Ptot_infected_drone = zeros(n+1,1);

H = zeros(n+1,nhive_levels);
H_infected = zeros(n+1,nhive_levels);
Htot = zeros(n+1,1);
Htot_infected = zeros(n+1,1);
transfer = zeros(nhive_levels,1);

Drone = zeros(n+1,ndrone_levels);
Drone_infected = zeros(n+1,ndrone_levels);
Dronetot = zeros(n+1,1);
Dronetot_infected = zeros(n+1,1);
transfer_drone = zeros(ndrone_levels,1);

F = zeros(n+1,nforager_levels);
F_infected = zeros(n+1,nforager_levels);
Ftot = zeros(n+1,1);
Ftot_infected = zeros(n+1,1);

Adult = zeros(n+1,1);
Brood = zeros(n+1,1);
Beetot = zeros(n+1,1);

Mite = zeros(n+1,nmite_levels);
Mite_tot = zeros(n+1,1);

Mite_daughter = zeros(n+1,npupae_levels);
Mite_daughter_tot = zeros(n+1);
Mite_foundress = zeros(n+1,npupae_levels);
Mite_foundress_tot= zeros(n+1);
Mite_daughter_drone = zeros(n+1,npupae_drone_levels);
Mite_daughter_drone_tot = zeros(n+1);
Mite_foundress_drone = zeros(n+1,npupae_drone_levels);
Mite_foundress_drone_tot = zeros(n+1);
Mite_brood = zeros(n+1,12);
Mites_in_brooda = zeros(n+1,1);
s_pia = zeros(n+ijump*npupae_levels,npupae_levels+1);
s_pia_drone = zeros(n+ijump*npupae_drone_levels,npupae_drone_levels+1);
multiple_worker_level = zeros(n+ijump*npupae_levels,npupae_levels+1);
multiple_drone_level = zeros(n+ijump*npupae_drone_levels,npupae_drone_levels+1);


% Survival rates of infested pupae
for i = 1:n+ijump*npupae_levels
    for j = 1:npupae_levels+1
        s_pia(i,j) = s_pi;
    end
end 

for i = 1:n+ijump*npupae_drone_levels
    for j = 1:npupae_drone_levels+1
        s_pia_drone(i,j) = s_pi;
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the initial populations on the first day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eggs
Etot(1) = 0.0;
for j = 1:negg_levels
   Egg(1,j) = initial_egg;
   Etot(1) = Etot(1) + Egg(1,j);
end

% Larvae
Btot(1) = 0.0;
for j = 1:nlarvae_levels 
    % Larvae 3 days old and younger eat royal jelly
    B(1,j) = initial_larvae;
    Btot(1) = Btot(1) + B(1,j);
end

% Pupae do not consume food
Ptot(1) = 0.0;
for j = 1:npupae_levels
    P(1,j) = initial_pupae;
    Ptot(1) = Ptot(1) + P(1,j);
end

% Drones
Dronetot(1) = 0.0;
for j = 1:ndrone_levels
    Drone(1,j) = initial_drone;
    Dronetot(1) = Dronetot(1) + Drone(1,j);
end

% Hive bees
Htot(1) = 0.0;
Nursetot = 0.0;
for j = 1:nhive_levels
    H(1,j) = initial_hive;
    if (j <= nnurse_levels)
        Nursetot = Nursetot + H(1,j);
    end
    Htot(1) = Htot(1) + H(1,j);
end


% Hive bees
Htot_infected(1) = 0.0;
Nursetot_infected = 0.0;
for j = 1:nhive_levels
    H_infected(1,j) = initial_hive_infected;
    if (j <= nnurse_levels)
        Nursetot_infected = Nursetot_infected + H_infected(1,j);
    end
    Htot_infected(1) = Htot_infected(1) + H_infected(1,j);
    Htot(1) = Htot(1) + H_infected(1,j);
end
Nursetot_infected = 0.0;

% Foraging bees
Ftot(1) = 0.0;
for j = 1:nforager_levels
    F(1,j) = initial_forager;
    Ftot(1) = Ftot(1) + F(1,j);
end

Mite(1,4) = initial_mite; % Initially all mites are 4 days old
Mite_tot(1) = Mite(1,4); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Adult(1) = Htot(1) + Ftot(1); % Hive and Forager bees
Brood(1) = Etot(1)+Btot(1)+Ptot(1); % Egg and Brood and Pupae
Beetot(1) = Etot(1)+Btot(1)+Ptot(1)+Htot(1)+Ftot(1)+Dronetot(1);


% Food required by colony in one day - pupae and eggs do not require food
% (units are grams/day)
foodreq = (Ftot(1)*gamma_forager + Htot(1)*gamma_hive + Dronetot(1)*gamma_drone + ...
           Btot(1)*gamma_larvae);    

food_inaccessible = 100.d0; % grams/day
food_accessible = max(food(1) - food_inaccessible,0.d0);
if (food_accessible < foodreq) 
   deficit = foodreq - food_accessible;
else
   deficit = 0.0;
end

mites_death_grooming = 0.0;

% Parameters used in computing seasonal egg laying rate
%x1 = 385.d0;
%x2 = 30.d0;
%x3 = 36.d0;
%x4 = 155.d0;
%x5 = 30.d0;

% pher: Variable which regulates age maturation of old hive bees due to
% pheromones
pheromone = 1.d0;

% Variables which modify the survival rate due to food shortage
fsb = 1.d0;
fsh = 1.d0;
fsd = 1.d0;
fsf = 1.d0;


if (time(1) >= summer_begin && time(1) <= summer_end)
    summer = 1;
else
    summer = 0;
end

hivemax = -1.d+20;
Mite_tot_dist = 0.0;
pupae_death = 0.0;

rate_sum1 = 0.0;
rate_sum2 = 0.0;
num_sum1 = 0.0;
num_sum2 = 0.0;

zsum = 0.0;
nzsum = 0.0;
% Loop over time
for i = 1:n

    ta(i) = time(i);
    timei = mod(time(i),365);
    
    
    % Calculating day length at a latitude
    pv =0.0;
    latitude = 35.0;
    theta = .2163108 + 2.0*atan(.9671396*tan(.00860*(timei-186)));
    phi = asin(.39795*cos(theta));
    arg = (sin(pv*pi/180) + sin(latitude*pi/180.)*sin(phi))/(cos(latitude*pi/180.)*cos(phi));
    day_length = 24.0 - (24.0/pi)*acos(arg);
    
    %season1 = 1.d0 - 1.d0/(1.d0 + x1*exp(-2.d0*timei/x2));
    %season2 = 1.d0/(1.d0 + x3*exp(-2.d0*(timei-x4)/x5));
    %season(i) = max(season1,season2);
    
    % Day at which egg laying rate peaks
    summer_mid = floor(0.5*(summer_begin+summer_end));
  
    % Parameters that insure egg-laying rate does not begin and end at zero
    % before ramping up to maximum rate
    incr_summer = .25*(summer_end - summer_begin);
    summer_early = floor(summer_begin - incr_summer);
    summer_late = floor(summer_end + .25*incr_summer);
    summer_egg = summer_begin - 10.;
    season(i) = 0.0;
    if (timei <= summer_mid && timei >= summer_early)
       season(i) = sin(0.5*pi*(timei - summer_early)/(summer_mid - summer_early))^2;
    end
    if (timei <= summer_late && timei > summer_mid)
       season(i) = cos(0.5*pi*(timei - summer_mid)/(summer_late - summer_mid))^2;
    end
    
    
    summerp = summer; % Previous day summer flag
    
    if (timei >= summer_begin && timei <= summer_end)
        % Summer
        summer = 1;
        
        c = forager_rate; % Grams of food gathered in one day by one forager (and processed by hive)
        ci = forager_rate/3.d0;  % Grams of food gathered in one day by one infected forager (and processed by hive)
        % Sumpter and Martin (2004) state that there are 50 drones per 1500
        % hive bees (3.3%) in the summer and 5 drones per 500 hive worker bees (1%) in
        % the autumn
        %fsh reduces the survival in the absence of sufficient food
        s_h = fsh*(1.d0 - hive_death_rate);  % Survival rate of hive bees in summer
        s_hi = fsh*(1.d0 - hive_death_rate_infected); % Survival rate of infected hive bees in summer
        
    else
        % Winter
        summer = 0;
        c = 0.; % Foraging rate set to zero
        ci = 0.;
        s_h = fsh*(1.d0 - hive_death_rate_winter);  % Survival rate of hive bees in winter
        s_hi = fsh*(1.d0 - hive_death_rate_winter_infected); % Survival rate of infected hive bees in winter
    end

    if (summerp == 0 && summer == 1) 
        switch_winter_to_summer = 1;
    else
        switch_winter_to_summer = 0;
    end
    
    if (summerp == 1 && summer == 0)
        switch_summer_to_winter = 1;
    else
        switch_summer_to_winter = 0;
    end
    
    if (Htot(i) > 1000. && timei >= summer_egg && timei <= summer_end)   
       %L = Lmax*(1.d0-season(i)); % Egg laying rate is affected by day in summer
       L = Lmax*season(i); % Egg laying rate is affected by day in summer
    else
       L = 0.d0; % Set egg laying rate to zero if there are 1000 or less hive bees
    end
 
    propz = max(log10(day_length*.1+1.d-16)*.284*log10(Htot(i)*.0006+1.d-16)*.797,0.0);
    propza(i) = (propz*L)/(L + eps);
    
    if (propza(i) > 1.d-3)
       zsum = zsum + propza(i);
       nzsum = nzsum + 1.0;
    end

    rexp = 1.d0/expl;
    ratio_H_L = Htot(i)/(Btot(i)+eps); % Ratio of total hive bees to total brood
    ratio_N_L = Nursetot/(Btot(i)+eps); % Ratio of total nurse bees to total brood
    
    Htot_healthy = Htot(i) - Htot_infected(i);
    %Nursetot_healthy = Nursetot - Nursetot_infected;
    %ratio_N_L_eff = (Nursetot_healthy*ratio_N_L_eq + Nursetot_infected*ratio_N_L_eq_infected)/Nursetot;
    ratio_H_L_eff = (Htot_healthy*ratio_H_L_eq + Htot_infected(i)*ratio_H_L_eq_infected)/(Htot(i)+eps);
   
    if (ratio_H_L  < ratio_H_L_eff)
       s_b = fsb*(1.d0 - larvae_death_rate)*(ratio_H_L/ratio_H_L_eff)^(rexp);
       %s_b = fsb*(1.d0 - larvae_death_rate)*(ratio_N_L/ratio_N_L_eff)^(rexp); % Brood survival rate is affected by number of nursing bees
    else
       s_b = fsb*(1.d0 - larvae_death_rate);
    end

    % Parameters for brood pheromone (affect maturation of hive bees)
    accel_brood = 1.d0 + 0.5d0*(ratio_H_L - ratio_H_L_eq)/ratio_H_L_eq;
    accel_brood = max(accel_brood,.5d0);
    accel_brood = min(accel_brood,2.d0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Solve differential equations
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing eggs
    
    % Eggs don't require food - After 3 days eggs become larvae
    Egg(i+1,1) = Egg(i,1)*(1.d0 - dt) + L*(1.0-propz)*dt;
    egg1 = Egg(i+1,1);
    Egg(i+1,1) = max(0.0,Egg(i+1,1));
    Egg(i+1,2) = Egg(i,2)*(1.d0 - dt) + Egg(i,1)*s_e*dt;
    egg2 = Egg(i+1,2);
    Egg(i+1,2) = max(0.0,Egg(i+1,2));
    Egg(i+1,3) = Egg(i,3)*(1.d0 - dt) + Egg(i,2)*s_e*dt;
    egg3 = Egg(i+1,3);
    Egg(i+1,3) = max(0.0,Egg(i+1,3));
    Etot(i+1) = Egg(i+1,1)+Egg(i+1,2)+Egg(i+1,3);
    
    % Eggs don't require food - After 3 days eggs become larvae
    Egg_drone(i+1,1) = Egg_drone(i,1)*(1.d0 - dt) + L*propz*dt;
    Egg_drone(i+1,1) = max(0.0,Egg_drone(i+1,1));
    Egg_drone(i+1,2) = Egg_drone(i,2)*(1.d0 - dt) + Egg_drone(i,1)*s_e*dt;
    Egg_drone(i+1,2) = max(0.0,Egg_drone(i+1,2));
    Egg_drone(i+1,3) = Egg_drone(i,3)*(1.d0 - dt) + Egg_drone(i,2)*s_e*dt;
    Egg_drone(i+1,3) = max(0.0,Egg_drone(i+1,3));
    Etot_drone(i+1) = Egg_drone(i+1,1)+Egg_drone(i+1,2)+Egg_drone(i+1,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing larvae
    % Larva - 5 days
    B(i+1,1) = B(i,1)*(1.d0 - dt) + Egg(i,3)*s_e*dt;
    B(i+1,1) = max(0.d0,B(i+1,1));
    Btot(i+1) = B(i+1,1);
    
    total_weight_deaths_larvae = 0.d0;

    for j = 2:nlarvae_levels
        B(i+1,j) = B(i,j)*(1.d0 - dt)+ B(i,j-1)*s_b*dt;
        total_weight_deaths_larvae =  total_weight_deaths_larvae + ...
        (1.d0 - s_b)*B(i,j-1)*weight(j-1);
        B(i+1,j) = max(0.0,B(i+1,j));
        Btot(i+1) = Btot(i+1) + B(i+1,j);
    end
    
    
    B_drone(i+1,1) = B_drone(i,1)*(1.d0 - dt) + Egg_drone(i,3)*s_e*dt;
    B_drone(i+1,1) = max(0.d0,B_drone(i+1,1));
    Btot_drone(i+1) = B_drone(i+1,1);
    for j = 2:nlarvae_drone_levels
        B_drone(i+1,j) = B_drone(i,j)*(1.d0 - dt)+ B_drone(i,j-1)*s_b*dt;
        %total_weight_deaths_larvae =  total_weight_deaths_larvae + ...
        %(1.d0 - s_b)*B_drone(i,j-1)*weight(j-1);
        B_drone(i+1,j) = max(0.0,B_drone(i+1,j));
        Btot_drone(i+1) = Btot_drone(i+1) + B_drone(i+1,j);
    end
    
    
    total_weight_deaths_larvae = total_weight_deaths_larvae + ...
    (1.d0 - s_b)*B(i,nlarvae_levels)*weight(nlarvae_levels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %total_weight_deaths_larvae = total_weight_deaths_larvae + ...
    %(1.d0 - s_b)*B_drone(i,nlarvae_levels)*weight(nlarvae_levels);
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switch to summer
    if (switch_winter_to_summer == 1) 
      % Triggered at beginning of summer
      
       Mitetot = 0.0;
       for j =1:nmite_levels
         Mitetot = Mitetot + Mite(i,j);
         Mite(i,j) = 0.0;
       end
       
       Hivetot = 0.0;
       for j = 1:nhive_levels
           Hivetot = Hivetot + H(i,j);
           H(i,j) = 0.0;
       end
       
       Hivetot_infected = 0.0;
       for j = 1:nhive_levels
           Hivetot_infected = Hivetot_infected + H_infected(i,j);
           H_infected(i,j) = 0.0;
       end
      
       % Hive bees that survive winter are made 1 day old to prevent
       % early maturation into foragers
       %H(i,10) = H(i,10) + H(i,nhive_levels);
       H(i,1) = Hivetot;
       H_infected(i,1) = Hivetot_infected;
       Mite(i,4) = Mitetot;
       
       Mite_tot(i) = Mitetot;
       Mite_tot_dist = 0.0;

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%   Distribute mites to soon to be capped brood cells    
    reduce_workera(i) = 1.0;
    Mites_in_brooda(i) = 0.0;  
    
    if (summer == 1) 
      
       %Mite_tot_dist are the number of fertile mites (nfertile_mites)
       Mite_theory_brood = .75*Mite_tot_dist;
       % Foundress mites in capped brood
       Mites_in_brood = Mite_foundress_tot(i) + Mite_foundress_drone_tot(i);

       Mites_in_brooda(i) = Mites_in_brood/(Mite_tot_dist+eps);
       
       % Mites that need to be distributed to available brood cells
       Mite_distribute = max(0.0,Mite_theory_brood - Mites_in_brood);
         
       if (Mite_distribute > 0) 
           
          Mite_worker_brood = .15*Mite_distribute;
          Mite_distribute = Mite_distribute - Mite_worker_brood;
 
          % First distribute mites to drones cells
          available_drone_cells = B_drone(i,nlarvae_drone_levels);
          
           for kji = 1:max_levels_drone
              multiple_drone(kji) = 0.0;
           end
          % Enough mites are available to multiply infest drone brood each with
          % multiple_infested_drone mites
          multiple_infested_drone = min(max_levels_drone,floor(Mite_distribute/available_drone_cells));
          
          for kji = 1:multiple_infested_drone
              multiple_drone(kji) = available_drone_cells; 
              % multiple_drone(kji) stores the number of drone cells that
              % have at least kji mites
          end
          
          % Additional mites may be available to spill over into worker
          % brood cells
          Mite_distribute = Mite_distribute - double(multiple_infested_drone)*available_drone_cells;
          if (multiple_infested_drone == max_levels_drone)
             rr_drone = 0.0;
          else
             multiple_drone(multiple_infested_drone+1) = Mite_distribute;
             rr_drone = double(Mite_distribute)/(double(available_drone_cells)+eps);
             Mite_distribute = 0;
          end
          
          % Distribute remaining mites and reserved mites (Mite_worker_brood) to worker brood
          Mite_distribute = Mite_distribute + Mite_worker_brood;

          available_worker_cells = B(i,nlarvae_levels);
          
          for kji = 1:max_levels_worker
              multiple_worker(kji) = 0.0;
              % multiple_worker(kji) stores the number of worker cells that
              % have at least kji mites
          end
          
          multiple_infested_worker = min(max_levels_worker,floor(Mite_distribute/available_worker_cells));
          for kji = 1:multiple_infested_worker
              multiple_worker(kji) = available_worker_cells;
          end
          
          Mite_distribute = Mite_distribute - double(multiple_infested_worker)*available_worker_cells;
          if (multiple_infested_worker == max_levels_worker)
             rr_worker = 0.0;
          else
             multiple_worker(multiple_infested_worker+1) = Mite_distribute;
             rr_worker = Mite_distribute/(double(available_worker_cells)+eps);
          end
                
          brood_drone_mites = 0.0;
          foundress_drone = 0.0;
          for kji = 1:max_levels_drone
              % brood_drone_mites are the number of daughter mites
              % generated in drone cells.  The survival rate of daughter
              % mites is reduced for multiple infestations
              brood_drone_mites = brood_drone_mites + multiple_drone(kji)*survival_multiple_drones(kji);
              % foundress_drone keeps track of foundress mites
              % in drone cells.  
              foundress_drone = foundress_drone + multiple_drone(kji);
          end

          brood_worker_mites = 0.0;
          foundress_worker = 0.0;
          for kji = 1:max_levels_worker
              % brood_worker_mites are the number of daughter mites
              % generated in worker cells.  The survival rate of daughter
              % mites is reduced for multiple infestations
              brood_worker_mites = brood_worker_mites + multiple_worker(kji)*survival_multiple_workers(kji);
              foundress_worker = foundress_worker + multiple_worker(kji);
          end
          
          % Fraction of soon to be capped drone brood infected in time step dt
          fraction_infected_pupae_drone = min(multiple_drone(1)/(B_drone(i,nlarvae_drone_levels)+eps),1.0);
          fraction_healthy_pupae_drone = 1.d0 - fraction_infected_pupae_drone;
          
          % Fraction of soon to be capped worker brood infected in time step dt
          fraction_infected_pupae = min(multiple_worker(1)/(B(i,nlarvae_levels)+eps),1.0);
          fraction_healthy_pupae = 1.d0 - fraction_infected_pupae;

          % Survival rate of drone is reduced if cells are multiply
          % infested
          reduce_drone = (1.0-rr_drone)*reduce_level_drone(multiple_infested_drone+1) + rr_drone*reduce_level_drone(multiple_infested_drone+2);
          
          % Survival rate of pupae is reduced if cells are multiply
          % infested
          reduce_worker = (1.0-rr_worker)*reduce_level(multiple_infested_worker+1) + rr_worker*reduce_level(multiple_infested_worker+2);
       
          reduce_workera(i) = reduce_worker^npupae_levels;
          
          
          % The survival rate of a dt cohort of drone pupae is accessed everytime
          % the cohort moves into the next day
          for kij = 1:npupae_drone_levels+1
              s_pia_drone(i+ijump*(kij-1),kij) = s_pi*reduce_drone;
          end
                   
          % The survival rate of a dt cohort of pupae is accessed everytime
          % the cohort moves into the next day
          for kij = 1:npupae_levels+1
              s_pia(i+ijump*(kij-1),kij) = s_pi*reduce_worker;
          end
          
          % The survival rate of a dt cohort of drone pupae is accessed everytime
          % the cohort moves into the next day
          for kij = 1:npupae_drone_levels+1
              if (multiple_infested_drone >= 1)
                 multiple_drone_level(i+ijump*(kij-1),kij) = double(multiple_infested_drone) + rr_drone;
              else
                 multiple_drone_level(i+ijump*(kij-1),kij) = 1.;
              end
          end
          
          % The survival rate of a dt cohort of pupae is accessed everytime
          % the cohort moves into the next day
          for kij = 1:npupae_levels+1
              if (multiple_infested_worker >= 1)
                 multiple_worker_level(i+ijump*(kij-1),kij) = double(multiple_infested_worker) + rr_worker;
              else
                 multiple_worker_level(i+ijump*(kij-1),kij) = 1.;
              end
          end

       else
          brood_worker_mites = 0.0;
          brood_drone_mites = 0.0;
          foundress_worker = 0.0;
          foundress_drone = 0.0;
          fraction_infected_pupae = 0.0;
          fraction_healthy_pupae = 1.0;
          fraction_infected_pupae_drone = 0.0;
          fraction_healthy_pupae_drone = 1.0;
          reduce_workera(i) = 1.;
       end
    
    else
       % Winter rates
       brood_worker_mites = 0.0;
       brood_drone_mites = 0.0;
       foundress_worker = 0.0;
       foundress_drone = 0.0;
       fraction_infected_pupae = 0.0;
       fraction_healthy_pupae = 1.0;
       fraction_infected_pupae_drone = 0.0;
       fraction_healthy_pupae_drone = 1.0;
       reduce_workera(i) = 1.0;
    end 
    
    if (summer == 1) 
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Daughter mite inside worker capped cells. Survival rates are
       % already accounted for in term brood_worker_mites. 
       Mite_daughter(i+1,1) = Mite_daughter(i,1)*(1.d0 - dt) + brood_worker_mites*dt;
       Mite_daughter(i+1,1) = max(0.d0,Mite_daughter(i+1,1));
       Mite_daughter_tot(i+1) = Mite_daughter(i+1,1);
       for j = 2:npupae_levels
           Mite_daughter(i+1,j) = Mite_daughter(i,j)*(1.0 - dt) + Mite_daughter(i,j-1)*dt;
           Mite_daughter(i+1,j) = max(0.d0,Mite_daughter(i+1,j));
           Mite_daughter_tot(i+1) = Mite_daughter_tot(i+1) + Mite_daughter(i+1,j);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Foundress mites in capped worker cells
       Mite_foundress(i+1,1) = Mite_foundress(i,1)*(1.d0 - dt) + foundress_worker*s_mite*dt;
       Mite_foundress(i+1,1) = max(0.d0,Mite_foundress(i+1,1));
       Mite_foundress_tot(i+1) = Mite_foundress(i+1,1);
       for j = 2:npupae_levels
           Mite_foundress(i+1,j) = Mite_foundress(i,j)*(1.0 - dt) + Mite_foundress(i,j-1)*s_mite*dt;
           Mite_foundress(i+1,j) = max(0.d0,Mite_foundress(i+1,j));
           Mite_foundress_tot(i+1) = Mite_foundress_tot(i+1) + Mite_foundress(i+1,j);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Daughter mites inside drone capped cells.  Survival rates are
       % already accounted for in term brood_drone_mites.
       Mite_daughter_drone(i+1,1) = Mite_daughter_drone(i,1)*(1.d0 - dt) + brood_drone_mites*dt;
       Mite_daughter_drone(i+1,1) = max(0.d0,Mite_daughter_drone(i+1,1));
       Mite_daughter_drone_tot(i+1) = Mite_daughter_drone(i+1,1);
       for j = 2:npupae_drone_levels
           Mite_daughter_drone(i+1,j) = Mite_daughter_drone(i,j)*(1.0 - dt) + Mite_daughter_drone(i,j-1)*dt;
           Mite_daughter_drone(i+1,j) = max(0.d0,Mite_daughter_drone(i+1,j));
           Mite_daughter_drone_tot(i+1) = Mite_daughter_drone_tot(i+1) + Mite_daughter_drone(i+1,j);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Foundress mites in capped drone cells
       Mite_foundress_drone(i+1,1) = Mite_foundress_drone(i,1)*(1.d0 - dt) + foundress_drone*s_mite*dt;
       Mite_foundress_drone(i+1,1) = max(0.d0,Mite_foundress_drone(i+1,1));
       Mite_foundress_drone_tot(i+1) = Mite_foundress_drone(i+1,1);
       for j = 2:npupae_drone_levels
           Mite_foundress_drone(i+1,j) = Mite_foundress_drone(i,j)*(1.0 - dt) + Mite_foundress_drone(i,j-1)*s_mite*dt;
           Mite_foundress_drone(i+1,j) = max(0.d0,Mite_foundress_drone(i+1,j));
           Mite_foundress_drone_tot(i+1) = Mite_foundress_drone_tot(i+1) + Mite_foundress_drone(i+1,j);
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    else
        Mite_daughter_tot(i+1) = 0.0;
        Mite_foundress_tot(i+1) = 0.0;
        Mite_daughter_drone_tot(i+1) = 0.0;
        Mite_foundress_drone_tot(i+1) = 0.0;
    end
    
   
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing mites

    Phoretic_mites = max(0.0,Mite_tot(i) - (Mite_foundress_tot(i) + Mite_foundress_drone_tot(i)));
    ratio_phoretic_worker = Phoretic_mites/(Htot_infected(i)+Dronetot_infected(i)+eps);
    

    if (ratio_phoretic_worker > xi)
       phoretic_mites_max = xi*(Htot_infected(i)+Dronetot_infected(i));
       reduce_mites_max = Phoretic_mites - phoretic_mites_max;
       ratio_phoretic_worker = xi;
    else
       reduce_mites_max = 0.0;
    end

    ratio_phoretic_workera(i) = ratio_phoretic_worker;
      
    mites_death_divide = (mites_death_grooming + dt*reduce_mites_max + dt*pupae_death)/double(nmite_levels);
    
    % Add mites emerging from worker pupae and drone pupae
    Mite(i+1,1) = Mite(i,1)*(1.d0 - dt) + Mite_daughter(i,npupae_levels)*s_mite*dt + Mite_daughter_drone(i,npupae_drone_levels)*s_mite*dt - mites_death_divide;
    Mite(i+1,1) = max(0.d0,Mite(i+1,1));
    Mite_tot(i+1) = Mite(i+1,1);
    for j =2:nmite_levels
        if (summer == 0 && j == nmite_levels)
           Mite(i+1,j) = Mite(i,j)*(1.0 - dt*mite_death_rate_winter) + Mite(i,j-1)*s_mite*dt - mites_death_divide;
        else
           Mite(i+1,j) = Mite(i,j)*(1.d0 - dt) + Mite(i,j-1)*s_mite*dt - mites_death_divide;
        end 
        Mite(i+1,j) = max(0.d0,Mite(i+1,j)); 
        Mite_tot(i+1) = Mite_tot(i+1) + Mite(i+1,j);
        if (j == nfertile_mites)
           Mite_tot_dist = Mite_tot(i+1);
        end
    end
    

    % Track the daily increase in mite population
    if (i >= ijump+1 && Mite_tot(i-ijump) > 1.0)
       rate_increase_mite(i) = max(100.*(Mite_tot(i)/(Mite_tot(i-ijump)+eps) - 1.0),0.0);
       if (rate_increase_mite(i) > 0 && double(i)*dt < 360)
          rate_sum1 = rate_sum1 + rate_increase_mite(i);
          num_sum1 = num_sum1 + 1.0;
       end
       if (rate_increase_mite(i) > 0 && double(i)*dt > 360)
          rate_sum2 = rate_sum2 + rate_increase_mite(i);
          num_sum2 = num_sum2 + 1.0;
       end
    else
       rate_increase_mite(i) = 0.0;
    end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing pupae       
    % Pupae stay pupae for 12 days.
    % Pupae don't need food.
    P(i+1,1) = P(i,1)*(1.d0 - dt) + fraction_healthy_pupae*B(i,nlarvae_levels)*s_b*dt;
    P(i+1,1) = max(0.0,P(i+1,1));
    Ptot(i+1) = P(i+1,1);
    for j = 2:npupae_levels
        P(i+1,j) = P(i,j)*(1.0 - dt) + P(i,j-1)*s_p*dt;
        P(i+1,j) = max(0.0,P(i+1,j));
        Ptot(i+1) = Ptot(i+1) + P(i+1,j);
    end
    
    %P_infected(i+1,1) = P_infected(i,1)*(1.d0 - dt) + fraction_infected_pupae*B(i,nlarvae_levels)*s_b*dt;
    P_infected(i+1,1) = P_infected(i,1)*(1.d0 - dt) + fraction_infected_pupae*B(i,nlarvae_levels)*s_pia(i,1)*dt;
    pupae_death = fraction_infected_pupae*B(i,nlarvae_levels)*(1.-s_pia(i,1))*multiple_worker_level(i,1);
    P_infected(i+1,1) = max(0.0,P_infected(i+1,1));
    Ptot(i+1) = Ptot(i+1) + P_infected(i+1,1);
    Ptot_infected(i+1) = P_infected(i+1,1);
    for j = 2:npupae_levels
        %ss = s_pia(i,j);
        %mw = multiple_worker_level(i,j);
        P_infected(i+1,j) = P_infected(i,j)*(1.0 - dt) + P_infected(i,j-1)*s_pia(i,j)*dt;
        pupae_death = pupae_death + P_infected(i,j-1)*(1. - s_pia(i,j))*multiple_worker_level(i,j);
        P_infected(i+1,j) = max(0.0,P_infected(i+1,j));
        Ptot_infected(i+1) = Ptot_infected(i+1) + P_infected(i+1,j);
        Ptot(i+1) = Ptot(i+1) + P_infected(i+1,j);
    end
    
    P_drone(i+1,1) = P_drone(i,1)*(1.d0 - dt) + fraction_healthy_pupae_drone*B_drone(i,nlarvae_drone_levels)*s_b*dt;
    P_drone(i+1,1) = max(0.0,P_drone(i+1,1));
    Ptot_drone(i+1) = P_drone(i+1,1);
    for j = 2:npupae_drone_levels
        P_drone(i+1,j) = P_drone(i,j)*(1.0 - dt) + P_drone(i,j-1)*s_p*dt;
        P_drone(i+1,j) = max(0.0,P_drone(i+1,j));
        Ptot_drone(i+1) = Ptot_drone(i+1) + P_drone(i+1,j);
    end
    
    P_infected_drone(i+1,1) = P_infected_drone(i,1)*(1.d0 - dt) + fraction_infected_pupae_drone*B_drone(i,nlarvae_drone_levels)*s_pia_drone(i,1)*dt;
    pupae_death = pupae_death + fraction_infected_pupae_drone*B_drone(i,nlarvae_drone_levels)*(1. - s_pia_drone(i,1))*multiple_drone_level(i,1);
    P_infected_drone(i+1,1) = max(0.0,P_infected_drone(i+1,1));
    Ptot_drone(i+1) = Ptot_drone(i+1) + P_infected_drone(i+1,1);
    Ptot_infected_drone(i+1) = P_infected_drone(i+1,1);
    for j = 2:npupae_drone_levels
        pupae_death = pupae_death + P_infected_drone(i,j-1)*(1. - s_pia_drone(i,j))*multiple_drone_level(i,j);
        P_infected_drone(i+1,j) = P_infected_drone(i,j)*(1.0 - dt) + P_infected_drone(i,j-1)*s_pia_drone(i,j)*dt;
        P_infected_drone(i+1,j) = max(0.0,P_infected_drone(i+1,j));
        Ptot_infected_drone(i+1) = Ptot_infected_drone(i+1) + P_infected_drone(i+1,j);
        Ptot_drone(i+1) = Ptot_drone(i+1) + P_infected_drone(i+1,j);
    end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % Ratio of hive bees to foraging bees
    ratio_H_F = Htot(i)/(Ftot(i)+eps);
    ratio_Hi_H = Htot_infected(i)/(Htot(i)+eps);
    
    rphor = transmission_rate*ratio_Hi_H;
    groom_transfer = H_infected(i,1)*grooming_rate*dt;
    mites_death_grooming = groom_transfer*ratio_phoretic_worker;    
  
    
    diff_Hi_mites =  Phoretic_mites - Htot_infected(i) - Dronetot_infected(i) - Ftot_infected(i);
    frachd = (Htot(i) + H_infected(i))/(Htot(i) + H_infected(i) + Dronetot(i)+Dronetot_infected(i) + eps);
    ofrachd = 1.d0 - frachd;
    
    if (diff_Hi_mites > 0) 
       transfer(1) = groom_transfer - H(i,1)*rphor*dt;
    else
       if (H_infected(i,1) > 0)
          transfer(1) = -diff_Hi_mites*dt*frachd/double(nhive_levels);
       else
          transfer(1) = 0.0;
       end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing hive bees     
    
    H(i+1,1) = H(i,1)*(1.0 - dt) + P(i,npupae_levels)*s_p*dt + transfer(1);
    H(i+1,1) = max(0.0,H(i+1,1));
    Nursetot = H(i+1,1);
    Htot(i+1) = H(i+1,1);
  
    for j = 2:nnurse_levels
        groom_transfer = H_infected(i,j)*grooming_rate*dt;
        mites_death_grooming = mites_death_grooming + groom_transfer*ratio_phoretic_worker;
        if (diff_Hi_mites > 0)
           transfer(j) = groom_transfer - H(i,j)*rphor*dt; 
        else
           if (H_infected(i,j) > 0)
              transfer(j) = -diff_Hi_mites*frachd*dt/double(nhive_levels);
           else
              transfer(j) = 0.0;
           end
        end

        H(i+1,j) = H(i,j)*(1.0 - dt) + H(i,j-1)*s_h*dt + transfer(j);
        H(i+1,j) = max(0.0,H(i+1,j));
        Htot(i+1) = Htot(i+1) + H(i+1,j);
        Nursetot = Nursetot + H(i+1,j);
    end
 
    for j = nnurse_levels+1:nhive_levels
        groom_transfer = H_infected(i,j)*grooming_rate*dt;
        mites_death_grooming = mites_death_grooming + groom_transfer*ratio_phoretic_worker;
        if (diff_Hi_mites > 0)
           transfer(j) = groom_transfer - H(i,j)*rphor*dt; 
        else     
           if (H_infected(i,j) > 0)
              transfer(j) = -diff_Hi_mites*frachd*dt/double(nhive_levels);
           else
              transfer(j) = 0.0;
           end
        end
        if (j == nhive_levels - ndays_transfer)
           H(i+1,j) = H(i,j)*(1.0 - dt*pheromone) + H(i,j-1)*s_h*dt + transfer(j);
        else
            if (j > nhive_levels - ndays_transfer)
                if (summer == 1 || j < nhive_levels ) 
                   % Summer OR j < nhive_levels (Winter and Summer)
                   % pheromone factor accelerates or decelerates development of
                   % hive bees
                   H(i+1,j) = H(i,j)*(1.0 - dt*pheromone) + H(i,j-1)*s_h*pheromone*dt + transfer(j);
                else
                   % Winter AND j = nhive_levels
                   % During the winter, all hive bees eventually accumulate in the oldest layer
                   % (nhive_levels) and die at the rate of
                   % hive_death_rate_winter per day
                   H(i+1,j) = H(i,j)*(1.0 - dt*hive_death_rate_winter) + H(i,j-1)*s_h*dt + transfer(j);
                end
            else
                % Winter or Summer j < nhive_levels - ndays_transfer
                H(i+1,j) = H(i,j)*(1.0 - dt) + H(i,j-1)*s_h*dt + transfer(j);
            end
        end
        H(i+1,j) = max(0.0,H(i+1,j));
        Htot(i+1) = Htot(i+1) + H(i+1,j);
    end
    
    
    
    H_infected(i+1,1) = H_infected(i,1)*(1.0 - dt) + P_infected(i,npupae_levels)*s_hi*dt ...
                      - transfer(1);
    H_infected(i+1,1) = max(0.0,H_infected(i+1,1));
    
    Nursetot_infected = H_infected(i+1,1);
    Htot_infected(i+1) = H_infected(i+1,1);
    Nursetot = Nursetot + H_infected(i+1,1);
    Htot(i+1) = Htot(i+1) + H_infected(i+1,1);
    
    for j = 2:nnurse_levels
        H_infected(i+1,j) = H_infected(i,j)*(1.0 - dt) + H_infected(i,j-1)*s_hi*dt - transfer(j);     
        H_infected(i+1,j) = max(0.0,H_infected(i+1,j));
        Htot_infected(i+1) = Htot_infected(i+1) + H_infected(i+1,j);
        Htot(i+1) = Htot(i+1) + H_infected(i+1,j);
        Nursetot_infected = Nursetot_infected + H_infected(i+1,j);
        Nursetot = Nursetot + H_infected(i+1,j);
    end
    
    for j = nnurse_levels+1:nhive_levels
        if (j == nhive_levels - ndays_transfer)
           H_infected(i+1,j) = H_infected(i,j)*(1.0 - dt*pheromone) + H_infected(i,j-1)*s_hi*dt - transfer(j);
        else
            if (j > nhive_levels - ndays_transfer)
                if (summer == 1 || j < nhive_levels ) 
                   % Summer
                   % Winter and j < nhive_levels
                   H_infected(i+1,j) = H_infected(i,j)*(1.0 - dt*pheromone) + H_infected(i,j-1)*s_hi*pheromone*dt  - transfer(j);
                else
                   % Winter and j = nhive_levels
                   H_infected(i+1,j) = H_infected(i,j)*(1.0 - dt*hive_death_rate_winter_infected) + ...
                   H_infected(i,j-1)*s_hi*dt - transfer(j);
                end
            else
                % Winter or Summer j < nhive_levels - ndays_transfer
                H_infected(i+1,j) = H_infected(i,j)*(1.0 - dt) + H_infected(i,j-1)*s_hi*dt ...
                - transfer(j);
            end
        end
        H_infected(i+1,j) = max(0.0,H_infected(i+1,j));
        Htot_infected(i+1) = Htot_infected(i+1) + H_infected(i+1,j);
        Htot(i+1) = Htot(i+1) + H_infected(i+1,j);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing drone bees    
    groom_transfer = Drone_infected(i,1)*grooming_rate*dt;
    mites_death_grooming = mites_death_grooming + groom_transfer*ratio_phoretic_worker;
    diff_Hi_mites =  Phoretic_mites - Htot_infected(i+1) - Dronetot_infected(i);
    
    if (diff_Hi_mites > 0)    
       transfer_drone(1) = groom_transfer - Drone(i,1)*rphor*dt;
    else
       if (Drone_infected(i,1) > 0.0)      
          transfer_drone(1) = -diff_Hi_mites*ofrachd*dt/double(ndrone_levels);
       else
          transfer_drone(1) = 0.0;
       end
    end
    Drone(i+1,1) = Drone(i,1)*(1.0 - dt) + P_drone(i,npupae_drone_levels)*s_p*dt + transfer_drone(1);
    Drone(i+1,1) = max(0.0,Drone(i+1,1));
    Dronetot(i+1) = Drone(i+1,1);
    for j = 2:ndrone_levels
        groom_transfer = Drone_infected(i,j)*grooming_rate*dt;
        mites_death_grooming = mites_death_grooming + groom_transfer*ratio_phoretic_worker;
        if (diff_Hi_mites > 0)  
           transfer_drone(j) = groom_transfer - Drone(i,j)*rphor*dt;
        else
           if (Drone_infected(i,j) > 0.0)       
              transfer_drone(j) = -diff_Hi_mites*ofrachd*dt/double(ndrone_levels);
           else
              transfer_drone(j) = 0.0;
           end
        end
        Drone(i+1,j) = Drone(i,j)*(1.0-dt) + Drone(i,j-1)*s_h*dt + transfer_drone(j);
        Drone(i+1,j) = max(0.0,Drone(i+1,j));
        Dronetot(i+1) = Dronetot(i+1) + Drone(i+1,j);
    end
    
    
    Drone_infected(i+1,1) = Drone_infected(i,1)*(1.0 - dt) + P_infected_drone(i,npupae_drone_levels)*s_hi*dt - transfer_drone(1);
    Drone_infected(i+1,1) = max(0.0,Drone_infected(i+1,1));
    Dronetot_infected(i+1) = Drone_infected(i+1,1);
    Dronetot(i+1) = Dronetot(i+1) + Drone_infected(i+1,1);
    for j = 2:ndrone_levels
        Drone_infected(i+1,j) = Drone_infected(i,j)*(1.0-dt) + Drone_infected(i,j-1)*s_hi*dt - transfer_drone(j);
        Drone_infected(i+1,j) = max(0.0,Drone_infected(i+1,j));
        Dronetot_infected(i+1) = Dronetot_infected(i+1) + Drone_infected(i+1,j);
        Dronetot(i+1) = Dronetot(i+1) + Drone_infected(i+1,j);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Differential equations governing forager bees   
    if (summer == 1)
       F(i+1,1) = F(i,1)*(1.0 - dt) + H(i,nhive_levels)*s_h*pheromone*dt;
    else
       F(i+1,1) = F(i,1)*(1.0 - dt);  
    end
    F(i+1,1) = max(0.0,F(i+1,1));
    Ftot(i+1) = F(i+1,1);
    for j = 2:nforager_levels-1
        if (j < ndays_transfer)
           F(i+1,j) =  F(i,j)*(1.0 - dt) + F(i,j-1)*s_f*fsf*dt;   
        else
            if (j == ndays_transfer)
                F(i+1,j) =  F(i,j)*(1.0 - dt) + F(i,j-1)*s_f*fsf*dt; 
            else
                F(i+1,j) =  F(i,j)*(1.0 - dt) + F(i,j-1)*s_f*fsf*dt;
            end
        end
        F(i+1,j) = max(0.0,F(i+1,j));
        Ftot(i+1) = Ftot(i+1) + F(i+1,j);
    end
    F(i+1,nforager_levels) = F(i,nforager_levels)*(1.d0-forager_death_rate_end*dt) + F(i,nforager_levels-1)*s_f*fsf*dt;
    F(i+1,nforager_levels) = max(0.0,F(i+1,nforager_levels));
    Ftot(i+1) = Ftot(i+1) + F(i+1,nforager_levels);
    
    if (summer == 1)
       F_infected(i+1,1) = F_infected(i,1)*(1.0 - dt) + H_infected(i,nhive_levels)*s_hi*pheromone*dt;
    else
       F_infected(i+1,1) = F_infected(i,1)*(1.0 - dt);  
    end
    F_infected(i+1,1) = max(0.0,F_infected(i+1,1));
    Ftot_infected(i+1) = F_infected(i+1,1);
    Ftot(i+1) = Ftot(i+1) + F_infected(i+1,1);

    for j = 2:nforager_levels-1
        if (j < ndays_transfer)
        F_infected(i+1,j) =  F_infected(i,j)*(1.0 - dt) + F_infected(i,j-1)*s_fi*fsf*dt;   
        else
            if (j == ndays_transfer)
                F_infected(i+1,j) =  F_infected(i,j)*(1.0 - dt) + F_infected(i,j-1)*s_fi*fsf*dt; 
            else
                F_infected(i+1,j) =  F_infected(i,j)*(1.0 - dt) + F_infected(i,j-1)*s_fi*fsf*dt;
            end
        end
        F_infected(i+1,j) = max(0.0,F_infected(i+1,j));
        Ftot_infected(i+1) = Ftot_infected(i+1) + F_infected(i+1,j);
        Ftot(i+1) = Ftot(i+1) + F_infected(i+1,j);
    end
    F_infected(i+1,nforager_levels) = F_infected(i,nforager_levels)*(1.d0-forager_death_rate_end*dt) + F_infected(i,nforager_levels-1)*s_fi*fsf*dt;
    F_infected(i+1,nforager_levels) = max(0.0,F_infected(i+1,nforager_levels));
    Ftot_infected(i+1) = Ftot_infected(i+1) + F_infected(i+1,nforager_levels);
    Ftot(i+1) = Ftot(i+1) + F_infected(i+1,nforager_levels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    mag_ethyl_oleate = 0.5d0;
    foodfac = mag_ethyl_oleate*deficit/(foodreq+eps);

    ratio_H_F_mod = ratio_H_F_eq - foodfac;

    diff = (ratio_H_F - ratio_H_F_mod)/ratio_H_F_mod;
    accel_ethyl = 1.d0 + diff;
    accel_ethyl = min(3.d0,accel_ethyl);
    accel_ethyl = max(1.d0/3.d0,accel_ethyl);
    accel = accel_brood*accel_ethyl;
    pheromone = min(3.d0,accel);
    pheromone = max(1.d0/3.d0,pheromone);
    if (c < .0000001)
        % Winter
       pheromone = 1.d0;
    end
 
    phera(i) = pheromone;
    
    Adult(i+1) = Htot(i+1)+Ftot(i+1);
    
    Brood(i+1) = Etot(i+1)+Btot(i+1)+Ptot(i+1);
    Beetot(i+1) = Etot(i+1)+Btot(i+1)+Ptot(i+1)+Htot(i+1)+Ftot(i+1)+Dronetot(i+1);

    food_need_forager = Ftot(i+1)*gamma_forager;
    food_need_hive = Htot(i+1)*gamma_hive;
    food_need_drone = Dronetot(i+1)*gamma_drone;
    food_need_larvae = Btot(i+1)*gamma_larvae;
    food_need_adult = food_need_forager+food_need_hive+food_need_drone;
    % Food required by colony in one day
    foodreq = food_need_adult + food_need_larvae;
    
    % Processors
    Q = Htot(i+1) - Nursetot;
    QFr = Q/(Ftot(i+1)+eps);
    fraction = double(nhive_levels - nnurse_levels + 1)/double(nhive_levels);
    QFrhealthy = .8d0*fraction*ratio_H_F_eq;
    QFrhealthy = 1.d0;
 
    if (QFr < QFrhealthy)
        %cred = ( QFr/QFrhealthy )^rexp;
        cred = ( QFr/QFrhealthy );
    else
        cred = 1.d0;
    end

    food(i+1) = food(i) +  c*cred*(Ftot(i)-Ftot_infected(i))*dt + ci*cred*Ftot_infected(i)*dt + foodreq*dt;
    food(i+1) = max(0.0,food(i+1));
    
    food_accessible = max(food(i+1) - food_inaccessible,0.d0);
    
    if (food_accessible < foodreq) 
       deficit = min(foodreq - food_accessible,0.d0);
    else
       deficit = 0.0;
    end
    
    if (food_accessible < foodreq) 
        % If the accessible food is greater than the daily food requirement,
        % there are no deaths due to food
        
        % Addition of food due to cannibalization
       
        food(i+1) = food(i+1) + total_weight_deaths_larvae*.5d0*dt;
	    %pause  % Changed from Github
      
        if (food_accessible < 0) 
           
           % Survival rates after one day
           fsf = nofood_surv_forager;
           fsh = nofood_surv_hive;
           fsd = nofood_surv_hive;
           fsb = nofood_surv_brood;
       
        else
            
            
            if (food_accessible > food_need_adult)
                
                food_left_larvae = food_accessible - food_need_adult;
                surv_larvae = food_left_larvae/(food_need_larvae+eps);
                fsf = 1.d0;
                fsh = 1.d0;
                fsd = 1.d0;
                fsb = max(surv_larvae,nofood_surv_brood);
                
            else
                
                'Not enough food'
                surv_adult = food_accessible/(food_need_adult+eps);
                fsf = max(surv_adult,nofood_surv_forager);
                fsh = max(surv_adult,nofood_surv_hive);
                fsd = max(surv_adult,nofood_surv_hive);
                fsb = nofood_surv_brood;
                
            end
                
                
        end
                     
   
    else
        
       fsf = 1.d0;
       fsh = 1.d0;
       fsd = 1.d0;
       fsb = 1.d0;
       
    end
    
    if (Adult(i+1) < 1000.) 
        fsf = 0.d0;
        fsh = 0.d0;
        fsd = 0.d0;
        fsb = 0.d0;
        %pause
    end
    time(i+1) = time(i) + dt;


end  %n

rr1 = rate_sum1/num_sum1;
rr2 = rate_sum2/num_sum2;

zavg = zsum/nzsum;
fprintf('Rate increase in first year %d \n',rr1);
fprintf('Rate increase in the second year %d \n',rr2);
fprintf('Average drone proportion %d \n',zavg);

if (iplot == 1)
        
   nsymbols = 50;
   r = n/nsymbols;
   iskip = int32(r);
   %iskip = 5;

   iii = 0;
   fileID = fopen('mite0','w');
   
   for i = 1:iskip:n+1
       iii =  iii + 1;
       Etot1a(iii) = Etot(i);
       EE = Etot1a(iii);
       Btot1a(iii) = Btot(i);
       BB = Btot1a(iii);
       Htot1a(iii) = Htot(i);
       HH = Htot1a(iii);
       Hitot1a(iii) = Htot_infected(i);
       HHi = Hitot1a(iii);
       Drone1a(iii) = Dronetot(i);
       DD = Drone1a(iii);
       Dronei1a(iii) = Drone_infected(i);
       DDi = Dronei1a(iii);
       Ftot1a(iii) = Ftot(i);
       FF = Ftot1a(iii);
       Fitot1a(iii) = Ftot_infected(i);
       FFi = Fitot1a(iii);
       Adult1a(iii) = Adult(i);
       Brood1a(iii) = Brood(i);
       Ptot1a(iii) = Ptot(i);
       PP = Ptot1a(iii);
       Pitot1a(iii) = Ptot_infected(i);
       PPi = Pitot1a(iii);
       Mite1a(iii) = Mite_tot(i);
       MM = Mite1a(iii);

       time1a(iii) = time(i);
       food1a(iii) = food(i);
       fprintf(fileID,'%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f \n',EE,BB,PP,PPi,DD,DDi,HH,HHi,FF,FFi,MM); 
   end
   fclose(fileID);
   
   plot(ta,100.*propza)
   xlabel('Days','FontSize',16);
   ylabel('Percentage of eggs that become drones','FontSize',12);
   axis([30 720 0 4])
   highresjpg('fig2.jpg')
   %set(h_l,'FontSize',14);
   legend boxoff
   pause(2)  
   
   plot(time,Mites_in_brooda*100)
   xlabel('Days','FontSize',18);
   ylabel('Percentage of mites in capped cells','FontSize',14);
   axis([30 720 0 85])
   highresjpg('mites_in_brood.jpg')
   %%set(h_l,'FontSize',14);
   %legend boxoff
   pause(2) 
   %rate_increase_avg = zeros(n+1-10,1);
   %tavg = zeros(n+1-10,1);
   
   for i = 50+1:n+1-50-1
       rate_increase_avg(i-50) = 0.0;
       for j = -50:50
           rate_increase_avg(i-50) = rate_increase_avg(i-50) + rate_increase_mite(i+j);
       end
       rate_increase_avg(i-50) = rate_increase_avg(i-50)/101.0;
       tavg(i-50) = ta(i);
   end
   
   %plot(ta,100.*(rate_increase_mite-1.))
   plot(ta,rate_increase_mite)
   %plot(tavg,rate_increase_avg)
   xlabel('Days','FontSize',18);
   ylabel('Daily percent rate of increase in mite population','FontSize',12);
   axis([30 720 0 8])
   highresjpg('fig3.jpg')
   %set(h_l,'FontSize',14);
   legend boxoff
   pause(2)  
   

   plot(time1a,Etot1a,'bo-',time1a,Btot1a,'kd-',time1a,Ptot1a,'bs-',time1a,Pitot1a,'b*-',time1a,Drone1a*10,'m<-',time1a,Dronei1a*10,'m*-','Markers',10);
   h_l = legend('Egg','Larvae','Pupae','Infested Pupae','Drones x 10','Infested Drones x 10','Location','North');
   axis([30 750 0 16000])
   xlabel('Days','FontSize',18);
   ylabel('Number of bees','FontSize',18);
   title('Simulation assumes that mites do not affect bee mortality','FontSize',12)
   curYTick = get(gca,'YTick');
   YLabel = cellstr( num2str( curYTick(:), '%5.0f') );
   set(gca, 'YTickLabel', YLabel);
   %ytickformat('%.2f')
   %h_l = legend('Egg','Larvae','Hive','Foragers','Drones x 10','Mites','Location','Northeast');
   set(h_l,'FontSize',9,'fontweight','bold');
   legend boxoff
   highresjpg('fig1a.jpg')
   pause(2)
   
   plot(time1a,Htot1a,'r+-',time1a,Hitot1a,'r*-',time1a,Ftot1a,'gv-',time1a,Fitot1a,'g*-',time1a,Mite1a,'kh-','Markers',10);
   h_l = legend('Hive bees','Infested Hive bees','Foragers','Infested Foragers','Mites','Location','North');

   axis([30 750 0 25000])
   xlabel('Days','FontSize',18);
   ylabel('Number of bees','FontSize',18);
   title('Simulation assumes that mites do not affect bee mortality','FontSize',12)
   curYTick = get(gca,'YTick');
   YLabel = cellstr( num2str( curYTick(:), '%5.0f') );
   set(gca, 'YTickLabel', YLabel);
   %ytickformat('%.2f')
   %h_l = legend('Egg','Larvae','Hive','Foragers','Drones x 10','Mites','Location','Northeast');
   set(h_l,'FontSize',9,'fontweight','bold');
   legend boxoff
   highresjpg('fig1b.jpg')
   
end  %plot


