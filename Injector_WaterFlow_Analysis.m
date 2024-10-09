clear 
clc 

% OBJECTIVES 
% Tank weights vs. time 
% 
% Characterize CdA values of fuel and oxidizer 
% Calculate average dp 
% Determine mass flow rate

% NOTES: 
   % For run 1, only the LOX tank weight decreases
   % For run 2, both tank weights decrease 

% Variables
densityWater = 997; %kg/m^3
junoCp= 3102750; %Pa 
OF = 1.3; 
deltaP = junoCp+(junoCp*.20); %Pa
mdot = 0.660; %kg/s 
AreaOx = 1.8322*10^-2;  %in^2
AreaF = 2.0414*10^-2;  % in^2

% DATA SETS
dracoDataRun2 = readtable("sep21_draco_waterflow_run2_r300.xlsx");
dracoDataRun1 = readtable("sep21_draco_waterflow_run1loxside_r300.xlsx");
 
%% Run 1
dracoDataRun1 = readtable("sep21_draco_waterflow_run1loxside_r300.xlsx");
% time = dracoDataRun1{:,1};

timeR1 = dracoDataRun1{399:441,1};                % Time (s) [65.01008s-71.96717s]       
FuelTankW1 = dracoDataRun1{399:441,2};            % Fuel tank weight (lbs) 
LOXTankW1 = dracoDataRun1{399:441,3};             % LOX Tank Weight 
manifoldPressureLOX = dracoDataRun1{399:441,12};  % LOX side manifold pressure (psi)

x = timeR1; 
z = -1.3388*x+144.1155;     % Curve Fit from curveFitter (mass flow = 1.3388 lb/s = 0.607269465 kg  kg/s)

figure(1)                                       % Tank Weights vs. Time 
plot(timeR1,FuelTankW1)
hold on 
plot(timeR1,LOXTankW1)

title('RUN 1: Tank Weights vs. Time')
xlabel('Time [s]')
ylabel('Tank Weight [lbf]')
legend('Fuel', 'Oxidizer')

figure(2)                                         % LOX Tank Weight vs. Time, Curve Fit (Mass Flow)
% Values from curveFitter 
% curveFitter                                     % [UNCOMMENT TO OPEN curveFitter]
x = timeR1; 
z = -1.3388*x+144.1155;                           % Curve Fit from curveFitter (mass flow = 1.3388 lb/s = 0.607269465 kg/s)
 
f_LOX1 = fit(timeR1,LOXTankW1,"rat23");   
plot(f_LOX1,timeR1,LOXTankW1, '-')
hold on 
plot(x,z)
title('LOX Tank Weight vs. Time [Mass Flow]')
xlabel('Time [s]')
ylabel('Tank Weight [lbf]')
legend('Fuel', 'Curve Fit')

figure(3)                                         % dp Downstream pressure (psi) - engine vs. Time
dp = dracoDataRun1{399:441,22};                   
Average_dp = mean(dracoDataRun1{399:441,22});
AvgdP = ones(size(timeR1))*Average_dp; 
plot(timeR1, dp)
hold on 
plot(timeR1, AvgdP) 
title('dP Pressure Engine vs. Time')
xlabel('Time [s]')
ylabel('dP Pressure Engine [psi]')
legend('dp', 'Average dp')

figure(4)                                         % LOX Manifold Pressure (psi) vs. Time 
LOXManifoldPressure = dracoDataRun1{399:441,12};
plot(timeR1, LOXManifoldPressure)
title('LOX Manifold Pressure vs. Time')
xlabel('Time [s]')
ylabel('LOX Manifold Pressure [psi]')
legend('LOX Manifold Pressure')


Average_dp_Pa = Average_dp*6894.75729;              % Pa 
mdotLOXRun1 = 0.607269465;                          % mass flow (kg/s) - from curveFitter
 
CdA = mdotLOXRun1/sqrt((2*densityWater*Average_dp_Pa));      % CdA (m^2)
CdA2 = CdA * 1550.0031;                                      % CdA (in^2)
Cd = CdA2/AreaOx;   


fprintf('\nRUN 1 DATA:\n\t')
fprintf('Mass Flow: %.4f kg/s\n\t',mdotLOXRun1)
fprintf("CdA: %.5f in^2\n\t", CdA2)
fprintf("Cd: %.5f \n\t", Cd)
%% Run 2 
dracoDataRun2 = readtable("sep21_draco_waterflow_run2_r300.xlsx");

timeR2 = dracoDataRun2{2786:3120,1};
FuelTankW2 = dracoDataRun2{2786:3120,2};
LOXTankW2 = dracoDataRun2{2786:3120,3};
manifoldPressureLOX2 = dracoDataRun2{2786:3120,22};

% timeR2 = dracoDataRun2{:,1};
% FuelTankW2 = dracoDataRun2{:,2};
% LOXTankW2 = dracoDataRun2{:,3};
% manifoldPressureLOX2 = dracoDataRun2{:,22};

x2 = timeR2; 
z_fuel = -0.3007*x+190.4960;      % from curveFitter (mass flow = 0.3007 lb/s = 0.136395226 kg/s)
z_ox = -0.3593*x+206.8065;        % from curveFitter (mass flow = 0.3593 lb/s = 0.162975739 kg/s)

figure(1)                               % Tank Weights vs. Time 
plot(timeR2, LOXTankW2)
hold on 
plot(timeR2, FuelTankW2)

title('RUN 2: Tank Weights vs. Time')
xlabel('Time [s]')
ylabel('Tank Weight [lbf]')
legend('Oxidizer', 'Fuel')

figure(2)                            % LOX Tank Weight vs. Time 
% curveFitter 
f_LOX2 = fit(timeR2,LOXTankW2,"rat23");
plot(f_LOX2,timeR2,LOXTankW2, '-')
hold on 
plot(x,z_ox)
title('LOX Tank Weight vs. Time [Mass Flow]')
xlabel('Time [s]')
ylabel('Tank Weight [lbf]')
legend('Fuel', 'Curve Fit')

figure(3)                         % Fuel Tank Weight vs. Time 
% curveFitter 
f_Fuel1 = fit(timeR2,FuelTankW2,"rat23");
plot(f_Fuel1,timeR2,FuelTankW2, '-')
hold on 
plot(x,z_fuel)
title('Fuel Tank Weight vs. Time [Mass Flow]')
xlabel('Time [s]')
ylabel('Tank Weight [lbf]')
legend('Fuel', 'Curve Fit')

figure(4)                                         % dp Downstream pressure (psi) - fuel vs. Time
dp_fuel = dracoDataRun2{2786:3120,22};                   
Average_dp_fuel = mean(dracoDataRun2{2786:3120,22});
AvgdP_fuel = ones(size(timeR2))*Average_dp_fuel; 
plot(timeR2, dp_fuel)
hold on 
plot(timeR2, AvgdP_fuel) 
title('dP Pressure Fuel vs. Time')
xlabel('Time [s]')
ylabel('dP Pressure Fuel [psi]')
legend('dp', 'Average dp')

figure(5)                                         % dp Downstream pressure (psi) - oxidizer vs. Time
dp_ox = dracoDataRun2{2786:3120,21};                   
Average_dp_ox = mean(dracoDataRun2{2786:3120,21});
AvgdP_ox = ones(size(timeR2))*Average_dp_ox; 
plot(timeR2, dp_ox)
hold on 
plot(timeR2, AvgdP_ox) 
title('dP Pressure Oxidizer vs. Time')
xlabel('Time [s]')
ylabel('dP Pressure Oxidizer [psi]')
legend('dp', 'Average dp')

figure(6)                                         % LOX Manifold Pressure (psi) vs. Time 
LOXManifoldPressure2 = dracoDataRun2{2786:3120,12};
plot(timeR2, LOXManifoldPressure2)
title('LOX Manifold Pressure vs. Time')
xlabel('Time [s]')
ylabel('LOX Manifold Pressure [psi]')
legend('LOX Manifold Pressure')

figure(7)                                         % Fuel Manifold Pressure (psi) vs. Time 
FuelManifoldPressure2 = dracoDataRun2{2786:3120,17};
plot(timeR2, FuelManifoldPressure2)
title('Fuel Manifold Pressure vs. Time')
xlabel('Time [s]')
ylabel('Fuel Manifold Pressure [psi]')
legend('Fuel Manifold Pressure')


mdotLOXRun2 = 0.162975739; % kg/s 
mdotFuelRun2 = 0.136395226; %kg/s 
 
Average_dp_ox_Pa = Average_dp_ox*6894.75729;                        % Pa 
Average_dp_fuel_Pa = Average_dp_fuel*6894.75729;                    % Pa 

CdA_LOX = mdotLOXRun2/sqrt((2*densityWater*Average_dp_ox_Pa));      % CdA (m^2)
CdA2_LOX = CdA_LOX * 1550.0031;                                     % CdA (in^2)
Cd_LOX = CdA2_LOX/AreaOx;   

CdA_Fuel = mdotFuelRun2/sqrt((2*densityWater*Average_dp_fuel_Pa));  % CdA (m^2)
CdA2_Fuel = CdA_Fuel * 1550.0031;                                   % CdA (in^2)
Cd_Fuel = CdA2_Fuel/AreaOx;   

fprintf('\nRUN 2 DATA:\n\t')
fprintf('Mass Flow: %.4f kg/s\n\t',mdotLOXRun2)
fprintf("CdA (Ox): %.5f in^2\n\t", CdA2_LOX)
fprintf("Cd (Ox): %.10f \n\t", Cd_LOX)
fprintf("CdA (Fuel): %.10f in^2\n\t", CdA2_Fuel)
fprintf("Cd (Fuel): %.10f \n\t", Cd_Fuel)
