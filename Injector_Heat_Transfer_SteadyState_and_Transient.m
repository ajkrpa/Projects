% PROPELLANT FLOW   
clear 
clc 

% m_dot = Cd * A = Cd*sqrt(2 * d * deltaP) 
% deltaP (pressure drop): 20% of chamber pressure 
% d (density) 
% m_dot (mass flow rate) 

% VARIABLES
cP =  2757902.92;                         % Pa (400 psi) 
deltaP = 0.2 * cP ;                       % Pa 
dF = .99*789.3 + .01* 999.7;              % kg/m^3 - density Ethanol (99% + 1% water) 
dO = 1141;                                % kg/m^3 - density LOX
ratioOF = 1.20;                           % O/F ratio
m_dot = 1.886;                            % kg/s

Cd_F = 0.57;                              % Fuel Discharge Coefficient - Juno's 
Cd_O = 0.47;                              % Oxidizer Discharge Coefficient - Juno's 
a_ox = 33;                               
a_f = 27;
numOrfices = 30;                          % Number of orfice pairs 
LD = 5;                                   % Optimal range for LD ratio is 5-7, but 4 is minimum [NASA]

% CALCULATIONS 
velF = (Cd_F*sqrt((deltaP*2)/dF));                     % m/s 
velO = (Cd_O*sqrt((deltaP*2)/dO));                     % m/s 

m_dotF = (m_dot)/(1+ratioOF);                          % kg/s
m_dotO = (ratioOF * m_dot)/(1+ratioOF);                % kg/s
InjAreaF_M = (m_dotF/(Cd_F*sqrt(2*dF*deltaP)));
InjAreaO_M = (m_dotO/(Cd_O*sqrt(2*dO*deltaP)));

% Conversion Values 
InchesToMeters = 0.0254; 
MMToM = 0.001;

% Heat Transfer Values
% visc - dynamic viscosity, cp- specific heat capacity, k - thermal
% conductivity
% NOTE: values for ethanol are taken for 298.15 K
n_cold = 0.3;
n_hot = 0.4;
visc_F = (0.00107*.099+0.0008891*0.01);        % kg/m*s, N*s/m^2 - fluid kinematic viscosity
visc_O = 0.00000693;                           % kg/m*s, N*s/m^2 - fluid viscosity
cp_F = (2570*0.95+75*0.05)*0.001;              % KJ/kg*K 
cp_O = 1700*0.001;                             % KJ/kg*K        
k_F = (0.17*0.99+0.60*0.01);                   % W/m*K
k_O = 0.152;                                   % W/m*K
DhManifold_F = 4.25*InchesToMeters;            % m - hydraulic diameter, fuel manifold                    
DhManifold_O = 2*InchesToMeters;               % m - hydraulic diameter, ox manifold 
DhChamber = 4*InchesToMeters;                  % m - hydraulic diameter, chamber 
DhOrfices_F = 28.484*InchesToMeters;
DhOrfices_O = 2.651*InchesToMeters;

% Heat Transfer Coefficient Calculations
% Pr - Prandtl's Number, Re- Reynold's Number, Nu - Nusselt Number, h -
% heat transfer coefficient
Pr_F = (visc_F*cp_F)/k_F;                      % Prandtl's Number - Ethanol    
Pr_O = (visc_O*cp_O)/k_O;                      % Prandtl's Number - LOX    

Re_F_Cold = (dF*velF*DhManifold_F )/visc_F;       
Re_O_Cold = (dO*velO*DhManifold_O)/visc_O;

Re_F_Orfices = (dF*velF*DhOrfices_F)/visc_F;   
Re_O_Orfices = (dO*velO*DhOrfices_O)/visc_O;   

% plot(Dh_Orfices, Re_F_Orfices)
% hold on
% plot(Dh_Orfices, Re_O_Orfices)
% xlabel("Orfice Diameter")
% ylabel("Reynold's Number")


Nu_F_Cold = 0.023*(Re_F_Cold^0.8)*(Pr_F^n_hot);
Nu_O_Cold = 0.023*(Re_O_Cold^0.8)*(Pr_O^n_hot);

h_cold_F = (Nu_F_Cold*k_F)/DhManifold_F;           % heat transfer coefficient, fuel manifold (h,cold)
h_cold_O = (Nu_O_Cold*k_O)/DhManifold_O;           % heat transfer coefficient, oxidizer manifold (h,cold)


%% Heat Transfer 
Tc = 2877;                                       % K 
Tm = 298.15;                                     % K 
h_cold = 3117;                                   % W/m^2K  - avg of fuel&ox h coeff                       
h_hot = 6000;                                    % W/m^2K - assumed
 
k_inj = 16.3;                                    % thermal conductivity, 316 stainless steel 
T_ip = 0.125*InchesToMeters;                      % meters, injector face plate thickness

R_total = (1/h_cold)+(T_ip/k_inj)+(1/h_hot);          % m^2K/W - total thermal resistance
q_total = ((Tc-Tm)/R_total);                          % W/m^2 - total heat transfer rate 
R_cold = 1/h_cold;                                    % m^2K/W 
T_wall_cold = q_total*R_cold+Tm;                      % K - temp of injector face (manifold) 
T_wall_hot = q_total*((1/h_hot)+(T_ip/k_inj))+Tm;   % K - temp of injector face (chamber)


fprintf("Injector Face Thickness: %d in\n", T_ip/InchesToMeters);
fprintf("Injector Face (Max Temp): %.2f", T_wall_hot)

%% Transient Heat Transfer
% Note: taken from 'The FTCS Method with MATLAB code (lect #02) 

%VARIABLES 
alpha = 3.95*10^-6; %m^2/s
time = [0:.5:5]; %s
t = repmat(time,11,1);
dist = [0:0.0003175:0.003175];
x = repmat(dist,11,1)';
h = 6000; %W/m^2K
k = 16.3; 
Ti = 293.15; %K
Tinf = 3000; %K 

[rowst, colst] = size(t);
[rowsx, colsx] = size(x);

Tinj = ones(11,11);
for i=1:rowst 
    for j=1:colst 
        Tinj(i,j)=(erf(x(i,j)/2*sqrt(alpha*t(i,j)))+(exp(h*x(i,j)/k+h^2*alpha*t(i,j)/k^2)*erfc(x(i,j)/2*sqrt(alpha*t(i,j))+h/k*sqrt(alpha*t(i,j)))))*(Ti-Tinf)+Tinf;
    end 
end

surf(x,t,Tinj)
xlabel("Distance")
ylabel("Time")
zlabel("Temperature")

% Tinj=(erf(x/2*sqrt(alpha*t))+(exp(h*x/k+h^2*alpha*t/k^2)*erfc(x/2*sqrt(alpha*t)+h/k*sqrt(alpha*t))))*(Ti-Tinf)+Tinf;


