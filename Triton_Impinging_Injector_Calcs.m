%% PROPELLANT FLOW   
clear 
clc 

% m_dot = Cd * A = Cd*sqrt(2 * d * deltaP) 
% deltaP (pressure drop): 20% of chamber pressure 
% d (density) 
% m_dot (mass flow rate) 

% VARIABLES
cP = 2757902.92;                          % Pa (400 psi) 
deltaP = 0.2 * cP ;                       % Pa 
dF = .99*789.3 + .01* 999.7;              % kg/m^3 - density Ethanol (99% + 1% water) 
dO = 1141;                                % kg/m^3 - density LOX
ratioOF = 1.20;                           % O/F ratio
m_dot = 1.886;                            % kg/s

Cd_F = 0.57;                              % Fuel Discharge Coefficient - Juno's 
Cd_O = 0.47;                              % Oxidizer Discharge Coefficient - Juno's 
a_ox = 33;                               
a_f = 27;
numOrfices = 30;                         % Number of orfice pairs 
LD = 5;                                  % Optimal range for LD ratio is 5-7, but 4 is minimum [NASA]
ID = 5;                                  % Impingement Distance:Orfice Diameter Ratio 

% CALCULATIONS 
m_dotF = (m_dot)/(1+ratioOF);                    % mass flow F 
m_dotO = (ratioOF * m_dot)/(1+ratioOF);          % mass flow O

InjAreaF = (m_dotF/(Cd_F*sqrt(2*dF*deltaP)))*1000000;  % mm^2 [RPE]
InjAreaO = (m_dotO/(Cd_O*sqrt(2*dO*deltaP)))*1000000;  % mm^2 [RPE]
InjAreaF_M = (m_dotF/(Cd_F*sqrt(2*dF*deltaP)));
InjAreaO_M = (m_dotO/(Cd_O*sqrt(2*dO*deltaP)));
   
velF = Cd_F*sqrt((deltaP*2)/dF);             % velocity of fuel [RPE]
velO = Cd_O*sqrt((deltaP*2)/dO);             % velocity of oxidizer [RPE]
 
tanBeta = (m_dotO*velO*sind(a_ox)-m_dotF*velF*sind(a_f))/(m_dotO*velO*cosd(a_ox)+m_dotF*velF*cosd(a_f));
resultantAngle = atan((m_dotO*velO*sind(a_ox)-m_dotF*velF*sind(a_f))/(m_dotO*velO*cosd(a_ox)+m_dotF*velF*cosd(a_f))); %resultant angle

OrficeAreaF = InjAreaF/numOrfices;                % mm^2
OrficeAreaO = InjAreaO/numOrfices;                % mm^2

diaFOrficeMM = 2*sqrt(OrficeAreaF/pi);              % mm
diaOOrficeMM = 2*sqrt(OrficeAreaO/pi);              % mm 

diaFOrficeInches = diaFOrficeMM * 0.0393700787;     % in
diaOOrficeInches = diaOOrficeMM * 0.0393700787;     % in

manifoldVolume_F = InjAreaF*4;                     % Manifold Volume 
manifoldVolume_O = InjAreaO*4; 

manifoldHeight = manifoldVolume_F/(pi*4^2);

lengthOrficeF = diaFOrficeInches*LD;           %inches
lengthOrficeO = diaOOrficeInches*LD;           %inches
avgOrficeLength = (lengthOrficeF+lengthOrficeO)/2;

impingementDistance = ((diaFOrficeInches+diaOOrficeInches)/2)*ID;  %impingment distance calculated from calculated avg orfice diameter 
orficeDistance = impingementDistance*(tand(a_f)+tand(a_ox));
facePlateThickness = lengthOrficeF*cosd(a_f); %consider heat tranfer...







% Command Window Outputs 
fprintf("Chamber Pressure: %.2f Pascals\n", cP)
fprintf("OF ratio: %.2f\n", ratioOF)
fprintf("Mdot: %.3f kg/s\n", m_dot)
fprintf("Cd(fuel): %.2f\n", Cd_F)
fprintf("Cd(oxidizer): %.2f\n", Cd_O)
fprintf("L/D ratio: %d\n", LD)
fprintf("Orfice pairs: %d\n", numOrfices)


 