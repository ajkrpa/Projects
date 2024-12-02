% Code Name: Launch Vehicle Database
% Author: AJ Giza
% Email: GIZAA @my.erau.edu 
% Class: EGR115 - Section 15
% Date: 04/14/24

clear 
clc 

launchVehicleData = readcell('LV_Data.xlsx');
[nRows, nCols] = size(launchVehicleData);

fprintf('The following code allows you to compare potential launch vehicles thet fit selected mission criteria (mission and payload).\n\n')

fprintf('Select your mission orbit, either Low Earth Orbit or Geostationary Transfer Orbit.\n')
fprintf('Orbit types:\n\n')
fprintf('\t (1) Low Earth Orbit is between 200 km and 1600 km above the Earth.\n')
fprintf('\t (2) Geostationary Orbit (GEO) or Geosynchronous Orbit is 35,786 km above the Earth.\n\n')

% User selects orbit - LEO or GTO
orbit = input('Input mission orbit (LEO or GTO): ', "s");
while isempty(orbit) || ~strcmpi(orbit,'GTO') && ~strcmpi(orbit, 'LEO') %<SM:WHILE:GIZA >
    disp('You must select a mission in either LEO or GTO')
    orbit = input('Input mission orbit (LEO or GTO): ', "s");
end
 
% User selects either a random payload or enters a payload 
selectPayload = input('Would you like to select a random payload? (Yes or No) ', "s");
while isempty(selectPayload) && ~strcmpi(selectPayload, 'Yes') &&  ~strcmpi(selectPayload, 'yes') &&  ~strcmpi(selectPayload, 'No') &&  ~strcmpi(selectPayload, 'no')
    disp('Please enter either "Yes" or "No".')
    selectPayload = input('Would you like to select a random payload? (Yes or No) ', "s");
end 
if strcmpi(selectPayload, 'Yes') || strcmpi(selectPayload, 'yes')
    payload = randi(1000); %<SM:RANDOM:GIZA >
elseif strcmpi(selectPayload, 'No') || strcmpi(selectPayload, 'no')
    payload = input('Input mission payload (kg): ');
    while isempty(payload) || payload < 0 
        disp('Invalid value inputted for payload.')
        payload = input('Input mission payload (kg): ');
    end
end

% Mission requirements entered by user are displayed 
fprintf('\nMission Requirements: \n')
fprintf('\t Mission Orbit: %s \n', orbit)
fprintf('\t Mission Payload: %.2f kg \n\n', payload)

% Launch vehicle database (cell array) is traversed 
% New cell array with launch vehicles that match criteria is created 
if strcmpi(orbit, 'LEO') %<SM:IF:GIZA >
    potentialLaunchVehicles = [];  %<SM:FILTER:GIZA >
    for row = 1:nRows
        if  launchVehicleData{row,2} >= payload %<SM:ROP:GIZA>
            fitLaunchVehicle = launchVehicleData{row,:};
            potentialLaunchVehicles = [potentialLaunchVehicles;launchVehicleData(row,:)]; 
        end
    end
    % A table with potential launch vehicles is displayed & information 
    fprintf('The following launch vehicles fit the mission criteria: \n')
    vehicleName = char(potentialLaunchVehicles(:,1));
    payloadLEO = cell2mat(potentialLaunchVehicles(:,2));
    payloadGTO = cell2mat(potentialLaunchVehicles(:,3));
    thrust = cell2mat(potentialLaunchVehicles(:,4));
    manufacturer = char(potentialLaunchVehicles(:,5));
    country = char(potentialLaunchVehicles(:,6));
    launchSite = char(potentialLaunchVehicles(:,7));
    numStages = cell2mat(potentialLaunchVehicles(:,8));
    fuel = char(potentialLaunchVehicles(:,9));
    oxidizer = char(potentialLaunchVehicles(:,10));
    costPerLaunch = cell2mat(potentialLaunchVehicles(:,11));
  
    tableHeaders = {'Launch Vehicle', 'Payload to LEO (kg)', 'Payload to GTO (kg)', 'Thrust (lbs)', 'Manufacturer(s)', 'Country/Origin', 'Launch Site', 'Number of Stages', 'Fuel', 'Oxidizer', 'Cost per Launch ($USD)'};
    fprintf('\n')
    potentialLaunchVehiclesTable = table(vehicleName, payloadLEO, payloadGTO, thrust, manufacturer, country, launchSite, numStages, fuel, oxidizer, costPerLaunch, 'VariableNames', tableHeaders);
    disp(potentialLaunchVehiclesTable) 

elseif strcmpi(orbit, 'GTO')
    potentialLaunchVehicles = []; %<SM:AUG:GIZA > 
    for row = 1:size(launchVehicleData,1)
        if  launchVehicleData{row,3} >= payload
            fitLaunchVehicle = launchVehicleData{row,:};
            potentialLaunchVehicles = [potentialLaunchVehicles;launchVehicleData(row,:)]; % create new cell array 
        end
    end
    % Display table with potential launch vehicles 
    fprintf('The following launch vehicles fit the mission criteria: \n')
    vehicleName = char(potentialLaunchVehicles(:,1));
    payloadLEO = cell2mat(potentialLaunchVehicles(:,2));
    payloadGTO = cell2mat(potentialLaunchVehicles(:,3));
    thrust = cell2mat(potentialLaunchVehicles(:,4));
    manufacturer = char(potentialLaunchVehicles(:,5));
    country = char(potentialLaunchVehicles(:,6));
    launchSite = char(potentialLaunchVehicles(:,7));
    numStages = cell2mat(potentialLaunchVehicles(:,8));
    fuel = char(potentialLaunchVehicles(:,9));
    oxidizer = char(potentialLaunchVehicles(:,10));
    costPerLaunch = cell2mat(potentialLaunchVehicles(:,11));

    tableHeaders = {'Launch Vehicle', 'Payload to LEO (kg)', 'Payload to GTO (kg)', 'Thrust (lbs)', 'Manufacturer(s)', 'Country/Origin', 'Launch Site', 'Number of Stages', 'Fuel', 'Oxidizer', 'Cost per Launch ($USD)'};
    fprintf('\n')
    potentialLaunchVehiclesTable = table(vehicleName, payloadLEO, payloadGTO, thrust, manufacturer, country, launchSite, numStages, fuel, oxidizer, costPerLaunch, 'VariableNames', tableHeaders);
    disp(potentialLaunchVehiclesTable)
end 

seeData = input('Would you liked to see the data comparing the potential launch vehicles? (Yes or No): ', "s");
if strcmpi(selectPayload, 'Yes') || strcmpi(selectPayload, 'yes')
    clf
    % Bar graph #1 - launch vehicles/cost per launch
    subplot(1,2,1) 
    launchVehiclePrices = potentialLaunchVehicles(1:end,11); %<SM:SLICE:GIZA >
    launchVehiclePrices = cell2mat(launchVehiclePrices);
    bar(launchVehiclePrices) %<SM:VIEW:GIZA >
    title('Launch Vehicles vs. Cost Per Launch')
    vehicleNames = char(potentialLaunchVehicles(1:end,1));
    xticklabels(vehicleNames)
    ylabel('Cost Per Launch (Millions)')
    xlabel('Launch Vehicle')
    % Bar graph #2 - launch vehicles/lbs thrust
    subplot(1,2,2)
    thrustComparison = potentialLaunchVehicles(1:end, 4);
    thrustComparison = cell2mat(thrustComparison);
    bar(thrustComparison)
    title('Launch Vehicles vs. Max Thrust')
    vehicleNames = char(potentialLaunchVehicles(1:end,1));
    xticklabels(vehicleNames)
    ylabel('Thrust (lbs)')
    xlabel('Launch Vehicle')
    selectVehicle = input('\nAfter reviewing the data and table, please select a launch vehicle for your mission: ', "s");
    fprintf('\nYour selected launch vehicle is the %s. \n', selectVehicle)
elseif strcmpi(selectPayload, 'No') || strcmpi(selectPayload, 'no')
    selectVehicle = input('\nOkay. Please select a launch vehicle for your mission: ', "s");
    fprintf('\nYour selected launch vehicle is the %s. \n', selectVehicle)

end
fprintf('\nThank you for using the Launch Vehicle Database to select a launch vehicle for your mission.')
