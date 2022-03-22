clear
close all
%%
data_directory = "./../../Data/";
number_of_samples = 1000; % Number of statistical samples

%% Loading data
d11B_data = importd11BData(data_directory+"/Rae_2021_Boron_Data.xlsx","d11Bdata_byStudy");
d11B_sw_data = importd11BswData(data_directory+"/Rae_2021_Boron_Data.xlsx","d11Bsw");

temperature = readtable(data_directory+"/temperature.csv");
salinity = readtable(data_directory+"/salinity.csv");
gmst = readtable(data_directory+"/gmst.csv");

% Use Anagnostou d11Bsw
d11B_sw = 38.5;

% Sensitivity tests...
%d11B_sw = 38;

d11B_sw_uncertainty = 0.2; % Assume uncertainty

%% Specifying data
% Get Mg from Zeebe 2019 for 56 Ma
age = 56; %Ma
magnesium = (52.82-35).*exp(-(age-0)./12) + 35;
calcium = magnesium./2.2; % Use 2.2 to match BayMAG setting

saturation_state = [5,8]; % From JRae

%% Make samplers
d11B_sw_sampler = Geochemistry_Helpers.Sampler(30:0.01:45,"Gaussian",[d11B_sw,d11B_sw_uncertainty],"latin_hypercube").normalise();
d11B_sw_sampler.getSamples(number_of_samples).shuffle();

saturation_state_sampler = Geochemistry_Helpers.Sampler(1:0.01:15,"Flat",saturation_state,"latin_hypercube").normalise();
saturation_state_sampler.getSamples(number_of_samples).shuffle();

% Age to Ma
d11B_data.age = d11B_data.age./1000;

% Filter data to PETM relevant by age
petm_data = d11B_data(d11B_data.age<=56.5 & d11B_data.age>=55.55 & ~d11B_data.exclude,:);
% Create a d11B sampler for each datapoint
for petm_index = 1:height(petm_data)
    petm_data.samplers(petm_index) = Geochemistry_Helpers.Sampler(10:0.01:20,"Gaussian",[petm_data.d11B(petm_index),petm_data.d11B_2SD(petm_index)/2],"latin_hypercube").normalise();
end
% Draw samples and shuffle
petm_data.samplers.getSamples(number_of_samples).shuffle();

% Split into two tables for two sites
core_401.data = petm_data(petm_data.site=="401",:);
core_1209.data = petm_data(petm_data.site=="1209" | petm_data.site=="1209B",:);

% Split each site into prePETM and PETM
% Currently using DeepMIP boundaries
core_401.prePETM.data = core_401.data(core_401.data.time=="LP",:);
core_401.PETM.data = core_401.data(core_401.data.time=="PETM",:);

core_1209.prePETM.data = core_1209.data(core_1209.data.time=="LP",:);
core_1209.PETM.data = core_1209.data(core_1209.data.time=="PETM",:);

% Create a new sampler using the mean of existing d11B samples
core_401.prePETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_401.prePETM.data.samplers.samples),"latin_hypercube");
core_401.PETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_401.PETM.data.samplers.samples),"latin_hypercube");

core_1209.prePETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_1209.prePETM.data.samplers.samples),"latin_hypercube");
core_1209.PETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_1209.PETM.data.samplers.samples),"latin_hypercube");

%% Temperature
% Create distributions
core_401.prePETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:42,temperature.DSDP401PrePETM).normalise();
core_401.PETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:42,temperature.DSDP401PETM).normalise();
core_1209.prePETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:42,temperature.ODP1209PrePETM).normalise();
core_1209.PETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:42,temperature.ODP1209PETM).normalise();

% Create samplers
core_401.prePETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_401.prePETM.temperature_distribution,'latin_hypercube');
core_401.PETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_401.PETM.temperature_distribution,'latin_hypercube');
core_1209.prePETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_1209.prePETM.temperature_distribution,'latin_hypercube');
core_1209.PETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_1209.PETM.temperature_distribution,'latin_hypercube');

% Draw samples
core_401.prePETM.temperature_sampler.getSamples(number_of_samples).shuffle();
core_401.PETM.temperature_sampler.getSamples(number_of_samples).shuffle();
core_1209.prePETM.temperature_sampler.getSamples(number_of_samples).shuffle();
core_1209.PETM.temperature_sampler.getSamples(number_of_samples).shuffle();

%% Salinity
% Create distributions
core_401.prePETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.DSDP401PrePETM).normalise();
core_401.PETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.DSDP401PETM).normalise();
core_1209.prePETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.ODP1209PrePETM).normalise();
core_1209.PETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.ODP1209PETM).normalise();

% Create samplers
core_401.prePETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_401.prePETM.salinity_distribution,'latin_hypercube');
core_401.PETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_401.PETM.salinity_distribution,'latin_hypercube');
core_1209.prePETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_1209.prePETM.salinity_distribution,'latin_hypercube');
core_1209.PETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_1209.PETM.salinity_distribution,'latin_hypercube');

% Draw samples
core_401.prePETM.salinity_sampler.getSamples(number_of_samples).shuffle();
core_401.PETM.salinity_sampler.getSamples(number_of_samples).shuffle();
core_1209.prePETM.salinity_sampler.getSamples(number_of_samples).shuffle();
core_1209.PETM.salinity_sampler.getSamples(number_of_samples).shuffle();

%% MyAMI
% Create a MyAMI object
myami = MyAMI.MyAMI("Precalculated",true);

%% 401 d11B-CO2
% prePETM
core_401.prePETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples); % Create a number of d11B-CO2 objects
core_401.prePETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_401.prePETM.data_distribution.samples); % Default species calibration of borate

core_401.prePETM.d11bco2.boron.assignToAll("epsilon",27.2); % Use Klochko? 27.2
core_401.prePETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("units"," mol/kg"); % Specify units as saturation state is unitless
core_401.prePETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_401.prePETM.temperature_sampler.samples);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_401.prePETM.salinity_sampler.samples);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_401.prePETM.d11bco2.carbonate_chemistry.assignToEach("saturation_state",saturation_state_sampler.samples);

core_401.prePETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_401.prePETM.d11bco2.calculate();

% PETM
core_401.PETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples); % Create a number of d11B-CO2 objects
core_401.PETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_401.PETM.data_distribution.samples); % Default species calibration of borate

core_401.PETM.d11bco2.boron.assignToAll("epsilon",27.2); % Use Klochko?
core_401.PETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_401.PETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_401.PETM.temperature_sampler.samples);
core_401.PETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_401.PETM.salinity_sampler.samples);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_401.PETM.d11bco2.carbonate_chemistry.assignToEach("alkalinity",core_401.prePETM.d11bco2.carbonate_chemistry.alkalinity); % Assume alkalinity the same as prePETM

core_401.PETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_401.PETM.d11bco2.calculate();

% Remove any complex value results
core_401.is_real = imag(core_401.prePETM.d11bco2.carbonate_chemistry.pH.pValue)==0 & imag(core_401.PETM.d11bco2.carbonate_chemistry.pH.pValue)==0;

core_401.prePETM.d11bco2 = core_401.prePETM.d11bco2(core_401.is_real);
core_401.PETM.d11bco2 = core_401.PETM.d11bco2(core_401.is_real);

% Output
core_401.deltaCO2 = (core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x - core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
core_401.deltaCO2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:8000,core_401.deltaCO2);

core_401.deltaCO2_doublings = log2(core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6) - log2(core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6);
core_401.deltaCO2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:2,core_401.deltaCO2_doublings);

%% Shuffle samples
% Shuffle saturation state because each site should have a different value
saturation_state_sampler.shuffle();
% No need to shuffle d11B as this should be the same at each site
%d11B_sw_sampler.shuffle();

%% 1209
% prePETM
core_1209.prePETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples); % Create a number of d11B-CO2 objects
core_1209.prePETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_1209.prePETM.data_distribution.samples);
core_1209.prePETM.d11bco2.species_calibration.assignToAll("coefficients",[1,0.7]); % Species calibration as in spreadsheet of d11Bc-0.7

core_1209.prePETM.d11bco2.boron.assignToAll("epsilon",27.2); % Use Klochko
core_1209.prePETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_1209.prePETM.d11bco2.carbonate_chemistry.assignToAll("units"," mol/kg");
core_1209.prePETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_1209.prePETM.temperature_sampler.samples);
core_1209.prePETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_1209.prePETM.salinity_sampler.samples);
core_1209.prePETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_1209.prePETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_1209.prePETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_1209.prePETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_1209.prePETM.d11bco2.carbonate_chemistry.assignToEach("saturation_state",saturation_state_sampler.samples);

core_1209.prePETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_1209.prePETM.d11bco2.calculate();

% PETM
core_1209.PETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples); % Create a number of d11B-CO2 objects
core_1209.PETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_1209.PETM.data_distribution.samples);
core_1209.PETM.d11bco2.species_calibration.assignToAll("coefficients",[1,0.7]); % Species calibration as in spreadsheet of d11Bc-0.7

core_1209.PETM.d11bco2.boron.assignToAll("epsilon",27.2); % Use Klochko?
core_1209.PETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_1209.PETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_1209.PETM.temperature_sampler.samples);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_1209.PETM.salinity_sampler.samples);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_1209.PETM.d11bco2.carbonate_chemistry.assignToEach("alkalinity",core_1209.prePETM.d11bco2.carbonate_chemistry.alkalinity); % Assume alkalinity the same as prePETM

core_1209.PETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_1209.PETM.d11bco2.calculate();

% Remove any complex value results
core_1209.is_real = imag(core_1209.prePETM.d11bco2.carbonate_chemistry.pH.pValue)==0 & imag(core_1209.PETM.d11bco2.carbonate_chemistry.pH.pValue)==0;

core_1209.prePETM.d11bco2 = core_1209.prePETM.d11bco2(core_1209.is_real);
core_1209.PETM.d11bco2 = core_1209.PETM.d11bco2(core_1209.is_real);

% Output
core_1209.deltaCO2 = (core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x - core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
core_1209.deltaCO2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:8000,core_1209.deltaCO2);

core_1209.deltaCO2_doublings = log2(core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6) - log2(core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6);
core_1209.deltaCO2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:2,core_1209.deltaCO2_doublings);

%% Combine estimates
combined_deltaCO2.sampler = Geochemistry_Helpers.Sampler(0:10:8000,"Manual",core_401.deltaCO2_distribution.probabilities.*core_1209.deltaCO2_distribution.probabilities,"latin_hypercube").normalise();
combined_deltaCO2.sampler.getSamples(number_of_samples).shuffle();

% Assume doublings is what is consistent with both sites
combined_doublings.distribution = Geochemistry_Helpers.Distribution(0:0.01:2,"Manual",core_401.deltaCO2_doublings_distribution.probabilities.*core_1209.deltaCO2_doublings_distribution.probabilities).normalise();
combined_doublings.sampler = Geochemistry_Helpers.Sampler(combined_doublings.distribution,"latin_hypercube");
combined_doublings.sampler.getSamples(number_of_samples).shuffle();

%% Calculate climate sensitivity
%Create GMST distribution + samples
global_deltaTemperature.distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:10,gmst.GMSTPETM-gmst.GMSTprePETM).normalise();
global_deltaTemperature.sampler = Geochemistry_Helpers.Sampler(global_deltaTemperature.distribution,'latin_hypercube');
global_deltaTemperature.sampler.getSamples(number_of_samples).shuffle();

%Calculate CS
climate_sensitivity.values = global_deltaTemperature.sampler.samples./combined_doublings.sampler.samples;
climate_sensitivity.distribution = Geochemistry_Helpers.Distribution.fromSamples(3:0.05:14,climate_sensitivity.values);

%% Save output - CO2, Alkalinity, Saturation State, pH
% CO2
output.core401.prePETMCO2 = core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;
output.core401.PETMCO2 = core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;

output.core1209.prePETMCO2 = core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;
output.core1209.PETMCO2 = core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;

% Alkalinity
output.core401.prePETMAlk  = core_401.prePETM.d11bco2.carbonate_chemistry.alkalinity.*1e6;
output.core401.PETMAlk = core_401.PETM.d11bco2.carbonate_chemistry.alkalinity.*1e6;

output.core1209.prePETMAlk = core_1209.prePETM.d11bco2.carbonate_chemistry.alkalinity.*1e6;
output.core1209.PETMAlk = core_1209.PETM.d11bco2.carbonate_chemistry.alkalinity.*1e6;

% Saturation State
output.core401.prePETMO = core_401.prePETM.d11bco2.carbonate_chemistry.saturation_state;
output.core401.PETMO = core_401.PETM.d11bco2.carbonate_chemistry.saturation_state;

output.core1209.prePETMO = core_1209.prePETM.d11bco2.carbonate_chemistry.saturation_state;
output.core1209.PETMO = core_1209.PETM.d11bco2.carbonate_chemistry.saturation_state;

% pH
output.core401.prePETMP = core_401.prePETM.d11bco2.carbonate_chemistry.pH.pValue;
output.core401.PETMP = core_401.PETM.d11bco2.carbonate_chemistry.pH.pValue;

output.core1209.prePETMP = core_1209.prePETM.d11bco2.carbonate_chemistry.pH.pValue;
output.core1209.PETMP = core_1209.PETM.d11bco2.carbonate_chemistry.pH.pValue;

% get combined distribution for CO2 values
output.core401.CO2_distPre = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,output.core401.prePETMCO2);
output.core401.CO2_distP = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,output.core401.PETMCO2);

output.core1209.CO2_distPre = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,output.core1209.prePETMCO2);
output.core1209.CO2_distP = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,output.core1209.PETMCO2);

%prePETM CO2
output.globalCPre.distribution = Geochemistry_Helpers.Distribution(0:50:10000,...
    "Manual",output.core1209.CO2_distPre.probabilities.*output.core401.CO2_distPre.probabilities).normalise();
output.globalCPre.sampler = Geochemistry_Helpers.Sampler(output.globalCPre.distribution,'latin_hypercube');
output.globalCPre.sampler.getSamples(number_of_samples).shuffle();

%PETM CO2
output.globalCP.distribution = Geochemistry_Helpers.Distribution(0:50:10000,"Manual",output.core1209.CO2_distP.probabilities.*output.core401.CO2_distP.probabilities).normalise();
output.globalCP.sampler = Geochemistry_Helpers.Sampler(output.globalCP.distribution,'latin_hypercube');
output.globalCP.sampler.getSamples(number_of_samples).shuffle();

% Save combined CO2 estimates
output.PETMCO2 = output.globalCP.sampler.samples;
output.prePETMCO2 = output.globalCPre.sampler.samples;

% Temperature
output.combined_doublings = combined_doublings;
output.globalT = global_deltaTemperature;
output.ClimateSens = climate_sensitivity;
output.combined_delCO2 = combined_deltaCO2;

save('petmCSMain.mat','-struct','output');

%% Make figures
base_colour = Geochemistry_Helpers.Colour.Colour("Goldenrod","ryb",0);
colourmap = base_colour.makePalette("triad");

figure(1);
clf
hold on
plot_handles(1) = core_401.deltaCO2_doublings_distribution.plot('Color',colourmap.colours(1).rgb,'LineWidth',1.2);
plot_handles(2) = core_1209.deltaCO2_doublings_distribution.plot('Color',colourmap.colours(2).rgb,'LineWidth',1.2);
plot_handles(3) = combined_doublings.distribution.plot('Color',colourmap.colours(3).rgb,'LineWidth',1.5);

xlabel("CO_2 doublings");
ylabel("Probability");

legend_handle = legend(["401","1209","Combined"],'Box','Off');

figure(2); 
clf
climate_sensitivity.distribution.plot('color','k','linewidth',1)
