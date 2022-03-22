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
% d11B_sw = d11B_sw_data(d11B_sw_data.age==53.2,:).d11Bsw;
d11B_sw = 38.5;
% sensitivity tests...
%d11B_sw = 38;
d11B_sw_uncertainty = 0.2; % Assume uncertainty

%% Specifying data
%get Mg from Zeebe 2019 for 56 Ma
yT = 56;
magnesium = (52.82-35).*exp(-(yT-0)./12) + 35;
calcium = magnesium./2.2; %use 2.2 to match BayMAG setting

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

% Output
core_401.deltaCO2 = (core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x - core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
%remove complex numbers
idx = core_401.deltaCO2 == real(core_401.deltaCO2);
core_401.deltaCO2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:8000,core_401.deltaCO2(idx));

core_401.deltaCO2_doublings = log2(core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6) - log2(core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6);
%remove complex numbers
idx = core_401.deltaCO2_doublings == real(core_401.deltaCO2_doublings);
core_401.deltaCO2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:2,core_401.deltaCO2_doublings(idx));

%% Shuffle samples
%shuffle saturation state b/c each site should have a different value
saturation_state_sampler.shuffle();
%no need to shuffle d11B as this should be the same at each site
%d11B_sw_sampler.shuffle();

%% 1209
% prePETM
core_1209.prePETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples); % Create a number of d11B-CO2 objects
core_1209.prePETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_1209.prePETM.data_distribution.samples);
core_1209.prePETM.d11bco2.species_calibration.assignToAll("coefficients",[1,0.7]); % Species calibration as in spreadsheet of d11Bc-0.7

core_1209.prePETM.d11bco2.boron.assignToAll("epsilon",27.2); % Use Klochko?
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

% Output
core_1209.deltaCO2 = (core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x - core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
% remove complex numbers
idx = core_1209.deltaCO2 == real(core_1209.deltaCO2);
core_1209.deltaCO2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:8000,core_1209.deltaCO2(idx));

core_1209.deltaCO2_doublings = log2(core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6) - log2(core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6);
% remove complex numbers
idx = core_1209.deltaCO2_doublings == real(core_1209.deltaCO2_doublings);
core_1209.deltaCO2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:2,core_1209.deltaCO2_doublings(idx));

%% Assume doublings is what is consistent with both sites
combined_doublings.distribution = Geochemistry_Helpers.Distribution(0:0.01:2,"Manual",core_401.deltaCO2_doublings_distribution.probabilities.*core_1209.deltaCO2_doublings_distribution.probabilities).normalise();
combined_doublings.sampler = Geochemistry_Helpers.Sampler(combined_doublings.distribution,'latin_hypercube');
combined_doublings.sampler.getSamples(number_of_samples).shuffle();

%but also save deltaCO2
core1209.deltaCO2 = core_1209.deltaCO2;
core401.deltaCO2 = core_401.deltaCO2;
%combined estimate
combined_delCO2.distribution = Geochemistry_Helpers.Distribution(0:10:8000,"Manual",core_401.deltaCO2_distribution.probabilities.*core_1209.deltaCO2_distribution.probabilities).normalise();
combined_delCO2.sampler = Geochemistry_Helpers.Sampler(combined_delCO2.distribution,'latin_hypercube');
combined_delCO2.sampler.getSamples(number_of_samples).shuffle();

%Create GMST distribution + samples
globalT.distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:10,gmst.GMSTPETM - gmst.GMSTprePETM).normalise();
globalT.sampler = Geochemistry_Helpers.Sampler(globalT.distribution,'latin_hypercube');
globalT.sampler.getSamples(number_of_samples).shuffle();

%Calculate CS
ClimateSens.values = globalT.sampler.samples ./ combined_doublings.sampler.samples;
ClimateSens.distribution = Geochemistry_Helpers.Distribution.fromSamples(3:.05:14,ClimateSens.values);

%save CO2, Alk, Saturation State, pH
core1209.PETMCO2 = core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;
core1209.prePETMCO2 = core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;
core401.PETMCO2 = core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;
core401.prePETMCO2 = core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6;
%alk
core1209.PETMAlk = core_1209.PETM.d11bco2.carbonate_chemistry.alkalinity.*1000000;
core1209.prePETMAlk = core_1209.prePETM.d11bco2.carbonate_chemistry.alkalinity.*1000000;
core401.PETMAlk = core_401.PETM.d11bco2.carbonate_chemistry.alkalinity.*1000000;
core401.prePETMAlk  = core_401.prePETM.d11bco2.carbonate_chemistry.alkalinity.*1000000;
%sat
core1209.PETMO = core_1209.PETM.d11bco2.carbonate_chemistry.saturation_state;
core1209.prePETMO = core_1209.prePETM.d11bco2.carbonate_chemistry.saturation_state;
core401.PETMO = core_401.PETM.d11bco2.carbonate_chemistry.saturation_state;
core401.prePETMO = core_401.prePETM.d11bco2.carbonate_chemistry.saturation_state;
%pH
core1209.PETMP = core_1209.PETM.d11bco2.carbonate_chemistry.pH.pValue;
core1209.prePETMP = core_1209.prePETM.d11bco2.carbonate_chemistry.pH.pValue;
core401.PETMP = core_401.PETM.d11bco2.carbonate_chemistry.pH.pValue;
core401.prePETMP = core_401.prePETM.d11bco2.carbonate_chemistry.pH.pValue;
%get combined distribution for CO2 values
idx = core1209.PETMCO2 == real(core1209.PETMCO2);
core_1209.CO2_distP = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,core1209.PETMCO2(idx));
core_1209.CO2_distPre = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,core1209.prePETMCO2(idx));
idx = core401.PETMCO2 == real(core401.PETMCO2);
core_401.CO2_distP = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,core401.PETMCO2(idx));
core_401.CO2_distPre = Geochemistry_Helpers.Distribution.fromSamples(0:50:10000,core401.prePETMCO2(idx));
%PETM CO2
globalCP.distribution = Geochemistry_Helpers.Distribution(0:50:10000,"Manual",core_1209.CO2_distP.probabilities.*core_401.CO2_distP.probabilities).normalise();
globalCP.sampler = Geochemistry_Helpers.Sampler(globalCP.distribution,'latin_hypercube');
globalCP.sampler.getSamples(number_of_samples).shuffle();
%prePETM CO2
globalCPre.distribution = Geochemistry_Helpers.Distribution(0:50:10000,...
    "Manual",core_1209.CO2_distPre.probabilities.*core_401.CO2_distPre.probabilities).normalise();
globalCPre.sampler = Geochemistry_Helpers.Sampler(globalCPre.distribution,'latin_hypercube');
globalCPre.sampler.getSamples(number_of_samples).shuffle();
%save combined CO2 estimates
PETMCO2 = globalCP.sampler.samples;
prePETMCO2 = globalCPre.sampler.samples;
%% Make a figure
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

figure(2); clf;
plot(ClimateSens.distribution.bin_midpoints,ClimateSens.distribution.probabilities,'color','k','linewidth',1);

save('petmCSMain.mat','combined_doublings','globalT','ClimateSens','core1209','core401','PETMCO2','prePETMCO2','combined_delCO2');