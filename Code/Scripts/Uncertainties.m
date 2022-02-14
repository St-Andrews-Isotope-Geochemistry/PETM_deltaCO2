clear

%%
number_of_samples = 1000;
data_directory = "./../../Data/";

%% Loading data
d11B_data = importd11BData(data_directory+"/Rae_2021_Boron_Data.xlsx","d11Bdata_byStudy");
d11B_sw_data = importd11BswData(data_directory+"/Rae_2021_Boron_Data.xlsx","d11Bsw");

temperature = readtable(data_directory+"/temperature.csv");
salinity = readtable(data_directory+"/salinity.csv");

% Get Anagnostou d11Bsw
d11B_sw = d11B_sw_data(d11B_sw_data.age==53.2,:).d11Bsw;
d11B_sw_uncertainty = 0.1;

%% Specifying data
calcium = 0.017566;
magnesium = 2.2*calcium;

saturation_state = [5,8];

%% Make samplers
d11B_sw_sampler = Geochemistry_Helpers.Sampler(30:0.01:45,"Gaussian",[d11B_sw,d11B_sw_uncertainty],"latin_hypercube").normalise();
d11B_sw_sampler.getSamples(number_of_samples).shuffle();

saturation_state_sampler = Geochemistry_Helpers.Sampler(1:0.01:15,"Flat",saturation_state,"latin_hypercube").normalise();
saturation_state_sampler.getSamples(number_of_samples).shuffle();

% Age to Ma
d11B_data.age = d11B_data.age/1000;

% Get PETM data
petm_data = d11B_data(d11B_data.age<=56.5 & d11B_data.age>=55.55 & ~d11B_data.exclude,:);
for petm_index = 1:height(petm_data)
    petm_data.samplers(petm_index) = Geochemistry_Helpers.Sampler(10:0.01:20,"Gaussian",[petm_data.d11B(petm_index),petm_data.d11B_2SD(petm_index)/2],"latin_hypercube").normalise();
end
petm_data.samplers.getSamples(number_of_samples).shuffle();

core_401.data = petm_data(petm_data.site=="401",:);
core_1209.data = petm_data(petm_data.site=="1209" | petm_data.site=="1209B",:);

core_401.prePETM.data = core_401.data(core_401.data.time=="LP",:);
core_401.PETM.data = core_401.data(core_401.data.time=="PETM",:);

core_1209.prePETM.data = core_1209.data(core_1209.data.time=="LP",:);
core_1209.PETM.data = core_1209.data(core_1209.data.time=="PETM",:);

core_401.prePETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_401.prePETM.data.samplers.samples),"latin_hypercube");
core_401.PETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_401.PETM.data.samplers.samples),"latin_hypercube");

core_1209.prePETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_1209.prePETM.data.samplers.samples),"latin_hypercube");
core_1209.PETM.data_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(core_1209.PETM.data.samplers.samples),"latin_hypercube");

%%
core_401.prePETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,temperature.DSDP401PrePETM).normalise();
core_401.PETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,temperature.DSDP401PETM).normalise();
core_1209.prePETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,temperature.ODP1209PrePETM).normalise();
core_1209.PETM.temperature_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,temperature.ODP1209PETM).normalise();

core_401.prePETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_401.prePETM.temperature_distribution,'latin_hypercube');
core_401.PETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_401.PETM.temperature_distribution,'latin_hypercube');
core_1209.prePETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_1209.prePETM.temperature_distribution,'latin_hypercube');
core_1209.PETM.temperature_sampler = Geochemistry_Helpers.Sampler(core_1209.PETM.temperature_distribution,'latin_hypercube');

core_401.prePETM.temperature_sampler.getSamples(number_of_samples).shuffle();
core_401.PETM.temperature_sampler.getSamples(number_of_samples).shuffle();
core_1209.prePETM.temperature_sampler.getSamples(number_of_samples).shuffle();
core_1209.PETM.temperature_sampler.getSamples(number_of_samples).shuffle();

%%
core_401.prePETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.DSDP401PrePETM).normalise();
core_401.PETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.DSDP401PETM).normalise();
core_1209.prePETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.ODP1209PrePETM).normalise();
core_1209.PETM.salinity_distribution = Geochemistry_Helpers.Distribution.fromSamples(10:0.1:40,salinity.ODP1209PETM).normalise();

core_401.prePETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_401.prePETM.salinity_distribution,'latin_hypercube');
core_401.PETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_401.PETM.salinity_distribution,'latin_hypercube');
core_1209.prePETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_1209.prePETM.salinity_distribution,'latin_hypercube');
core_1209.PETM.salinity_sampler = Geochemistry_Helpers.Sampler(core_1209.PETM.salinity_distribution,'latin_hypercube');

core_401.prePETM.salinity_sampler.getSamples(number_of_samples).shuffle();
core_401.PETM.salinity_sampler.getSamples(number_of_samples).shuffle();
core_1209.prePETM.salinity_sampler.getSamples(number_of_samples).shuffle();
core_1209.PETM.salinity_sampler.getSamples(number_of_samples).shuffle();

%%
myami = MyAMI.MyAMI("Precalculated",true);

%%
core_401.prePETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples);
core_401.prePETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_401.prePETM.data_distribution.samples);

core_401.prePETM.d11bco2.boron.assignToAll("epsilon",27.2);
core_401.prePETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("units"," mol/kg");
core_401.prePETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_401.prePETM.temperature_sampler.samples);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_401.prePETM.salinity_sampler.samples);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_401.prePETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_401.prePETM.d11bco2.carbonate_chemistry.assignToEach("saturation_state",saturation_state_sampler.samples);

core_401.prePETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_401.prePETM.d11bco2.calculate();

% 
core_401.PETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples);
core_401.PETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_401.PETM.data_distribution.samples);

core_401.PETM.d11bco2.boron.assignToAll("epsilon",27.2);
core_401.PETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_401.PETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_401.PETM.temperature_sampler.samples);
core_401.PETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_401.PETM.salinity_sampler.samples);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_401.PETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_401.PETM.d11bco2.carbonate_chemistry.assignToEach("alkalinity",core_401.prePETM.d11bco2.carbonate_chemistry.alkalinity);

core_401.PETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_401.PETM.d11bco2.calculate();

% Output
core_401.deltaCO2 = (core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x - core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
core_401.deltaCO2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:5000,core_401.deltaCO2);

core_401.deltaCO2_doublings = log2(core_401.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6) - log2(core_401.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6);
core_401.deltaCO2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:2,core_401.deltaCO2_doublings);

%%

% Shuffle samples
saturation_state_sampler.shuffle();
d11B_sw_sampler.shuffle();

core_1209.prePETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples);
core_1209.prePETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_1209.prePETM.data_distribution.samples);
core_1209.prePETM.d11bco2.species_calibration.assignToAll("coefficients",[1,0.7]);

core_1209.prePETM.d11bco2.boron.assignToAll("epsilon",27.2);
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

% 
core_1209.PETM.d11bco2 = BuCC.d11BCO2().create(number_of_samples);
core_1209.PETM.d11bco2.species_calibration.d11B_measured.assignToEach("value",core_1209.PETM.data_distribution.samples);
core_1209.PETM.d11bco2.species_calibration.assignToAll("coefficients",[1,0.7]);

core_1209.PETM.d11bco2.boron.assignToAll("epsilon",27.2);
core_1209.PETM.d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

core_1209.PETM.d11bco2.carbonate_chemistry.assignToEach("temperature",core_1209.PETM.temperature_sampler.samples);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToEach("salinity",core_1209.PETM.salinity_sampler.samples);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("calcium",calcium);
core_1209.PETM.d11bco2.carbonate_chemistry.assignToAll("magnesium",magnesium);

core_1209.PETM.d11bco2.carbonate_chemistry.assignToEach("alkalinity",core_1209.prePETM.d11bco2.carbonate_chemistry.alkalinity);

core_1209.PETM.d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

core_1209.PETM.d11bco2.calculate();

core_1209.deltaCO2 = (core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x - core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
core_1209.deltaCO2_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:5000,core_1209.deltaCO2);

core_1209.deltaCO2_doublings = log2(core_1209.PETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6) - log2(core_1209.prePETM.d11bco2.carbonate_chemistry.atmospheric_co2.x*1e6);
core_1209.deltaCO2_doublings_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:0.01:2,core_1209.deltaCO2_doublings);

%%
combined_doublings = Geochemistry_Helpers.Distribution(0:0.01:2,"Manual",core_401.deltaCO2_doublings_distribution.probabilities.*core_1209.deltaCO2_doublings_distribution.probabilities).normalise();

%%
figure(1);
clf
hold on
core_401.deltaCO2_doublings_distribution.plot();
core_1209.deltaCO2_doublings_distribution.plot();
combined_doublings.plot();

xlabel("CO_2 doublings");
ylabel("Probability");

