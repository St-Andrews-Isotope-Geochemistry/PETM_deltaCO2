clear

%%
number_of_samples = 1000;

d11B_data = importd11BData("./../../Data/Rae_2021_Boron_Data.xlsx","d11Bdata_byStudy");
d11B_sw_data = importd11BswData("./../../Data/Rae_2021_Boron_Data.xlsx","d11Bsw");

% Get Anagnostou d11Bsw
d11B_sw = d11B_sw_data(d11B_sw_data.age==53.2,:).d11Bsw;
d11B_sw_uncertainty = 0.1;

saturation_state = [5,12];
saturation_state_sampler = Geochemistry_Helpers.Sampler(1:0.01:15,"Flat",saturation_state,"latin_hypercube").normalise();
saturation_state_sampler.getSamples(number_of_samples).shuffle();

d11B_sw_sampler = Geochemistry_Helpers.Sampler(30:0.01:45,"Gaussian",[d11B_sw,d11B_sw_uncertainty],"latin_hypercube").normalise();
d11B_sw_sampler.getSamples(number_of_samples).shuffle();

% Age to Ma
d11B_data.age = d11B_data.age/1000;

% Get PETM data
petm_data = d11B_data(d11B_data.age<=56.5 & d11B_data.age>=55.55 & ~d11B_data.exclude,:);
for petm_index = 1:height(petm_data)
    petm_data.samplers(petm_index) = Geochemistry_Helpers.Sampler(10:0.01:20,"Gaussian",[petm_data.d11B(petm_index),petm_data.d11B_2SD(petm_index)/2],"latin_hypercube").normalise();
end
petm_data.samplers.getSamples(number_of_samples).shuffle();

core_401 = petm_data(petm_data.site=="401",:);
core_1209 = petm_data(petm_data.site=="1209" | petm_data.site=="1209B",:);

prePETM_401 = core_401(core_401.time=="LP",:);
PETM_401 = core_401(core_401.time=="PETM",:);

prePETM_1209 = core_1209(core_1209.time=="LP",:);
PETM_1209 = core_1209(core_1209.time=="PETM",:);

prePETM_401_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(prePETM_401.samplers.samples),"latin_hypercube");
PETM_401_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(PETM_401.samplers.samples),"latin_hypercube");

prePETM_1209_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(prePETM_1209.samplers.samples),"latin_hypercube");
PETM_1209_distribution = Geochemistry_Helpers.Sampler.fromSamples(10:0.01:20,mean(PETM_1209.samplers.samples),"latin_hypercube");

%%
temperature = [28.7,33.3;
               34.1,38.5];

prePETM_401_d11bco2 = BuCC.d11BCO2().create(number_of_samples);
prePETM_401_d11bco2.species_calibration.d11B_measured.assignToEach("value",prePETM_401_distribution.samples);

prePETM_401_d11bco2.boron.assignToAll("epsilon",27.2);
prePETM_401_d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

prePETM_401_d11bco2.carbonate_chemistry.assignToAll("units"," mol/kg");
prePETM_401_d11bco2.carbonate_chemistry.assignToAll("temperature",temperature(1,1));
prePETM_401_d11bco2.carbonate_chemistry.assignToAll("salinity",35);
prePETM_401_d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
prePETM_401_d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
prePETM_401_d11bco2.carbonate_chemistry.assignToAll("calcium",20);
prePETM_401_d11bco2.carbonate_chemistry.assignToAll("magnesium",30);

prePETM_401_d11bco2.carbonate_chemistry.assignToEach("saturation_state",saturation_state_sampler.samples);

myami = MyAMI.MyAMI("Precalculated",true);
prePETM_401_d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

prePETM_401_d11bco2.calculate();

% 
PETM_401_d11bco2 = BuCC.d11BCO2().create(number_of_samples);
PETM_401_d11bco2.species_calibration.d11B_measured.assignToEach("value",PETM_401_distribution.samples);

PETM_401_d11bco2.boron.assignToAll("epsilon",27.2);
PETM_401_d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

PETM_401_d11bco2.carbonate_chemistry.assignToAll("temperature",temperature(1,2));
PETM_401_d11bco2.carbonate_chemistry.assignToAll("salinity",35);
PETM_401_d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
PETM_401_d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
PETM_401_d11bco2.carbonate_chemistry.assignToAll("calcium",20);
PETM_401_d11bco2.carbonate_chemistry.assignToAll("magnesium",30);

PETM_401_d11bco2.carbonate_chemistry.assignToEach("alkalinity",prePETM_401_d11bco2.carbonate_chemistry.alkalinity);

% myami = MyAMI.MyAMI("Precalculated",true);
PETM_401_d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

PETM_401_d11bco2.calculate();

prePETM_401_bad = imag(prePETM_401_d11bco2.carbonate_chemistry.pH.pValue)~=0;
prePETM_401_d11bco2 = prePETM_401_d11bco2(~prePETM_401_bad);

prePETM_401_pH = Geochemistry_Helpers.Distribution.fromSamples([],prePETM_401_d11bco2.carbonate_chemistry.pH.pValue);

PETM_401_bad = imag(PETM_401_d11bco2.carbonate_chemistry.pH.pValue)~=0;
PETM_401_d11bco2 = PETM_401_d11bco2(~PETM_401_bad);
PETM_401_pH = Geochemistry_Helpers.Distribution.fromSamples([],PETM_401_d11bco2.carbonate_chemistry.pH.pValue);

pH_change_401 = Geochemistry_Helpers.Distribution.fromSamples([],PETM_401_d11bco2.carbonate_chemistry.pH.pValue-prePETM_401_d11bco2.carbonate_chemistry.pH.pValue);

change_401 = (PETM_401_d11bco2.carbonate_chemistry.atmospheric_co2.x - prePETM_401_d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
change_401_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:5000,change_401);

%%
saturation_state_sampler.shuffle();
d11B_sw_sampler.shuffle();

prePETM_1209_d11bco2 = BuCC.d11BCO2().create(number_of_samples);
prePETM_1209_d11bco2.species_calibration.d11B_measured.assignToEach("value",prePETM_1209_distribution.samples);

prePETM_1209_d11bco2.boron.assignToAll("epsilon",27.2);
prePETM_1209_d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("units"," mol/kg");
prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("temperature",temperature(2,1));
prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("salinity",35);
prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("calcium",20);
prePETM_1209_d11bco2.carbonate_chemistry.assignToAll("magnesium",30);

prePETM_1209_d11bco2.carbonate_chemistry.assignToEach("saturation_state",saturation_state_sampler.samples);

myami = MyAMI.MyAMI("Precalculated",true);
prePETM_1209_d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

prePETM_1209_d11bco2.calculate();

% 
PETM_1209_d11bco2 = BuCC.d11BCO2().create(number_of_samples);
PETM_1209_d11bco2.species_calibration.d11B_measured.assignToEach("value",PETM_1209_distribution.samples);

PETM_1209_d11bco2.boron.assignToAll("epsilon",27.2);
PETM_1209_d11bco2.boron.d11B_sw.assignToEach("value",d11B_sw_sampler.samples);

PETM_1209_d11bco2.carbonate_chemistry.assignToAll("temperature",temperature(2,2));
PETM_1209_d11bco2.carbonate_chemistry.assignToAll("salinity",35);
PETM_1209_d11bco2.carbonate_chemistry.assignToAll("oceanic_pressure",0);
PETM_1209_d11bco2.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
PETM_1209_d11bco2.carbonate_chemistry.assignToAll("calcium",20);
PETM_1209_d11bco2.carbonate_chemistry.assignToAll("magnesium",30);

PETM_1209_d11bco2.carbonate_chemistry.assignToEach("alkalinity",prePETM_1209_d11bco2.carbonate_chemistry.alkalinity);

% myami = MyAMI.MyAMI("Precalculated",true);
PETM_1209_d11bco2.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

PETM_1209_d11bco2.calculate();

change_1209 = (PETM_1209_d11bco2.carbonate_chemistry.atmospheric_co2.x - prePETM_1209_d11bco2.carbonate_chemistry.atmospheric_co2.x)*1e6;
change_1209_distribution = Geochemistry_Helpers.Distribution.fromSamples(0:10:5000,change_1209);

%%
figure(1);
clf
hold on
change_401_distribution.plot();
change_1209_distribution.plot();

