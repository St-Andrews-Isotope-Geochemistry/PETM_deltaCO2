clear

%%
d11B_data = importd11BData("./../../Data/Rae_2021_Boron_Data.xlsx","d11Bdata_byStudy");
d11B_sw_data = importd11BswData("./../../Data/Rae_2021_Boron_Data.xlsx","d11Bsw");

% Get Anagnostou d11Bsw
d11B_sw = d11B_sw_data(d11B_sw_data.age==53.2,:).d11Bsw;

% Age to Ma
d11B_data.age = d11B_data.age/1000;

% Get PETM data
petm_data = d11B_data(d11B_data.age<=56.5 & d11B_data.age>=55.55 & ~d11B_data.exclude,:);
% potential_datasets = string(categories(removecats(petm_data.ref)));

core_1209 = petm_data(petm_data.site=="1209" | petm_data.site=="1209B",:);
core_401 = petm_data(petm_data.site=="401",:);

prePETM_1209 = core_1209(core_1209.time=="LP",:);
PETM_1209 = core_1209(core_1209.time=="PETM",:);

prePETM_401 = core_401(core_401.time=="LP",:);
PETM_401 = core_401(core_401.time=="PETM",:);

%%
d11B = [mean(prePETM_401.d11B),mean(PETM_401.d11B);
        mean(prePETM_1209.d11B),mean(PETM_1209.d11B)];

temperature = [28.7,33.3;
               34.1,38.5];



%%
d11b_co2_401 = BuCC.d11BCO2().create(2);
d11b_co2_401.species_calibration.d11B_measured.assignToEach("value",d11B(1,:));

d11b_co2_401.boron.assignToAll("epsilon",27.2);
d11b_co2_401.boron.d11B_sw.assignToAll("value",d11B_sw);

d11b_co2_401.carbonate_chemistry.assignToEach("temperature",temperature(1,:));
d11b_co2_401.carbonate_chemistry.assignToAll("salinity",35);
d11b_co2_401.carbonate_chemistry.assignToAll("oceanic_pressure",0);
d11b_co2_401.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
d11b_co2_401.carbonate_chemistry.assignToAll("calcium",20);
d11b_co2_401.carbonate_chemistry.assignToAll("magnesium",30);

d11b_co2_401.carbonate_chemistry.assignToAll("alkalinity",2300);

myami = MyAMI.MyAMI("Precalculated",true);
d11b_co2_401.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

d11b_co2_401.calculate();



d11b_co2_1209 = BuCC.d11BCO2().create(2);
d11b_co2_1209.species_calibration.d11B_measured.assignToEach("value",d11B(2,:));

d11b_co2_1209.boron.assignToAll("epsilon",27.2);
d11b_co2_1209.boron.d11B_sw.assignToAll("value",d11B_sw);

d11b_co2_1209.carbonate_chemistry.assignToEach("temperature",temperature(2,:));
d11b_co2_1209.carbonate_chemistry.assignToAll("salinity",35);
d11b_co2_1209.carbonate_chemistry.assignToAll("oceanic_pressure",0);
d11b_co2_1209.carbonate_chemistry.assignToAll("atmospheric_pressure",1);
d11b_co2_1209.carbonate_chemistry.assignToAll("calcium",20);
d11b_co2_1209.carbonate_chemistry.assignToAll("magnesium",30);

d11b_co2_1209.carbonate_chemistry.assignToAll("alkalinity",2300);

myami = MyAMI.MyAMI("Precalculated",true);
d11b_co2_1209.carbonate_chemistry.equilibrium_coefficients.assignToAll("MyAMI",myami);

d11b_co2_1209.calculate();

co2_change = [diff(d11b_co2_401.carbonate_chemistry.atmospheric_co2.x);diff(d11b_co2_1209.carbonate_chemistry.atmospheric_co2.x)];

%%
figure(1);
clf
hold on
plot([mean(prePETM_401.age),mean(PETM_401.age)],d11b_co2_401.boron.pH.pValue,'x');
plot([mean(prePETM_1209.age),mean(PETM_1209.age)],d11b_co2_1209.boron.pH.pValue,'x');

set(gca,'XDir','Reverse');

xlabel("Age (Ma)");
ylabel("pH");
