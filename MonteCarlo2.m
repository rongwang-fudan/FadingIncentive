% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.31
% Assessing the contribution of different parameters to the uncertainty
% 1	Economic damage of 1C warming	dcoef	-
% 2	Equilibrium sensitivity of climate	ESC	K/(W/m2)
% 3	Inertia of atmospheric surface temperature	INT	yr
% 4	Industrial CO2 emissions	indemi	-
% 5	LUC emissions	lucemi	-
% 6	C flux from air to terrestrial biosphere	Flux14	GtC/yr
% 7	C flux from air to surface ocean	Flux15	GtC/yr
% 8	C flux from terrestrial biosphere to soil	Flux42	GtC/yr
% 9	Turnover time of C in soil	LTsoil	yr
% 10	Turnover time of C in deep ocean	LTocean	yr
% 11	Radiative forcing of methan	rf_ch4	-
% 12	Radiative forcing of nitrogen oxide	rf_n2o	-
% 13	Radiative forcing of CFCs	rf_cfc	-
% 14	Radiative forcing of aerosols	rf_aer	-
% 15	Elasticity of substitution between energy and non-energy inputs	elas	-
% 16	Sensitivity of induced EUE change rate to the share of energy expenditure in total costs	ss_eue	-
% 17	Sensitivity of induced ENE change rate to the share of energy expenditure in total costs	ss_ene	-
% 18	Learning rate	LR	-

function [ para_range2 ] = MonteCarlo2( FFlux, L, dpo, dcoef, xy_damage, LR, covidyear, abtfrac, abtlen, unc, modes )

% modes: 0 for monte carlo + sensitivty; 1-18 for sensitivty to one parameter; 19 for sensitivty to all parameters

global  alpha elas theta2 econo0 clim0 EFco2 Egreen Fex Eland realtime carbonbudget20

num1=21; % 21 for 0, 5, 10, ..., 100 percentiles; 11 for 0, 10, 20, ..., 100 percentiles
num2=18; % number of parameters
num3=num1*num2;
T = size(realtime,1);

EFco20=EFco2;
carbonbudget20s=carbonbudget20;
Eland0 = Eland;
Fex0=Fex;
dcoef0=dcoef;
FFlux0=FFlux;
LR0=LR;
elas0=max(0.02,elas);

para_range2 = zeros(num2,num1); % range of parameters

for mc=1:(num3+10000)
    
% if mod(mc,500)==1
    display(mc);
% end

pset=ones(1,num2);
if mc>num3
    % 1, generate a random value
    if modes~=0
        continue;
    end
    for j=1:num2
        pset(j) = min(2, max(0.1,randvar(unc(2:5,j),1,0)));
    end
else
    % 2, select the percentile i for parameter j
    i=mod(mc-1,num1) * 100 / (num1-1); % percentile
    j=floor((mc-1)/num1)+1; % parameter
    if modes~=j && modes~=(num2+1) && modes~=0
        continue;
    end
    pset(j) = randvar(unc(2:5,j),2,i);
end
pset=pset.*unc(1,:);
pset(16)=min(2,max(0,pset(16)));
pset(17)=min(2,max(0,pset(17)));

%Coefficient in the damage function
dcoef = dcoef0 * pset(1);
%Elasticity of substitution (avoid a zero-like elas)
elas = elas0 * max(0.04/elas0, pset(15));
%Learning rate on the cost curve
LR = LR0 * pset(18);

%Equilibrium sensitivity of climate
FFlux(1) = FFlux0(1) * pset(2);
%Time inertia of climate system to reach equilibirum (year)
FFlux(2) = FFlux0(2) * pset(3);
%air to land biosphere GtC/yr
FFlux(4) = FFlux0(4) * pset(6);
%air to surface ocean GtC/yr
FFlux(5) = FFlux0(5) * pset(7);
%land biosphere to soil GtC/yr
FFlux(6) = FFlux0(6) * pset(8);
%surface soil to deep soil GtC/yr
FFlux(7) = FFlux0(7) / pset(9);
%surface ocean to deep ocean GtC/yr
FFlux(8) = FFlux0(8) / pset(10);

%CO2 emission factors for fossil fuel only tCO2 / MJ
EFco2 = EFco20 * pset(4);
carbonbudget20(:,2) = carbonbudget20s(:,2) * pset(4);
%CO2 emissions from land use change
Eland = Eland0 * pset(5);
%Radiative forcing by 1 CH4, 2 N2O, 3 CFCs, 4 aerosol
for i=1:4
    Fex(:,i) = Fex0(:,i) * pset(i+10);
end

%Rates of induced efficiency changes
[iec, output_iec, xy_iec] = Calibration_IEC( L );
eue0 = mean(xy_iec(:,1),1);
epe0 = mean(xy_iec(:,2),1);
ene0 = mean(xy_iec(:,3),1);
for i=1:71
    iec(1,i) = eue0 + (iec(1,i)-eue0) * pset(16);
    iec(2,i) = epe0 + (iec(2,i)-epe0) * pset(16);
    iec(4,i) = ene0 + (iec(4,i)-ene0) * pset(17);
end
for i=1:26
    xy_iec(i,6) = eue0 + (xy_iec(i,6)-eue0) * pset(16);
    xy_iec(i,7) = epe0 + (xy_iec(i,7)-epe0) * pset(16);
    xy_iec(i,8) = ene0 + (xy_iec(i,8)-ene0) * pset(17);
end

%Initial population (millions)
L0 = L(1,1);
%Initial level of total factor productivity
A0 = econo0(13) / (econo0(10)^alpha) / (L0/1000)^(1-alpha);
%Energy use efficiency $ / KJ
econo0(1) = econo0(15)^(elas/(elas-1)) / (econo0(12)/econo0(13));
%Energy production efficiency PJ / (trillion $)^0.3 / (billion cap)^0.7
econo0(2) = (A0^(elas-1) * econo0(15))^(1/(elas-1)) / econo0(1);
%Non-energy efficiency (trillion $)^0.7 / (billion cap)^0.7
econo0(3) = (A0^(elas-1) * (1-econo0(15)))^(1/(elas-1));
%Abatement cost as a percentage of GDP
econo0(5) = econo0(4) / theta2 * econo0(18)^theta2  * econo0(12) / econo0(13) * EFco2(1) / 1000;
%Industrial emissions (Gt CO2 per year)
econo0(20) = econo0(12) * EFco2(1,1) * (1-Egreen(1,8));

% switcher for  C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
switcher = ones(1,10);

%Calibration of climate damage function
[output_dam] = Calibration_DAM( dcoef, dpo, xy_damage );
% save('..\output\output_dam.dat','output_dam');

%Calibration of equilibrium sensitivity of climate
[output_esc] = Calibration_ESC( FFlux, 0 );

%Calibration of savings rate by capital, energy and ouput
[calrsav, output_cap] = Calibration_CAP( L, iec, dpo, dcoef, LR, switcher, 0 );

%Calibration of ENE reduction by COVID-19
[output_covid, deffs] = Calibration_COVID( FFlux, L, iec, calrsav, dpo, dcoef, LR, switcher, covidyear, 0 );

%Factor contributions
S8 = zeros(T,34,10,14);
for sce=1:14
    S = Factor_Attri( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, 2020+sce*5, abtfrac, abtlen, 3 );
    S8(:,:,:,sce) = S(:,:,4:13);
end

output_util=zeros(51,14,10); % 51 for growth rate of consumption which will be used to compute discounting
elasmu=1.45; % elasticity of marginal utility of consumption
U0 = ((1-calrsav(5)).* output_cap(45,29)./L(45,1)*1000)^(1-elasmu)/(1-elasmu); % utility in 2015

j=121; % 2025
while realtime(j,1)<=2300
    for p=1:10
        for s=1:14
            Us = (((1-calrsav(5)).* S8(j,7,p,s)./L(j,1)*1000)^(1-elasmu)/(1-elasmu)-U0)*L(j,1)/1000; % utility of action in 2025, 2030, 2035 ...
            for r=1:50
                disc=(1-r*0.001)^(realtime(j,1)-2020);                
                output_util(r,s,p) = output_util(r,s,p) + Us * disc;  % Sum of NPV of utility
            end
            if s==1
                output_util(51,s,p) = log(S8(126,7,p,s) / S8(61,7,p,s) )/10; % growth rate in consumption in 2025
            else
                output_util(51,s,p) = log(S8(121+s*5,7,p,s) / S8(111+s*5,7,p,s) )/10; % growth rate in consumption
            end
        end
    end
    j=j+1;
end

%Output parameters
parameters=zeros(56+71*3,1);
parameters(1,1) = dcoef;
parameters(2,1) = FFlux(1); % ESC
parameters(3,1) = FFlux(2); % Inertia time
parameters(4,1) = output_cap(45,20); % total emissions for 2015
parameters(5,1) = Eland(45,1); % total emissions for 2015
parameters(6:8,1) = FFlux(4:6);
parameters(9,1) = clim0(2)/FFlux(7);
parameters(10,1) = clim0(5)/FFlux(8);
parameters(11:14,1) = Fex(45,1:4);
parameters(15,1) =  elas;
parameters(16,1) = (xy_iec(1,6)-xy_iec(4,6))/(log(xy_iec(1,4))-log(xy_iec(4,4))); % (dEUE/EUE)/(dOmega/Omega)
parameters(17,1) = (xy_iec(1,8)-xy_iec(4,8))/(log(xy_iec(1,5))-log(xy_iec(4,5))); % (dENE/ENE)/(dOmega/Omega)
parameters(18,1) = LR;
parameters(19:23,1) = calrsav(1,1:5);
parameters(24:34,1) = deffs(1,1:11);
parameters(35:45,1) = deffs(2,1:11);
parameters(46:56,1) = deffs(3,1:11);
parameters(57:127,1) = iec(1,1:71); % EUE rate
parameters(128:198,1) = iec(2,1:71); % EPE rate
parameters(199:269,1) = iec(4,1:71); % ENE rate

output_obsm = zeros(45,15*2);
output_obsm(1:35,1)=output_dam(:,3); % damage obs
output_obsm(1:35,2)=output_dam(:,2); % damage model
output_obsm(1:45,3)=output_esc(:,15); % temperature obs
output_obsm(1:45,4)=output_esc(:,8); % temperature model
output_obsm(1:45,5)=output_esc(:,14); % atmco2 obs
output_obsm(1:45,6)=output_esc(:,12); % atmco2 model
output_obsm(1:45,7)=output_esc(:,16); % LandC obs
output_obsm(1:45,8)=output_esc(:,10); % LandC model
output_obsm(1:45,9)=output_esc(:,17); % OceanC obs
output_obsm(1:45,10)=output_esc(:,11); % OceanC model
output_obsm(1:45,11)=output_cap(:,27); % Capital obs
output_obsm(1:45,12)=output_cap(:,10); % Capital model
output_obsm(1:45,13)=output_cap(:,29); % GDP obs
output_obsm(1:45,14)=output_cap(:,13); % GDP model
output_obsm(1:45,15)=output_cap(:,28); % Energy obs
output_obsm(1:45,16)=output_cap(:,12); % Energy model
output_obsm(1:45,17)=output_cap(:,31); % Energy price obs
output_obsm(1:45,18)=output_cap(:,14); % Energy price model    
output_obsm(1:26,19:21)=xy_iec(1:26,1:3); % EUE/EPE/ENE rates
output_obsm(1:26,22:24)=xy_iec(1:26,6:8); % EUE/EPE/ENE rates
output_obsm(1:28,25)=output_covid(1:28,22); % Energy obs
output_obsm(1:28,26)=output_covid(1:28,12); % Energy model
output_obsm(1:28,27)=output_covid(1:28,23); % GDP obs
output_obsm(1:28,28)=output_covid(1:28,13); % GDP model
output_obsm(1:28,29)=output_covid(1:28,24); % Energy price obs
output_obsm(1:28,30)=output_covid(1:28,14); % Energy price model

save(strcat('D:\monte carlo2\parameters-mc',num2str(mc),'.dat'),'parameters');
save(strcat('D:\monte carlo2\output_obsm-mc',num2str(mc),'.dat'),'output_obsm');
save(strcat('D:\monte carlo2\output_util-mc',num2str(mc),'.dat'),'output_util');
save(strcat('D:\monte carlo2\output_dam-mc',num2str(mc),'.dat'),'output_dam');
save(strcat('D:\monte carlo2\output_esc-mc',num2str(mc),'.dat'),'output_esc');
save(strcat('D:\monte carlo2\output_cap-mc',num2str(mc),'.dat'),'output_cap');
save(strcat('D:\monte carlo2\output_covid-mc',num2str(mc),'.dat'),'output_covid');
save(strcat('D:\monte carlo2\S8-mc',num2str(mc),'.dat'),'S8');

if mc<=num3
    i=mod(mc-1,num1)+1; % percentile
    j=floor((mc-1)/num1)+1; % parameter
    para_range2(j,i) = parameters(j,1);
end

if mc==num3
    save('D:\monte carlo2\para_range2.dat','para_range2');
    save('D:\monte carlo2\population.dat','L');
end

end

end
