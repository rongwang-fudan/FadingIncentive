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

function [ duds ] = MIC( FFlux, L, dpo, dcoef, xy_damage, LR, covidyear, abtfrac, abtlen )

global  alpha elas theta2 econo0 EFco2 Egreen Fex Eland realtime carbonbudget20

T = size(realtime,1);

EFco20=EFco2;
carbonbudget20s=carbonbudget20;
Eland0 = Eland;
Fex0=Fex;
dcoef0=dcoef;
FFlux0=FFlux;
LR0=LR;
elas0=max(0.02,elas);

pset=ones(18,1);

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

duds=zeros(10,50);
U1971=(((1-0.2).* 23.3236./3.7585e+03*1000)^(1-elasmu)/(1-elasmu))*7.3339e+03/1000;
U2015=(((1-0.3).* 105.0582./7.3339e+03*1000)^(1-elasmu)/(1-elasmu))*7.3339e+03/1000;
deltaU=U2015-U1971;

t=[5:5:25];
for s=1:10
    for r=1:50
        A=zeros(1,5);
        A(1:5)=output_util(r,s:(s+4),10);
        [rs,ms,bs] = regression(t,A);
        duds(s,r)=-ms/(1-r*0.001)^(s*5); % discounting to the mitigation year
    end
end

duds=duds./deltaU;

end
