% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
tic
clear;
clear global;

global  elas ESC INT FFlux

%Elasticity of substitution
elas = 0.4;
%Equilibrium sensitivity of climate K Calderia, Contribution of Sea Ice Response to Climate Sensitivity, 2014
ESC = 1.05;
%Time inertia of climate system to reach equilibirum (year)
INT = 53;
%Climate damage function: dpo for power coefficient on temperature, dcoef for damage as a percentage of GDP for 1 degree warming
dpo = 2;
%Data of climate damage: 1 (high damage, unconfirmed) 1999-Climate change policy_ quantifying uncertainties for damages and optimal carbon taxes; 2 (moderate, used by DICE) 2017-A Survey of Global Impacts of Climate Change: Replication, Survey Methods, and a Statistical Analysis
damagedata = 2;
%Learning rate on the cost curve
LR = 0.2;
%Peak population
LA = 11500;
%Year of COVID-19 outbreak
covidyear = 2020;
%year to initiate mitigation
abtyear = 2025;
%Ultimate fraction of CO2 emission abatement (1 for zero emissions)
abtfrac = 1;
%Time (years) to abate CO2 emissions
abtlen = 10;
% switcher for  C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
switcher = ones(1,10);

%Historical data of economy
Initialset;

%Historical data of climate
InitialsetC;

%Land-use change emissions GtC and radiative forcing for non-CO2
AerosolsLUC;

%Time series of population: L
L = population( LA );

%Calibration of climate damage function
[dcoef, xy_damage] = damage( dpo, damagedata );

%Calibration of induced efficiency change
[iec, output_iec, xy_iec] = Calibration_IEC( L );
% save('..\output\xy_iec.dat','xy_iec');

%Calibration of equilibrium sensitivity of climate
[output_esc] = Calibration_ESC( FFlux, 1 );

%Calibration of savings rate by capital, energy and ouput
[calrsav, output_cap] = Calibration_CAP( L, iec, dpo, dcoef, LR, switcher, 1 );

%Calibration of energy price during COVID-19
[oilprice1971_2019, oilprice2020_normalize, energyprice2020, SE, xy_price, co2emi_2019_2021] = energyprice( covidyear );

%Calibration of ENE reduction by COVID-19
[output_covid, deffs] = Calibration_COVID( FFlux, L, iec, calrsav, dpo, dcoef, LR, switcher, covidyear, 1 );

%Monte Carlo without observational constraints
[ para_range ] = MonteCarlo( FFlux, L, dpo, dcoef, xy_damage, LR, covidyear, abtfrac, abtlen, 0 );

%Constraining the parameters by observations
[ uncerts2 ] = Calibration_MC( [1	0.1	2	0.4	0.45	0.1	0.1	0.1	0.01	1	1	1	0.04	0.01	0.04] );

%Output results
save('D:\monte carlo1\uncerts2.dat','uncerts2');

clear
toc
