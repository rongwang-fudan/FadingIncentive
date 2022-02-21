% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function [ S ] = Abatement( FFlux, L, iec, calrsav, dpower, dcoef, LR, covidyear, eff_covid, yabt, fabt, alen, switcher )

%   Damage function: dpower for power coefficient on temperature, dcoef for damage as a percentage of GDP for 1 degree warming
%   data from eia (available by country) https://www.eia.gov/international/data/world
%   iec:  regression coeffients of efficiencies against omega
%   1-2 for slope/offset; 1-2 for eue/epe - omega; 3 for ene - (1-omega); 4 for ene - omega
%   LR: learning rate
%   L: labor available
%   eff_covid: efficiency change in COVID-19
%   yabt: year to initiate mitigation
%   fabt: fraction of CO2 emission abatement
%   alen: time to reach the fraction of CO2 emission abatement
%   S0:  initital variables
%   switcher:  C1	C2	S1	S2	S3	S4	T1	T2

global  realtime Egreen S0
%   realtime:  time
%   ESC: equilibrium sensitivity of climate
%   calrsav: 1, EUE0; 2, ENE0; 3, saving rate 1971-2003; 4, saving rate 2004-2008; 5, saving rate 2009-2015
%   Egreen: renewable energy data quad Btu/yr: 1 total; 2 coal; 3 natural gas; 4 oil; 5 nuclear+renewable; 6 nuclear; 7 renewable

T = size(realtime,1);

% inertia of the adjustment of labor / investment allocation
inertia=2;
% month after the outbreak of COVID-19
covid=1;

% Run the simulation
t=1;
S=zeros(T,34);
S(1,1:end) = S0(1,1:end);
while realtime(t,1)<2301
    % using the simulation results before abatements
    if realtime(t,1)<min(yabt,2025)
        econ1 = S0(t,1:23);
        clim1 = S0(t,24:34);
        t=t+1;
        S(t,1:end) = S0(t,1:end);
        continue;
    end
    %
    if t<33
        rsav=calrsav(3);
    elseif t<38
        rsav=calrsav(4);
    else
        rsav=calrsav(5);
    end
    %Fraction of investment allocated to carbon-emission-free energy: S transition
    if floor(realtime(t+1,1))>=yabt
        Nabt=realtime(t+1,1)-yabt;
        fracinv=(Egreen(end,8)-fabt)*exp(-(Nabt^2)/2/alen/alen)+fabt;
    else
        fracinv=Egreen(end,8);
    end
    %Energy cost share (Omega) in the past 20 years
    omega = 0; tt=0;
    for i=1:t
        if (realtime(t,1)-realtime(i,1))<20
            omega=omega+S(i,15)*realtime(i,2);
            tt=tt+realtime(i,2);
        end
    end
    omega = omega/tt;
    %Investment for the previous calender year
    investment=0; tt2=0;
    nextyear=floor(realtime(t+1,1)+realtime(t+1,2)/2);
    for i=1:t
        if nextyear==(floor(realtime(i,1)+realtime(i,2)/2)+1)
            investment = investment + rsav * S(i,7) * realtime(i,2);
            tt2 = tt2 + realtime(i,2);
        end
    end
    investment=investment/tt2;
    %Change of efficiencies in COVID-19
    if realtime(t,1)>covidyear && covid<=11 && realtime(t,1)<(covidyear+1)
        deff=[1+eff_covid(1,covid),1+eff_covid(2,covid),1+eff_covid(3,covid)];
        covid=covid+1;
    else
        deff=[1,1,1];
    end
    %Evolution of state from time t to time (t+1)
    econ2 = econdyn(t, L(t), econ1, fracinv, iec, omega, investment, LR, clim1(8), dpower, dcoef, inertia, deff, switcher);
    clim2 = climdyn(t, clim1, FFlux, econ2(20));
    t=t+1;
    econ1 = econ2;
    clim1 = clim2;
    S(t,1:23)=econ2(1,1:23);
    S(t,24:34)=clim2(1,1:11);
end

S(:,1) = S(:,1) * 3600; % EUE $/KJ -> $/kWh
S(:,2) = S(:,2) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
S(:,12) = S(:,12) / 3600; % energy PJ -> PWh
S(:,21) = S(:,21) / 3600; % cumulative green energy PJ -> PWh

end

