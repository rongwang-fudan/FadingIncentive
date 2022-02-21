% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.31
% Assessing the contribution of different factors to the mitigation of climate change
% Scenarios	Y0	s (y)	C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
% Zero mitigation	
% z0	2300	10	×	×	×	×	×	×	×	×	×	×	×
% zC1	2300	10	√	×	×	×	×	×	×	×	×	×	×
% zC2	2300	10	√	√	×	×	×	×	×	×	×	×	×
% Immediate mitigation (The result of i0 is identical to zC2)
% i0	2025	10	√	√	×	×	×	×	×	×	×	×	×
% iS1	2025	10	√	√	√	×	×	×	×	×	×	×	×
% iS2	2025	10	√	√	√	√	×	×	×	×	×	×	×
% iS3	2025	10	√	√	√	√	√	×	×	×	×	×	×
% iS4	2025	10	√	√	√	√	√	√	×	×	×	×	×
% iS5	2025	10	√	√	√	√	√	√	√	×	×	×	×
% iT1	2025	10	√	√	√	√	√	√	√	√	×	×	×
% iT2	2025	10	√	√	√	√	√	√	√	√	√	×	×
% iT3	2025	10	√	√	√	√	√	√	√	√	√	×	×
% iT4	2025	10	√	√	√	√	√	√	√	√	√	√	√

function [ S ] = Factor_Attri( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, att )
% att: 1 calculate output_att; 2 not calculate output_att

global  realtime
T = size(realtime,1);

%outputs
S = zeros(T,34,13);

% switcher for  C1	C2	S1	S2	S3	S4	S5	T1	T2	T3	T4
switcher = zeros(1,10);
if att~=3
    S(:,:,1) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, 2300, abtfrac, abtlen, switcher );
end

switcher(1) = 1; % C1: direct economic damage of climate change
if att~=3
    S(:,:,2) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, 2300, abtfrac, abtlen, switcher );
end

switcher(2) = 1; % C2: indirect economic damage due to reduced non-energy services by climate change
if att~=3
    S(:,:,3) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, 2300, abtfrac, abtlen, switcher );
end

% Shift from no mitigation to mitigation
S(:,:,4) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

switcher(3) = 1; % S1: reduced energy use due to structural change under a carbon price
S(:,:,5) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

switcher(4) = 1; % S2: reduced energy-use-efficiency improvements due to structural change under a carbon price
S(:,:,6) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

switcher(5) = 1; % S3: reduced non-energy-efficiency improvements due to structural change under a carbon price
S(:,:,7) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

switcher(6) = 1; % S4: increased energy by using carbon-emission-free energy under a carbon price
S(:,:,8) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

switcher(7) = 1; switcher(8) = 1; % S5 + T1 + T2: reduced direct economic damage of climate change by using carbon-emission-free energy
S(:,:,11) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

% S5: reduced non-energy loss due to climate change by using carbon-emission-free energy
S(:,:,9) = S(:,:,11);
S(:,7,9) = S(:,7,9)./(1-S(:,5,9)-S(:,6,9)).*(1-S(:,5,8)-S(:,6,8));

% T1: direct costs of technological change from fossil fuel to carbon-emission-free energy
S(:,:,10) = S(:,:,11);
S(:,7,10) = S(:,7,10)./(1-S(:,5,10)-S(:,6,10)).*(1-S(:,5,10)-S(:,6,8));

switcher(9) = 1; % T3: reduced costs of carbon-emission-free energy through learning-by-doing
S(:,:,12) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

switcher(10) = 1; % T4: increased rate of learning by deployment programs in technological change 
S(:,:,13) = Abatement( FFlux, L, iec, calrsav, dpo, dcoef, LR, covidyear, deffs, abtyear, abtfrac, abtlen, switcher );

end
