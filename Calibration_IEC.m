% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.6.18

function [iec, output_iec, xy_iec] = Calibration_IEC ( L )
%   L(time,1); labor mill
%   iec:  rates of induced efficiency changes
%   output_iec: output of iec variables for model calibration
%   xy_iec: the rate of efficiency change agains omega

global alpha elas inputs econo0 realtime
%   inputs 45x6; 1 energy PWh; 2 capital trill $; 3 GDP trill $; 4 population mill; 5 energy price ($/kWh); 6 omega
%   alpha:  elasticity of output to capital
%   elas:   elasticity of substitution between energy and non-energy in output
%   inputs 45x6: 1 energy PWh; 2 capital trill $; 3 GDP trill $; 4 population mill; 5 energy price ($/kWh); 6 omega
%   econo:  economic variables over time
%     1 EUE; 2 EPE; 3 ENE; 4 backstop price $/tCO2
%     5 abatement cost as a percentage of GDP; 6 climate damage as a percentage of GDP; 7 net output
%     8 fraction of labor allocation to energy; 9 fraction of investment allocation to energy
%     10 total capital t$; 11 energy capital (trill $); 12 energy (PJ); 13 output (t$); 14 energy price $/kWh; 15 omega
%     16 fraction of energy investment allocated to carbon-emission-free energy
%     17 energy capital carbon-emission-free t$
%     18 fraction to abate CO2 emission; 19 carbon price $/tCO2; 20 CO2 emissions Gt CO2; 21 cumulative green energy PJ

lra=zeros(26,7);

for i=1:26
    x=[(1970+i):(1989+i)]; y=inputs(i:(i+19),5); [sR,lr_pe0,bb0] = regression(x,log(y'));
    y=inputs(i:(i+19),1)./inputs(i:(i+19),4); [sR,lr_e0,bb0] = regression(x,log(y')); % e=E/L
    y=inputs(i:(i+19),3)./inputs(i:(i+19),4); [sR,lr_y0,bb0] = regression(x,log(y')); % y=Y/L
    y=inputs(i:(i+19),6); avef0=mean(y,1); startf0=inputs(i,6); startpe0=inputs(i,5); % omega
    y=inputs(i:(i+19),7); [sR,lr_k0,bb0] = regression(x,log(y')); % k=K/L
    lr_se0 = lr_pe0 + lr_e0 - lr_y0;
    lr_B = lr_e0 - alpha*lr_k0;
    lr_A = lr_y0 - alpha*lr_k0;
    lr_taue0 = elas/(elas-1)*lr_se0 - lr_e0 + lr_y0;
    lr_we0 = lr_B - lr_se0;
    lr_b0 = (lr_A - avef0*(lr_we0 + lr_taue0))/(1-avef0);
    lra(i,1) = startf0;
    lra(i,2) = lr_taue0;
    lra(i,3) = lr_we0;
    lra(i,4) = lr_b0;
    lra(i,5) = 0; % damage of climate change
    lra(i,6) = lra(i,1);
    lra(i,7) = 1-lra(i,1);
    lra(i,8) = lr_se0;
end
xy_iec=zeros(26,9);
xy_iec(:,1)=lra(:,2); % eue rate
xy_iec(:,2)=lra(:,3); % epe rate
xy_iec(:,3)=lra(:,4); % ene rate
xy_iec(:,4)=lra(:,6); % omega
xy_iec(:,5)=lra(:,7); % 1 - omega
xy_iec(:,9)=lra(:,8); % omega rate

omegas=(-1.52:0.01:-0.82); sn=size(omegas,2);
omegas2=(-0.071:0.001:-0.013); sn2=size(omegas2,2);
iec=zeros(6,sn);
idx1=find(lra(:,7)<0.93);
idx2=find(lra(:,7)>0.92);

% EUE
x=log10(lra(:,6)); y=lra(:,2);
a=polyfit(x,y,1);
iec(1,1:sn)=polyval(a,omegas);
xy_iec(:,6)=polyval(a,x);
% [b2,bint2,r2,rint2,stats2]=regress(y,[ones(size(x,1),1) x]);
% euestd = (bint2(2,2)-bint2(2,1))/1.96/2/b2(2);

% EPE
y=lra(:,3);
a=polyfit(x,y,1);
iec(2,:)=polyval(a,omegas);
xy_iec(:,7)=polyval(a,x);
% [b2,bint2,r2,rint2,stats2]=regress(y,[ones(size(x,1),1) x]);
% epestd = (bint2(2,2)-bint2(2,1))/1.96/2/b2(2);

% ENE - omega
x=log10(lra(idx1,6)); y=lra(idx1,4);
a=polyfit(x,y,1);
a1=polyval(a,omegas);
x=log10(lra(idx2,6)); y=lra(idx2,4);
a=polyfit(x,y,1);
a2=max(a1,polyval(a,omegas));
iec(3,1:sn)=max(a1,a2);

% ENE - 1-omega
x=log10(lra(idx1,7)); y=lra(idx1,4);
a=polyfit(x,y,1);
a1=polyval(a,omegas2);
xy_iec(:,8)=polyval(a,log10(lra(:,7)));

x=log10(lra(idx2,7)); y=lra(idx2,4);
a=polyfit(x,y,1);
a2=polyval(a,omegas2);
iec(4,1:sn2)=max(a1,a2);
xy_iec(:,8)=max(xy_iec(:,8),polyval(a,log10(lra(:,7))));

% x=log10(lra(:,7)); y=lra(:,4);
% [b2,bint2,r2,rint2,stats2]=regress(y,[ones(size(x,1),1) x]);
% enestd = abs((bint2(2,2)-bint2(2,1))/1.96/2/b2(2));

iec(5,:)=omegas;
iec(6,1:sn2)=omegas2;

for j=1:sn2
    iec(4,j)=max(iec(4,j),0);
end

%Initial states
eue=econo0(1); %Energy use efficiency $ / KJ
epe=econo0(2); %Energy production efficiency PJ / (trillion $)^0.3 / (billion cap)^0.7
ene=econo0(3); %Non-energy efficiency (trillion $)^0.7 / (billion cap)^0.7
se=econo0(15); %Share of energy expenditure in GDP
pe=econo0(14); %Energy price $/kWh

%Calibration of the model
output_iec=zeros(45,10); % 1-5 for model; 6-10 for observations
output_iec(1,1:10) = [eue, epe, ene, se, pe, eue, epe, ene, se, pe];

i=1;
while realtime(i,1)<2015
    if i<=10
        omega = inputs(1,6)+(inputs(2,6)-inputs(1,6))*(i-11);
    else
        omega = inputs(i-10,6);
    end
    output_iec(i+1,1)= output_iec(i,1)*(1+iec(1,min(71,max(1,floor((log10(omega)+1.52)*100)+1))))^realtime(i,2);
    output_iec(i+1,2)= output_iec(i,2)*(1+iec(2,min(71,max(1,floor((log10(omega)+1.52)*100)+1))))^realtime(i,2);
    output_iec(i+1,3)= output_iec(i,3)*(1+iec(4,min(59,max(1,floor((log10(1-omega)+0.071)*1000)+1))))^realtime(i,2);
    output_iec(i+1,4)= 1/(1+(output_iec(i+1,3)/output_iec(i+1,1)/output_iec(i+1,2))^(elas-1)); % omega
    output_iec(i+1,5)= output_iec(i+1,1) / output_iec(i+1,4)^(1/(elas-1)) *3600; % energy price $/KJ -> $/kWh    
    %
    pe = inputs(i+1,5); % $/kWh
    se = inputs(i+1,6);
    K = inputs(i+1,2); % capital stock ($ trillion 2010 USD)
    E = inputs(i+1,1)*3600; % energy PJ
    Y = inputs(i+1,3); % gross output (trill 2010 USD)
    A = Y / (K^alpha) / (L(i+1,1)/1000)^(1-alpha); %Initial level of total factor productivity
    eue = se^(elas/(elas-1)) / (E/Y); %Energy use efficiency $ / KJ
    epe = (A^(elas-1) * se)^(1/(elas-1)) / eue; %Energy production efficiency PJ / (trillion $)^0.3 / (billion cap)^0.7
    ene = (A^(elas-1) * (1-se))^(1/(elas-1)); %Non-energy efficiency (t$)^0.7/(billion cap)^0.7        
    %
    output_iec(i+1,6)= eue;
    output_iec(i+1,7)= epe;
    output_iec(i+1,8)= ene;
    output_iec(i+1,9)= se;
    output_iec(i+1,10)= pe;    
    i=i+1;
end

output_iec(:,1) = output_iec(:,1) * 3600; % EUE $/KJ -> $/kWh
output_iec(:,6) = output_iec(:,6) * 3600; % EUE $/KJ -> $/kWh
output_iec(:,2) = output_iec(:,2) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7
output_iec(:,7) = output_iec(:,7) / 3600; % EPE PJ/(t$)^0.3 / (billion cap)^0.7 -> PWh/(t$)^0.3 / (billion cap)^0.7

end

%Regression coefficient
% k_EUE = 0.1357; % coefficient for EUE (slope to energy cost share)
% b_EUE = 0.1717; % coefficient for EUE (offset)
% k_EPE = 0.1985; % coefficient for EPE (slope to energy cost share)
% b_EPE = 0.2256; % coefficient for EPE (offset)
% k_ENE = -0.013; % coefficient for ENE (slope to energy cost share)
% b_ENE = - k_ENE * log10(0.1) + ENErate(Ts1,1); % coefficient for ENE (offset)


