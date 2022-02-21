% Author: Rong Wang // contact rongwang@fudan.edu.cn //
% Date: 2021.7.5
function [ duyy ] = plot_dUdt( elasmu )

num1=21; % 21 for 0, 5, 10, ..., 100 percentiles; 11 for 0, 10, 20, ..., 100 percentiles
num2=18; % number of parameters
num3=num1*num2;
t=[5:5:25]; 
r=1; % discounting 0.015 0.02 0.025
duyy = zeros(18,12);
duds_decom=zeros(10*6,18*2);
duyy2=zeros(10,18*2);
duyy3=zeros(10,6*2);
duyy4=zeros(19,16);

load('D:\monte carlo1\output_obsm-mc1.dat','-mat');
load('D:\monte carlo1\parameters-mc1.dat','-mat');
load('output\population.dat','-mat');
U1971=(((1-parameters(21)).* output_obsm(1,13)./L(1,1)*1000)^(1-elasmu)/(1-elasmu))*L(45,1)/1000;
U2015=(((1-parameters(23)).* output_obsm(45,13)./L(45,1)*1000)^(1-elasmu)/(1-elasmu))*L(45,1)/1000;
deltaU=U2015-U1971;

% different discounting
load('output\montecarlo1_duds_rou010.dat','-mat');
duds=duds./deltaU;
duyy(15,1:10) = prctile(duds(num3:(num3+100000),1:10,10),[50],1); % du/ds by year
duyy(15,11) = prctile(mean(duds(num3:(num3+100000),1:3,10),2),[50],1); % du/ds in 2035-2045
duyy(15,12) = prctile(mean(duds(num3:(num3+100000),5:7,10),2),[50],1); % du/ds in 2055-2065
load('output\montecarlo1_duds_rou025.dat','-mat');
duds=duds./deltaU;
duyy(16,1:10) = prctile(duds(num3:(num3+100000),1:10,10),[50],1); % du/ds by year
duyy(16,11) = prctile(mean(duds(num3:(num3+100000),1:3,10),2),[50],1); % du/ds in 2035-2045
duyy(16,12) = prctile(mean(duds(num3:(num3+100000),5:7,10),2),[50],1); % du/ds in 2055-2065
load('output\montecarlo2_duds_rou010.dat','-mat');
duds=duds./deltaU;
duyy(17,1:10) = prctile(duds(num3:(num3+100000),1:10,10),[50],1); % du/ds by year
duyy(17,11) = prctile(mean(duds(num3:(num3+100000),1:3,10),2),[50],1); % du/ds in 2035-2045
duyy(17,12) = prctile(mean(duds(num3:(num3+100000),5:7,10),2),[50],1); % du/ds in 2055-2065
load('output\montecarlo2_duds_rou025.dat','-mat');
duds=duds./deltaU;
duyy(18,1:10) = prctile(duds(num3:(num3+100000),1:10,10),[50],1); % du/ds by year
duyy(18,11) = prctile(mean(duds(num3:(num3+100000),1:3,10),2),[50],1); % du/ds in 2035-2045
duyy(18,12) = prctile(mean(duds(num3:(num3+100000),5:7,10),2),[50],1); % du/ds in 2055-2065

% without being constrained by observations
load('output\montecarlo1_duds_rou015.dat','-mat');
duds=duds./deltaU;
duyy(1:5,1:10) = prctile(duds(num3:(num3+100000),1:10,10),[5 25 50 75 95],1); % du/ds by year
duyy(1:5,11) = prctile(mean(duds(num3:(num3+100000),1:3,10),2),[5 25 50 75 95],1); % du/ds in 2035-2045
duyy(1:5,12) = prctile(mean(duds(num3:(num3+100000),5:7,10),2),[5 25 50 75 95],1); % du/ds in 2055-2065
duyy(6,1:10) = mean(duds(num3:(num3+100000),1:10,10),1); % du/ds by year
duyy(6,11) = mean(duyy(6,1:3),2); % du/ds in 2035-2045
duyy(6,12) = mean(duyy(6,5:7),2); % du/ds in 2055-2065
duyy(7,1:10) = duds(11,1:10,10); % du/ds by year
duyy(7,11) = mean(duds(11,1:3,10),2); % du/ds in 2035-2045
duyy(7,12) = mean(duds(11,5:7,10),2); % du/ds in 2055-2065
for p=1:10
    if p<=9
        A2=prctile(mean(duds(num3:(num3+100000),1:3,p+1),2),[25 50 75],1)-prctile(mean(duds(num3:(num3+100000),1:3,p),2),[25 50 75],1);
        B2=prctile(mean(duds(num3:(num3+100000),5:7,p+1),2),[25 50 75],1)-prctile(mean(duds(num3:(num3+100000),5:7,p),2),[25 50 75],1);
    else
        A2=prctile(mean(duds(num3:(num3+100000),1:3,p),2),[25 50 75],1);
        B2=prctile(mean(duds(num3:(num3+100000),5:7,p),2),[25 50 75],1);
    end
    duyy3(p,1)=A2(2);
    duyy3(p,2)=A2(3)-A2(2);
    duyy3(p,3)=A2(2)-A2(1);
    duyy3(p,4)=B2(2);
    duyy3(p,5)=B2(3)-B2(2);
    duyy3(p,6)=B2(2)-B2(1);
end
duds_t = zeros(100,10); n2=0;
duds2035=zeros(10,18,21);
duds2060=zeros(10,18,21);
for mc=1:num3
    i=mod(mc-1,num1)+1; % percentile 1 - 21
    j=floor((mc-1)/num1)+1; % parameter 1 - 18
    n2=n2+1; duds_t(n2,1:10)=duds(mc,1:10,10);
    for p=1:9
        duds2035(p,j,i)=mean(duds(mc,1:3,p+1),2)-mean(duds(mc,1:3,p),2);
        duds2060(p,j,i)=mean(duds(mc,5:7,p+1),2)-mean(duds(mc,5:7,p),2);
        if p==9
            duds2035(10,j,i)=mean(duds(mc,1:3,p+1),2);
            duds2060(10,j,i)=mean(duds(mc,5:7,p+1),2);
        end
    end
end
A=prctile(duds2035,[5 50 95],3);
B=prctile(duds2060,[5 50 95],3);
for i=1:3
    duds_decom((i*10-9):(i*10),1:18)=A(:,:,i);
    duds_decom((i*10-9):(i*10),19:36)=B(:,:,i);
end
A2=mean(duds2035,3); duyy4(1:18,1)=A2(10,1:18);
duyy4(1:18,2)=A(10,1:18,2);
duyy4(1:18,3)=A(10,1:18,1);
duyy4(1:18,4)=A(10,1:18,3);
A2=mean(duds2060,3); duyy4(1:18,5)=A2(10,1:18);
duyy4(1:18,6)=B(10,1:18,2);
duyy4(1:18,7)=B(10,1:18,1);
duyy4(1:18,8)=B(10,1:18,3);
duyy4(19,1:8)=[duyy(6,11) duyy(3,11) duyy(1,11) duyy(5,11) duyy(6,12) duyy(3,12) duyy(1,12) duyy(5,12)];
for i=1:10
    for j=1:18
        duyy2(i,j)=(duds_decom(i+20,j)-duds_decom(i,j))/(duyy(5,11)-duyy(1,11));
        duyy2(i,j+18)=(duds_decom(i+20,j+18)-duds_decom(i,j+18))/(duyy(5,12)-duyy(1,12));
    end
end

% Constrained by observations
load('output\montecarlo2_duds_rou015.dat','-mat');
duds=duds./deltaU;
load('output\montecarlo2_mcbias.dat','-mat');
A=duds((num3+1):(num3+100000),:,:);
B=sum(mcbias((num3+1):(num3+100000),16:30),2);
idx=find(B<3);
duyy(8:12,1:10) = prctile(A(idx,1:10,10),[5 25 50 75 95],1); % du/ds by year
C=zeros(size(idx,1),10,2);
for i=1:size(idx,1)
    for p=1:9
        C(i,p,1) = mean(A(idx(i),1:3,p+1),2)-mean(A(idx(i),1:3,p),2);
        C(i,p,2) = mean(A(idx(i),5:7,p+1),2)-mean(A(idx(i),5:7,p),2);
        if p==9
            C(i,10,1)=mean(duds(i,1:3,p+1),2);
            C(i,10,2)=mean(duds(i,5:7,p+1),2);
        end
    end
end
duyy(8:12,11) = prctile(C(:,10,1),[5 25 50 75 95],1); % du/ds in 2035-2045
duyy(8:12,12) = prctile(C(:,10,2),[5 25 50 75 95],1); % for actions in different years
duyy(13,1:10) = mean(A(idx,1:10,10),1); % du/ds by year
duyy(13,11) = mean(duyy(13,1:3),2); % du/ds in 2035-2045
duyy(13,12) = mean(duyy(13,5:7),2); % du/ds in 2055-2065
duyy(14,1:10) = duds(11,1:10,10);
duyy(14,11) = mean(duds(11,1:3,10),2); % for actions in different years
duyy(14,12) = mean(duds(11,5:7,10),2); % for actions in different years
for p=1:10
    if p<=9
        A2=prctile(mean(duds(num3:(num3+100000),1:3,p+1),2),[25 50 75],1)-prctile(mean(duds(num3:(num3+100000),1:3,p),2),[25 50 75],1);
        B2=prctile(mean(duds(num3:(num3+100000),5:7,p+1),2),[25 50 75],1)-prctile(mean(duds(num3:(num3+100000),5:7,p),2),[25 50 75],1);
    else
        A2=prctile(mean(duds(num3:(num3+100000),1:3,p),2),[25 50 75],1);
        B2=prctile(mean(duds(num3:(num3+100000),5:7,p),2),[25 50 75],1);
    end
    duyy3(p,7)=A2(2);
    duyy3(p,8)=A2(3)-A2(2);
    duyy3(p,9)=A2(2)-A2(1);
    duyy3(p,10)=B2(2);
    duyy3(p,11)=B2(3)-B2(2);
    duyy3(p,12)=B2(2)-B2(1);
end
duds_t2=zeros(100,10); n4=0;
duds2035=zeros(10,18,21);
duds2060=zeros(10,18,21);
for mc=1:num3
    i=mod(mc-1,num1)+1; % percentile 1 - 21
    j=floor((mc-1)/num1)+1; % parameter 1 - 18
    n4=n4+1; duds_t2(n4,1:10)=duds(mc,1:10,10);
    for p=1:9
        duds2035(p,j,i)=mean(duds(mc,1:3,p+1),2)-mean(duds(mc,1:3,p),2);
        duds2060(p,j,i)=mean(duds(mc,5:7,p+1),2)-mean(duds(mc,5:7,p),2);
        if p==9
            duds2035(10,j,i)=mean(duds(mc,1:3,p+1),2);
            duds2060(10,j,i)=mean(duds(mc,5:7,p+1),2);
        end
    end
end
A=prctile(duds2035,[5 50 95],3);
B=prctile(duds2060,[5 50 95],3);
for i=1:3
    duds_decom((i*10+21):(i*10+30),1:18)=A(:,:,i);
    duds_decom((i*10+21):(i*10+30),19:36)=B(:,:,i);
end
A2=mean(duds2035,3); duyy4(1:18,9)=A2(10,1:18);
duyy4(1:18,10)=A(10,1:18,2);
duyy4(1:18,11)=A(10,1:18,1);
duyy4(1:18,12)=A(10,1:18,3);
A2=mean(duds2060,3); duyy4(1:18,13)=A2(10,1:18);
duyy4(1:18,14)=B(10,1:18,2);
duyy4(1:18,15)=B(10,1:18,1);
duyy4(1:18,16)=B(10,1:18,3);
duyy4(19,9:16)=[duyy(13,11) duyy(10,11) duyy(8,11) duyy(12,11) duyy(13,12) duyy(10,12) duyy(8,12) duyy(12,12)];

subplot(1,2,1);
for i=1:n2
    plot(duds_t(i,:),'LineStyle','-','LineWidth',0.1,'Color',[0.39 0.97 0.74]); hold on; % green - mitigation is good before 2025
end
plot(duyy(1,1:10),'LineStyle','--','LineWidth',0.1,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(2,1:10),'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(3,1:10),'LineStyle','--','LineWidth',3  ,'Color',[0.2 0.2 0.2]); hold on; % black
plot(duyy(4,1:10),'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(5,1:10),'LineStyle','--','LineWidth',0.1,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(6,1:10),'LineStyle','-','LineWidth',3,'Color',[0.2 0.2 0.2]); hold on; % black
axis([1 10 -0.7 0.3]);
% set(gca,'ytick',0:75:75); set(gca,'xtick',0:100:100);
subplot(1,2,2);
for i=1:n4
    plot(duds_t2(i,:),'LineStyle','-','LineWidth',0.1,'Color',[0.39 0.97 0.74]); hold on; % green - mitigation is good before 2025
end
plot(duyy(7,1:10),'LineStyle','--','LineWidth',0.1,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(8,1:10),'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(9,1:10),'LineStyle','--','LineWidth',3  ,'Color',[0.2 0.2 0.2]); hold on; % black
plot(duyy(10,1:10),'LineStyle','--','LineWidth',3,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(11,1:10),'LineStyle','--','LineWidth',0.1,'Color',[0.7 0.7 0.7]); hold on; % grey
plot(duyy(12,1:10),'LineStyle','-','LineWidth',3,'Color',[0.2 0.2 0.2]); hold on; % black
axis([1 10 -0.7 0.3]);
% set(gca,'ytick',0:75:75); set(gca,'xtick',0:100:100);

end